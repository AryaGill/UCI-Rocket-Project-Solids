/**

******************************************************************************

* @file : main.c

* @brief : LPS22HHTR Barometer Test with Live Expressions

******************************************************************************

*/

/* USER CODE END Header */


/* Includes ------------------------------------------------------------------*/

#include "main.h"

#include "fatfs.h"


/* Private includes ----------------------------------------------------------*/

/* USER CODE BEGIN Includes */

#include "sensors.h"

#include "telemetry.h"

#include <string.h>

#include <math.h>

/* USER CODE END Includes */


/* Private variables ---------------------------------------------------------*/

SPI_HandleTypeDef hspi1;

UART_HandleTypeDef huart3;


/* USER CODE BEGIN PV */


// Telemetry structure

Telemetry_t telemetry = {0};


// Debug variables for Live Expressions

volatile float debug_pressure = 0.0f; // hPa (mbar)

volatile float debug_temperature = 0.0f; // °C

volatile float debug_altitude = 0.0f; // meters (calculated)


// NEW: IMU Debug Variables

volatile float lsm_accel_x, lsm_accel_y, lsm_accel_z;

volatile float lsm_gyro_x, lsm_gyro_y, lsm_gyro_z;


volatile float icm_accel_x, icm_accel_y, icm_accel_z;

volatile float icm_gyro_x, icm_gyro_y, icm_gyro_z;


// WHO_AM_I and status

volatile uint8_t whoami = 0; // Should be 0xB3 for LPS22HH

volatile uint8_t lsm_whoami = 0; // Should be 0x6A

volatile uint8_t icm_whoami = 0;

volatile uint8_t status_reg = 0;

volatile uint8_t ctrl_reg1 = 0;


// Counter

volatile uint32_t read_counter = 0;


/* USER CODE END PV */


/* Private function prototypes -----------------------------------------------*/

void SystemClock_Config(void);

static void MPU_Config(void);

static void MX_GPIO_Init(void);

static void MX_SPI1_Init(void);

static void MX_USART3_UART_Init(void);

/* USER CODE BEGIN PFP */

float Calculate_Altitude(float pressure_hPa);

/* USER CODE END PFP */


/* USER CODE BEGIN 0 */


// Calculate altitude from pressure (standard atmosphere model)

float Calculate_Altitude(float pressure_hPa)

{

const float sea_level_pressure = 1013.25f; // hPa at sea level


// Barometric formula: h = 44330 * (1 - (P/P0)^(1/5.255))

float ratio = pressure_hPa / sea_level_pressure;

float altitude = 44330.0f * (1.0f - powf(ratio, 0.1903f));


return altitude;

}


/* USER CODE END 0 */


int main(void)

{

/* USER CODE BEGIN 1 */

/* USER CODE END 1 */


/* MPU Configuration--------------------------------------------------------*/

MPU_Config();


/* MCU Configuration--------------------------------------------------------*/

HAL_Init();


/* Configure the system clock */

SystemClock_Config();


/* USER CODE BEGIN SysInit */

/* USER CODE END SysInit */


/* Initialize all configured peripherals */

MX_GPIO_Init();

MX_SPI1_Init();

MX_USART3_UART_Init();

MX_FATFS_Init();


/* USER CODE BEGIN 2 */


// Startup LED blink

for(int i = 0; i < 3; i++) {

HAL_GPIO_WritePin(LED_GPIO_Port, LED_Pin, GPIO_PIN_SET);

HAL_Delay(100);

HAL_GPIO_WritePin(LED_GPIO_Port, LED_Pin, GPIO_PIN_RESET);

HAL_Delay(100);

}


// Initialize LPS22HHTR barometer

init_sensors(&hspi1);

HAL_Delay(100);


// Verify sensor communication

whoami = LPS22HH_WhoAmI(); // Should be 0xB3 (179 decimal)

lsm_whoami = LSM6DSL_WhoAmI(); // 0x6A

HAL_Delay(10);


// Read control register to verify configuration

ctrl_reg1 = LPS22HH_ReadReg(LPS22HH_CTRL_REG1);


/* USER CODE END 2 */


/* Infinite loop */

/* USER CODE BEGIN WHILE */

while (1)

{

// Read sensor data

read_sensors(&telemetry);


// Update debug variables for Live Expressions

debug_pressure = telemetry.pressure;

debug_temperature = telemetry.temperature;


// 3. Map IMU values (LSM6DSL)

lsm_accel_x = telemetry.lsm_accel_r; // mapped to roll axis

lsm_accel_y = telemetry.lsm_accel_p; // mapped to pitch axis

lsm_accel_z = telemetry.lsm_accel_y; // mapped to yaw axis


lsm_gyro_x = telemetry.lsm_gyro_r;

lsm_gyro_y = telemetry.lsm_gyro_p;

lsm_gyro_z = telemetry.lsm_gyro_y;


icm_accel_x = telemetry.icm_accel_r;

icm_accel_y = telemetry.icm_accel_p;

icm_accel_z = telemetry.icm_accel_y;


icm_gyro_x = telemetry.icm_gyro_r;

icm_gyro_y = telemetry.icm_gyro_p;

icm_gyro_z = telemetry.icm_gyro_y;


// Calculate altitude from pressure

debug_altitude = Calculate_Altitude(debug_pressure);


// Read status register periodically

if(read_counter % 10 == 0) {

status_reg = LPS22HH_ReadReg(LPS22HH_STATUS_REG);

whoami = LPS22HH_WhoAmI();

lsm_whoami = LSM6DSL_WhoAmI();

icm_whoami = ICM45686_WhoAmI();

}


read_counter++;


// Blink LED to show loop is running

HAL_GPIO_TogglePin(LED_GPIO_Port, LED_Pin);


HAL_Delay(100); // 10 Hz update rate


/* USER CODE END WHILE */


/* USER CODE BEGIN 3 */

}

/* USER CODE END 3 */

}


/**

* @brief System Clock Configuration

* @retval None

*/

void SystemClock_Config(void)

{

RCC_OscInitTypeDef RCC_OscInitStruct = {0};

RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};


HAL_PWREx_ConfigSupply(PWR_LDO_SUPPLY);

__HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE3);


while(!__HAL_PWR_GET_FLAG(PWR_FLAG_VOSRDY)) {}


RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI|RCC_OSCILLATORTYPE_HSE;

RCC_OscInitStruct.HSEState = RCC_HSE_ON;

RCC_OscInitStruct.HSIState = RCC_HSI_DIV1;

RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;

RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;

RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSE;

RCC_OscInitStruct.PLL.PLLM = 2;

RCC_OscInitStruct.PLL.PLLN = 12;

RCC_OscInitStruct.PLL.PLLP = 2;

RCC_OscInitStruct.PLL.PLLQ = 3;

RCC_OscInitStruct.PLL.PLLR = 2;

RCC_OscInitStruct.PLL.PLLRGE = RCC_PLL1VCIRANGE_3;

RCC_OscInitStruct.PLL.PLLVCOSEL = RCC_PLL1VCOMEDIUM;

RCC_OscInitStruct.PLL.PLLFRACN = 0;

if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)

{

Error_Handler();

}


RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK

|RCC_CLOCKTYPE_PCLK1|RCC_CLOCKTYPE_PCLK2

|RCC_CLOCKTYPE_D3PCLK1|RCC_CLOCKTYPE_D1PCLK1;

RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_HSI;

RCC_ClkInitStruct.SYSCLKDivider = RCC_SYSCLK_DIV1;

RCC_ClkInitStruct.AHBCLKDivider = RCC_HCLK_DIV1;

RCC_ClkInitStruct.APB3CLKDivider = RCC_APB3_DIV1;

RCC_ClkInitStruct.APB1CLKDivider = RCC_APB1_DIV2;

RCC_ClkInitStruct.APB2CLKDivider = RCC_APB2_DIV1;

RCC_ClkInitStruct.APB4CLKDivider = RCC_APB4_DIV1;


if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_1) != HAL_OK)

{

Error_Handler();

}

}


/**

* @brief SPI1 Initialization Function

* @param None

* @retval None

*/

static void MX_SPI1_Init(void)

{

hspi1.Instance = SPI1;

hspi1.Init.Mode = SPI_MODE_MASTER;

hspi1.Init.Direction = SPI_DIRECTION_2LINES;

hspi1.Init.DataSize = SPI_DATASIZE_8BIT;

hspi1.Init.CLKPolarity = SPI_POLARITY_HIGH;

hspi1.Init.CLKPhase = SPI_PHASE_2EDGE;

hspi1.Init.NSS = SPI_NSS_SOFT;

hspi1.Init.BaudRatePrescaler = SPI_BAUDRATEPRESCALER_16;

hspi1.Init.FirstBit = SPI_FIRSTBIT_MSB;

hspi1.Init.TIMode = SPI_TIMODE_DISABLE;

hspi1.Init.CRCCalculation = SPI_CRCCALCULATION_DISABLE;

hspi1.Init.CRCPolynomial = 0x0;

hspi1.Init.NSSPMode = SPI_NSS_PULSE_ENABLE;

hspi1.Init.NSSPolarity = SPI_NSS_POLARITY_LOW;

hspi1.Init.FifoThreshold = SPI_FIFO_THRESHOLD_01DATA;

hspi1.Init.TxCRCInitializationPattern = SPI_CRC_INITIALIZATION_ALL_ZERO_PATTERN;

hspi1.Init.RxCRCInitializationPattern = SPI_CRC_INITIALIZATION_ALL_ZERO_PATTERN;

hspi1.Init.MasterSSIdleness = SPI_MASTER_SS_IDLENESS_00CYCLE;

hspi1.Init.MasterInterDataIdleness = SPI_MASTER_INTERDATA_IDLENESS_00CYCLE;

hspi1.Init.MasterReceiverAutoSusp = SPI_MASTER_RX_AUTOSUSP_DISABLE;

hspi1.Init.MasterKeepIOState = SPI_MASTER_KEEP_IO_STATE_DISABLE;

hspi1.Init.IOSwap = SPI_IO_SWAP_DISABLE;


if (HAL_SPI_Init(&hspi1) != HAL_OK)

{

Error_Handler();

}

}


/**

* @brief USART3 Initialization Function

* @param None

* @retval None

*/

static void MX_USART3_UART_Init(void)

{

huart3.Instance = USART3;

huart3.Init.BaudRate = 115200;

huart3.Init.WordLength = UART_WORDLENGTH_8B;

huart3.Init.StopBits = UART_STOPBITS_1;

huart3.Init.Parity = UART_PARITY_NONE;

huart3.Init.Mode = UART_MODE_TX_RX;

huart3.Init.HwFlowCtl = UART_HWCONTROL_NONE;

huart3.Init.OverSampling = UART_OVERSAMPLING_16;

huart3.Init.OneBitSampling = UART_ONE_BIT_SAMPLE_DISABLE;

huart3.Init.ClockPrescaler = UART_PRESCALER_DIV1;

huart3.AdvancedInit.AdvFeatureInit = UART_ADVFEATURE_NO_INIT;


if (HAL_UART_Init(&huart3) != HAL_OK)

{

Error_Handler();

}

if (HAL_UARTEx_SetTxFifoThreshold(&huart3, UART_TXFIFO_THRESHOLD_1_8) != HAL_OK)

{

Error_Handler();

}

if (HAL_UARTEx_SetRxFifoThreshold(&huart3, UART_RXFIFO_THRESHOLD_1_8) != HAL_OK)

{

Error_Handler();

}

if (HAL_UARTEx_DisableFifoMode(&huart3) != HAL_OK)

{

Error_Handler();

}

}


/**

* @brief GPIO Initialization Function

* @param None

* @retval None

*/

static void MX_GPIO_Init(void)

{

GPIO_InitTypeDef GPIO_InitStruct = {0};


__HAL_RCC_GPIOB_CLK_ENABLE(); // Port B for both IMUs

__HAL_RCC_GPIOE_CLK_ENABLE(); // Port E for Baro

__HAL_RCC_GPIOC_CLK_ENABLE(); // Port C for LED


/* Set all CS pins HIGH before configuring to avoid bus noise */

HAL_GPIO_WritePin(GPIOE, Baro_CS_Pin, GPIO_PIN_SET); // PE11

HAL_GPIO_WritePin(GPIOB, IMU_CS_Pin, GPIO_PIN_SET); // PB10

HAL_GPIO_WritePin(GPIOB, IMU_2_CS_Pin, GPIO_PIN_SET); // PB2


/* Configure IMU_CS_Pin (LSM6DSL) - PB10 */

GPIO_InitStruct.Pin = IMU_CS_Pin;

GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;

GPIO_InitStruct.Pull = GPIO_NOPULL;

GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_VERY_HIGH;

HAL_GPIO_Init(GPIOB, &GPIO_InitStruct);


/* Configure IMU_2_CS_Pin (ICM45686) - PB2 */

GPIO_InitStruct.Pin = IMU_2_CS_Pin;

HAL_GPIO_Init(GPIOB, &GPIO_InitStruct);


/* Configure Baro_CS_Pin - PE11 */

GPIO_InitStruct.Pin = Baro_CS_Pin;

HAL_GPIO_Init(GPIOE, &GPIO_InitStruct);

}


/* USER CODE BEGIN 4 */

/* USER CODE END 4 */


void MPU_Config(void)

{

MPU_Region_InitTypeDef MPU_InitStruct = {0};


HAL_MPU_Disable();


MPU_InitStruct.Enable = MPU_REGION_ENABLE;

MPU_InitStruct.Number = MPU_REGION_NUMBER0;

MPU_InitStruct.BaseAddress = 0x0;

MPU_InitStruct.Size = MPU_REGION_SIZE_4GB;

MPU_InitStruct.SubRegionDisable = 0x87;

MPU_InitStruct.TypeExtField = MPU_TEX_LEVEL0;

MPU_InitStruct.AccessPermission = MPU_REGION_NO_ACCESS;

MPU_InitStruct.DisableExec = MPU_INSTRUCTION_ACCESS_DISABLE;

MPU_InitStruct.IsShareable = MPU_ACCESS_SHAREABLE;

MPU_InitStruct.IsCacheable = MPU_ACCESS_NOT_CACHEABLE;

MPU_InitStruct.IsBufferable = MPU_ACCESS_NOT_BUFFERABLE;


HAL_MPU_ConfigRegion(&MPU_InitStruct);

HAL_MPU_Enable(MPU_PRIVILEGED_DEFAULT);

}


void Error_Handler(void)

{

__disable_irq();

while (1)

{

}

}


#ifdef USE_FULL_ASSERT

void assert_failed(uint8_t *file, uint32_t line)

{

}

#endif
