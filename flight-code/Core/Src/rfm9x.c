#include "rfm9x.h"

/* SX1276 Register Map */
#define REG_FIFO                0x00
#define REG_OP_MODE             0x01
#define REG_FRF_MSB             0x06
#define REG_FRF_MID             0x07
#define REG_FRF_LSB             0x08
#define REG_PA_CONFIG           0x09
#define REG_FIFO_ADDR_PTR       0x0D
#define REG_IRQ_FLAGS           0x12
#define REG_MODEM_CONFIG1       0x1D
#define REG_MODEM_CONFIG2       0x1E
#define REG_MODEM_CONFIG3		0x26
#define REG_PAYLOAD_LENGTH      0x22
#define REG_SYNC_WORD           0x39
#define REG_FIFO_RX_BASE_ADDR   0x0F
#define REG_FIFO_TX_BASE_ADDR	0x0E
#define REG_FIFO_RX_CURRENT     0x10
#define REG_RX_NB_BYTES         0x13
#define REG_PREAMBLE_MSB		0x20
#define REG_PREAMBLE_LSB		0x21
#define REG_LNA					0x0C

#define IRQ_RX_DONE             0x40
#define IRQ_PAYLOAD_CRC_ERROR   0x20
#define IRQ_TX_DONE             0x08

#define MODE_SLEEP              0x80
#define MODE_STDBY              0x81
#define MODE_TX                 0x83
#define MODE_RX_CONTINUOUS		0x85

// Globals
SPI_HandleTypeDef *hspi;
GPIO_TypeDef *cs_port;
uint16_t cs_pin;
GPIO_TypeDef *reset_port;
uint16_t reset_pin;
GPIO_TypeDef *en_port;
uint16_t en_pin;
uint8_t tx_busy;

/* -------------------------------------------------- */
/* Low-level SPI helpers */
/* -------------------------------------------------- */

static void cs_select_rfm9x()
{
    HAL_GPIO_WritePin(cs_port, cs_pin, GPIO_PIN_RESET);
}

static void cs_deselect_rfm9x()
{
    HAL_GPIO_WritePin(cs_port, cs_pin, GPIO_PIN_SET);
}

static void write_reg(uint8_t addr, uint8_t value)
{
    uint8_t buf[2];
    buf[0] = addr | 0x80;
    buf[1] = value;

    cs_select_rfm9x();
    HAL_SPI_Transmit(hspi, buf, 2, HAL_MAX_DELAY);
    cs_deselect_rfm9x();
}

static uint8_t read_reg(uint8_t addr)
{
    uint8_t tx[2] = {addr & 0x7F, 0x00};
    uint8_t rx[2];

    cs_select_rfm9x();
    HAL_SPI_TransmitReceive(hspi, tx, rx, 2, HAL_MAX_DELAY);
    cs_deselect_rfm9x();

    return rx[1];
}

/* -------------------------------------------------- */

void RFM9X_Reset()
{
    HAL_GPIO_WritePin(reset_port, reset_pin, GPIO_PIN_RESET);
    HAL_Delay(1);
    HAL_GPIO_WritePin(reset_port, reset_pin, GPIO_PIN_SET);
    HAL_Delay(10);
}

/* -------------------------------------------------- */

void RFM9X_Init(SPI_HandleTypeDef *hspi_p, GPIO_TypeDef *cs_port_p, uint16_t cs_pin_p, GPIO_TypeDef *reset_port_p, uint16_t reset_pin_p, GPIO_TypeDef *en_port_p, uint16_t en_pin_p)
{
	hspi = hspi_p;
	cs_port = cs_port_p;
	cs_pin = cs_pin_p;
	reset_port = reset_port_p;
	reset_pin = reset_pin_p;
	en_port = en_port_p;
	en_pin = en_pin_p;
    tx_busy = 0;

    /* Ensure CS idle */
	cs_deselect_rfm9x();

	/* Enable module */
    HAL_GPIO_WritePin(en_port, en_pin, GPIO_PIN_SET);
    HAL_Delay(10);

    /* Hardware reset */
	RFM9X_Reset();

	/* Enter sleep mode with LoRa enabled */
	write_reg(REG_OP_MODE, MODE_SLEEP);
	HAL_Delay(1);

	// Verify LoRa mode enabled
	if (read_reg(REG_OP_MODE) != MODE_SLEEP){
		//error
		return;
	}

	/* Frequency */
	RFM9X_SetFrequency(433000000);

	/* Power */
	RFM9X_SetTxPower(17);

	/* Configure FIFO */
	write_reg(REG_FIFO_TX_BASE_ADDR, 0x00);
	write_reg(REG_FIFO_RX_BASE_ADDR, 0x00);
	write_reg(REG_FIFO_ADDR_PTR,     0x00);

	/* Modem config */
	write_reg(REG_MODEM_CONFIG1, 0x72); // BW 125kHz, CR 4/5
	write_reg(REG_MODEM_CONFIG2, 0x74); // SF7
	write_reg(REG_MODEM_CONFIG3, 0x04); // AGC ON

	write_reg(REG_LNA, 0x23);
	write_reg(REG_PREAMBLE_MSB, 0x00);
	write_reg(REG_PREAMBLE_LSB, 0x08);

	/* Sync word */
	write_reg(REG_SYNC_WORD, 0x12);

	/* Clear IRQ flags */
	write_reg(REG_IRQ_FLAGS, 0xFF);

	/* Enter continuous receive mode */
	write_reg(REG_OP_MODE, MODE_STDBY);
	HAL_Delay(1);

	write_reg(REG_OP_MODE, MODE_RX_CONTINUOUS);
	HAL_Delay(1);
}

/* -------------------------------------------------- */

void RFM9X_SetFrequency(uint32_t freq_hz)
{
	uint64_t frf = ((uint64_t)freq_hz << 19) / 32000000;

    write_reg(REG_FRF_MSB, (uint8_t)(frf >> 16));
    write_reg(REG_FRF_MID, (uint8_t)(frf >> 8));
    write_reg(REG_FRF_LSB, (uint8_t)(frf));
}

/* -------------------------------------------------- */

void RFM9X_SetTxPower(uint8_t power)
{
	if(power > 17) power = 17;
    write_reg(REG_PA_CONFIG, 0x80 | (power - 2));
}

/* -------------------------------------------------- */

void RFM9X_Send(uint8_t *data, uint8_t len)
{
    if(tx_busy) return;

    /* MUST enter standby first */
    write_reg(REG_OP_MODE, MODE_STDBY);
    HAL_Delay(1);

    /* Set FIFO pointer */
    write_reg(REG_FIFO_ADDR_PTR, 0x00);

    /* Write payload */
    for(uint8_t i = 0; i < len; i++)
        write_reg(REG_FIFO, data[i]);

    /* Set payload length */
    write_reg(REG_PAYLOAD_LENGTH, len);

    /* Clear IRQ flags */
    write_reg(REG_IRQ_FLAGS, 0xFF);

    /* Start TX */
    write_reg(REG_OP_MODE, MODE_TX);

    tx_busy = 1;
}

/* -------------------------------------------------- */

void RFM9X_Poll()
{
	uint8_t irq = read_reg(REG_IRQ_FLAGS);

	/* TX finished */
	if(tx_busy && (irq & IRQ_TX_DONE))
	{
		write_reg(REG_IRQ_FLAGS, IRQ_TX_DONE);
		tx_busy = 0;

		write_reg(REG_FIFO_RX_BASE_ADDR, 0x00);
		write_reg(REG_FIFO_ADDR_PTR, 0x00);

		/* Immediately return to RX */
		write_reg(REG_OP_MODE, MODE_RX_CONTINUOUS);
	}
}

/* -------------------------------------------------- */

uint8_t RFM9X_IsTxBusy()
{
    return tx_busy;
}

/* -------------------------------------------------- */

uint8_t RFM9X_Receive(uint8_t *buf, uint8_t max_len)
{
	uint8_t irq = read_reg(REG_IRQ_FLAGS);

	if(!(irq & IRQ_RX_DONE))
		return 0;

	if(irq & IRQ_PAYLOAD_CRC_ERROR)
	{
		write_reg(REG_IRQ_FLAGS, 0xFF);
		return 0;
	}

	uint8_t len = read_reg(REG_RX_NB_BYTES);
	if(len > max_len) len = max_len;

	uint8_t fifo_addr = read_reg(REG_FIFO_RX_CURRENT);
	write_reg(REG_FIFO_ADDR_PTR, fifo_addr);

	for(uint8_t i = 0; i < len; i++)
		buf[i] = read_reg(REG_FIFO);

	if(len >= max_len)
	    len = max_len - 1;
	buf[len] = '\0';

	write_reg(REG_IRQ_FLAGS, 0xFF);

	return len;
}
