################################################################################
# Automatically-generated file. Do not edit!
# Toolchain: GNU Tools for STM32 (14.3.rel1)
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../Core/Src/FATFS_SD.c \
../Core/Src/airbrakes.c \
../Core/Src/buzzer.c \
../Core/Src/cameras.c \
../Core/Src/commands.c \
../Core/Src/complementary_filter.c \
../Core/Src/fsm.c \
../Core/Src/kalman_filter.c \
../Core/Src/launch_buffer.c \
../Core/Src/madgwick.c \
../Core/Src/main.c \
../Core/Src/parachute.c \
../Core/Src/parachutes.c \
../Core/Src/rfm9x.c \
../Core/Src/sd_card.c \
../Core/Src/sensors.c \
../Core/Src/stm32h7xx_hal_msp.c \
../Core/Src/stm32h7xx_it.c \
../Core/Src/syscalls.c \
../Core/Src/sysmem.c \
../Core/Src/system_stm32h7xx.c \
../Core/Src/telemetry.c 

OBJS += \
./Core/Src/FATFS_SD.o \
./Core/Src/airbrakes.o \
./Core/Src/buzzer.o \
./Core/Src/cameras.o \
./Core/Src/commands.o \
./Core/Src/complementary_filter.o \
./Core/Src/fsm.o \
./Core/Src/kalman_filter.o \
./Core/Src/launch_buffer.o \
./Core/Src/madgwick.o \
./Core/Src/main.o \
./Core/Src/parachute.o \
./Core/Src/parachutes.o \
./Core/Src/rfm9x.o \
./Core/Src/sd_card.o \
./Core/Src/sensors.o \
./Core/Src/stm32h7xx_hal_msp.o \
./Core/Src/stm32h7xx_it.o \
./Core/Src/syscalls.o \
./Core/Src/sysmem.o \
./Core/Src/system_stm32h7xx.o \
./Core/Src/telemetry.o 

C_DEPS += \
./Core/Src/FATFS_SD.d \
./Core/Src/airbrakes.d \
./Core/Src/buzzer.d \
./Core/Src/cameras.d \
./Core/Src/commands.d \
./Core/Src/complementary_filter.d \
./Core/Src/fsm.d \
./Core/Src/kalman_filter.d \
./Core/Src/launch_buffer.d \
./Core/Src/madgwick.d \
./Core/Src/main.d \
./Core/Src/parachute.d \
./Core/Src/parachutes.d \
./Core/Src/rfm9x.d \
./Core/Src/sd_card.d \
./Core/Src/sensors.d \
./Core/Src/stm32h7xx_hal_msp.d \
./Core/Src/stm32h7xx_it.d \
./Core/Src/syscalls.d \
./Core/Src/sysmem.d \
./Core/Src/system_stm32h7xx.d \
./Core/Src/telemetry.d 


# Each subdirectory must supply rules for building sources it contributes
Core/Src/%.o Core/Src/%.su Core/Src/%.cyclo: ../Core/Src/%.c Core/Src/subdir.mk
	arm-none-eabi-gcc "$<" -mcpu=cortex-m7 -std=gnu11 -g3 -DDEBUG -DUSE_HAL_DRIVER -DSTM32H743xx -DUSE_PWR_LDO_SUPPLY -c -I../Core/Inc -I../Drivers/STM32H7xx_HAL_Driver/Inc -I../Drivers/STM32H7xx_HAL_Driver/Inc/Legacy -I../Drivers/CMSIS/Device/ST/STM32H7xx/Include -I../Drivers/CMSIS/Include -I../FATFS/Target -I../FATFS/App -I../Middlewares/Third_Party/FatFs/src -O0 -ffunction-sections -fdata-sections -Wall -fstack-usage -fcyclomatic-complexity -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" --specs=nano.specs -mfpu=fpv5-d16 -mfloat-abi=hard -mthumb -o "$@"

clean: clean-Core-2f-Src

clean-Core-2f-Src:
	-$(RM) ./Core/Src/FATFS_SD.cyclo ./Core/Src/FATFS_SD.d ./Core/Src/FATFS_SD.o ./Core/Src/FATFS_SD.su ./Core/Src/airbrakes.cyclo ./Core/Src/airbrakes.d ./Core/Src/airbrakes.o ./Core/Src/airbrakes.su ./Core/Src/buzzer.cyclo ./Core/Src/buzzer.d ./Core/Src/buzzer.o ./Core/Src/buzzer.su ./Core/Src/cameras.cyclo ./Core/Src/cameras.d ./Core/Src/cameras.o ./Core/Src/cameras.su ./Core/Src/commands.cyclo ./Core/Src/commands.d ./Core/Src/commands.o ./Core/Src/commands.su ./Core/Src/complementary_filter.cyclo ./Core/Src/complementary_filter.d ./Core/Src/complementary_filter.o ./Core/Src/complementary_filter.su ./Core/Src/fsm.cyclo ./Core/Src/fsm.d ./Core/Src/fsm.o ./Core/Src/fsm.su ./Core/Src/kalman_filter.cyclo ./Core/Src/kalman_filter.d ./Core/Src/kalman_filter.o ./Core/Src/kalman_filter.su ./Core/Src/launch_buffer.cyclo ./Core/Src/launch_buffer.d ./Core/Src/launch_buffer.o ./Core/Src/launch_buffer.su ./Core/Src/madgwick.cyclo ./Core/Src/madgwick.d ./Core/Src/madgwick.o ./Core/Src/madgwick.su ./Core/Src/main.cyclo ./Core/Src/main.d ./Core/Src/main.o ./Core/Src/main.su ./Core/Src/parachute.cyclo ./Core/Src/parachute.d ./Core/Src/parachute.o ./Core/Src/parachute.su ./Core/Src/parachutes.cyclo ./Core/Src/parachutes.d ./Core/Src/parachutes.o ./Core/Src/parachutes.su ./Core/Src/rfm9x.cyclo ./Core/Src/rfm9x.d ./Core/Src/rfm9x.o ./Core/Src/rfm9x.su ./Core/Src/sd_card.cyclo ./Core/Src/sd_card.d ./Core/Src/sd_card.o ./Core/Src/sd_card.su ./Core/Src/sensors.cyclo ./Core/Src/sensors.d ./Core/Src/sensors.o ./Core/Src/sensors.su ./Core/Src/stm32h7xx_hal_msp.cyclo ./Core/Src/stm32h7xx_hal_msp.d ./Core/Src/stm32h7xx_hal_msp.o ./Core/Src/stm32h7xx_hal_msp.su ./Core/Src/stm32h7xx_it.cyclo ./Core/Src/stm32h7xx_it.d ./Core/Src/stm32h7xx_it.o ./Core/Src/stm32h7xx_it.su ./Core/Src/syscalls.cyclo ./Core/Src/syscalls.d ./Core/Src/syscalls.o ./Core/Src/syscalls.su ./Core/Src/sysmem.cyclo ./Core/Src/sysmem.d ./Core/Src/sysmem.o ./Core/Src/sysmem.su ./Core/Src/system_stm32h7xx.cyclo ./Core/Src/system_stm32h7xx.d ./Core/Src/system_stm32h7xx.o ./Core/Src/system_stm32h7xx.su ./Core/Src/telemetry.cyclo ./Core/Src/telemetry.d ./Core/Src/telemetry.o ./Core/Src/telemetry.su

.PHONY: clean-Core-2f-Src

