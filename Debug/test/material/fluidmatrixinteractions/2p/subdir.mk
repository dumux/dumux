################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/material/fluidmatrixinteractions/2p/test_thermalconductivityjohansen.cc \
../test/material/fluidmatrixinteractions/2p/test_thermalconductivitysomerton.cc 

CC_DEPS += \
./test/material/fluidmatrixinteractions/2p/test_thermalconductivityjohansen.d \
./test/material/fluidmatrixinteractions/2p/test_thermalconductivitysomerton.d 

OBJS += \
./test/material/fluidmatrixinteractions/2p/test_thermalconductivityjohansen.o \
./test/material/fluidmatrixinteractions/2p/test_thermalconductivitysomerton.o 


# Each subdirectory must supply rules for building sources it contributes
test/material/fluidmatrixinteractions/2p/%.o: ../test/material/fluidmatrixinteractions/2p/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


