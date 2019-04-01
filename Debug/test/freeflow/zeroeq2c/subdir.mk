################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/freeflow/zeroeq2c/test_zeroeq2c.cc 

CC_DEPS += \
./test/freeflow/zeroeq2c/test_zeroeq2c.d 

OBJS += \
./test/freeflow/zeroeq2c/test_zeroeq2c.o 


# Each subdirectory must supply rules for building sources it contributes
test/freeflow/zeroeq2c/%.o: ../test/freeflow/zeroeq2c/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


