################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/freeflow/zeroeq2cni/test_zeroeq2cni.cc 

CC_DEPS += \
./test/freeflow/zeroeq2cni/test_zeroeq2cni.d 

OBJS += \
./test/freeflow/zeroeq2cni/test_zeroeq2cni.o 


# Each subdirectory must supply rules for building sources it contributes
test/freeflow/zeroeq2cni/%.o: ../test/freeflow/zeroeq2cni/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


