################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/multidomain/2cstokes2p2c/test_2cstokes2p2c.cc 

CC_DEPS += \
./test/multidomain/2cstokes2p2c/test_2cstokes2p2c.d 

OBJS += \
./test/multidomain/2cstokes2p2c/test_2cstokes2p2c.o 


# Each subdirectory must supply rules for building sources it contributes
test/multidomain/2cstokes2p2c/%.o: ../test/multidomain/2cstokes2p2c/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


