################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/material/fluidsystems/test_fluidsystems.cc 

CC_DEPS += \
./test/material/fluidsystems/test_fluidsystems.d 

OBJS += \
./test/material/fluidsystems/test_fluidsystems.o 


# Each subdirectory must supply rules for building sources it contributes
test/material/fluidsystems/%.o: ../test/material/fluidsystems/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


