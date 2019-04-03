################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/common/propertysystem/test_propertysystem.cc 

CC_DEPS += \
./test/common/propertysystem/test_propertysystem.d 

OBJS += \
./test/common/propertysystem/test_propertysystem.o 


# Each subdirectory must supply rules for building sources it contributes
test/common/propertysystem/%.o: ../test/common/propertysystem/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


