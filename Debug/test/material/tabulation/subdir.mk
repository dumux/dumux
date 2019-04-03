################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/material/tabulation/test_tabulation.cc 

CC_DEPS += \
./test/material/tabulation/test_tabulation.d 

OBJS += \
./test/material/tabulation/test_tabulation.o 


# Each subdirectory must supply rules for building sources it contributes
test/material/tabulation/%.o: ../test/material/tabulation/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


