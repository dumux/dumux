################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/material/immiscibleflash/test_immiscibleflash.cc 

CC_DEPS += \
./test/material/immiscibleflash/test_immiscibleflash.d 

OBJS += \
./test/material/immiscibleflash/test_immiscibleflash.o 


# Each subdirectory must supply rules for building sources it contributes
test/material/immiscibleflash/%.o: ../test/material/immiscibleflash/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


