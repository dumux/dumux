################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../tutorial/tutorial_sequential/tutorial_sequential.cc 

CC_DEPS += \
./tutorial/tutorial_sequential/tutorial_sequential.d 

OBJS += \
./tutorial/tutorial_sequential/tutorial_sequential.o 


# Each subdirectory must supply rules for building sources it contributes
tutorial/tutorial_sequential/%.o: ../tutorial/tutorial_sequential/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


