################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../tutorial/ex2/exercise2.cc 

CC_DEPS += \
./tutorial/ex2/exercise2.d 

OBJS += \
./tutorial/ex2/exercise2.o 


# Each subdirectory must supply rules for building sources it contributes
tutorial/ex2/%.o: ../tutorial/ex2/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


