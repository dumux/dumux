################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../tutorial/ex3/ex3_a.cc \
../tutorial/ex3/ex3_b.cc 

CC_DEPS += \
./tutorial/ex3/ex3_a.d \
./tutorial/ex3/ex3_b.d 

OBJS += \
./tutorial/ex3/ex3_a.o \
./tutorial/ex3/ex3_b.o 


# Each subdirectory must supply rules for building sources it contributes
tutorial/ex3/%.o: ../tutorial/ex3/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


