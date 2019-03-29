################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../tutorial/solution/ex1/exercise1_2pni.cc 

CC_DEPS += \
./tutorial/solution/ex1/exercise1_2pni.d 

OBJS += \
./tutorial/solution/ex1/exercise1_2pni.o 


# Each subdirectory must supply rules for building sources it contributes
tutorial/solution/ex1/%.o: ../tutorial/solution/ex1/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


