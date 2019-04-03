################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../tutorial/ex1/exercise1_2p.cc \
../tutorial/ex1/exercise1_2p2c.cc 

CC_DEPS += \
./tutorial/ex1/exercise1_2p.d \
./tutorial/ex1/exercise1_2p2c.d 

OBJS += \
./tutorial/ex1/exercise1_2p.o \
./tutorial/ex1/exercise1_2p2c.o 


# Each subdirectory must supply rules for building sources it contributes
tutorial/ex1/%.o: ../tutorial/ex1/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


