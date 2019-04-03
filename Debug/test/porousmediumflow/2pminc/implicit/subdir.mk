################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/2pminc/implicit/test_box2pminc.cc 

CC_DEPS += \
./test/porousmediumflow/2pminc/implicit/test_box2pminc.d 

OBJS += \
./test/porousmediumflow/2pminc/implicit/test_box2pminc.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/2pminc/implicit/%.o: ../test/porousmediumflow/2pminc/implicit/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


