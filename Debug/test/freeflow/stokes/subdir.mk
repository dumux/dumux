################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/freeflow/stokes/test_stokes.cc 

CC_DEPS += \
./test/freeflow/stokes/test_stokes.d 

OBJS += \
./test/freeflow/stokes/test_stokes.o 


# Each subdirectory must supply rules for building sources it contributes
test/freeflow/stokes/%.o: ../test/freeflow/stokes/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


