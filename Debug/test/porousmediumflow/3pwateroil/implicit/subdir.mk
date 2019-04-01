################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/3pwateroil/implicit/test_box3pwateroil.cc 

CC_DEPS += \
./test/porousmediumflow/3pwateroil/implicit/test_box3pwateroil.d 

OBJS += \
./test/porousmediumflow/3pwateroil/implicit/test_box3pwateroil.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/3pwateroil/implicit/%.o: ../test/porousmediumflow/3pwateroil/implicit/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


