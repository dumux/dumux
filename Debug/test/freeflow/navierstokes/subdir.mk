################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/freeflow/navierstokes/test_navierstokes.cc 

CC_DEPS += \
./test/freeflow/navierstokes/test_navierstokes.d 

OBJS += \
./test/freeflow/navierstokes/test_navierstokes.o 


# Each subdirectory must supply rules for building sources it contributes
test/freeflow/navierstokes/%.o: ../test/freeflow/navierstokes/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


