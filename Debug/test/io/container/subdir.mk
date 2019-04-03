################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/io/container/test_container_io.cc 

CC_DEPS += \
./test/io/container/test_container_io.d 

OBJS += \
./test/io/container/test_container_io.o 


# Each subdirectory must supply rules for building sources it contributes
test/io/container/%.o: ../test/io/container/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


