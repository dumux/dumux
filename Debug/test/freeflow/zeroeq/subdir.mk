################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/freeflow/zeroeq/test_zeroeq.cc \
../test/freeflow/zeroeq/test_zeroeq_channel.cc 

CC_DEPS += \
./test/freeflow/zeroeq/test_zeroeq.d \
./test/freeflow/zeroeq/test_zeroeq_channel.d 

OBJS += \
./test/freeflow/zeroeq/test_zeroeq.o \
./test/freeflow/zeroeq/test_zeroeq_channel.o 


# Each subdirectory must supply rules for building sources it contributes
test/freeflow/zeroeq/%.o: ../test/freeflow/zeroeq/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


