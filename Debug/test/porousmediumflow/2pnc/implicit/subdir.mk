################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/2pnc/implicit/test_box2pnc.cc \
../test/porousmediumflow/2pnc/implicit/test_cc2pnc.cc 

CC_DEPS += \
./test/porousmediumflow/2pnc/implicit/test_box2pnc.d \
./test/porousmediumflow/2pnc/implicit/test_cc2pnc.d 

OBJS += \
./test/porousmediumflow/2pnc/implicit/test_box2pnc.o \
./test/porousmediumflow/2pnc/implicit/test_cc2pnc.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/2pnc/implicit/%.o: ../test/porousmediumflow/2pnc/implicit/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


