################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/2p2c/implicit/test_box2p2c.cc \
../test/porousmediumflow/2p2c/implicit/test_box2p2cni.cc \
../test/porousmediumflow/2p2c/implicit/test_cc2p2c.cc \
../test/porousmediumflow/2p2c/implicit/test_cc2p2cni.cc 

CC_DEPS += \
./test/porousmediumflow/2p2c/implicit/test_box2p2c.d \
./test/porousmediumflow/2p2c/implicit/test_box2p2cni.d \
./test/porousmediumflow/2p2c/implicit/test_cc2p2c.d \
./test/porousmediumflow/2p2c/implicit/test_cc2p2cni.d 

OBJS += \
./test/porousmediumflow/2p2c/implicit/test_box2p2c.o \
./test/porousmediumflow/2p2c/implicit/test_box2p2cni.o \
./test/porousmediumflow/2p2c/implicit/test_cc2p2c.o \
./test/porousmediumflow/2p2c/implicit/test_cc2p2cni.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/2p2c/implicit/%.o: ../test/porousmediumflow/2p2c/implicit/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


