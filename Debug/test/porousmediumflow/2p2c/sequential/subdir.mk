################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/2p2c/sequential/test_adaptive2p2c2d.cc \
../test/porousmediumflow/2p2c/sequential/test_adaptive2p2c3d.cc \
../test/porousmediumflow/2p2c/sequential/test_dec2p2c.cc \
../test/porousmediumflow/2p2c/sequential/test_multiphysics2p2c.cc 

CC_DEPS += \
./test/porousmediumflow/2p2c/sequential/test_adaptive2p2c2d.d \
./test/porousmediumflow/2p2c/sequential/test_adaptive2p2c3d.d \
./test/porousmediumflow/2p2c/sequential/test_dec2p2c.d \
./test/porousmediumflow/2p2c/sequential/test_multiphysics2p2c.d 

OBJS += \
./test/porousmediumflow/2p2c/sequential/test_adaptive2p2c2d.o \
./test/porousmediumflow/2p2c/sequential/test_adaptive2p2c3d.o \
./test/porousmediumflow/2p2c/sequential/test_dec2p2c.o \
./test/porousmediumflow/2p2c/sequential/test_multiphysics2p2c.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/2p2c/sequential/%.o: ../test/porousmediumflow/2p2c/sequential/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


