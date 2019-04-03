################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/1p2c/implicit/test_box1p2c.cc \
../test/porousmediumflow/1p2c/implicit/test_box1p2cniconduction.cc \
../test/porousmediumflow/1p2c/implicit/test_box1p2cniconvection.cc \
../test/porousmediumflow/1p2c/implicit/test_cc1p2c.cc \
../test/porousmediumflow/1p2c/implicit/test_cc1p2cniconduction.cc \
../test/porousmediumflow/1p2c/implicit/test_cc1p2cniconvection.cc 

CC_DEPS += \
./test/porousmediumflow/1p2c/implicit/test_box1p2c.d \
./test/porousmediumflow/1p2c/implicit/test_box1p2cniconduction.d \
./test/porousmediumflow/1p2c/implicit/test_box1p2cniconvection.d \
./test/porousmediumflow/1p2c/implicit/test_cc1p2c.d \
./test/porousmediumflow/1p2c/implicit/test_cc1p2cniconduction.d \
./test/porousmediumflow/1p2c/implicit/test_cc1p2cniconvection.d 

OBJS += \
./test/porousmediumflow/1p2c/implicit/test_box1p2c.o \
./test/porousmediumflow/1p2c/implicit/test_box1p2cniconduction.o \
./test/porousmediumflow/1p2c/implicit/test_box1p2cniconvection.o \
./test/porousmediumflow/1p2c/implicit/test_cc1p2c.o \
./test/porousmediumflow/1p2c/implicit/test_cc1p2cniconduction.o \
./test/porousmediumflow/1p2c/implicit/test_cc1p2cniconvection.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/1p2c/implicit/%.o: ../test/porousmediumflow/1p2c/implicit/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


