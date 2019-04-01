################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/co2/implicit/test_boxco2.cc \
../test/porousmediumflow/co2/implicit/test_boxco2ni.cc \
../test/porousmediumflow/co2/implicit/test_ccco2.cc \
../test/porousmediumflow/co2/implicit/test_ccco2ni.cc 

CC_DEPS += \
./test/porousmediumflow/co2/implicit/test_boxco2.d \
./test/porousmediumflow/co2/implicit/test_boxco2ni.d \
./test/porousmediumflow/co2/implicit/test_ccco2.d \
./test/porousmediumflow/co2/implicit/test_ccco2ni.d 

OBJS += \
./test/porousmediumflow/co2/implicit/test_boxco2.o \
./test/porousmediumflow/co2/implicit/test_boxco2ni.o \
./test/porousmediumflow/co2/implicit/test_ccco2.o \
./test/porousmediumflow/co2/implicit/test_ccco2ni.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/co2/implicit/%.o: ../test/porousmediumflow/co2/implicit/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


