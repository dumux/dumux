################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/mpnc/implicit/test_boxmpnc.cc \
../test/porousmediumflow/mpnc/implicit/test_boxmpnckinetic.cc \
../test/porousmediumflow/mpnc/implicit/test_boxmpncthermalnonequil.cc \
../test/porousmediumflow/mpnc/implicit/test_ccmpnc.cc \
../test/porousmediumflow/mpnc/implicit/test_forchheimer1p.cc \
../test/porousmediumflow/mpnc/implicit/test_forchheimer2p.cc 

CC_DEPS += \
./test/porousmediumflow/mpnc/implicit/test_boxmpnc.d \
./test/porousmediumflow/mpnc/implicit/test_boxmpnckinetic.d \
./test/porousmediumflow/mpnc/implicit/test_boxmpncthermalnonequil.d \
./test/porousmediumflow/mpnc/implicit/test_ccmpnc.d \
./test/porousmediumflow/mpnc/implicit/test_forchheimer1p.d \
./test/porousmediumflow/mpnc/implicit/test_forchheimer2p.d 

OBJS += \
./test/porousmediumflow/mpnc/implicit/test_boxmpnc.o \
./test/porousmediumflow/mpnc/implicit/test_boxmpnckinetic.o \
./test/porousmediumflow/mpnc/implicit/test_boxmpncthermalnonequil.o \
./test/porousmediumflow/mpnc/implicit/test_ccmpnc.o \
./test/porousmediumflow/mpnc/implicit/test_forchheimer1p.o \
./test/porousmediumflow/mpnc/implicit/test_forchheimer2p.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/mpnc/implicit/%.o: ../test/porousmediumflow/mpnc/implicit/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


