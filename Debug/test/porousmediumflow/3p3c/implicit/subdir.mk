################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/3p3c/implicit/test_box3p3c.cc \
../test/porousmediumflow/3p3c/implicit/test_box3p3cnicolumnxylol.cc \
../test/porousmediumflow/3p3c/implicit/test_box3p3cnikuevette.cc \
../test/porousmediumflow/3p3c/implicit/test_cc3p3c.cc \
../test/porousmediumflow/3p3c/implicit/test_cc3p3cnicolumnxylol.cc \
../test/porousmediumflow/3p3c/implicit/test_cc3p3cnikuevette.cc 

CC_DEPS += \
./test/porousmediumflow/3p3c/implicit/test_box3p3c.d \
./test/porousmediumflow/3p3c/implicit/test_box3p3cnicolumnxylol.d \
./test/porousmediumflow/3p3c/implicit/test_box3p3cnikuevette.d \
./test/porousmediumflow/3p3c/implicit/test_cc3p3c.d \
./test/porousmediumflow/3p3c/implicit/test_cc3p3cnicolumnxylol.d \
./test/porousmediumflow/3p3c/implicit/test_cc3p3cnikuevette.d 

OBJS += \
./test/porousmediumflow/3p3c/implicit/test_box3p3c.o \
./test/porousmediumflow/3p3c/implicit/test_box3p3cnicolumnxylol.o \
./test/porousmediumflow/3p3c/implicit/test_box3p3cnikuevette.o \
./test/porousmediumflow/3p3c/implicit/test_cc3p3c.o \
./test/porousmediumflow/3p3c/implicit/test_cc3p3cnicolumnxylol.o \
./test/porousmediumflow/3p3c/implicit/test_cc3p3cnikuevette.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/3p3c/implicit/%.o: ../test/porousmediumflow/3p3c/implicit/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


