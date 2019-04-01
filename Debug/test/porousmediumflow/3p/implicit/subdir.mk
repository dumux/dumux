################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/3p/implicit/test_box3p.cc \
../test/porousmediumflow/3p/implicit/test_box3pniconduction.cc \
../test/porousmediumflow/3p/implicit/test_box3pniconvection.cc \
../test/porousmediumflow/3p/implicit/test_cc3p.cc \
../test/porousmediumflow/3p/implicit/test_cc3pniconduction.cc \
../test/porousmediumflow/3p/implicit/test_cc3pniconvection.cc 

CC_DEPS += \
./test/porousmediumflow/3p/implicit/test_box3p.d \
./test/porousmediumflow/3p/implicit/test_box3pniconduction.d \
./test/porousmediumflow/3p/implicit/test_box3pniconvection.d \
./test/porousmediumflow/3p/implicit/test_cc3p.d \
./test/porousmediumflow/3p/implicit/test_cc3pniconduction.d \
./test/porousmediumflow/3p/implicit/test_cc3pniconvection.d 

OBJS += \
./test/porousmediumflow/3p/implicit/test_box3p.o \
./test/porousmediumflow/3p/implicit/test_box3pniconduction.o \
./test/porousmediumflow/3p/implicit/test_box3pniconvection.o \
./test/porousmediumflow/3p/implicit/test_cc3p.o \
./test/porousmediumflow/3p/implicit/test_cc3pniconduction.o \
./test/porousmediumflow/3p/implicit/test_cc3pniconvection.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/3p/implicit/%.o: ../test/porousmediumflow/3p/implicit/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


