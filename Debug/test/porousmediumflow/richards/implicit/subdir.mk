################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/richards/implicit/test_boxrichards.cc \
../test/porousmediumflow/richards/implicit/test_boxrichardsniconduction.cc \
../test/porousmediumflow/richards/implicit/test_boxrichardsniconvection.cc \
../test/porousmediumflow/richards/implicit/test_ccrichards.cc \
../test/porousmediumflow/richards/implicit/test_ccrichardsanalytical.cc \
../test/porousmediumflow/richards/implicit/test_ccrichardsniconduction.cc \
../test/porousmediumflow/richards/implicit/test_ccrichardsniconvection.cc 

CC_DEPS += \
./test/porousmediumflow/richards/implicit/test_boxrichards.d \
./test/porousmediumflow/richards/implicit/test_boxrichardsniconduction.d \
./test/porousmediumflow/richards/implicit/test_boxrichardsniconvection.d \
./test/porousmediumflow/richards/implicit/test_ccrichards.d \
./test/porousmediumflow/richards/implicit/test_ccrichardsanalytical.d \
./test/porousmediumflow/richards/implicit/test_ccrichardsniconduction.d \
./test/porousmediumflow/richards/implicit/test_ccrichardsniconvection.d 

OBJS += \
./test/porousmediumflow/richards/implicit/test_boxrichards.o \
./test/porousmediumflow/richards/implicit/test_boxrichardsniconduction.o \
./test/porousmediumflow/richards/implicit/test_boxrichardsniconvection.o \
./test/porousmediumflow/richards/implicit/test_ccrichards.o \
./test/porousmediumflow/richards/implicit/test_ccrichardsanalytical.o \
./test/porousmediumflow/richards/implicit/test_ccrichardsniconduction.o \
./test/porousmediumflow/richards/implicit/test_ccrichardsniconvection.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/richards/implicit/%.o: ../test/porousmediumflow/richards/implicit/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


