################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/geomechanics/elastic/test_elastic.cc 

CC_DEPS += \
./test/geomechanics/elastic/test_elastic.d 

OBJS += \
./test/geomechanics/elastic/test_elastic.o 


# Each subdirectory must supply rules for building sources it contributes
test/geomechanics/elastic/%.o: ../test/geomechanics/elastic/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


