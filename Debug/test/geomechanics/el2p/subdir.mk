################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/geomechanics/el2p/test_el2p.cc 

CC_DEPS += \
./test/geomechanics/el2p/test_el2p.d 

OBJS += \
./test/geomechanics/el2p/test_el2p.o 


# Each subdirectory must supply rules for building sources it contributes
test/geomechanics/el2p/%.o: ../test/geomechanics/el2p/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


