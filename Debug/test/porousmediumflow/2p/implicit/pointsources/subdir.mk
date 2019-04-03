################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/2p/implicit/pointsources/test_ccadaptive2ppointsource.cc 

CC_DEPS += \
./test/porousmediumflow/2p/implicit/pointsources/test_ccadaptive2ppointsource.d 

OBJS += \
./test/porousmediumflow/2p/implicit/pointsources/test_ccadaptive2ppointsource.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/2p/implicit/pointsources/%.o: ../test/porousmediumflow/2p/implicit/pointsources/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


