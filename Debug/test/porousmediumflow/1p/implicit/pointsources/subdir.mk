################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/1p/implicit/pointsources/test_box1p_pointsources.cc \
../test/porousmediumflow/1p/implicit/pointsources/test_cc1p_pointsources.cc \
../test/porousmediumflow/1p/implicit/pointsources/test_cc1p_pointsources_timedependent.cc 

CC_DEPS += \
./test/porousmediumflow/1p/implicit/pointsources/test_box1p_pointsources.d \
./test/porousmediumflow/1p/implicit/pointsources/test_cc1p_pointsources.d \
./test/porousmediumflow/1p/implicit/pointsources/test_cc1p_pointsources_timedependent.d 

OBJS += \
./test/porousmediumflow/1p/implicit/pointsources/test_box1p_pointsources.o \
./test/porousmediumflow/1p/implicit/pointsources/test_cc1p_pointsources.o \
./test/porousmediumflow/1p/implicit/pointsources/test_cc1p_pointsources_timedependent.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/1p/implicit/pointsources/%.o: ../test/porousmediumflow/1p/implicit/pointsources/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


