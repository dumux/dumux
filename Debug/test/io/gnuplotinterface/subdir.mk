################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/io/gnuplotinterface/test_gnuplotinterface.cc 

CC_DEPS += \
./test/io/gnuplotinterface/test_gnuplotinterface.d 

OBJS += \
./test/io/gnuplotinterface/test_gnuplotinterface.o 


# Each subdirectory must supply rules for building sources it contributes
test/io/gnuplotinterface/%.o: ../test/io/gnuplotinterface/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


