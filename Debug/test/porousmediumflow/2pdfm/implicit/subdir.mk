################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/2pdfm/implicit/test_2pdfm.cc 

CC_DEPS += \
./test/porousmediumflow/2pdfm/implicit/test_2pdfm.d 

OBJS += \
./test/porousmediumflow/2pdfm/implicit/test_2pdfm.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/2pdfm/implicit/%.o: ../test/porousmediumflow/2pdfm/implicit/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


