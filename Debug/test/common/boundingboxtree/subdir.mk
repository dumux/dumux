################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/common/boundingboxtree/test_bboxtree.cc 

CC_DEPS += \
./test/common/boundingboxtree/test_bboxtree.d 

OBJS += \
./test/common/boundingboxtree/test_bboxtree.o 


# Each subdirectory must supply rules for building sources it contributes
test/common/boundingboxtree/%.o: ../test/common/boundingboxtree/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


