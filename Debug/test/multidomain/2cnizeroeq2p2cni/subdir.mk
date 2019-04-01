################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/multidomain/2cnizeroeq2p2cni/test_2cnizeroeq2p2cni.cc 

CC_DEPS += \
./test/multidomain/2cnizeroeq2p2cni/test_2cnizeroeq2p2cni.d 

OBJS += \
./test/multidomain/2cnizeroeq2p2cni/test_2cnizeroeq2p2cni.o 


# Each subdirectory must supply rules for building sources it contributes
test/multidomain/2cnizeroeq2p2cni/%.o: ../test/multidomain/2cnizeroeq2p2cni/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


