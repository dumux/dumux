################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/1p/sequential/test_1p.cc \
../test/porousmediumflow/1p/sequential/test_diffusion.cc \
../test/porousmediumflow/1p/sequential/test_diffusion3d.cc 

CC_DEPS += \
./test/porousmediumflow/1p/sequential/test_1p.d \
./test/porousmediumflow/1p/sequential/test_diffusion.d \
./test/porousmediumflow/1p/sequential/test_diffusion3d.d 

OBJS += \
./test/porousmediumflow/1p/sequential/test_1p.o \
./test/porousmediumflow/1p/sequential/test_diffusion.o \
./test/porousmediumflow/1p/sequential/test_diffusion3d.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/1p/sequential/%.o: ../test/porousmediumflow/1p/sequential/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


