################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/material/fluidmatrixinteractions/test_effectivediffusivityconstant.cc \
../test/material/fluidmatrixinteractions/test_effectivediffusivityconstanttau.cc \
../test/material/fluidmatrixinteractions/test_effectivediffusivitymillingtonquirk.cc 

CC_DEPS += \
./test/material/fluidmatrixinteractions/test_effectivediffusivityconstant.d \
./test/material/fluidmatrixinteractions/test_effectivediffusivityconstanttau.d \
./test/material/fluidmatrixinteractions/test_effectivediffusivitymillingtonquirk.d 

OBJS += \
./test/material/fluidmatrixinteractions/test_effectivediffusivityconstant.o \
./test/material/fluidmatrixinteractions/test_effectivediffusivityconstanttau.o \
./test/material/fluidmatrixinteractions/test_effectivediffusivitymillingtonquirk.o 


# Each subdirectory must supply rules for building sources it contributes
test/material/fluidmatrixinteractions/%.o: ../test/material/fluidmatrixinteractions/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


