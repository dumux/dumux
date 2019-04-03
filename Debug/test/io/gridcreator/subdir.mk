################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/io/gridcreator/test_gridcreator_cake.cc \
../test/io/gridcreator/test_gridcreator_gmsh.cc 

CC_DEPS += \
./test/io/gridcreator/test_gridcreator_cake.d \
./test/io/gridcreator/test_gridcreator_gmsh.d 

OBJS += \
./test/io/gridcreator/test_gridcreator_cake.o \
./test/io/gridcreator/test_gridcreator_gmsh.o 


# Each subdirectory must supply rules for building sources it contributes
test/io/gridcreator/%.o: ../test/io/gridcreator/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


