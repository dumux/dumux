################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/1p/implicit/test_box1p.cc \
../test/porousmediumflow/1p/implicit/test_box1p1d3d.cc \
../test/porousmediumflow/1p/implicit/test_box1p2d3d.cc \
../test/porousmediumflow/1p/implicit/test_box1pniconduction.cc \
../test/porousmediumflow/1p/implicit/test_box1pniconvection.cc \
../test/porousmediumflow/1p/implicit/test_box1pwithamg.cc \
../test/porousmediumflow/1p/implicit/test_cc1p.cc \
../test/porousmediumflow/1p/implicit/test_cc1p1d3d.cc \
../test/porousmediumflow/1p/implicit/test_cc1p2d3d.cc \
../test/porousmediumflow/1p/implicit/test_cc1pniconduction.cc \
../test/porousmediumflow/1p/implicit/test_cc1pniconvection.cc \
../test/porousmediumflow/1p/implicit/test_cc1pwithamg.cc \
../test/porousmediumflow/1p/implicit/test_cc1pwithgstat.cc 

CC_DEPS += \
./test/porousmediumflow/1p/implicit/test_box1p.d \
./test/porousmediumflow/1p/implicit/test_box1p1d3d.d \
./test/porousmediumflow/1p/implicit/test_box1p2d3d.d \
./test/porousmediumflow/1p/implicit/test_box1pniconduction.d \
./test/porousmediumflow/1p/implicit/test_box1pniconvection.d \
./test/porousmediumflow/1p/implicit/test_box1pwithamg.d \
./test/porousmediumflow/1p/implicit/test_cc1p.d \
./test/porousmediumflow/1p/implicit/test_cc1p1d3d.d \
./test/porousmediumflow/1p/implicit/test_cc1p2d3d.d \
./test/porousmediumflow/1p/implicit/test_cc1pniconduction.d \
./test/porousmediumflow/1p/implicit/test_cc1pniconvection.d \
./test/porousmediumflow/1p/implicit/test_cc1pwithamg.d \
./test/porousmediumflow/1p/implicit/test_cc1pwithgstat.d 

OBJS += \
./test/porousmediumflow/1p/implicit/test_box1p.o \
./test/porousmediumflow/1p/implicit/test_box1p1d3d.o \
./test/porousmediumflow/1p/implicit/test_box1p2d3d.o \
./test/porousmediumflow/1p/implicit/test_box1pniconduction.o \
./test/porousmediumflow/1p/implicit/test_box1pniconvection.o \
./test/porousmediumflow/1p/implicit/test_box1pwithamg.o \
./test/porousmediumflow/1p/implicit/test_cc1p.o \
./test/porousmediumflow/1p/implicit/test_cc1p1d3d.o \
./test/porousmediumflow/1p/implicit/test_cc1p2d3d.o \
./test/porousmediumflow/1p/implicit/test_cc1pniconduction.o \
./test/porousmediumflow/1p/implicit/test_cc1pniconvection.o \
./test/porousmediumflow/1p/implicit/test_cc1pwithamg.o \
./test/porousmediumflow/1p/implicit/test_cc1pwithgstat.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/1p/implicit/%.o: ../test/porousmediumflow/1p/implicit/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


