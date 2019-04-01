################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/2p/sequential/test_3d2p.cc \
../test/porousmediumflow/2p/sequential/test_impes.cc \
../test/porousmediumflow/2p/sequential/test_impesadaptive.cc \
../test/porousmediumflow/2p/sequential/test_impesadaptiverestart.cc \
../test/porousmediumflow/2p/sequential/test_impeswithamg.cc \
../test/porousmediumflow/2p/sequential/test_mpfa2p.cc \
../test/porousmediumflow/2p/sequential/test_transport.cc 

CC_DEPS += \
./test/porousmediumflow/2p/sequential/test_3d2p.d \
./test/porousmediumflow/2p/sequential/test_impes.d \
./test/porousmediumflow/2p/sequential/test_impesadaptive.d \
./test/porousmediumflow/2p/sequential/test_impesadaptiverestart.d \
./test/porousmediumflow/2p/sequential/test_impeswithamg.d \
./test/porousmediumflow/2p/sequential/test_mpfa2p.d \
./test/porousmediumflow/2p/sequential/test_transport.d 

OBJS += \
./test/porousmediumflow/2p/sequential/test_3d2p.o \
./test/porousmediumflow/2p/sequential/test_impes.o \
./test/porousmediumflow/2p/sequential/test_impesadaptive.o \
./test/porousmediumflow/2p/sequential/test_impesadaptiverestart.o \
./test/porousmediumflow/2p/sequential/test_impeswithamg.o \
./test/porousmediumflow/2p/sequential/test_mpfa2p.o \
./test/porousmediumflow/2p/sequential/test_transport.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/2p/sequential/%.o: ../test/porousmediumflow/2p/sequential/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


