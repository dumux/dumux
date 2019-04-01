################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../test/porousmediumflow/2p/implicit/test_box2p.cc \
../test/porousmediumflow/2p/implicit/test_box2pni.cc \
../test/porousmediumflow/2p/implicit/test_boxadaptive2p.cc \
../test/porousmediumflow/2p/implicit/test_cc2p.cc \
../test/porousmediumflow/2p/implicit/test_cc2pcornerpoint.cc \
../test/porousmediumflow/2p/implicit/test_cc2pni.cc \
../test/porousmediumflow/2p/implicit/test_ccadaptive2p.cc \
../test/porousmediumflow/2p/implicit/test_fracture_box2p.cc \
../test/porousmediumflow/2p/implicit/test_generalizeddirichlet.cc 

CC_DEPS += \
./test/porousmediumflow/2p/implicit/test_box2p.d \
./test/porousmediumflow/2p/implicit/test_box2pni.d \
./test/porousmediumflow/2p/implicit/test_boxadaptive2p.d \
./test/porousmediumflow/2p/implicit/test_cc2p.d \
./test/porousmediumflow/2p/implicit/test_cc2pcornerpoint.d \
./test/porousmediumflow/2p/implicit/test_cc2pni.d \
./test/porousmediumflow/2p/implicit/test_ccadaptive2p.d \
./test/porousmediumflow/2p/implicit/test_fracture_box2p.d \
./test/porousmediumflow/2p/implicit/test_generalizeddirichlet.d 

OBJS += \
./test/porousmediumflow/2p/implicit/test_box2p.o \
./test/porousmediumflow/2p/implicit/test_box2pni.o \
./test/porousmediumflow/2p/implicit/test_boxadaptive2p.o \
./test/porousmediumflow/2p/implicit/test_cc2p.o \
./test/porousmediumflow/2p/implicit/test_cc2pcornerpoint.o \
./test/porousmediumflow/2p/implicit/test_cc2pni.o \
./test/porousmediumflow/2p/implicit/test_ccadaptive2p.o \
./test/porousmediumflow/2p/implicit/test_fracture_box2p.o \
./test/porousmediumflow/2p/implicit/test_generalizeddirichlet.o 


# Each subdirectory must supply rules for building sources it contributes
test/porousmediumflow/2p/implicit/%.o: ../test/porousmediumflow/2p/implicit/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


