################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include tutorial/tutorial_sequential/subdir.mk
-include tutorial/tutorial_implicit/subdir.mk
-include tutorial/solution/ex1/subdir.mk
-include tutorial/ex3/subdir.mk
-include tutorial/ex2/subdir.mk
-include tutorial/ex1/subdir.mk
-include test/porousmediumflow/richards/implicit/subdir.mk
-include test/porousmediumflow/mpnc/implicit/subdir.mk
-include test/porousmediumflow/co2/implicit/subdir.mk
-include test/porousmediumflow/3pwateroil/implicit/subdir.mk
-include test/porousmediumflow/3p3c/implicit/subdir.mk
-include test/porousmediumflow/3p/implicit/subdir.mk
-include test/porousmediumflow/2pncmin/implicit/subdir.mk
-include test/porousmediumflow/2pnc/implicit/subdir.mk
-include test/porousmediumflow/2pminc/implicit/subdir.mk
-include test/porousmediumflow/2pdfm/implicit/subdir.mk
-include test/porousmediumflow/2p2c/sequential/subdir.mk
-include test/porousmediumflow/2p2c/implicit/subdir.mk
-include test/porousmediumflow/2p1c/implicit/subdir.mk
-include test/porousmediumflow/2p/sequential/subdir.mk
-include test/porousmediumflow/2p/implicit/pointsources/subdir.mk
-include test/porousmediumflow/2p/implicit/subdir.mk
-include test/porousmediumflow/1p2c/implicit/subdir.mk
-include test/porousmediumflow/1p/sequential/subdir.mk
-include test/porousmediumflow/1p/implicit/pointsources/subdir.mk
-include test/porousmediumflow/1p/implicit/subdir.mk
-include test/multidomain/2czeroeq2p2c/subdir.mk
-include test/multidomain/2cstokes2p2c/subdir.mk
-include test/multidomain/2cnizeroeq2p2cni/subdir.mk
-include test/multidomain/2cnistokes2p2cni/subdir.mk
-include test/material/tabulation/subdir.mk
-include test/material/pengrobinson/subdir.mk
-include test/material/ncpflash/subdir.mk
-include test/material/immiscibleflash/subdir.mk
-include test/material/fluidsystems/subdir.mk
-include test/material/fluidmatrixinteractions/2p/subdir.mk
-include test/material/fluidmatrixinteractions/subdir.mk
-include test/io/gridcreator/subdir.mk
-include test/io/gnuplotinterface/subdir.mk
-include test/io/container/subdir.mk
-include test/geomechanics/elastic/subdir.mk
-include test/geomechanics/el2p/subdir.mk
-include test/geomechanics/el1p2c/subdir.mk
-include test/freeflow/zeroeq2cni/subdir.mk
-include test/freeflow/zeroeq2c/subdir.mk
-include test/freeflow/zeroeq/subdir.mk
-include test/freeflow/stokes2cni/subdir.mk
-include test/freeflow/stokes2c/subdir.mk
-include test/freeflow/stokes/subdir.mk
-include test/freeflow/navierstokes/subdir.mk
-include test/common/spline/subdir.mk
-include test/common/propertysystem/subdir.mk
-include test/common/generalproblem/subdir.mk
-include test/common/boundingboxtree/subdir.mk
-include build-cmake/test/geomechanics/el2p/CMakeFiles/test_el2p.dir/subdir.mk
-include build-cmake/CMakeFiles/FortranCInterface/symbols/subdir.mk
-include build-cmake/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/VerifyFortranC.dir/subdir.mk
-include build-cmake/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/VerifyFortran.dir/subdir.mk
-include build-cmake/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface/symbols/subdir.mk
-include build-cmake/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface/CMakeFiles/symbols.dir/symbols/subdir.mk
-include build-cmake/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface/CMakeFiles/symbols.dir/subdir.mk
-include build-cmake/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface/CMakeFiles/myfort.dir/subdir.mk
-include build-cmake/CMakeFiles/FortranCInterface/VerifyCXX/CMakeFiles/FortranCInterface/CMakeFiles/FortranCInterface.dir/subdir.mk
-include build-cmake/CMakeFiles/FortranCInterface/CMakeFiles/symbols.dir/symbols/subdir.mk
-include build-cmake/CMakeFiles/FortranCInterface/CMakeFiles/symbols.dir/subdir.mk
-include build-cmake/CMakeFiles/FortranCInterface/CMakeFiles/myfort.dir/subdir.mk
-include build-cmake/CMakeFiles/FortranCInterface/CMakeFiles/FortranCInterface.dir/subdir.mk
-include build-cmake/CMakeFiles/3.5.1/CompilerIdCXX/subdir.mk
-include build-cmake/CMakeFiles/3.5.1/CompilerIdC/subdir.mk
-include build-cmake/CMakeFiles/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: DuMux2.12C

# Tool invocations
DuMux2.12C: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++  -o "DuMux2.12C" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CPP_DEPS)$(C_DEPS) DuMux2.12C
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
