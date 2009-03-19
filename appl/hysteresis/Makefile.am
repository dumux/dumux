# tests where program to build and program to run are equal
NORMALTESTS = test_pc_hysteresis_lenhard1 test_pc_hysteresis_lenhard2 test_pc_hysteresis_lens

# list of tests to run
#TESTS = $(NORMALTESTS)

# programs just to build when "make check" is used
check_PROGRAMS = $(NORMALTESTS)

dist_noinst_DATA = 	lenhardexp1main.cc \
	           	lenhardexp2main.cc \
	           	lenhardmain.cc \
	           	lensmain.cc \
			lenhardconfig.hh \
			lenharddomain.hh \
			lenhardnewtoncontroller.hh \
			lenhardproblem.hh \
			lenhardstate.hh \
			lensconfig.hh \
			lensdomain.hh \
			lensnewtoncontroller.hh \
			lensproblem.hh \
			lensstate.hh

test_pc_hysteresis_lenhard1_SOURCES = lenhardexp1main.cc
test_pc_hysteresis_lenhard2_SOURCES = lenhardexp2main.cc
test_pc_hysteresis_lens_SOURCES = lensmain.cc

include $(top_srcdir)/am/global-rules
