SUBDIRS = 2p

check_PROGRAMS = test_effectivediffusivityconstanttau \
                 test_effectivediffusivitymillingtonquirk

test_effectivediffusivityconstanttau_SOURCES = test_effectivediffusivityconstanttau.cc
test_effectivediffusivitymillingtonquirk_SOURCES = test_effectivediffusivitymillingtonquirk.cc

EXTRA_DIST=CMakeLists.txt
noinst_HEADERS:=$(wildcard *.hh)
include $(top_srcdir)/am/global-rules