2cstokes2p2cdir = $(includedir)/dumux/multidomain/2cstokes2p2c
2cstokes2p2c_HEADERS = *.hh

if !HAVE_DUNE_MULTIDOMAIN
headercheck_IGNORE = \
  2cstokes2p2clocaloperator.hh \
  2cstokes2p2cnewtoncontroller.hh
endif

include $(top_srcdir)/am/global-rules
