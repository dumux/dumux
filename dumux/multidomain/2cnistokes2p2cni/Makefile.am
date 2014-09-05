2cnistokes2p2cnidir = $(includedir)/dumux/multidomain/2cnistokes2p2cni
2cnistokes2p2cni_HEADERS := $(wildcard *.hh)

if !HAVE_DUNE_MULTIDOMAIN
headercheck_IGNORE = 2cnistokes2p2cnilocaloperator.hh
endif

include $(top_srcdir)/am/global-rules
