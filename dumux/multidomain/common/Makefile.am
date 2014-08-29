commondir = $(includedir)/dumux/multidomain/common
common_HEADERS := $(wildcard *.hh)

if !HAVE_DUNE_MULTIDOMAIN
headercheck_IGNORE = \
  multidomainassembler.hh \
  multidomainmodel.hh \
  multidomainproperties.hh \
  splitandmerge.hh \
  multidomainconvergencewriter.hh \
  multidomainnewtoncontroller.hh \
  multidomainpropertydefaults.hh \
  subdomainproperties.hh \
  multidomainlocaloperator.hh \
  multidomainproblem.hh \
  pdelablocaloperator.hh \
  subdomainpropertydefaults.hh
endif

include $(top_srcdir)/am/global-rules
