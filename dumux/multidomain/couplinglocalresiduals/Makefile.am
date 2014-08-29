couplinglocalresidualsdir = $(includedir)/dumux/multidomain/couplinglocalresiduals
couplinglocalresiduals_HEADERS := $(wildcard *.hh)

if !HAVE_DUNE_MULTIDOMAIN
headercheck_IGNORE = \
  2p2ccouplinglocalresidual.hh \
  boxcouplinglocalresidual.hh \
  stokesncnicouplinglocalresidual.hh \
  2p2cnicouplinglocalresidual.hh \
  stokesnccouplinglocalresidual.hh
endif

include $(top_srcdir)/am/global-rules
