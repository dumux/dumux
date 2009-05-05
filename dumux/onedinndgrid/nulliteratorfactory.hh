#ifndef DUNE_ONEDINNDGRID_NULL_ITERATORS_HH
#define DUNE_ONEDINNDGRID_NULL_ITERATORS_HH

#include "onedinndgridlist.hh"

namespace Dune {

  template <int mydim, int dimworld>
  class OneDInNDEntityImp;

  template <int mydim, int dimworld>
    class OneDInNDGridNullIteratorFactory {};

    template <int dimworld>
    class OneDInNDGridNullIteratorFactory<0, dimworld> {

    public:

      static typename OneDInNDGridList<OneDInNDEntityImp<0, dimworld> >::iterator null() {
            return emptyList_.end();
        }

    private:
      static OneDInNDGridList<OneDInNDEntityImp<0, dimworld> > emptyList_;
    };

    template <int dimworld>
    class OneDInNDGridNullIteratorFactory<1, dimworld> {

    public:

      static typename OneDInNDGridList<OneDInNDEntityImp<1, dimworld> >::iterator null() {
            return emptyList_.end();
        }

    private:
      static OneDInNDGridList<OneDInNDEntityImp<1, dimworld> > emptyList_;
    };

}

#endif
