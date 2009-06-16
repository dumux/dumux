#ifndef DUNE_ONEDINNDGRID_NULL_ITERATORS_HH
#define DUNE_ONEDINNDGRID_NULL_ITERATORS_HH

#include "onedinndgridlist.hh"

namespace Dune {

  template <int mydim, int dimworld>
  class OneDInNDEntityImp;

  template <int dim, int dimworld>
    class OneDInNDGridNullIteratorFactory {};

    template <>
    class OneDInNDGridNullIteratorFactory<0, 1> {

    public:

      static OneDInNDGridList<OneDInNDEntityImp<0, 1> >::iterator null() {
            return emptyList_.end();
        }

    private:
      static OneDInNDGridList<OneDInNDEntityImp<0, 1> > emptyList_;
    };

    template <>
    class OneDInNDGridNullIteratorFactory<0, 2> {

    public:

      static OneDInNDGridList<OneDInNDEntityImp<0, 2> >::iterator null() {
            return emptyList_.end();
        }

    private:
      static OneDInNDGridList<OneDInNDEntityImp<0, 2> > emptyList_;
    };

    template <>
    class OneDInNDGridNullIteratorFactory<0, 3> {

    public:

      static OneDInNDGridList<OneDInNDEntityImp<0, 3> >::iterator null() {
            return emptyList_.end();
        }

    private:
      static OneDInNDGridList<OneDInNDEntityImp<0, 3> > emptyList_;
    };

    template <>
    class OneDInNDGridNullIteratorFactory<1, 1> {

    public:

      static OneDInNDGridList<OneDInNDEntityImp<1, 1> >::iterator null() {
            return emptyList_.end();
        }

    private:
      static OneDInNDGridList<OneDInNDEntityImp<1, 1> > emptyList_;
    };

    template <>
    class OneDInNDGridNullIteratorFactory<1, 2> {

    public:

      static OneDInNDGridList<OneDInNDEntityImp<1, 2> >::iterator null() {
            return emptyList_.end();
        }

    private:
      static OneDInNDGridList<OneDInNDEntityImp<1, 2> > emptyList_;
    };

    template <>
    class OneDInNDGridNullIteratorFactory<1, 3> {

    public:

      static OneDInNDGridList<OneDInNDEntityImp<1, 3> >::iterator null() {
            return emptyList_.end();
        }

    private:
      static OneDInNDGridList<OneDInNDEntityImp<1, 3> > emptyList_;
    };

}

#endif
