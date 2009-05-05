#include "nulliteratorfactory.hh"

template<>
Dune::OneDInNDGridList<Dune::OneDInNDEntityImp<0, 1> > Dune::OneDInNDGridNullIteratorFactory<0, 1>::emptyList_;

template<>
Dune::OneDInNDGridList<Dune::OneDInNDEntityImp<1, 1> > Dune::OneDInNDGridNullIteratorFactory<1, 1>::emptyList_;

template<>
Dune::OneDInNDGridList<Dune::OneDInNDEntityImp<0, 2> > Dune::OneDInNDGridNullIteratorFactory<0, 2>::emptyList_;

template<>
Dune::OneDInNDGridList<Dune::OneDInNDEntityImp<1, 2> > Dune::OneDInNDGridNullIteratorFactory<1, 2>::emptyList_;

template<>
Dune::OneDInNDGridList<Dune::OneDInNDEntityImp<0, 3> > Dune::OneDInNDGridNullIteratorFactory<0, 3>::emptyList_;

template<>
Dune::OneDInNDGridList<Dune::OneDInNDEntityImp<1, 3> > Dune::OneDInNDGridNullIteratorFactory<1, 3>::emptyList_;

template<int dimworld>
Dune::OneDInNDGridList<Dune::OneDInNDEntityImp<0, dimworld> > Dune::OneDInNDGridNullIteratorFactory<0, dimworld>::emptyList_;

template<int dimworld>
Dune::OneDInNDGridList<Dune::OneDInNDEntityImp<1, dimworld> > Dune::OneDInNDGridNullIteratorFactory<1, dimworld>::emptyList_;
