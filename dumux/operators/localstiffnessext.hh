#ifndef DUNE_LOCALSTIFFNESSEXT_HH
#define DUNE_LOCALSTIFFNESSEXT_HH

#include <dune/disc/operators/localstiffness.hh>

namespace Dune
{
  template<class Imp, class G, class RT, int m>
  class LocalStiffnessExt : public LocalStiffness<Imp, G, RT, m> 
  {
    // grid types
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    enum {n=G::dimension};
    
  public:
    template<class TypeTag>
    void assembleElementMatrices(const Entity& e, Dune::FieldVector<DT,2*n>& faceVol, 
				 Dune::FieldMatrix<RT,2*n,2*n>& W, Dune::FieldVector<DT,2*n>& c, 
				 Dune::FieldMatrix<RT,2*n,2*n>& Pi, RT& dinv, Dune::FieldVector<DT,2*n>& F, RT& qmean) {
      this->getImp().template assembleElementMatrices<TypeTag>(e, faceVol, W, c, Pi, dinv, F, qmean);
    }
  };
}
#endif
