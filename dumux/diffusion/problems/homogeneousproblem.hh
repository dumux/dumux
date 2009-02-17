// $Id$

#ifndef HOMOGENEOUSPROBLEM_HH
#define HOMOGENEOUSPROBLEM_HH

#include "dumux/diffusion/diffusionproblem_deprecated.hh"

namespace Dune
{
  //! \ingroup diffusionProblems
  //! example class for diffusion problems
  template<class G, class RT, class VC>
  class HomogeneousProblem : public DiffusionProblem<G,RT,VC>
  {
    typedef typename G::ctype DT;
    enum {n=G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;

  public:
    HomogeneousProblem(VC& variableobj, TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false,
               RT K=1e-7)
      : DiffusionProblem<G,RT,VC>(variableobj,law, cap)
    {
        if (n == 1) {
            K_=K;
        }
        else if (n == 2) {
            K_[0][0]=K_[1][1]=K;
            K_[0][1]=K_[1][0]=0;
    }
    }

    virtual FieldMatrix<DT,n,n>& K(const FieldVector<DT,n>& x,
            const Entity& e, const FieldVector<DT,n>& xi) {
        return K_;
    }

    virtual RT source   (const FieldVector<DT,n>& x, const Entity& e,
            const FieldVector<DT,n>& xi)
    {
      return 0;
    }

    virtual typename BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e,
                               const FieldVector<DT,n>& xi) const
    {
      if (x[0] > 300 - 1e-6 || x[0] < 1e6)
    return BoundaryConditions::dirichlet;
      // all other boundaries
      return BoundaryConditions::neumann;
    }

    virtual RT dirichletPress (const FieldVector<DT,n>& x, const Entity& e,
          const FieldVector<DT,n>& xi) const
    {
      return (x[0] < 1e-6) ? 2e5 : 1.9e5;
    }

    virtual RT dirichletSat (const FieldVector<DT,n>& x, const Entity& e,
          const FieldVector<DT,n>& xi) const
    {
      return (x[0] < 1e-6) ? 1 : 0;
    }

    virtual RT neumannPress (const FieldVector<DT,n>& x, const Entity& e,
          const FieldVector<DT,n>& xi) const
    {
      return 0;
    }


  private:
    FieldMatrix<DT,n,n> K_;
//    G& grid;
    };
}

#endif
