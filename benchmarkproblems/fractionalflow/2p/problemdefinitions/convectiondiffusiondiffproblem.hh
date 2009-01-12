#ifndef CONVECTIONDIFFUSIONDIFFPROBLEM_HH
#define CONVECTIONDIFFUSIONDIFFPROBLEM_HH

#include "dumux/diffusion/problems/homogeneousproblem.hh"

namespace Dune
{
  //! \ingroup diffusionProblems
  //! example class for diffusion problems
  template<class G,class RT, class VC>
  class ConvectionDiffusionDiffProblem : public HomogeneousProblem<G,RT,VC>
  {
    typedef typename G::ctype DT;
    enum {n=G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;

  public:
    ConvectionDiffusionDiffProblem(VC& variableobj, TwoPhaseRelations& law = *(new LinearLaw),
                   const FieldVector<DT,n> Left = 0, const FieldVector<DT,n> Right = 0,
                   RT pleftbc=200005.3, RT prightbc=1.999986e5, const bool cap = false)
      : HomogeneousProblem<G,RT,VC>(variableobj, law, cap),
        Left_(Left[0]), Right_(Right[0]),
        eps_(1e-8),
        pleftbc_(pleftbc),prightbc_(prightbc)
    {}

                   RT source (const FieldVector<DT,n>& x, const Entity& e,
                     const FieldVector<DT,n>& xi)
                   {
                       return 0;
                   }

    typename BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e,
                           const FieldVector<DT,n>& xi) const
    {
      if (x[0] < eps_)//  || x[0] > (Right_ - eps_))
    return BoundaryConditions::dirichlet;
      // all other boundaries
      return BoundaryConditions::neumann;
    }

    RT dirichletPress (const FieldVector<DT,n>& x, const Entity& e,
      const FieldVector<DT,n>& xi) const
    {
      if (x[0] <  eps_)
    return pleftbc_;
      // all other boundaries
      return prightbc_;
    }

    RT dirichletSat(const FieldVector<DT,n>& x, const Entity& e,
            const FieldVector<DT,n>& xi) const {
        if (x[0] < eps_)
            return 0.9;
        // all other boundaries
        return 0;
    }

    RT neumannPress (const FieldVector<DT,n>& x, const Entity& e,
      const FieldVector<DT,n>& xi) const
    {
      if (x[0] > Right_ - eps_) return 3e-7;
      //if (x[0] < eps_) return -3e-7;
      return 0;
    }
  private:
    DT Left_;
    DT Right_;

    RT eps_;
    RT bcf_;
    RT pleftbc_, prightbc_;
  };
}

#endif
