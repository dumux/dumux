#ifndef DUNE_MCWHORTERTRANSPORTPROBLEM_HH
#define DUNE_MCWHORTERTRANSPORTPROBLEM_HH

#include "dumux/transport/transportproblem_deprecated.hh"

namespace Dune
{
  //! \ingroup transportProblems
  //! @brief example class for a transport problem
  template<class G, class RT, class VC>
  class McWhorterTransportProblem
    : public DeprecatedTransportProblem<G, RT, VC> {

    typedef typename G::ctype DT;
    enum {n=G::dimension, m=1};
    bool analytical_;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef Dune::FieldVector<double, n> R1;

  public:
    BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e,
                      const FieldVector<DT,n>& xi) const
    {
      if (x[0] < eps_ || x[0] > right - eps_)
    return Dune::BoundaryConditions::dirichlet;
      else
    return Dune::BoundaryConditions::neumann;
    }

    RT dirichlet (const FieldVector<DT,n>& x, const Entity& e,
      const FieldVector<DT,n>& xi) const
    {
      if (x[0] < eps_)
    return 1;
      else
    return 0;
    }

    RT initSat (const FieldVector<DT,n>& x, const Entity& e,
       const FieldVector<DT,n>& xi) const
    {
        const RT initlength = 2.6/this->variables.grid.size(0);
        if (x[0]< initlength)
            return (1-Sinit_-(1-2*Sinit_)/initlength*x[0]);

        return Sinit_;
    }


    RT porosity () const {
      return poro_;
    }

    McWhorterTransportProblem(VC& variableobj, DeprecatedTwoPhaseRelations& law = *(new DeprecatedLinearLaw),
             const FieldVector<DT,n> Left = 0, const FieldVector<DT,n> Right = 0,
             bool exsol = false, const bool cap = false, const int level = 0, RT poro=0.3,RT Si=0.0)
      : DeprecatedTransportProblem<G, RT, VC>(variableobj,law, cap, exsol), left(Left[0]), right(Right[0]),
    eps_(1e-8),
    poro_(poro),
    Sinit_(Si)
    {}

  private:
    DT left;
    DT right;
    RT eps_;
    RT poro_;
  protected:
    RT Sinit_;
  };
}
#endif
