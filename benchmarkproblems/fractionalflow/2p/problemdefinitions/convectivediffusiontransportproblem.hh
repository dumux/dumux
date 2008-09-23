#ifndef DUNE_CONVECTIONDIFFUSIONTRANSPORTPROBLEM_HH
#define DUNE_CONVECTIONDIFFUSIONTRANSPORTPROBLEM_HH

#include "dumux/transport/transportproblem.hh"

namespace Dune
{
  //! \ingroup transportProblems
  //! @brief example class for a transport problem
  template<class G, class RT, class VC>
  class ConvectionDiffusionTransportProblem
    : public TransportProblem<G, RT, VC> {

    typedef typename G::ctype DT;
    enum {n=G::dimension, m=1, blocksize=2*G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef Dune::FieldVector<double, n> R1;

  private:
    DT left;
    DT right;
    RT eps_;
    RT poro_;

  public:
    BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e,
				      const FieldVector<DT,n>& xi) const
    {
      if (x[0] < eps_)// (x[0] > (right - eps_) || x[0] < eps_)
	return Dune::BoundaryConditions::dirichlet;
      else
	return Dune::BoundaryConditions::neumann;
    }

    RT g (const FieldVector<DT,n>& x, const Entity& e,
	  const FieldVector<DT,n>& xi) const
    {
      if (x[0] < eps_)
	return 0.9;
      else
	return 0;
    }

    RT initSat (const FieldVector<DT,n>& x, const Entity& e,
	   const FieldVector<DT,n>& xi) const
    {
      //if (x[0] < eps_)
      //	return 0.9;
      //if (x[0] < 20)
      //	return (0.9)/10*x[0];
      //else
	return 0;
    }


    RT porosity () const {
      return poro_;
    }

    ConvectionDiffusionTransportProblem(VC& variableobj, TwoPhaseRelations& law = *(new LinearLaw),
			   const FieldVector<DT,n> Left = 0, const FieldVector<DT,n> Right = 0,
			   const int level = 0, const bool cap =
			   false,RT poro=0.2)
      : TransportProblem<G, RT, VC>(variableobj,law, cap), left(Left[0]), right(Right[0]),
	eps_(1e-8),
	poro_(poro)
    {}
  };

}
#endif
