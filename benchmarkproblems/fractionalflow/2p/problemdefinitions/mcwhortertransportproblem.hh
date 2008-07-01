#ifndef DUNE_MCWHORTERTRANSPORTPROBLEM_HH
#define DUNE_MCWHORTERTRANSPORTPROBLEM_HH

#include "dumux/transport/transportproblem.hh"

namespace Dune
{
  //! \ingroup transportProblems
  //! @brief example class for a transport problem
  template<class G, class RT, class VC>
  class McWhorterTransportProblem 
    : public TransportProblem<G, RT, VC> {
		  
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

    RT g (const FieldVector<DT,n>& x, const Entity& e, 
	  const FieldVector<DT,n>& xi) const 
    {
      if (x[0] < eps_) 
	return 1;
      else
	return 0;
    }
	  
    RT S0 (const FieldVector<DT,n>& x, const Entity& e, 
	   const FieldVector<DT,n>& xi) const 
    {
      if (x[0]< 0.1)
      return (1-Sinit_-(1-2*Sinit_)/0.1*x[0]);
       
      return Sinit_;
    }

  
    RT porosity () const {
      return poro_;
    }

    McWhorterTransportProblem(VC& variableobj, TwoPhaseRelations& law = *(new LinearLaw),
		     const FieldVector<DT,n> Left = 0, const FieldVector<DT,n> Right = 0,
		     bool analytic = false, const bool cap = false, const int level = 0, RT poro=0.3,RT Si=0.0) 
      : TransportProblem<G, RT, VC>(variableobj,law, cap), left(Left[0]), right(Right[0]),
	eps_(1e-8),
	poro_(poro),
	Sinit_(Si)
    {}
		     
  private:
    DT left;
    DT right;
    RT eps_;
    RT poro_;
    RT Sinit_;
  };
}
#endif
