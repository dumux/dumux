#ifndef DUNE_INJECTIONTRANSPORTPROBLEM_HH
#define DUNE_INJECTIONTRANSPORTPROBLEM_HH

#include "dumux/transport/transportproblem.hh"

namespace Dune
{
  //! \ingroup transportProblems
  //! @brief example class for a transport problem
  template<class G, class RT, class VC>
  class InjectionTransportProblem
    : public TransportProblem<G, RT, VC> {

    typedef typename G::ctype DT;
    enum {n=G::dimension, m=1, blocksize=2*G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef Dune::FieldVector<double, n> R1;

  private:
		FieldVector<DT,n> outerLowerLeft_;
		FieldVector<DT,n> outerUpperRight_;
		FieldVector<DT,n> innerLowerLeft_;
		FieldVector<DT,n> innerUpperRight_;
    RT eps_;
	RT outerPorosity_, innerPorosity_;

  public:
    BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e,
				      const FieldVector<DT,n>& xi) const
    {
		if (x[0] < outerLowerLeft_[0] + eps_)
			return BoundaryConditions::dirichlet;
		// all other boundaries
		return BoundaryConditions::neumann;
    }

    RT g (const FieldVector<DT,n>& x, const Entity& e,
	  const FieldVector<DT,n>& xi) const
    {
   	if (x[0] < outerLowerLeft_[0] + eps_)
    		return 1;
  	return 0;
    }

    RT initSat (const FieldVector<DT,n>& x, const Entity& e,
	   const FieldVector<DT,n>& xi) const
    {
    	return 1;
    }


    RT porosity () const {
//		if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
//		    && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
			return innerPorosity_;
//		else
//			return outerPorosity_;
    }

    InjectionTransportProblem(VC& variableobj, TwoPhaseRelations& law = *(new LinearLaw),
			const FieldVector<DT,n> outerLowerLeft = 0., const FieldVector<DT,n> outerUpperRight = 0.,
			const FieldVector<DT,n> innerLowerLeft = 0., const FieldVector<DT,n> innerUpperRight = 0.,
			   const int level = 0, const bool cap =
			   false,RT outerPorosity = 0.4, RT innerPorosity = 0.4)
      : TransportProblem<G, RT, VC>(variableobj,law, cap),
        outerLowerLeft_(outerLowerLeft), outerUpperRight_(outerUpperRight),
        innerLowerLeft_(innerLowerLeft), innerUpperRight_(innerUpperRight),
      	eps_(1e-8), outerPorosity_(outerPorosity), innerPorosity_(innerPorosity)
          {}
  };

}
#endif
