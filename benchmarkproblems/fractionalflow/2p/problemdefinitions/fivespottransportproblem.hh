#ifndef DUNE_FIVESPOTTRANSPORTPROBLEM_HH
#define DUNE_FIVESPOTTRANSPORTPROBLEM_HH

#include "dumux/transport/transportproblem.hh"

namespace Dune
{
  //! \ingroup transportProblems
  //! @brief example class for a transport problem
  template<class G, class RT, class VC>
  class Fivespotcase1TransportProblem
    : public TransportProblem<G, RT,VC> {

    typedef typename G::ctype DT;
    enum {n=G::dimension, m=1, blocksize=2*G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef Dune::FieldVector<double, n> R1;

  private:
    FieldVector<DT,n>  LowerLeft_;
    FieldVector<DT,n> UpperRight_;
    RT eps_;
    RT bcf_;
    RT poro_;


  public:
    BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e,
				      const FieldVector<DT,n>& xi) const
    {
      if ((x[0] < LowerLeft_[0] + eps_ && x[1] < LowerLeft_[1] + bcf_) ||
	  (x[1] < LowerLeft_[1] + eps_ && x[0] < LowerLeft_[0] + bcf_) ||
	  (x[0] > UpperRight_[0] - eps_ && x[1] > UpperRight_[1] - bcf_) ||
	  (x[1] > UpperRight_[1] - eps_ && x[0] > UpperRight_[0] - bcf_))
	return BoundaryConditions::dirichlet;
      // all other boundaries
      return BoundaryConditions::neumann;
    }

    RT g (const FieldVector<DT,n>& x, const Entity& e,
	  const FieldVector<DT,n>& xi) const
    {
      if ((x[0] < LowerLeft_[0] + eps_ && x[1] < LowerLeft_[1] + bcf_) ||
	  (x[1] < LowerLeft_[1] + eps_ && x[0] < LowerLeft_[0] + bcf_))
	return 0.8;
      else
	return 0.2;
    }

    RT initSat (const FieldVector<DT,n>& x, const Entity& e,
	   const FieldVector<DT,n>& xi) const
    {
      return 0.2;
    }

    RT porosity () const
    {
      return poro_;
    }


    Fivespotcase1TransportProblem(VC& variableobj, TwoPhaseRelations& law = *(new LinearLaw), RT bcf = 11,
			 const int level = 0, const bool cap =
			 false)
      : TransportProblem<G, RT, VC>(variableobj,law, cap), LowerLeft_(variableobj.grid.lowerLeft()), UpperRight_(variableobj.grid.upperRight()),
	eps_(1e-8*(variableobj.grid.upperRight())[0]),poro_(0.2),bcf_(bcf)
    {}
  };
  template<class G, class RT,class VC>
  class Fivespotcase2TransportProblem
    : public TransportProblem<G, RT, VC> {

    typedef typename G::ctype DT;
    enum {n=G::dimension, m=1, blocksize=2*G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef Dune::FieldVector<double, n> R1;

  private:
    FieldVector<DT,n>  LowerLeft_;
    FieldVector<DT,n> UpperRight_;
    RT eps_;
    RT bcf_;
    RT poro_;


  public:
    BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e,
				      const FieldVector<DT,n>& xi) const
    {
      if ((x[0] < LowerLeft_[0] + eps_ && x[1] < LowerLeft_[1] + bcf_) || //lower left
	  (x[1] <  LowerLeft_[1] + eps_ && x[0] < LowerLeft_[0] + bcf_) || //lower left
	  (x[0] > UpperRight_[0] - eps_ && x[1] > UpperRight_[1] - bcf_) || //upper right
	  (x[1] > UpperRight_[1] - eps_ && x[0] > UpperRight_[0] - bcf_) || //upper right
	  (x[0] < LowerLeft_[0] + eps_ && x[1] > UpperRight_[1] - bcf_) || //upper left
	  (x[1] >  UpperRight_[1] - eps_ && x[0] < LowerLeft_[0] + bcf_) || //upper left
	  (x[0] > UpperRight_[0] - eps_ && x[1] < LowerLeft_[1] + bcf_) ||  //lower right
	  (x[1] < LowerLeft_[1] + eps_ && x[0] > UpperRight_[0] - bcf_))   //lower right
	return BoundaryConditions::dirichlet;
      // all other boundaries
      return BoundaryConditions::neumann;
    }

    RT g (const FieldVector<DT,n>& x, const Entity& e,
	  const FieldVector<DT,n>& xi) const
    {
      if ((x[0] < LowerLeft_[0] + eps_ && x[1] < LowerLeft_[1] + bcf_) || //lower left
	  (x[1] <  LowerLeft_[1] + eps_ && x[0] < LowerLeft_[0] + bcf_) || //lower left
	  (x[0] > UpperRight_[0] - eps_ && x[1] > UpperRight_[1] - bcf_) || //upper right
	  (x[1] > UpperRight_[1] - eps_ && x[0] > UpperRight_[0] - bcf_))  //upper right
	return 0.8;
      // all other boundaries
      return 0.2;
    }

    RT initSat (const FieldVector<DT,n>& x, const Entity& e,
	   const FieldVector<DT,n>& xi) const
    {
      return 0.2;
    }

    RT porosity () const
    {
      return poro_;
    }


    Fivespotcase2TransportProblem(VC& variableobj, TwoPhaseRelations& law = *(new LinearLaw), RT bcf = 11,
			 const int level = 0, const bool cap =
			 false)
      : TransportProblem<G, RT, VC>(variableobj,law, cap), LowerLeft_(variableobj.grid.lowerLeft()), UpperRight_(variableobj.grid.upperRight()),
	eps_(1e-8*(variableobj.grid.upperRight())[0]),poro_(0.2),bcf_(bcf)
    {}
  };
}
#endif
