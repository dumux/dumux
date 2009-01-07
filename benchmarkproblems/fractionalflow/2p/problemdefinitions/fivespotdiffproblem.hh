#ifndef FIVESPOTDIFFPROBLEM_HH
#define FIVESPOTDIFFPROBLEM_HH

#include "dumux/diffusion/problems/homogeneousproblem.hh"

namespace Dune
{
  //! \ingroup diffusionProblems
  //! example class for diffusion problems
  template<class G,class RT,class VC>
  class Fivespotcase1DiffProblem : public HomogeneousProblem<G,RT,VC>
  {
    typedef typename G::ctype DT;
    enum {n=G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;

  public:
    Fivespotcase1DiffProblem(VC& variableobj, TwoPhaseRelations& law = *(new LinearLaw), RT bcf = 11,
			     RT pleftbc=2e5, RT prightbc=1.999986e5, const bool cap = false)
      : HomogeneousProblem<G,RT, VC>(variableobj, law, cap),
        LowerLeft_(variableobj.grid.lowerLeft()), UpperRight_(variableobj.grid.upperRight()),
        eps_(1e-8*UpperRight_[0]),bcf_(bcf),
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
      if ((x[0] < LowerLeft_[0] + eps_ && x[1] < LowerLeft_[1] + bcf_) || //lower left
	  (x[1] < LowerLeft_[1] + eps_ && x[0] < LowerLeft_[0] + bcf_))    //lower left
	return BoundaryConditions::dirichlet;

      // all other boundaries
      return BoundaryConditions::neumann;
    }

    RT dirichletPress (const FieldVector<DT,n>& x, const Entity& e,
	  const FieldVector<DT,n>& xi) const
    {
      if ((x[0] < LowerLeft_[0] + eps_ && x[1] < LowerLeft_[1] + bcf_) ||  //lower left
	  (x[1] < LowerLeft_[1] + eps_ && x[0] < LowerLeft_[0] + bcf_))     //lower left
	return pleftbc_;
      // all other boundaries
      return prightbc_;
    }


    RT dirichletSat (const FieldVector<DT,n>& x, const Entity& e,
	  const FieldVector<DT,n>& xi) const
    {
      if ((x[0] < LowerLeft_[0] + eps_ && x[1] < LowerLeft_[1] + bcf_) ||
	  (x[1] < LowerLeft_[1] + eps_ && x[0] < LowerLeft_[0] + bcf_))
	return 0.8;
      else
	return 0.2;
    }

    RT neumannPress (const FieldVector<DT,n>& x, const Entity& e,
	  const FieldVector<DT,n>& xi) const
    {
      if ((x[0] > UpperRight_[0] - eps_ && x[1] > UpperRight_[1] - bcf_) ||
	  (x[1] > UpperRight_[1] - eps_ && x[0] > UpperRight_[0] - bcf_))
	//return 1e-6;//15x15 cells: bcf 21
	return 2e-6; //30x30 cells: bcf 11
	//return 4e-6;// 60x60 cells: bcf 6
      //all other boundaries
      return 0;
    }
  private:
    FieldVector<DT,n> LowerLeft_;
    FieldVector<DT,n> UpperRight_;

    RT eps_;
    RT bcf_;
    RT pleftbc_, prightbc_;
  };
  template<class G,class RT, class VC>
  class Fivespotcase2DiffProblem : public HomogeneousProblem<G,RT,VC>
  {
    typedef typename G::ctype DT;
    enum {n=G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;

  public:
    Fivespotcase2DiffProblem(VC& variableobj, TwoPhaseRelations& law = *(new LinearLaw), RT bcf = 11,
			     RT pleftbc=2e5, RT prightbc=0, const bool cap = false)
      : HomogeneousProblem<G,RT,VC>(variableobj, law, cap),
        LowerLeft_(variableobj.grid.lowerLeft()), UpperRight_(variableobj.grid.upperRight()),
        eps_(1e-8*UpperRight_[0]),bcf_(bcf),
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
      if ((x[0] < LowerLeft_[0] + eps_ && x[1] < LowerLeft_[1] + bcf_) || //lower left
	  (x[1] <  LowerLeft_[1] + eps_ && x[0] < LowerLeft_[0] + bcf_) || //lower left
	  (x[0] > UpperRight_[0] - eps_ && x[1] > UpperRight_[1] - bcf_) || //upper right
	  (x[1] > UpperRight_[1] - eps_ && x[0] > UpperRight_[0] - bcf_))  //upper right
	return BoundaryConditions::dirichlet;
      // all other boundaries
      return BoundaryConditions::neumann;
    }

    RT dirichletPress (const FieldVector<DT,n>& x, const Entity& e,
	  const FieldVector<DT,n>& xi) const
    {
      if ((x[0] < LowerLeft_[0] + eps_ && x[1] < LowerLeft_[1] + bcf_) || //lower left
	  (x[1] <  LowerLeft_[1] + eps_ && x[0] < LowerLeft_[0] + bcf_) || //lower left
	  (x[0] > UpperRight_[0] - eps_ && x[1] > UpperRight_[1] - bcf_) || //upper right
	  (x[1] > UpperRight_[1] - eps_ && x[0] > UpperRight_[0] - bcf_))  //upper right
	return pleftbc_;
      // all other boundaries
      return prightbc_;
    }

    RT dirichletSat (const FieldVector<DT,n>& x, const Entity& e,
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

    RT neumannPress (const FieldVector<DT,n>& x, const Entity& e,
	  const FieldVector<DT,n>& xi) const
    {
      if ((x[0] < LowerLeft_[0] + eps_ && x[1] > UpperRight_[1] - bcf_) || //upper left
	  (x[1] >  UpperRight_[1] - eps_ && x[0] < LowerLeft_[0] + bcf_) || //upper left
	  (x[0] > UpperRight_[0] - eps_ && x[1] < LowerLeft_[1] + bcf_) ||  //lower right
	  (x[1] < LowerLeft_[1] + eps_ && x[0] > UpperRight_[0] - bcf_))   //lower right
	//return 1e-6;//15x15 cells: bcf 21
	return 2e-6; //30x30 cells: bcf 11
	//return 4e-6;// 60x60 cells: bcf 6;
      // all other boundaries
      return 0;
    }
  private:
    FieldVector<DT,n> LowerLeft_;
    FieldVector<DT,n> UpperRight_;

    RT eps_;
    RT bcf_;
    RT pleftbc_, prightbc_;
  };

}

#endif
