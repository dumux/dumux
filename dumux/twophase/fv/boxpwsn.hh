#ifndef DUNE_BOXPWSN_HH
#define DUNE_BOXPWSN_HH

#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>

//#include <dune/common/array.hh>        // defines simple array class
#include <dune/common/fixedarray.hh>   // defines simple array classes
#include <dune/common/geometrytype.hh>
#include <dune/grid/sgrid.hh>          // a complete structured grid
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/common/universalmapper.hh>
#include <dune/grid/common/quadraturerules.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/vbvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>
#include <dune/istl/gsetc.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/disc/functions/functions.hh>
#include <dune/disc/functions/p0function.hh>
#include <dune/disc/functions/p1function.hh>
#include <dune/disc/operators/p1operator.hh>
#include <dune/disc/operators/boundaryconditions.hh>
#include <dune/disc/groundwater/groundwater.hh>
#include <dune/disc/groundwater/p1groundwater.hh>
#include <dune/disc/groundwater/p1groundwaterestimator.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/paamg/amg.hh>
#include "dumux/pardiso/pardiso.hh"
#include "dumux/pardiso/identity.hh"
#include "dumux/nonlinear/newtonmethod.hh"
#include "dumux/twophase/twophasemodel.hh"
#include "dumux/twophase/twophaseproblem.hh"
#include "dumux/twophase/fv/boxpwsnjacobian.hh"

namespace Dune
{

/**
 \brief Two phase model with Pw and Sn as primary unknowns
 
 This implements a two phase model with Pw and Sn as primary unknowns.
 */
  template<class G, class RT>
  class BoxPwSn 
  : public LeafP1TwoPhaseModel<G, RT, TwoPhaseProblem<G, RT>, 
                                 BoxPwSnJacobian<G, RT> >
  {
		
		
		 

 

  public:
	// define the problem type (also change the template argument above)
	typedef TwoPhaseProblem<G, RT> ProblemType;
	
	// define the local Jacobian (also change the template argument above)
	typedef BoxPwSnJacobian<G, RT> LocalJacobian;
	
	typedef LeafP1TwoPhaseModel<G, RT, ProblemType, LocalJacobian> LeafP1TwoPhaseModel;
	
    typedef typename LeafP1TwoPhaseModel::FunctionType FunctionType;

    typedef typename G::Traits::LeafIndexSet IS;

    enum{m = 2};

	typedef BoxPwSn<G, RT> ThisType;
		typedef typename LeafP1TwoPhaseModel::FunctionType::RepresentationType VectorType;
		typedef typename LeafP1TwoPhaseModel::OperatorAssembler::RepresentationType MatrixType;
		typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator; 
#ifdef HAVE_PARDISO
	SeqPardiso<MatrixType,VectorType,VectorType> pardiso; 
#endif
	
	BoxPwSn(const G& g, ProblemType& prob) 
	: LeafP1TwoPhaseModel(g, prob)
	{ 	}

	void solve() 
	{
		
		
		 
		
		Operator op(*(this->A));  // make operator out of matrix
		double red=1E-8;

#ifdef HAVE_PARDISO 
//	SeqPardiso<MatrixType,VectorType,VectorType> ilu0(*(this->A)); 
		pardiso.factorize(*(this->A));
		BiCGSTABSolver<VectorType> solver(op,pardiso,red,100,2);         // an inverse operator 
	//	SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
		//LoopSolver<VectorType> solver(op, ilu0, red, 10, 2);
#else
		SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner

		//SeqIdentity<MatrixType,VectorType,VectorType> ilu0(*(this->A));// a precondtioner
		BiCGSTABSolver<VectorType> solver(op,ilu0,red,10000,1);         // an inverse operator 
#endif
		InverseOperatorResult r;
		solver.apply(*(this->u), *(this->f), r);
		
		return;		
	}

	void update (double& dt)
	{
		this->localJacobian.setDt(dt);
		this->localJacobian.setOldSolution(this->uOldTimeStep);
		NewtonMethod<G, ThisType> newtonMethod(this->grid, *this);
		newtonMethod.execute();
		dt = this->localJacobian.getDt();
		double upperMass, oldUpperMass;
		double totalMass = this->injected(upperMass, oldUpperMass);
		std::cout << totalMass << "\t" << upperMass 
			  << "\t" << oldUpperMass << "\t# totalMass, upperMass, oldUpperMass" << std::endl;
		
		*(this->uOldTimeStep) = *(this->u);

		return;
	}

    void globalDefect(FunctionType& defectGlobal)
    {   
      typedef typename G::Traits::template Codim<0>::Entity Entity;
      typedef typename G::ctype DT;
      typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator Iterator;
      enum{dim = G::dimension};
      typedef array<BoundaryConditions::Flags, m> BCBlockType;     
      
      const IS& indexset(this->grid.leafIndexSet());
      (*defectGlobal)=0;
           
      // allocate flag vector to hold flags for essential boundary conditions
      std::vector<BCBlockType> essential(this->vertexmapper.size());
      for (typename std::vector<BCBlockType>::size_type i=0; i<essential.size(); i++)
	essential[i].assign(BoundaryConditions::neumann);

      // iterate through leaf grid 
      Iterator eendit = indexset.template end<0, All_Partition>();
      for (Iterator it = indexset.template begin<0, All_Partition>(); it != eendit; ++it)
	{
	  // get geometry type
	  Dune::GeometryType gt = it->geometry().type();
	  
	  const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
	    sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt, 1);
	  int size = sfs.size();
	  
	  // get entity 
	  const Entity& entity = *it;
	  
	  this->localJacobian.fvGeom.update(entity);
	  
	  this->localJacobian.setLocalSolution(entity);
	  this->localJacobian.template localDefect<LeafTag>(entity,this->localJacobian.u);

	  // begin loop over vertices
	  for(int i=0; i < size; i++)
	    {
	      int globalId = this->vertexmapper.template map<dim>(entity, sfs[i].entity());
	    
	      if (this->localJacobian.bc(i)[0] == BoundaryConditions::neumann)
		(*defectGlobal)[globalId] += this->localJacobian.def[i];
	      else 
		essential[globalId].assign(BoundaryConditions::dirichlet);
	    }
	}

      for (typename std::vector<BCBlockType>::size_type i=0; i<essential.size(); i++)
	if (essential[i][0] == BoundaryConditions::dirichlet)
	  (*defectGlobal)[i] = 0;
    }
  };

}
#endif
