#ifndef DUNE_BOX2P2C_HH
#define DUNE_BOX2P2C_HH

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
#include "dumux/operators/p1operatorextended.hh"
#include <dune/disc/operators/boundaryconditions.hh>
#include <dune/disc/groundwater/groundwater.hh>
#include <dune/disc/groundwater/p1groundwater.hh>
#include <dune/disc/groundwater/p1groundwaterestimator.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/paamg/amg.hh>
#include "dumux/nonlinear/newtonmethod.hh"
#include "dumux/twophase/twophasemodel.hh"
#include "dumux/2p2c/2p2cproblem.hh"
#include "dumux/2p2c/fv/box2p2cjacobian.hh"

namespace Dune
{
  template<class G, class RT>
  class Box2P2C 
  : public LeafP1TwoPhaseModel<G, RT, TwoPTwoCProblem<G, RT>, Box2P2CJacobian<G, RT> >
  {
  public:
	// define the problem type (also change the template argument above)
	typedef TwoPTwoCProblem<G, RT> ProblemType;

	// define the local Jacobian (also change the template argument above)
	typedef Box2P2CJacobian<G, RT> LocalJacobian;
	typedef LeafP1TwoPhaseModel<G, RT, ProblemType, LocalJacobian> LeafP1TwoPhaseModel;
	typedef Box2P2C<G, RT> ThisType;

	typedef typename LeafP1TwoPhaseModel::FunctionType FunctionType;

   typedef typename G::Traits::LeafIndexSet IS;

    enum{m = 2};

		typedef typename LeafP1TwoPhaseModel::FunctionType::RepresentationType VectorType;
		typedef typename LeafP1TwoPhaseModel::OperatorAssembler::RepresentationType MatrixType;
		typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator; 

	
	Box2P2C(const G& g, ProblemType& prob) 
	: LeafP1TwoPhaseModel(g, prob), Xwn(this->size), Xaw(this->size) // (this->size) vectors
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
     
     const IS& indexset(this->grid.leafIndexSet());
     (*defectGlobal)=0;
          
     // iterate through leaf grid 
     Iterator eendit = indexset.template end<0, All_Partition>();
     for (Iterator it = indexset.template begin<0, All_Partition>(); it != eendit; ++it)
	{
	  // get geometry type
	  Dune::GeometryType gt = it->geometry().type();
	  
	  const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
	    sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt, 1);
	  int size = sfs.size();
	  Dune::FieldVector<RT,2> defhelp[size];
	  
	  // get entity 
	  const Entity& entity = *it;
	  
	  this->localJacobian.fvGeom.update(entity);

	  this->localJacobian.getLocalDefect(entity,defhelp);
	  //std::cout<<" defhelp: "<<*defhelp<<std::endl;
	  // begin loop over vertices
	  for(int i=0; i < size; i++)
	    {
	      int globalId = this->vertexmapper.template map<dim>(entity, sfs[i].entity());
	    
	      if (this->localJacobian.bc(i)[0] == BoundaryConditions::neumann)
		(*defectGlobal)[globalId] += defhelp[i];
	      else 
		(*defectGlobal)[globalId] = 0;
	    }
	}
   }
 

	
	void vtkout (const char* name, int k) 
	{
		VTKWriter<G> vtkwriter(this->grid);
		char fname[128];	
		sprintf(fname,"%s-%05d",name,k);
		for (int i = 0; i < this->size; i++) {
			this->pW[i] = (*(this->u))[i][0];
			this->satN[i] = (*(this->u))[i][1];
			this->satW[i] = 1 - this->satN[i];
			//const FieldVector<RT, 4> parameters(this->problem.materialLawParameters
			//	 		 (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal));
			//			parameters = problem.materialLawParameters
			//			 		 (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
			//			this->pC[i] = this->problem.materialLaw().pC(this->satW[i], parameters);			
			Xwn[i] = this->problem.solu().Xwn(this->pW[i], 283.15); //Achtung!! pW instead of pN!!!
			Xaw[i] = this->problem.solu().Xaw(this->pW[i], 283.15); //Achtung!! pW instead of pN!!!
		}
		vtkwriter.addVertexData(this->pW,"wetting phase pressure");
		vtkwriter.addVertexData(this->satW,"wetting phase saturation");
		vtkwriter.addVertexData(this->satN,"nonwetting phase saturation");
		vtkwriter.addVertexData(Xwn, "water in air");
		vtkwriter.addVertexData(Xaw, "dissolved air");
		vtkwriter.write(fname, VTKOptions::ascii);		
	}

  protected:
	  BlockVector<FieldVector<RT, 1> > Xwn;
	  BlockVector<FieldVector<RT, 1> > Xaw;
  };

}
#endif
