#ifndef DUNE_BOXPWSNTE_HH
#define DUNE_BOXPWSNTE_HH

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
#include "dumux/pardiso/pardiso.hh"
#include "dumux/pardiso/identity.hh"
#include "dumux/nonlinear/newtonmethod.hh"
#include "dumux/2penergy/2penergymodel.hh"
#include "dumux/2penergy/2penergyproblem.hh"
#include "dumux/2penergy/fv/boxpwsntejacobian.hh"
namespace Dune
{

/**
 \brief Two phase model with Pw and Sn as primary unknowns
 
 This implements a two phase model with Pw and Sn as primary unknowns.
 */
  template<class G, class RT>
  class BoxPwSnTe 
  : public LeafP1TwoPhaseModel<G, RT, TwoPhaseHeatProblem<G, RT>, BoxPwSnTeJacobian<G, RT> >
  {
		
		
		 

 

  public:
	// define the problem type (also change the template argument above)
	typedef TwoPhaseHeatProblem<G, RT> ProblemType;
	
	// define the local Jacobian (also change the template argument above)
	typedef BoxPwSnTeJacobian<G, RT> LocalJacobian;
	
	typedef LeafP1TwoPhaseModel<G, RT, ProblemType, LocalJacobian> LeafP1TwoPhaseModel;
	
    typedef typename LeafP1TwoPhaseModel::FunctionType FunctionType;

    typedef typename G::Traits::LeafIndexSet IS;

    enum{m = 3};

	typedef BoxPwSnTe<G, RT> ThisType;
		typedef typename LeafP1TwoPhaseModel::FunctionType::RepresentationType VectorType;
		typedef typename LeafP1TwoPhaseModel::OperatorAssembler::RepresentationType MatrixType;
		typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator; 
#ifdef HAVE_PARDISO
	SeqPardiso<MatrixType,VectorType,VectorType> pardiso; 
#endif
	
	BoxPwSnTe(const G& g, ProblemType& prob) 
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
		NewtonMethod<G, ThisType> newtonMethod(this->grid, *this, 1.e-8, 1.e+5);
		newtonMethod.execute();
		dt = this->localJacobian.getDt();
		double upperMass, oldUpperMass;
		double totalMass = this->injected(upperMass, oldUpperMass);
		std::cout << "total CO2 Mass: "<<totalMass << "\t" << std::endl;
//		double MassFlux = 0;
//		double upperValue[2], lowerValue[2];
//		upperValue[0] = upperValue[1] = 0.;
//		lowerValue[0] = lowerValue[1] = 0.;
//		MassFlux = this->ComputeFlux(upperValue,lowerValue);
//		std::cout << MassFlux <<" "<< upperValue[0] << " "<< upperValue[1] << " "<< 
//		lowerValue[0] <<" "<< lowerValue[1] <<" ";
		*(this->uOldTimeStep) = *(this->u);

		return;
	}

           void globalDefect(FunctionType& defectGlobal) {
                typedef typename G::Traits::template Codim<0>::Entity Entity;
                typedef typename G::ctype DT;
                typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator
                                Iterator;
                enum {dim = G::dimension};
                typedef array<BoundaryConditions::Flags, m> BCBlockType;

                const IS& indexset(this->grid.leafIndexSet());
                (*defectGlobal)=0;

                // allocate flag vector to hold flags for essential boundary conditions
                std::vector<BCBlockType> essential(this->vertexmapper.size());
                for (typename std::vector<BCBlockType>::size_type i=0; i
                                <essential.size(); i++)
                        essential[i].assign(BoundaryConditions::neumann);

                // iterate through leaf grid
                Iterator eendit = indexset.template end<0, All_Partition>();
                for (Iterator it = indexset.template begin<0, All_Partition>(); it
                                != eendit; ++it) {
                        // get geometry type
                        Dune::GeometryType gt = it->geometry().type();

                        const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
                                        &sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt,
                                                        1);
                        int size = sfs.size();

                        // get entity
                        const Entity& entity = *it;

                        this->localJacobian.fvGeom.update(entity);

                        this->localJacobian.setLocalSolution(entity);
                        this->localJacobian.template localDefect<LeafTag>(entity,
                                        this->localJacobian.u);

                        // begin loop over vertices
                        for (int i=0; i < size; i++) {
                                int globalId = this->vertexmapper.template map<dim>(entity,
                                                sfs[i].entity());
                                for (int equationnumber = 0; equationnumber < m; equationnumber++) {
                                        if (this->localJacobian.bc(i)[equationnumber] == BoundaryConditions::neumann)
                                                (*defectGlobal)[globalId][equationnumber]
                                                                += this->localJacobian.def[i][equationnumber];
                                        else
                                                essential[globalId].assign(BoundaryConditions::dirichlet);
                                }
                        }
                }

                for (typename std::vector<BCBlockType>::size_type i=0; i
                                <essential.size(); i++)
                        for (int equationnumber = 0; equationnumber < m; equationnumber++) {
                        if (essential[i][equationnumber] == BoundaryConditions::dirichlet)
                                (*defectGlobal)[i][equationnumber] = 0;
                        }
        }
 
  };

}
#endif
