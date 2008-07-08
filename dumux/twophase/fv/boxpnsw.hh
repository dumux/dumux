#ifndef DUNE_BOXPNSW_HH
#define DUNE_BOXPNSW_HH

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
#include "dumux/twophase/twophasemodel.hh"
#include "dumux/twophase/twophaseproblem.hh"
#include "dumux/twophase/fv/boxpnswjacobian.hh"

namespace Dune {

/**
 \brief Two phase model with Pn and Sw as primary unknowns

 This implements a two phase model with Pn and Sw as primary unknowns.
 */
template<class G, class RT> class BoxPnSw :
public LeafP1TwoPhaseModel<G, RT, TwoPhaseProblem<G, RT>,
BoxPnSwJacobian<G, RT> > {

public:
	// define the problem type (also change the template argument above)
	typedef TwoPhaseProblem<G, RT> ProblemType;

	// define the local Jacobian (also change the template argument above)
	typedef BoxPnSwJacobian<G, RT> LocalJacobian;

	typedef LeafP1TwoPhaseModel<G, RT, ProblemType, LocalJacobian>
	LeafP1TwoPhaseModel;

	typedef typename LeafP1TwoPhaseModel::FunctionType FunctionType;

	typedef typename G::Traits::LeafIndexSet IS;

	enum {m = 2};

	typedef BoxPnSw<G, RT> ThisType;
	typedef typename LeafP1TwoPhaseModel::FunctionType::RepresentationType VectorType;
	typedef typename LeafP1TwoPhaseModel::OperatorAssembler::RepresentationType MatrixType;
	typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
#ifdef HAVE_PARDISO
	SeqPardiso<MatrixType,VectorType,VectorType> pardiso;
#endif

	BoxPnSw(const G& g, ProblemType& prob) :
		LeafP1TwoPhaseModel(g, prob) {
	}

	void solve() {

		Operator op(*(this->A)); // make operator out of matrix
		double red=1E-8;

#ifdef HAVE_PARDISO 
		//	SeqPardiso<MatrixType,VectorType,VectorType> ilu0(*(this->A)); 
		pardiso.factorize(*(this->A));
		BiCGSTABSolver<VectorType> solver(op,pardiso,red,100,2); // an inverse operator 
		//	SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
		//LoopSolver<VectorType> solver(op, ilu0, red, 10, 2);
#else
		SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A), 1.0);// a precondtioner

		//SeqIdentity<MatrixType,VectorType,VectorType> ilu0(*(this->A));// a precondtioner
		BiCGSTABSolver<VectorType> solver(op, ilu0, red, 10000, 1); // an inverse operator 
#endif
		InverseOperatorResult r;
		solver.apply(*(this->u), *(this->f), r);

		return;
	}

	void update(double& dt) {
		this->localJacobian.setDt(dt);
		this->localJacobian.setOldSolution(this->uOldTimeStep);
		NewtonMethod<G, ThisType> newtonMethod(this->grid, *this, 1e-6, 1e-4);
		newtonMethod.execute();
		dt = this->localJacobian.getDt();
		double upperMass, oldUpperMass;
		double totalMass = this->injected(upperMass, oldUpperMass);
		std::cout << totalMass << "\t"<< upperMass<< "\t"<< oldUpperMass
		<< "\t# totalMass, upperMass, oldUpperMass"<< std::endl;

		*(this->uOldTimeStep) = *(this->u);

		if (this->problem.exsolution)
			this->problem.updateExSol(dt, *(this->u));

		return;
	}

	//overwrite vtkout for pnSw formulation
	virtual void vtkout(const char* name, int k) {
		int size=this->vertexmapper.size();
		VTKWriter<G> vtkwriter(this->grid);
		char fname[128];
		sprintf(fname, "%s-%05d", name, k);
		double minSat = 1e100;
		double maxSat = -1e100;
		if (this->problem.exsolution){
			this->satEx.resize(size);
			this->satError.resize(size);
		}
		for (int i = 0; i < size; i++) {
			this->pN[i] = (*(this->u))[i][0];
			this->satW[i] = (*(this->u))[i][1];
			this->satN[i] = 1 - this->satW[i];
			double satWI = this->satW[i];
			minSat = std::min(minSat, satWI);
			maxSat = std::max(maxSat, satWI);
			if (this->problem.exsolution){
				this->satEx[i]=this->problem.uExOutVertex(i, 1);
				this->satError[i]=this->problem.uExOutVertex(i, 2);
			}
		}
		vtkwriter.addVertexData(this->pN, "nonwetting phase pressure");
		vtkwriter.addVertexData(this->satW, "wetting phase saturation");
		vtkwriter.addVertexData(this->satN, "nonwetting phase saturation");
		if (this->problem.exsolution){
			vtkwriter.addVertexData(this->satEx, "saturation, exact solution");
			vtkwriter.addVertexData(this->satError, "saturation error");
		}
		vtkwriter.write(fname, VTKOptions::ascii);
		std::cout << "nonwetting phase saturation: min = "<< minSat
		<< ", max = "<< maxSat << std::endl;
		if (minSat< -0.5 || maxSat > 1.5)DUNE_THROW(MathError, "Saturation exceeds range.");
	}

};

}
#endif
