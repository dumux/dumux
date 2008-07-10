#ifndef DUNE_BOX2P2CNI_HH
#define DUNE_BOX2P2CNI_HH

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
#include "dumux/nonlinear/newtonmethod.hh"
#include "dumux/2p2cni/2p2cnimodel.hh"
#include "dumux/2p2cni/2p2cniproblem.hh"
#include "dumux/2p2cni/fv/box2p2cnijacobian.hh"

namespace Dune
{
/**
 \brief Non-isothermal two phase two component model with Pw, Sn/X and Temp as primary unknowns
 
 This implements a non-isothermal two phase two component model with Pw, Sn/X and Temp as primary unknowns
 */  
 template<class G, class RT, class VtkMultiWriter>
  class Box2P2CNI 
  : public LeafP1TwoPhaseModel<G, RT, TwoPTwoCNIProblem<G, RT>, Box2P2CNIJacobian<G, RT> >
  {
  public:
	// define the problem type (also change the template argument above)
	typedef TwoPTwoCNIProblem<G, RT> ProblemType;

	// define the local Jacobian (also change the template argument above)
	typedef Box2P2CNIJacobian<G, RT> LocalJacobian;
	typedef LeafP1TwoPhaseModel<G, RT, ProblemType, LocalJacobian> LeafP1TwoPhaseModel;
	typedef Box2P2CNI<G, RT, VtkMultiWriter> ThisType;

	typedef typename LeafP1TwoPhaseModel::FunctionType FunctionType;

   typedef typename G::Traits::LeafIndexSet IS;

    enum{m = 3};

		typedef typename LeafP1TwoPhaseModel::FunctionType::RepresentationType VectorType;
		typedef typename LeafP1TwoPhaseModel::OperatorAssembler::RepresentationType MatrixType;
		typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator; 
#ifdef HAVE_PARDISO
	SeqPardiso<MatrixType,VectorType,VectorType> pardiso;
#endif

	
	Box2P2CNI(const G& g, ProblemType& prob) 
	: LeafP1TwoPhaseModel(g, prob)// (this->size) vectors
	{ }

	void initial() {
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::ctype DT;
		typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator
				Iterator;
		typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

		enum {dim = G::dimension};
		enum {dimworld = G::dimensionworld};

//		this->localJacobian.hackySaturationN = hackyVtkWriter->template createField<RT, 1>(this->size);
//		this->localJacobian.hackyMassFracAir = hackyVtkWriter->template createField<RT, 1>(this->size);
//		this->localJacobian.hackyMassFracWater = hackyVtkWriter->template createField<RT, 1>(this->size);
		
		const IS& indexset(this->grid.leafIndexSet());

		// iterate through leaf grid an evaluate c0 at cell center
		Iterator eendit = indexset.template end<0, All_Partition>();
		for (Iterator it = indexset.template begin<0, All_Partition>(); it
				!= eendit; ++it) {
			// get geometry type
			Dune::GeometryType gt = it->geometry().type();

			// get entity 
			const Entity& entity = *it;

			const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
					&sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt,1);
			int size = sfs.size();

			for (int i = 0; i < size; i++) {
				// get cell center in reference element
				const Dune::FieldVector<DT,dim>&local = sfs[i].position();

				// get global coordinate of cell center
				Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

				int globalId = this->vertexmapper.template map<dim>(entity,
						sfs[i].entity());

				// initialize cell concentration
				(*(this->u))[globalId] = this->problem.initial(
						global, entity, local);
				
				// initialize phase state
				this->localJacobian.sNDat[globalId].phaseState = 
					this->problem.initialPhaseState(global, entity, local);
			}
		}

		// set Dirichlet boundary conditions
		for (Iterator it = indexset.template begin<0, All_Partition>(); it
				!= eendit; ++it) {
			// get geometry type
			Dune::GeometryType gt = it->geometry().type();

			// get entity 
			const Entity& entity = *it;

			const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
					&sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt,
							1);
			int size = sfs.size();

			// set type of boundary conditions 
			this->localJacobian.template assembleBC<LeafTag>(entity);

			IntersectionIterator
					endit = IntersectionIteratorGetter<G, LeafTag>::end(entity);
			for (IntersectionIterator is = IntersectionIteratorGetter<G,
					LeafTag>::begin(entity); is!=endit; ++is)
				if (is->boundary()) {
					for (int i = 0; i < size; i++)
						// handle subentities of this face
						for (int j = 0; j < ReferenceElements<DT,dim>::general(gt).size(is->numberInSelf(), 1, sfs[i].codim()); j++)
							if (sfs[i].entity()
									== ReferenceElements<DT,dim>::general(gt).subEntity(is->numberInSelf(), 1,
											j, sfs[i].codim())) {
								for (int equationNumber = 0; equationNumber<m; equationNumber++) {
									if (this->localJacobian.bc(i)[equationNumber]
											== BoundaryConditions::dirichlet) {
										// get cell center in reference element
										Dune::FieldVector<DT,dim>
												local = sfs[i].position();

										// get global coordinate of cell center
										Dune::FieldVector<DT,dimworld>
												global = it->geometry().global(local);

										int
												globalId = this->vertexmapper.template map<dim>(
														entity, sfs[i].entity());
										FieldVector<int,m> dirichletIndex;
										FieldVector<BoundaryConditions::Flags, m>
												bctype = this->problem.bctype(
														global, entity, is,
														local);
												this->problem.dirichletIndex(global, entity, is,
														local, dirichletIndex);	

										if (bctype[equationNumber]
												== BoundaryConditions::dirichlet) {
											FieldVector<RT,m>
													ghelp = this->problem.g(
															global, entity, is,
															local);
											(*(this->u))[globalId][dirichletIndex[equationNumber]]
													= ghelp[dirichletIndex[equationNumber]];
										}
									}
								}
							}
				}
		}

		*(this->uOldTimeStep) = *(this->u);
		
		return;
	}

	virtual void globalDefect(FunctionType& defectGlobal) 
	{
		LeafP1TwoPhaseModel::globalDefect(defectGlobal);
	}
	
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
//		this->localJacobian.hackySaturationN = hackyVtkWriter->template createField<RT, 1>(this->size);
//		this->localJacobian.hackyMassFracAir = hackyVtkWriter->template createField<RT, 1>(this->size);
//		this->localJacobian.hackyMassFracWater = hackyVtkWriter->template createField<RT, 1>(this->size);

		this->localJacobian.setDt(dt);
		this->localJacobian.setOldSolution(this->uOldTimeStep);

		///////////////////////////////////
		// define solver tolerances here
		///////////////////////////////////
////////////////
//		typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator
//		Iterator;
//		typedef typename G::Traits::template Codim<0>::Entity Entity;
//		typedef typename G::ctype DT;
//		enum {dim = G::dimension};
//
//		const IS& indexset(this->grid.leafIndexSet());
//
//		Iterator eendit = indexset.template end<0, All_Partition>();
//		for (Iterator it = indexset.template begin<0, All_Partition>(); 
//				it != eendit; ++it) 
//		{
//			const Entity& e = *it;
//			Dune::GeometryType gt = it->geometry().type();
//
//			const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
//			&sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt,	1);
//
//			int size = sfs.size();
//			this->localJacobian.updateVariableData(e, this->localJacobian.u);
//			for (int i=0; i < size; i++)
//			{	
//				int globalId = this->vertexmapper.template map<dim>(e,
//						sfs[i].entity());
//				this->localJacobian.primaryVarSwitch(e, globalId, this->localJacobian.u, i);
//			}
//		}
/////////////////////		
		RT absTol = 2e+7;
		RT relTol = 1e-5;
		NewtonMethod<G, ThisType> newtonMethod(this->grid, *this, relTol, absTol);
		newtonMethod.execute();
		dt = this->localJacobian.getDt();
//		double upperMass, oldUpperMass;
//		double totalMass = this->injected(upperMass, oldUpperMass);
//		std::cout << totalMass << "\t" << upperMass 
//			  << "\t" << oldUpperMass << "\t# totalMass, upperMass, oldUpperMass" << std::endl;
		
		*(this->uOldTimeStep) = *(this->u);

		return;
	}

	const G& getGrid() const
	{
		return this->grid;
	}
	
	template<class MultiWriter>
	void addvtkfields (MultiWriter& writer) 
	{
//		BlockVector<FieldVector<RT, 1> > &xWN = *writer.template createField<RT, 1>(this->size);
//		BlockVector<FieldVector<RT, 1> > &xAW = *writer.template createField<RT, 1>(this->size);
//		BlockVector<FieldVector<RT, 1> > &satW = *writer.template createField<RT, 1>(this->size);

		for (int i = 0; i < this->size; i++) {
//			RT pW = (*(this->u))[i][0];
//			RT satN = (*(this->u))[i][1];
//			satW[i] = 1 - satN;
//			xWN[i] = this->problem.multicomp().xWN(pW, 283.15); //Achtung!! pW instead of pN!!!
//			xAW[i] = this->problem.multicomp().xAW(pW, 283.15); //Achtung!! pW instead of pN!!!
		}

//		writer.addScalarVertexFunction("nonwetting phase saturation", 
//										this->u, 
//										1);
		writer.addScalarVertexFunction("wetting phase pressure", 
										this->u, 
										0);
//		writer.addVertexData(&satW,"wetting phase saturation");
		writer.addVertexData(this->localJacobian.hackySaturationN,"nonwetting phase saturation");
		writer.addVertexData(this->localJacobian.hackyMassFracAir,"air in water");
		writer.addVertexData(this->localJacobian.hackyMassFracWater,"water in Gasphase");
//		writer.addVertexData(&xWN, "water in air");
//		writer.addVertexData(&xAW, "dissolved air");
	}

	
//	void vtkout (const char* name, int k) 
//	{
//		
//	}
//
//	void vtkout (const char* name, int k) 
//	{
//		VTKWriter<G> vtkwriter(this->grid);
//		char fname[128];	
//		sprintf(fname,"%s-%05d",name,k);
//	  BlockVector<FieldVector<RT, 1> > xWN(this->size);
//	  BlockVector<FieldVector<RT, 1> > xAW(this->size);
//		for (int i = 0; i < this->size; i++) {
//			this->pW[i] = (*(this->u))[i][0];
//			this->satN[i] = (*(this->u))[i][1];
//			this->satW[i] = 1 - this->satN[i];
//			//const FieldVector<RT, 4> parameters(this->problem.materialLawParameters
//			//	 		 (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal));
//			//			this->pC[i] = this->problem.materialLaw().pC(this->satW[i], parameters);			
//			xWN[i] = this->problem.multicomp().xWN(this->pW[i], 283.15); //Achtung!! pW instead of pN!!!
//			xAW[i] = this->problem.multicomp().xAW(this->pW[i], 283.15); //Achtung!! pW instead of pN!!!
//		}
//		vtkwriter.addVertexData(this->pW,"wetting phase pressure");
//		vtkwriter.addVertexData(this->satW,"wetting phase saturation");
//		vtkwriter.addVertexData(this->satN,"nonwetting phase saturation");
//		vtkwriter.addVertexData(xWN, "water in air");
//		vtkwriter.addVertexData(xAW, "dissolved air");
//		vtkwriter.write(fname, VTKOptions::ascii);		
//	}
//
	
	void setHackyVtkMultiWriter(VtkMultiWriter *writer)
	{ hackyVtkWriter = writer; }
  protected:
    VtkMultiWriter *hackyVtkWriter;
   
  };
}
#endif
