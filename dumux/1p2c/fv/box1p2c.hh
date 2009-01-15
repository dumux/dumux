// $Id$

#ifndef DUNE_BOX1P2C_HH
#define DUNE_BOX1P2C_HH

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
#include "dumux/operators/p1operatorextended.hh"
#include <dune/disc/operators/boundaryconditions.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/paamg/amg.hh>
#include "dumux/pardiso/pardiso.hh"
#include "dumux/nonlinear/newtonmethod.hh"
#include "dumux/1p2c/onephasemodel.hh"
#include "dumux/1p2c/1p2cproblem.hh"
#include "dumux/1p2c/fv/box1p2cjacobian.hh"

#include "dumux/nonlinear/new_newtonmethod.hh"
#include "dumux/nonlinear/new_newtoncontroller.hh"

namespace Dune
{
/**
 \brief Isothermal one phase two component model with P and X as primary unknowns

 This implements an isothermal one phase two component model with P and X as primary unknowns
 */  template<class G, class RT, class VtkMultiWriter>
  class Box1P2C
  : public LeafP1OnePhaseModel<G, RT, OnePTwoCProblem<G, RT>, Box1P2CJacobian<G, RT> >
  {
  public:
      // define the problem type (also change the template argument above)
	typedef OnePTwoCProblem<G, RT> ProblemType;

	// define the local Jacobian (also change the template argument above)
	typedef Box1P2CJacobian<G, RT> LocalJacobian;
	typedef LeafP1OnePhaseModel<G, RT, ProblemType, LocalJacobian> ThisLeafP1OnePhaseModel;
	typedef Box1P2C<G, RT, VtkMultiWriter> ThisType;

	typedef typename ThisLeafP1OnePhaseModel::FunctionType FunctionType;
	typedef typename ThisLeafP1OnePhaseModel::FunctionType::RepresentationType VectorType;
	typedef typename ThisLeafP1OnePhaseModel::OperatorAssembler::RepresentationType MatrixType;
    typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
    
	typedef typename G::LeafGridView GV;
	

	enum{m = 2};

      //////////////////////
      // Stuff required for the new newton method

        //! The traits class for the new newton method.
        struct NewtonTraits {
            typedef RT                                                  Scalar;
            typedef typename ThisLeafP1OnePhaseModel::FunctionType      Function;
            typedef typename ThisType::LocalJacobian                    LocalJacobian;
            typedef typename ThisLeafP1OnePhaseModel::OperatorAssembler JacobianAssembler;
        };

        // HACK: traits for the domain of the problem. this is incomplete...
        struct DomainTraits {
            typedef RT   Scalar;
        };

        typedef NewNewtonMethod<ThisType> NewtonMethod;
        typedef Dune::NewtonController<NewtonMethod> NewtonController;

        typedef typename NewtonTraits::Function Function;
        Function &currentSolution()
          { return this->u; };

        LocalJacobian &getLocalJacobian()
          { return this->localJacobian; }

        typedef typename NewtonTraits::JacobianAssembler JacobianAssembler;
        JacobianAssembler &jacobianAssembler()
          { return this->A; }
        // End of stuff for new newton method
        //////////////////////

	Box1P2C(const G& g, ProblemType& prob)
	: ThisLeafP1OnePhaseModel(g, prob)// (this->size) vectors
	{ }

	void initial()
	{
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::ctype DT;
		typedef typename GV::template Codim<0>::Iterator Iterator;
		typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

		enum {dim = G::dimension};
		enum {dimworld = G::dimensionworld};

		const GV& gridview(this->grid().leafView());


		this->localJacobian().clearVisited();
	

		// iterate through leaf grid an evaluate c0 at cell center
		Iterator eendit = gridview.template end<0>();
		for (Iterator it = gridview.template begin<0>(); it
				!= eendit; ++it)
		{
			// get geometry type
			Dune::GeometryType gt = it->geometry().type();

			// get entity
			const Entity& entity = *it;

			this->localJacobian().fvGeom.update(entity);

			const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
					&sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt,1);

			int size = sfs.size();

			for (int i = 0; i < size; i++)
			{
				// get cell center in reference element
				const Dune::FieldVector<DT,dim>&local = sfs[i].position();

				// get global coordinate of cell center
				Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

				int globalId = this->vertexmapper.template map<dim>(entity,
						sfs[i].entity());

				// initialize cell concentration
				(*(this->u))[globalId] = this->problem.initial(
						global, entity, local);
			}
			this->localJacobian().clearVisited();
			this->localJacobian().setLocalSolution(entity);
			this->localJacobian().updateStaticData(entity, this->localJacobian().u);
		}
		
		// set Dirichlet boundary conditions
		for (Iterator it = gridview.template begin<0>();
				it != eendit; ++it)
		{
			// get geometry type
			Dune::GeometryType gt = it->geometry().type();

			// get entity
			const Entity& entity = *it;

			const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
					&sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt,1);
			int size = sfs.size();

			// set type of boundary conditions
			this->localJacobian().template assembleBC<LeafTag>(entity);

			IntersectionIterator
					endit = IntersectionIteratorGetter<G, LeafTag>::end(entity);

			for (IntersectionIterator is = IntersectionIteratorGetter<G,
					LeafTag>::begin(entity); is!=endit; ++is)
				if (is->boundary())
				{
					for (int i = 0; i < size; i++)
						// handle subentities of this face
						for (int j = 0; j < ReferenceElements<DT,dim>::general(gt).size(is->numberInSelf(), 1, sfs[i].codim()); j++)
							if (sfs[i].entity()
									== ReferenceElements<DT,dim>::general(gt).subEntity(is->numberInSelf(), 1,
											j, sfs[i].codim()))
							{
								for (int equationNumber = 0; equationNumber<m; equationNumber++)
								{
									if (this->localJacobian().bc(i)[equationNumber]
											== BoundaryConditions::dirichlet)
									{
										// get cell center in reference element
										Dune::FieldVector<DT,dim>
												local = sfs[i].position();

										// get global coordinate of cell center
										Dune::FieldVector<DT,dimworld>
												global = it->geometry().global(local);

										int globalId = this->vertexmapper.template map<dim>(entity, sfs[i].entity());
										FieldVector<int,m> dirichletIndex;
										FieldVector<BoundaryConditions::Flags, m>
												bctype = this->problem.bctype(global, entity, is, local);
												this->problem.dirichletIndex(global, entity, is,
														local, dirichletIndex);

										if (bctype[equationNumber] == BoundaryConditions::dirichlet)
										{
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
		this->localJacobian().setLocalSolution(entity);
		for (int i = 0; i < size; i++)
			this->localJacobian().updateVariableData(entity, this->localJacobian().u, i, false);

		}

		*(this->uOldTimeStep) = *(this->u);

		return;
	}

	void updateModel (double& dt, double& nextDt)
	{
		
		this->localJacobian.setDt(dt);
		this->localJacobian.setOldSolution(this->uOldTimeStep);

		// execute newton method
		typedef typename GV::template Codim<0>::Iterator Iterator;
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::ctype DT;
		enum {dim = G::dimension};

		bool newtonLoop = false;
		while(!newtonLoop)
		{
			    nextDt = this->localJacobian.getDt();
                NewtonMethod newton(*this); // *this means object itself (box1p2c)
                NewtonController newtonCtl(1e-7, 6);
                newtonLoop = newton.execute(*this, newtonCtl);
                nextDt = newtonCtl.suggestTimeStepSize(nextDt);
                this->localJacobian.setDt(nextDt);
                if(!newtonLoop){
                	*this->u = *this->uOldTimeStep;
                	this->localJacobian.resetPhaseState();
                }
                std::cout<<"timeStep resized to: "<<nextDt<<std::endl;
		}

		
		this->localJacobian.clearVisited();

		*(this->uOldTimeStep) = *(this->u);

		return;
	}

//	void restart()
//	{
//		typedef typename G::Traits::template Codim<0>::Entity Entity;
//		typedef typename G::ctype DT;
//		typedef typename GV::template Codim<0>::Iterator Iterator;
//		typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;
//
//		enum {dim = G::dimension};
//		enum {dimworld = G::dimensionworld};
//
//		const GV& gridview(this->grid().leafView());
//
//		this->localJacobian.outPressureN = vtkMultiWriter->template createField<RT, 1>(this->size);
//		this->localJacobian.outCapillaryP = vtkMultiWriter->template createField<RT, 1>(this->size);
//		this->localJacobian.outSaturationW = vtkMultiWriter->template createField<RT, 1>(this->size);
//		this->localJacobian.outSaturationN = vtkMultiWriter->template createField<RT, 1>(this->size);
//		this->localJacobian.outMassFracAir = vtkMultiWriter->template createField<RT, 1>(this->size);
//		this->localJacobian.outMassFracWater = vtkMultiWriter->template createField<RT, 1>(this->size);
//		this->localJacobian.outDensityW = vtkMultiWriter->template createField<RT, 1>(this->size);
//		this->localJacobian.outDensityN = vtkMultiWriter->template createField<RT, 1>(this->size);
//		this->localJacobian.outMobilityW = vtkMultiWriter->template createField<RT, 1>(this->size);
//		this->localJacobian.outMobilityN = vtkMultiWriter->template createField<RT, 1>(this->size);
//		this->localJacobian.outPhaseState = vtkMultiWriter->template createField<RT, 1>(this->size);
//
//		int size = this->vertexmapper.size();
//		Dune::BlockVector<FieldVector<double, m+1> > data(size);
//		data=0;
//
////		importFromDGF<GV>(data, "data", false);
//
//		for (int i=0;i<size;i++)
//		{
//			for (int j=0;j<m;j++)
//			{
//				(*(this->u))[i][j]=data[i][j];
//			}
//		}
//
//		// iterate through leaf grid an evaluate c0 at cell center
//		Iterator eendit = gridview.template end<0>();
//		for (Iterator it = gridview.template begin<0>();
//			it != eendit; ++it)
//		{
//			// get geometry type
//			Dune::GeometryType gt = it->geometry().type();
//
//			// get entity
//			const Entity& entity = *it;
//
//			this->localJacobian.fvGeom.update(entity);
//
//			const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
//					&sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt,1);
//			int size = sfs.size();
//
//			for (int i = 0; i < size; i++)
//			{
//				// get cell center in reference element
//				const Dune::FieldVector<DT,dim>&local = sfs[i].position();
//
//				// get global coordinate of cell center
//				Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);
//
//				int globalId = this->vertexmapper.template map<dim>(entity,
//						sfs[i].entity());
//
//				// initialize cell concentration
//
//				// initialize variable phaseState
//				this->localJacobian.sNDat[globalId].phaseState = data[globalId][m];
//				// initialize variable oldPhaseState
//				this->localJacobian.sNDat[globalId].oldPhaseState = data[globalId][m];
//
//			}
//			this->localJacobian.clearVisited();
//			this->localJacobian.setLocalSolution(entity);
//			this->localJacobian.updateStaticData(entity, this->localJacobian.u);
//		}
//	}

	virtual void globalDefect(FunctionType& defectGlobal)
	{
		ThisLeafP1OnePhaseModel::globalDefect(defectGlobal);
	}

	void solve()
		{
		Operator op(*(this->A));  // make operator out of matrix
		double red=1E-11;

#ifdef HAVE_PARDISO
//	SeqPardiso<MatrixType,VectorType,VectorType> ilu0(*(this->A));
		pardiso.factorize(*(this->A));
		BiCGSTABSolver<VectorType> solver(op,pardiso,red,100,2);         // an inverse operator
	//	SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
		//LoopSolver<VectorType> solver(op, ilu0, red, 10, 2);
#else
		SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner

		//SeqIdentity<MatrixType,VectorType,VectorType> ilu0(*(this->A));// a precondtioner
		BiCGSTABSolver<VectorType> solver(op,ilu0,red,1000,1);         // an inverse operator
		//iteration numbers were 1e4 before
#endif
		InverseOperatorResult r;
		solver.apply(*(this->u), *(this->f), r);

		return;
	}

	void updateState()
	{
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::ctype DT;
		typedef typename GV::template Codim<0>::Iterator Iterator;
		typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

		enum {dim = G::dimension};
		enum {dimworld = G::dimensionworld};

		const GV& gridview(this->grid_.leafView());
		// iterate through leaf grid and evaluate c0 at cell center
		Iterator eendit = gridview.template end<0>();
		for (Iterator it = gridview.template begin<0>();
			it != eendit; ++it) {

			const Entity& entity = *it;
			this->localJacobian.fvGeom.update(entity);
			this->localJacobian.setLocalSolution(entity);
			this->localJacobian.computeElementData(entity);
			this->localJacobian.updateStaticData(entity, this->localJacobian.u);
		}

		return;
	}



	template<class MultiWriter>
	void addvtkfields (MultiWriter& writer)
	{

	}

	void vtkout (const char* name, int k)
	{
		VTKWriter<typename G::LeafGridView> vtkwriter(this->grid_.leafView());
		char fname[128];
		sprintf(fname,"%s-%05d",name,k);
		BlockVector<FieldVector<RT, 1> > p(this->size);
		BlockVector<FieldVector<RT, 1> > x(this->size);
	  	for (int i = 0; i < this->size; i++) {
			p[i] = (*(this->u))[i][0];  
			x[i] = (*(this->u))[i][1];
						//const FieldVector<RT, 4> parameters(this->problem.materialLawParameters
			//	 		 (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal));
			//			this->pC[i] = this->problem.materialLaw().pC(this->satW[i], parameters);
			}
		
		vtkwriter.addVertexData(p, "pressure");
		vtkwriter.addVertexData(x, "mole fraction");
		vtkwriter.write(fname, VTKOptions::ascii);
	}


	void setVtkMultiWriter(VtkMultiWriter *writer)
	{ vtkMultiWriter = writer; }

  protected:
    VtkMultiWriter *vtkMultiWriter;
    
  };
}

#endif
