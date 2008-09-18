// $Id$ 

#ifndef DUNE_TWOPHASEHEATMODEL_HH
#define DUNE_TWOPHASEHEATMODEL_HH

#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include "dumux/operators/p1operatorextended.hh"
#include "dumux/nonlinear/nonlinearmodel.hh"
#include "dumux/fvgeometry/fvelementgeometry.hh"

namespace Dune {
template<class G, class RT, class ProblemType, class LocalJacobian,
		class FunctionType, class OperatorAssembler> class TwoPhaseHeatModel :
	public NonlinearModel<G, RT, ProblemType, LocalJacobian, FunctionType, OperatorAssembler> {
public:
	typedef NonlinearModel<G, RT, ProblemType, LocalJacobian,
	FunctionType, OperatorAssembler> ThisNonlinearModel;

	TwoPhaseHeatModel(const G& g, ProblemType& prob) :
		ThisNonlinearModel(g, prob), uOldTimeStep(g) {
	}

	TwoPhaseHeatModel(const G& g, ProblemType& prob, int level) :
		ThisNonlinearModel(g, prob, level), uOldTimeStep(g, level) {
	}

	virtual void initial() = 0;

	virtual void update(double& dt) = 0;

	virtual void solve() = 0;

	FunctionType uOldTimeStep;
};

template<class G, class RT, class ProblemType, class LocalJac, int m=3> class LeafP1TwoPhaseModel :
	public TwoPhaseHeatModel<G, RT, ProblemType, LocalJac,
		LeafP1FunctionExtended<G, RT, m>, LeafP1OperatorAssembler<G, RT, m> > {
public:
	// define the function type:
	typedef LeafP1FunctionExtended<G, RT, m> FunctionType;

	// define the operator assembler type:
	typedef LeafP1OperatorAssembler<G, RT, m> OperatorAssembler;

	typedef TwoPhaseHeatModel<G, RT, ProblemType, LocalJac,
	FunctionType, OperatorAssembler> ThisTwoPhaseHeatModel;

	typedef LeafP1TwoPhaseModel<G, RT, ProblemType, LocalJac, m> ThisType;

	typedef LocalJac LocalJacobian;

	// mapper: one data element per vertex
	template<int dim> struct P1Layout {
		bool contains(Dune::GeometryType gt) {
			return gt.dim() == 0;
		}
	};

	typedef typename G::LeafGridView GV;
    typedef typename GV::IndexSet IS;
	typedef MultipleCodimMultipleGeomTypeMapper<G,IS,P1Layout> VertexMapper;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator
			IntersectionIterator;

	LeafP1TwoPhaseModel(const G& g, ProblemType& prob) :
		ThisTwoPhaseHeatModel(g, prob), problem(prob), grid(g), vertexmapper(g,	g.leafIndexSet()), size((*(this->u)).size())
		{
	}

	virtual void initial() {
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::ctype DT;
		typedef typename GV::template Codim<0>::Iterator Iterator;
		enum {dim = G::dimension};
		enum {dimworld = G::dimensionworld};

		const GV& gridview(this->grid.leafView());

		// iterate through leaf grid an evaluate c0 at cell center
		Iterator eendit = gridview.template end<0>();
		for (Iterator it = gridview.template begin<0>(); it
				!= eendit; ++it) {
			// get geometry type
			Dune::GeometryType gt = it->geometry().type();

			// get entity
			const Entity& entity = *it;

			const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
					&sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt,
							1);
			int size = sfs.size();

			for (int i = 0; i < size; i++) {
				// get cell center in reference element
				const Dune::FieldVector<DT,dim>&local = sfs[i].position();

				// get global coordinate of cell center
				Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

				int globalId = vertexmapper.template map<dim>(entity,
						sfs[i].entity());

				// initialize cell concentration
				(*(this->u))[globalId] = this->problem.initial(
						global, entity, local);
			}
			this->localJacobian.clearVisited();
			this->localJacobian.initiateStaticData(entity);
		}

		// set Dirichlet boundary conditions
		for (Iterator it = gridview.template begin<0>(); it
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
												globalId = vertexmapper.template map<dim>(
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


    virtual double computeFlux ()
     {
		  typedef typename G::Traits::template Codim<0>::Entity Entity;
		  typedef typename G::ctype DT;
		  typedef typename GV::template Codim<0>::Iterator Iterator;
		  enum{dim = G::dimension};
		  enum{dimworld = G::dimensionworld};
	   	  double sign;
	   	  const GV& gridview(grid.leafView());
		  Iterator eendit = gridview.template end<0>();
		  FieldVector<RT,m> flux(0);
		  double Flux(0);

		  for (Iterator it = gridview.template begin<0>(); it != eendit; ++it) // loop over all entities
		  {

			  	// get geometry type
			  	Dune::GeometryType gt = it->geometry().type();

			  	// get entity
			  	const Entity& entity = *it;

				FVElementGeometry<G> fvGeom;
			    fvGeom.update(entity);

				for (int k = 0; k < fvGeom.nEdges; k++)
			     {
				    int i = fvGeom.subContVolFace[k].i;

				    int j = fvGeom.subContVolFace[k].j;

				    int flag_i, flag_j;

				    // 2D case: give y or x value of the line over which flux is to be
				    //			calculated.
				    // up to now only flux calculation to lines or planes (3D) parallel to
				    // x, y and z axis possible

				    // Flux across plane with z = 80 m
				     if(fvGeom.subContVol[i].global[2] < 80.)
			   			  flag_i = 1;
			   		  else flag_i = -1;

			   		  if(fvGeom.subContVol[j].global[2] < 80.)
			   			  flag_j = 1;
			   		  else flag_j = -1;

			   		  if(flag_i == flag_j)
			   		   {
			   			  sign = 0;
			   		   }
			   		  else {
			   			  	if(flag_i > 0)
			   			  		sign = -1;
			   			  	else sign = 1; }

			   			  // get variables

			   		  if(flag_i != flag_j)
			   		  {
						this->localJacobian.setLocalSolution(entity);
						this->localJacobian.computeElementData(entity);
						this->localJacobian.updateVariableData(entity, this->localJacobian.u);


						flux = this->localJacobian.computeA(entity, this->localJacobian.u, k);
						Flux += sign*flux[1];
			   		  }
			     }

		  }
		  return Flux; // Co2 flux
     }


	virtual double totalCO2Mass() {
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::ctype DT;
		typedef typename GV::template Codim<0>::Iterator Iterator;
		enum {dim = G::dimension};
		enum {dimworld = G::dimensionworld};
		enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};	// Phase state
		const GV& gridview(grid.leafView());
		double totalMass = 0;
		double minSat = 1e100;
		double maxSat = -1e100;
		double minP  = 1e100;
		double maxP = -1e100;
		double minTe = 1e100;
		double maxTe = -1e100;
		double minX = 1e100;
		double maxX = -1e100;

		// iterate through leaf grid an evaluate c0 at cell center
		Iterator eendit = gridview.template end<0>();
		for (Iterator it = gridview.template begin<0>(); it
				!= eendit; ++it) {
			// get geometry type
			Dune::GeometryType gt = it->geometry().type();

			// get entity
			const Entity& entity = *it;

			FVElementGeometry<G> fvGeom;
			fvGeom.update(entity);

			const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
					&sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt,
							1);
			int size = sfs.size();

			for (int i = 0; i < size; i++) {
				// get cell center in reference element
				const Dune::FieldVector<DT,dim>&local = sfs[i].position();

				// get global coordinate of cell center
				Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

				int globalId = vertexmapper.template map<dim>(entity,
						sfs[i].entity());

				int state;
				state = this->localJacobian.sNDat[globalId].phaseState;
				RT vol = fvGeom.subContVol[i].volume;
				RT poro = this->problem.soil().porosity(global, entity, local);

				RT rhoN = (*(this->localJacobian.outDensityN))[globalId];
				RT rhoW = (*(this->localJacobian.outDensityW))[globalId];
				RT satN = (*(this->localJacobian.outSaturationN))[globalId];
				RT satW = (*(this->localJacobian.outSaturationW))[globalId];
				RT xAW = (*(this->localJacobian.outMassFracAir))[globalId];
				RT xWN = (*(this->localJacobian.outMassFracWater))[globalId];
				RT xAN = 1 - xWN;
				RT pW = (*(this->u))[globalId][0];
				RT Te = (*(this->u))[globalId][2];
				RT mass = vol * poro * (satN * rhoN * xAN + satW * rhoW * xAW);



				minSat = std::min(minSat, satN);
				maxSat = std::max(maxSat, satN);
				minP = std::min(minP, pW);
				maxP = std::max(maxP, pW);
				minX = std::min(minX, xAW);
				maxX = std::max(maxX, xAW);
				minTe = std::min(minTe, Te);
				maxTe = std::max(maxTe, Te);

				totalMass += mass;
			}

		}

		// print minimum and maximum values
//		std::cout << "nonwetting phase saturation: min = "<< minSat
//				<< ", max = "<< maxSat << std::endl;
//		std::cout << "wetting phase pressure: min = "<< minP
//				<< ", max = "<< maxP << std::endl;
//		std::cout << "mole fraction CO2: min = "<< minX
//				<< ", max = "<< maxX << std::endl;
//		std::cout << "temperature: min = "<< minTe
//				<< ", max = "<< maxTe << std::endl;

		return totalMass;
	}

	virtual void globalDefect(FunctionType& defectGlobal) {
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::ctype DT;
		typedef typename GV::template Codim<0>::Iterator Iterator;
		enum {dim = G::dimension};
		typedef array<BoundaryConditions::Flags, m> BCBlockType;

		const GV& gridview(grid.leafView());
		(*defectGlobal)=0;
		// allocate flag vector to hold flags for essential boundary conditions
		std::vector<BCBlockType> essential(this->vertexmapper.size());
		for (typename std::vector<BCBlockType>::size_type i=0; i
				<essential.size(); i++)
			essential[i].assign(BoundaryConditions::neumann);
		// iterate through leaf grid
		Iterator eendit = gridview.template end<0>();
		for (Iterator it = gridview.template begin<0>(); it
				!= eendit; ++it) {
			// get geometry type
			Dune::GeometryType gt = it->geometry().type();

			// get entity
			const Entity& entity = *it;
			this->localJacobian.fvGeom.update(entity);
			int size = this->localJacobian.fvGeom.nNodes;
			this->localJacobian.setLocalSolution(entity);
			this->localJacobian.computeElementData(entity);
			bool old = true;
			this->localJacobian.updateVariableData(entity, this->localJacobian.uold, old);
			this->localJacobian.updateVariableData(entity, this->localJacobian.u);
			this->localJacobian.template localDefect<LeafTag>(entity, this->localJacobian.u);

			// begin loop over vertices
			for (int i=0; i < size; i++) {
				int globalId = this->vertexmapper.template map<dim>(entity,i);
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


	virtual void vtkout(const char* name, int k) {
	}

protected:
  ProblemType& problem;
  const G& grid;
  VertexMapper vertexmapper;
  int size;
};

}
#endif
