#ifndef DUNE_TWOPHASEMODEL_HH
#define DUNE_TWOPHASEMODEL_HH

#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include <dune/disc/functions/p1function.hh>
#include "dumux/operators/p1operatorextended.hh"
#include "dumux/nonlinear/nonlinearmodel.hh"
#include "dumux/fvgeometry/fvelementgeometry.hh"

namespace Dune {
template<class G, class RT, class ProblemType, class LocalJacobian,
		class FunctionType, class OperatorAssembler> class TwoPhaseModel :
	public NonlinearModel<G, RT, ProblemType, LocalJacobian, FunctionType, OperatorAssembler> {
public:
	typedef NonlinearModel<G, RT, ProblemType, LocalJacobian,
	FunctionType, OperatorAssembler> NonlinearModel;

	TwoPhaseModel(const G& g, ProblemType& prob) :
		NonlinearModel(g, prob), uOldTimeStep(g) {
	}

	TwoPhaseModel(const G& g, ProblemType& prob, int level) :
		NonlinearModel(g, prob, level), uOldTimeStep(g, level) {
	}

	virtual void initial() = 0;

	virtual void update(double& dt) = 0;

	virtual void solve() = 0;

	FunctionType uOldTimeStep;
};

template<class G, class RT, class ProblemType, class LocalJac, int m=2> class LeafP1TwoPhaseModel :
	public TwoPhaseModel<G, RT, ProblemType, LocalJac,
		LeafP1Function<G, RT, m>, LeafP1OperatorAssembler<G, RT, m> > {
public:
	// define the function type:
	typedef LeafP1Function<G, RT, m> FunctionType;

	// define the operator assembler type:
	typedef LeafP1OperatorAssembler<G, RT, m> OperatorAssembler;

	typedef TwoPhaseModel<G, RT, ProblemType, LocalJac,
	FunctionType, OperatorAssembler> TwoPhaseModel;

	typedef LeafP1TwoPhaseModel<G, RT, ProblemType, LocalJac, m> ThisType;

	typedef LocalJac LocalJacobian;

	// mapper: one data element per vertex
	template<int dim> struct P1Layout {
		bool contains(Dune::GeometryType gt) {
			return gt.dim() == 0;
		}
	};

	typedef typename G::Traits::LeafIndexSet IS;
	typedef MultipleCodimMultipleGeomTypeMapper<G,IS,P1Layout> VertexMapper;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator
			IntersectionIterator;

	LeafP1TwoPhaseModel(const G& g, ProblemType& prob) :
		TwoPhaseModel(g, prob), problem(prob), grid(g), vertexmapper(g,
				g.leafIndexSet()), size((*(this->u)).size()), pW(size), pN(size), pC(size),
				satW(size), satN(size), satEx(0), pEx(0), satError(0) {
	}

	virtual void initial() {
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::ctype DT;
		typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator
				Iterator;
		enum {dim = G::dimension};
		enum {dimworld = G::dimensionworld};

		const IS& indexset(grid.leafIndexSet());

		// iterate through leaf grid an evaluate c0 at cell center
		Iterator eendit = indexset.template end<0, All_Partition>();
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

	virtual double injected(double& upperMass, double& oldUpperMass) {
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::ctype DT;
		typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator
				Iterator;
		enum {dim = G::dimension};
		enum {dimworld = G::dimensionworld};

		const IS& indexset(grid.leafIndexSet());
		double totalMass = 0;
		upperMass = 0;
		oldUpperMass = 0;
		// iterate through leaf grid an evaluate c0 at cell center
		Iterator eendit = indexset.template end<0, All_Partition>();
		for (Iterator it = indexset.template begin<0, All_Partition>(); it
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

				double volume = fvGeom.subContVol[i].volume;

				double porosity = this->problem.porosity(global, entity, local);

				double density = this->problem.materialLaw().nonwettingPhase.density();

				double mass = volume*porosity*density*((*(this->u))[globalId][1]);

				totalMass += mass;

				if (global[2] > 80.0) {
					upperMass += mass;
					oldUpperMass += volume*porosity*density*((*(this->uOldTimeStep))[globalId][1]);
				}
			}
		}

		return totalMass;
	}

	virtual void vtkout(const char* name, int k) {
		VTKWriter<G> vtkwriter(this->grid);
		char fname[128];
		sprintf(fname, "%s-%05d", name, k);
		double minSat = 1e100;
		double maxSat = -1e100;
		if (problem.exsolution) {
			satEx.resize(size);
			satError.resize(size);
		}
		for (int i = 0; i < size; i++) {
			pW[i] = (*(this->u))[i][0];
			satN[i] = (*(this->u))[i][1];
			satW[i] = 1 - satN[i];
			double satNI = satN[i];
			minSat = std::min(minSat, satNI);
			maxSat = std::max(maxSat, satNI);
			//pN[i] = (*(this->u))[i][1];
			//pC[i] = pN[i] - pW[i];
			//satW[i] = this->problem.materialLaw().saturationW(pC[i]);
			//satN[i] = 1 - satW[i];
			if (problem.exsolution) {
				satEx[i]=problem.uExOutVertex(i, 1);
				satError[i]=problem.uExOutVertex(i, 2);
			}
		}
		vtkwriter.addVertexData(pW, "wetting phase pressure");
		vtkwriter.addVertexData(satW, "wetting phase saturation");
		vtkwriter.addVertexData(satN, "nonwetting phase saturation");
		if (problem.exsolution) {
			vtkwriter.addVertexData(satEx, "saturation, exact solution");
			vtkwriter.addVertexData(satError, "saturation error");
		}
		vtkwriter.write(fname, VTKOptions::ascii);
		std::cout << "nonwetting phase saturation: min = "<< minSat
				<< ", max = "<< maxSat << std::endl;
		if (minSat< -0.5 || maxSat > 1.5)DUNE_THROW(MathError, "Saturation exceeds range.");
	}

protected:
	ProblemType& problem;
	const G& grid;
	VertexMapper vertexmapper;
	int size;
	BlockVector<FieldVector<RT, 1> > pW;
	BlockVector<FieldVector<RT, 1> > pN;
	BlockVector<FieldVector<RT, 1> > pC;
	BlockVector<FieldVector<RT, 1> > satW;
	BlockVector<FieldVector<RT, 1> > satN;
	BlockVector<FieldVector<RT, 1> > satEx;
	BlockVector<FieldVector<RT, 1> > pEx;
	BlockVector<FieldVector<RT, 1> > satError;
};

}
#endif
