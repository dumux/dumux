// $Id$

#ifndef DUNE_TWOPHASEHEATMODEL_HH
#define DUNE_TWOPHASEHEATMODEL_HH

#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include "dumux/operators/p1operatorextended.hh"
#include "dumux/nonlinear/nonlinearmodel.hh"
#include "dumux/fvgeometry/fvelementgeometry.hh"
#include "dumux/io/exporttodgf.hh"
#include <boost/format.hpp>

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

	virtual void restart() {}

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
		ThisTwoPhaseHeatModel(g, prob), problem(prob), _grid(g), vertexmapper(g,	g.leafIndexSet()), size((*(this->u)).size())
		{
	}

        virtual void update(double &dt) {
            DUNE_THROW(NotImplemented, "This method is obsolete. Use updateModel()!");
        }

	virtual void initial() {}

	virtual void restart(int restartNum=0) {}

	virtual void globalDefect(FunctionType& defectGlobal) {
		typedef typename G::Traits::template Codim<0>::Entity Entity;
		typedef typename G::ctype DT;
		typedef typename GV::template Codim<0>::Iterator Iterator;
		enum {dim = G::dimension};
		typedef array<BoundaryConditions::Flags, m> BCBlockType;

		const GV& gridview(this->_grid.leafView());
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

	void writerestartfile(int restartNum=0)
	{
		enum {dim = G::dimension};
		typedef typename GV::template Codim<dim>::Iterator Iterator;

//		exportToDGF(_grid.leafView(), *(this->u), m, "primvar", false);

		const int size = vertexmapper.size();
		BlockVector<FieldVector<double, m> > data(size);
		data=0;

		Iterator endIt = _grid.leafView().template end<dim>();
		for (Iterator it = _grid.leafView().template begin<dim>(); it != endIt;	++it)
		{
			int index = vertexmapper.map(*it);
			for (int i = 0; i < m;i++)
			{
				data[index][i]=(*(this->u))[index][i];
			}
		}

		restartFileName = (boost::format("data-%05d")
                           %restartNum).str();
		exportToDGF(_grid.leafView(), data, (m), restartFileName, false);
	}

	virtual void vtkout(const char* name, int k) {}

    const G &grid() const
        { return _grid; }



protected:
  ProblemType& problem;
  const G& _grid;
  VertexMapper vertexmapper;
  int size;
  std::string restartFileName;

};

}
#endif
