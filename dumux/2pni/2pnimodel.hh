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

/** \todo Please doc me! */

template<class Grid, class Scalar, class ProblemType, class LocalJacobian,
        class FunctionType, class OperatorAssembler> class TwoPhaseHeatModel :
    public NonlinearModel<Grid, Scalar, ProblemType, LocalJacobian, FunctionType, OperatorAssembler> {
public:
    typedef NonlinearModel<Grid, Scalar, ProblemType, LocalJacobian,
    FunctionType, OperatorAssembler> ThisNonlinearModel;

    TwoPhaseHeatModel(const Grid& grid, ProblemType& prob) :
        ThisNonlinearModel(grid, prob), uOldTimeStep(grid) {
    }

    TwoPhaseHeatModel(const Grid& grid, ProblemType& prob, int level) :
        ThisNonlinearModel(grid, prob, level), uOldTimeStep(grid, level) {
    }

    virtual void initial() = 0;

    virtual void restart() {}

    virtual void update(double& dt) = 0;

    virtual void solve() = 0;

    FunctionType uOldTimeStep;
};

/** \todo Please doc me! */

template<class Grid, class Scalar, class ProblemType, class LocalJac, int numEq=3> class LeafP1TwoPhaseModel :
    public TwoPhaseHeatModel<Grid, Scalar, ProblemType, LocalJac,
        LeafP1Function<Grid, Scalar, numEq>, LeafP1OperatorAssembler<Grid, Scalar, numEq> > {
public:
    // define the function type:
    typedef LeafP1Function<Grid, Scalar, numEq> FunctionType;

    // define the operator assembler type:
    typedef LeafP1OperatorAssembler<Grid, Scalar, numEq> OperatorAssembler;

    typedef TwoPhaseHeatModel<Grid, Scalar, ProblemType, LocalJac,
    FunctionType, OperatorAssembler> ThisTwoPhaseHeatModel;

    typedef LeafP1TwoPhaseModel<Grid, Scalar, ProblemType, LocalJac, numEq> ThisType;

    typedef LocalJac LocalJacobian;

    // mapper: one data element per vertex
    template<int dim> struct P1Layout {
        bool contains(Dune::GeometryType gt) {
            return gt.dim() == 0;
        }
    };

    typedef typename Grid::LeafGridView GridView;
    typedef typename GridView::IndexSet IS;
    typedef MultipleCodimMultipleGeomTypeMapper<Grid,IS,P1Layout> VertexMapper;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator
            IntersectionIterator;

    LeafP1TwoPhaseModel(const Grid& grid, ProblemType& prob) :
        ThisTwoPhaseHeatModel(grid, prob), problem(prob), _grid(grid), vertexmapper(grid,    grid.leafIndexSet()), size((*(this->u)).size())
        {
    }

        virtual void update(double &dt) {
            DUNE_THROW(NotImplemented, "This method is obsolete. Use updateModel()!");
        }

    virtual void initial() {}

    virtual void restart(int restartNum=0) {}

    virtual void globalDefect(FunctionType& defectGlobal) {
        typedef typename Grid::Traits::template Codim<0>::Entity Element;
        typedef typename Grid::ctype CoordScalar;
        typedef typename GridView::template Codim<0>::Iterator ElementIterator;
        enum {dim = Grid::dimension};
        typedef array<BoundaryConditions::Flags, numEq> BCBlockType;

        const GridView& gridview(this->_grid.leafView());
        (*defectGlobal)=0;

        // allocate flag vector to hold flags for essential boundary conditions
        std::vector<BCBlockType> essential(this->vertexmapper.size());
        for (typename std::vector<BCBlockType>::size_type globalIdx=0; globalIdx
                <essential.size(); globalIdx++)
            essential[globalIdx].assign(BoundaryConditions::neumann);

        // iterate through leaf grid
        ElementIterator eendit = gridview.template end<0>();
        for (ElementIterator eIt = gridview.template begin<0>(); eIt
                != eendit; ++eIt) {
            // get geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // get element
            const Element& element = *eIt;
            this->localJacobian().fvGeom.update(element);
            int size = this->localJacobian().fvGeom.numVertices;

            this->localJacobian().setLocalSolution(element);
            this->localJacobian().computeElementData(element);
            bool old = true;
            this->localJacobian().updateVariableData(element, this->localJacobian().uold, old);
            this->localJacobian().updateVariableData(element, this->localJacobian().u);
            this->localJacobian().template localDefect<LeafTag>(element, this->localJacobian().u);

            // begin loop over vertices
            for (int idx=0; idx < size; idx++) {
                int globalId = this->vertexmapper.template map<dim>(element,idx);
                for (int equationnumber = 0; equationnumber < numEq; equationnumber++) {
                    if (this->localJacobian().bc(idx)[equationnumber] == BoundaryConditions::neumann)
                        (*defectGlobal)[globalId][equationnumber]
                                += this->localJacobian().def[idx][equationnumber];
                    else
                        essential[globalId].assign(BoundaryConditions::dirichlet);
                }
            }
        }

        for (typename std::vector<BCBlockType>::size_type globalIdx=0; globalIdx
                <essential.size(); globalIdx++)
            for (int equationnumber = 0; equationnumber < numEq; equationnumber++) {
            if (essential[globalIdx][equationnumber] == BoundaryConditions::dirichlet)
                (*defectGlobal)[globalIdx][equationnumber] = 0;
            }
    }

    void writerestartfile(int restartNum=0)
    {
        enum {dim = Grid::dimension};
        typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

//        exportToDGF(_grid.leafView(), *(this->u), numEq, "primvar", false);

        const int size = vertexmapper.size();
        BlockVector<FieldVector<double, numEq> > data(size);
        data=0;

        VertexIterator endIt = _grid.leafView().template end<dim>();
        for (VertexIterator vIt = _grid.leafView().template begin<dim>(); vIt != endIt;    ++vIt)
        {
            int globalIdx = vertexmapper.map(*vIt);
            for (int equationnumber = 0; equationnumber < numEq;equationnumber++)
            {
                data[globalIdx][equationnumber]=(*(this->u))[globalIdx][equationnumber];
            }
        }

        restartFileName = (boost::format("data-%05d")
                           %restartNum).str();
        exportToDGF(_grid.leafView(), data, (numEq), restartFileName, false);
    }

    virtual void vtkout(const char* name, int k) {}

    const Grid &grid() const
        { return _grid; }



protected:
  ProblemType& problem;
  const Grid& _grid;
  VertexMapper vertexmapper;
  int size;
  std::string restartFileName;

};

}
#endif
