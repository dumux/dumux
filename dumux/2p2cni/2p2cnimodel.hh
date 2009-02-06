// $Id$

#ifndef DUNE_TWOPHASEHEATMODEL_HH
#define DUNE_TWOPHASEHEATMODEL_HH

#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include "dumux/operators/p1operatorextended.hh"
#include "dumux/nonlinear/nonlinearmodel.hh"
#include "dumux/fvgeometry/fvelementgeometry.hh"
#include "dumux/io/exporttodgf.hh"
#include <boost/format.hpp>

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar, class ProblemType, class LocalJacobian,
        class FunctionType, class OperatorAssembler> class TwoPhaseHeatModel: public NonlinearModel<
        Grid, Scalar, ProblemType, LocalJacobian, FunctionType,
        OperatorAssembler>
{
public:
    typedef NonlinearModel<Grid, Scalar, ProblemType, LocalJacobian,
            FunctionType, OperatorAssembler> ThisNonlinearModel;

    TwoPhaseHeatModel(const Grid& grid, ProblemType& prob) :
        ThisNonlinearModel(grid, prob), uOldTimeStep(grid)
    {
    }

    TwoPhaseHeatModel(const Grid& grid, ProblemType& prob, int level) :
        ThisNonlinearModel(grid, prob, level), uOldTimeStep(grid, level)
    {
    }

    virtual void initial() = 0;

    virtual void restart()
    {
    }

    virtual void update(double& dt) = 0;

    virtual void solve() = 0;

    FunctionType uOldTimeStep;
};

template<class Grid, class Scalar, class ProblemType, class LocalJac,
        int numEq = 3> class LeafP1TwoPhaseModel: public TwoPhaseHeatModel<
        Grid, Scalar, ProblemType, LocalJac, LeafP1FunctionExtended<Grid,
                Scalar, numEq> , LeafP1OperatorAssembler<Grid, Scalar, numEq> >
{
public:
    // define the function type:
    typedef LeafP1FunctionExtended<Grid, Scalar, numEq> FunctionType;

    // define the operator assembler type:
    typedef LeafP1OperatorAssembler<Grid, Scalar, numEq> OperatorAssembler;

    typedef TwoPhaseHeatModel<Grid, Scalar, ProblemType, LocalJac,
            FunctionType, OperatorAssembler> ThisTwoPhaseHeatModel;

    typedef LeafP1TwoPhaseModel<Grid, Scalar, ProblemType, LocalJac, numEq>
            ThisType;

    typedef LocalJac LocalJacobian;

    // mapper: one data element per vertex
    template<int dim> struct P1Layout
    {
        bool contains(Dune::GeometryType gt)
        {
            return gt.dim() == 0;
        }
    };

typedef    typename Grid::LeafGridView GridView;
    typedef typename GridView::IndexSet IS;
    typedef MultipleCodimMultipleGeomTypeMapper<Grid,IS,P1Layout> VertexMapper;
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator
    IntersectionIterator;

    LeafP1TwoPhaseModel(const Grid& grid, ProblemType& prob) :
    ThisTwoPhaseHeatModel(grid, prob), problem(prob), _grid(grid), vertexmapper(grid, grid.leafIndexSet()), size((*(this->u)).size())
    {
    }

    virtual void update(double &dt)
    {
        DUNE_THROW(NotImplemented, "This method is obsolete. Use updateModel()!");
    }

    virtual void initial()
    {}

    virtual void restart(int restartNum=0)
    {}

    virtual double computeFlux ()
    {
        typedef typename Grid::Traits::template Codim<0>::Entity Element;
        typedef typename Grid::ctype CoordScalar;
        typedef typename GridView::template Codim<0>::Iterator Iterator;
        enum
        {   dim = Grid::dimension};
        enum
        {   dimworld = Grid::dimensionworld};
        double sign;
        const GridView& gridview(_grid.leafView());
        Iterator eendit = gridview.template end<0>();
        FieldVector<Scalar,numEq> flux(0);
        double Flux(0);

        for (Iterator eIt = gridview.template begin<0>(); eIt != eendit; ++eIt) // loop over all entities

        {

            // get geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // get element
            const Element& element = *eIt;

            FVElementGeometry<Grid> fvGeom;
            fvGeom.update(element);

            for (int k = 0; k < fvGeom.numEdges; k++)
            {
                int idx_i = fvGeom.subContVolFace[k].i;

                int idx_j = fvGeom.subContVolFace[k].j;

                int flag_i, flag_j;

                // 2D case: give y or x value of the line over which flux is to be
                //            calculated.
                // up to now only flux calculation to lines or planes (3D) parallel to
                // x, y and z axis possible

                // Flux across plane with z = 80 numEq
                if(fvGeom.subContVol[idx_i].global[2] < 80.)
                flag_i = 1;
                else flag_i = -1;

                if(fvGeom.subContVol[idx_j].global[2] < 80.)
                flag_j = 1;
                else flag_j = -1;

                if(flag_i == flag_j)
                {
                    sign = 0;
                }
                else
                {
                    if(flag_i> 0)
                    sign = -1;
                    else sign = 1;}

                // get variables

                if(flag_i != flag_j)
                {
                    this->localJacobian().setLocalSolution(element);
                    this->localJacobian().computeElementData(element);
                    this->localJacobian().updateVariableData(element, this->localJacobian().u);

                    flux = this->localJacobian().computeA(element, this->localJacobian().u, k);
                    Flux += sign*flux[1];
                }
            }

        }
        return Flux; // Co2 flux
    }

    virtual double totalCO2Mass()
    {
        typedef typename Grid::Traits::template Codim<0>::Entity Element;
        typedef typename Grid::ctype CoordScalar;
        typedef typename GridView::template Codim<0>::Iterator Iterator;
        enum
        {   dim = Grid::dimension};
        enum
        {   dimworld = Grid::dimensionworld};
        enum
        {   gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase state
        const GridView& gridview(_grid.leafView());
        double totalMass = 0;
        double minSat = 1e100;
        double maxSat = -1e100;
        double minP = 1e100;
        double maxP = -1e100;
        double minTe = 1e100;
        double maxTe = -1e100;
        double minX = 1e100;
        double maxX = -1e100;

        // iterate through leaf grid an evaluate c0 at element center
        Iterator eendit = gridview.template end<0>();
        for (Iterator eIt = gridview.template begin<0>(); eIt
                != eendit; ++eIt)
        {
            // get geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // get element
            const Element& element = *eIt;

            FVElementGeometry<Grid> fvGeom;
            fvGeom.update(element);

            const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,Scalar,dim>::value_type
            &sfs=Dune::LagrangeShapeFunctions<CoordScalar, Scalar, dim>::general(gt,
                    1);
            int size = sfs.size();

            for (int idx = 0; idx < size; idx++)
            {
                // get element center in reference element
                const Dune::FieldVector<CoordScalar,dim>&localPos = sfs[idx].position();

                // get globalPos coordinate of element center
                Dune::FieldVector<CoordScalar,dimworld> globalPos = eIt->geometry().global(localPos);

                int globalIdx = vertexmapper.template map<dim>(element,
                        sfs[idx].entity());

                int state;
                state = this->localJacobian().sNDat[globalIdx].phaseState;
                Scalar vol = fvGeom.subContVol[idx].volume;
                Scalar poro = this->problem.soil().porosity(globalPos, element, localPos);

                Scalar rhoN = (*(this->localJacobian().outDensityN))[globalIdx];
                Scalar rhoW = (*(this->localJacobian().outDensityW))[globalIdx];
                Scalar satN = (*(this->localJacobian().outSaturationN))[globalIdx];
                Scalar satW = (*(this->localJacobian().outSaturationW))[globalIdx];
                Scalar xAW = (*(this->localJacobian().outMassFracAir))[globalIdx];
                Scalar xWN = (*(this->localJacobian().outMassFracWater))[globalIdx];
                Scalar xAN = 1 - xWN;
                Scalar pW = (*(this->u))[globalIdx][0];
                Scalar Te = (*(this->u))[globalIdx][2];
                Scalar mass = vol * poro * (satN * rhoN * xAN + satW * rhoW * xAW);

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
        std::cout << "nonwetting phase saturation: min = "<< minSat
        << ", max = "<< maxSat << std::endl;
        std::cout << "wetting phase pressure: min = "<< minP
        << ", max = "<< maxP << std::endl;
        std::cout << "mass fraction CO2: min = "<< minX
        << ", max = "<< maxX << std::endl;
        std::cout << "temperature: min = "<< minTe
        << ", max = "<< maxTe << std::endl;

        return totalMass;
    }

    virtual void globalDefect(FunctionType& defectGlobal)
    {
        typedef typename Grid::Traits::template Codim<0>::Entity Element;
        typedef typename Grid::ctype CoordScalar;
        typedef typename GridView::template Codim<0>::Iterator Iterator;
        enum
        {   dim = Grid::dimension};
        typedef array<BoundaryConditions::Flags, numEq> BCBlockType;

        const GridView& gridview(_grid.leafView());
        (*defectGlobal)=0;
        // allocate flag vector to hold flags for essential boundary conditions
        std::vector<BCBlockType> essential(this->vertexmapper.size());
        for (typename std::vector<BCBlockType>::size_type globalIdx=0; globalIdx
                <essential.size(); globalIdx++)
        essential[globalIdx].assign(BoundaryConditions::neumann);
        // iterate through leaf grid
        Iterator eendit = gridview.template end<0>();
        for (Iterator eIt = gridview.template begin<0>(); eIt
                != eendit; ++eIt)
        {
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
            for (int idx=0; idx < size; idx++)
            {
                int globalIdx = this->vertexmapper.template map<dim>(element,idx);
                for (int equationnumber = 0; equationnumber < numEq; equationnumber++)
                {
                    if (this->localJacobian().bc(idx)[equationnumber] == BoundaryConditions::neumann)
                    (*defectGlobal)[globalIdx][equationnumber]
                    += this->localJacobian().def[idx][equationnumber];
                    else
                    essential[globalIdx].assign(BoundaryConditions::dirichlet);
                }
            }
        }

        for (typename std::vector<BCBlockType>::size_type globalIdx=0; globalIdx
                <essential.size(); globalIdx++)
        for (int equationnumber = 0; equationnumber < numEq; equationnumber++)
        {
            if (essential[globalIdx][equationnumber] == BoundaryConditions::dirichlet)
            (*defectGlobal)[globalIdx][equationnumber] = 0;
        }
    }

    virtual void vtkout(const char* name, int k)
    {
    }

    void writerestartfile(int restartNum=0)
    {
        enum
        {   dim = Grid::dimension};
        typedef typename GridView::template Codim<dim>::Iterator Iterator;

        //        exportToDGF(_grid.leafView(), *(this->u), numEq, "primvar", false);

        const int size = vertexmapper.size();
        BlockVector<FieldVector<double, numEq+1> > data(size);
        data=0;

        Iterator endIt = _grid.leafView().template end<dim>();
        for (Iterator vIt = _grid.leafView().template begin<dim>(); vIt != endIt; ++vIt)
        {
            int globalIdx = vertexmapper.map(*vIt);
            for (int eq = 0; eq < numEq; eq++)
            {
                data[globalIdx][eq]=(*(this->u))[globalIdx][eq];
            }
            data[globalIdx][numEq]=this->localJacobian().sNDat[globalIdx].phaseState;
        }
        restartFileName = (boost::format("data-%05d")
                %restartNum).str();
        exportToDGF(_grid.leafView(), data, (numEq+1), restartFileName, false);
    }
    const Grid &grid() const
    {   return _grid;}

protected:
    ProblemType& problem;
    const Grid& _grid;
    VertexMapper vertexmapper;
    int size;
    std::string restartFileName;
};

}
#endif
