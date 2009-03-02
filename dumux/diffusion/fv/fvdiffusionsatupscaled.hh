// $Id$
#ifndef DUNE_FVDIFFUSION_HH
#define DUNE_FVDIFFUSION_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include "dumux/diffusion/diffusion.hh"
#include "dumux/pardiso/pardiso.hh"

namespace Dune
{
//! \ingroup diffusion
//! Finite Volume Diffusion Model
/*! Provides a Finite Volume implementation for the evaluation
 * of equations of the form
 * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = q, \f$,
 * \f$p = g\f$ on \f$\Gamma_1\f$, and
 * \f$\lambda K \text{grad}\, p \cdot \mathbf{n} = J\f$
 * on \f$\Gamma_2\f$. Here,
 * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability,
 * and \f$\lambda\f$ the total mobility, possibly depending on the
 * saturation, \f$q\f$ the source term.
 Template parameters are:

 - G         a DUNE grid type
 - RT        type used for return values
*/
template<class Grid, class Scalar, class VC> class FVDiffusion: public Diffusion<
    Grid, Scalar, VC>
{
    template<int dim> struct ElementLayout
    {
        bool contains(GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    enum
        {
            dim = Grid::dimension
        };
    enum
        {
            dimWorld = Grid::dimensionworld
        };

    typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LevelGridView GridView;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename Grid::template Codim<0>::HierarchicIterator
    HierarchicIterator;

    typedef MultipleCodimMultipleGeomTypeMapper<Grid,IndexSet,ElementLayout> ElementMapper;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    typedef Dune::FieldMatrix<Scalar, 1, 1> MB;
    typedef BCRSMatrix<MB> MatrixType;
    typedef BlockVector<FieldVector<Scalar, 1> > Vector;

public:

    void assemble(const Scalar t);

    void solve();

    void pressure(const Scalar t=0)
    {
        assemble(t);
        solve();

        return;
    }

    void initializeMatrix();

    FVDiffusion(Grid& grid, FractionalFlowProblem<Grid, Scalar, VC>& prob) :
        Diffusion<Grid, Scalar, VC>(grid, prob),
        elementMapper(grid, grid.levelIndexSet(this->level())),
        globalMatrix_(grid.size(this->level(), 0), grid.size(this->level(), 0), (2*dim+1)*grid.size(this->level(), 0), BCRSMatrix<MB>::random),
        rightHandSide_(grid.size(this->level(), 0)), solverName_("Loop"),
        preconditionerName_("SeqPardiso")
    {
        initializeMatrix();
    }

    FVDiffusion(Grid& grid, FractionalFlowProblem<Grid, Scalar, VC>& prob, std::string solverName,
                std::string preconditionerName) :
        Diffusion<Grid, Scalar, VC>(grid, prob),
        elementMapper(grid, grid.levelIndexSet(this->level())),globalMatrix_(grid.size(
                                                                                       this->level(), 0), grid.size(this->level(), 0), (2*dim+1)
                                                                             *grid.size(this->level(), 0), BCRSMatrix<MB>::random),
        rightHandSide_(grid.size(this->level(), 0)), solverName_(solverName),
        preconditionerName_(preconditionerName)
    {
        initializeMatrix();
    }

    ElementMapper elementMapper;

private:
    MatrixType globalMatrix_;
    BlockVector< FieldVector<Scalar,1> > rightHandSide_;
    std::string solverName_;
    std::string preconditionerName_;
};

template<class Grid, class Scalar, class VC> void FVDiffusion<Grid, Scalar, VC>::initializeMatrix()
{

    const GridView& gridView = this->grid.levelView(this->level());

    // determine matrix row sizes
    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {

        int globalIdxI = elementMapper.map(*eIt);

        // initialize row size
        int rowSize = 1;

        // run through all intersections with neighbors
        IntersectionIterator
            isItEnd = gridView.template iend(*eIt);
        for (IntersectionIterator
                 isIt = gridView.template ibegin(*eIt); isIt
                 !=isItEnd; ++isIt)
            if (isIt->neighbor())
                rowSize++;
        globalMatrix_.setrowsize(globalIdxI, rowSize);
    }
    globalMatrix_.endrowsizes();

    // determine position of matrix entries
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = elementMapper.map(*eIt);

        // add diagonal index
        globalMatrix_.addindex(globalIdxI, globalIdxI);

        // run through all intersections with neighbors
        IntersectionIterator
            isItEnd = gridView.template iend(*eIt);
        for (IntersectionIterator
                 isIt = gridView.template ibegin(*eIt); isIt
                 !=isItEnd; ++isIt)
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer outside = isIt->outside();
                int globalIdxJ = elementMapper.map(*outside);

                // add off diagonal index
                globalMatrix_.addindex(globalIdxI, globalIdxJ);
            }
    }
    globalMatrix_.endindices();

    return;
}

template<class Grid, class Scalar, class VC> void FVDiffusion<Grid, Scalar, VC>::assemble(const Scalar t=0)
{
    // initialization: set matrix globalMatrix_ to zero
    globalMatrix_ = 0;

    const GridView& gridView = this->grid.levelView(this->level());

    // find out whether gravity effects are relevant
    bool hasGravity = false;
    const FieldVector<Scalar,dim>& gravity = this->diffProblem.gravity();
    for (int k = 0; k < dim; k++)
        if (gravity[k] != 0)
            hasGravity = true;

    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell geometry type
        GeometryType gt = eIt->geometry().type();

        // cell center in reference element
        const LocalPosition& localPos = ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        const GlobalPosition& globalPos = eIt->geometry().global(localPos);

        // cell index
        int globalIdxI = elementMapper.map(*eIt);

        Element& entityI = *eIt;
        ElementPointer fatherPointerI = entityI.father();
        int fatherLevelI = fatherPointerI.level();
        while (this->diffProblem.variables.transLevel != fatherLevelI)
        {
            Element& fatherErentityI = *fatherPointerI;
            fatherPointerI = fatherErentityI.father();
            fatherLevelI = fatherPointerI.level();
        }

        int globalIdxCoarseI = this->diffProblem.variables.transMapper.map(*fatherPointerI);

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().integrationElement(localPos)*ReferenceElements<Scalar,dim>::general(gt).volume();

        // set right side to zero
        rightHandSide_[globalIdxI] = volume*this->diffProblem.sourcePress(globalPos, *eIt, localPos);

        // get absolute permeability
        FieldMatrix kI(this->diffProblem.soil.K(globalPos, *eIt, localPos));

        //compute total mobility
        Scalar lambdaI, fractionalWI;
        Scalar satI = this->diffProblem.variables.saturation[globalIdxCoarseI];
        lambdaI = this->diffProblem.materialLaw.mobTotal(satI,globalPos, *eIt, localPos);
        if (hasGravity)
            fractionalWI = this->diffProblem.materialLaw.fractionalW(satI,globalPos, *eIt, localPos);

        IntersectionIterator
            isItEnd = gridView.template iend(*eIt);
        for (IntersectionIterator
                 isIt = gridView.template ibegin(*eIt); isIt
                 !=isItEnd; ++isIt)
        {
            // get geometry type of face
            GeometryType faceGT = isIt->intersectionSelfLocal().type();

            // center in face's reference element
            const FieldVector<Scalar,dim-1>&
                faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

            // center of face inside volume reference element
            const LocalPosition& localPosFace = ReferenceElements<Scalar,dim>::general(faceGT).position(isIt->numberInSelf(),1);

            // get normal vector
            FieldVector<Scalar,dimWorld> unitOuterNormal
                = isIt->unitOuterNormal(faceLocal);

            // get normal vector scaled with volume
            FieldVector<Scalar,dimWorld> integrationOuterNormal
                = isIt->integrationOuterNormal(faceLocal);
            integrationOuterNormal
                *= ReferenceElements<Scalar,dim-1>::general(faceGT).volume();

            // get face volume
            Scalar faceVol = 1;
            switch (Grid::dimension)
            {
            case 1: break;
            default: faceVol = isIt->intersectionGlobal().volume();
                break;
            }

            // compute directed permeability vector kI.n
            FieldVector<Scalar,dim> kNI(0);
            kI.umv(unitOuterNormal, kNI);

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = elementMapper.map(*neighborPointer);

                Element& entityJ = *neighborPointer;
                ElementPointer fatherPointerJ = entityJ.father();
                int fatherLevelJ = fatherPointerJ.level();
                while (this->diffProblem.variables.transLevel != fatherLevelJ)
                {
                    Element& fatherEntityJ = *fatherPointerJ;
                    fatherPointerJ = fatherEntityJ.father();
                    fatherLevelJ = fatherPointerJ.level();
                }

                int globalIdxCoarseJ = this->diffProblem.variables.transMapper.map(*fatherPointerJ);

                // compute factor in neighbor
                GeometryType neighborGT = neighborPointer->geometry().type();
                const LocalPosition&
                    localPosNeighbor = ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

                // neighbor cell center in globalPos coordinates
                const GlobalPosition&
                    globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

                // distance vector between barycenters
                FieldVector<Scalar,dimWorld>
                    distVec = globalPos - globalPosNeighbor;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                // get absolute permeability
                FieldMatrix kJ(this->diffProblem.soil.K(globalPosNeighbor, *neighborPointer, localPosNeighbor));

                // compute vectorized permeabilities
                FieldVector<Scalar,dim> kNJ(0);
                kJ.umv(unitOuterNormal, kNJ);
                Scalar k_N_I = kNI * unitOuterNormal;
                Scalar k_N_J = kNJ * unitOuterNormal;
                Scalar Kn = 2 * k_N_I * k_N_J / (k_N_I + k_N_J);
                // compute permeability tangential to intersection and take arithmetic mean
                FieldVector<Scalar,dim> uON = unitOuterNormal;
                FieldVector<Scalar,dim> k_T_I = kNI - (uON *= k_N_I);
                uON = unitOuterNormal;
                FieldVector<Scalar,dim> k_T_J = kNJ - (uON *= k_N_J);
                FieldVector<Scalar,dim> kT = (k_T_I += k_T_J);
                kT *= 0.5;
                // Build vectorized averaged permeability
                uON = unitOuterNormal;
                FieldVector<Scalar,dim> permeability = (kT += (uON *=Kn));

                //compute total mobility
                Scalar lambdaJ, fractionalWJ;
                Scalar satJ = this->diffProblem.variables.saturation[globalIdxCoarseJ];

                lambdaJ = this->diffProblem.materialLaw.mobTotal(satJ,globalPosNeighbor, *neighborPointer, localPosNeighbor);
                if (hasGravity)
                {
                    fractionalWJ = this->diffProblem.materialLaw.fractionalW(satJ,globalPosNeighbor, *neighborPointer, localPosNeighbor);
                }

                // compute averaged total mobility
                // CAREFUL: Harmonic weightig can generate zero matrix entries,
                // use arithmetic weighting instead:
                Scalar fractionalW;
                Scalar lambda = 0.5*(lambdaI + lambdaJ);
                if (hasGravity)
                {
                    fractionalW = 0.5*(fractionalWI + fractionalWJ);
                }

                // update diagonal entry
                Scalar entry = fabs(lambda*faceVol*(permeability*distVec)/(dist*dist));
                globalMatrix_[globalIdxI][globalIdxI] += entry;

                // set off-diagonal entry
                globalMatrix_[globalIdxI][globalIdxJ] = -entry;

                if (hasGravity)
                {
                    Scalar factor = fractionalW*(this->diffProblem.wettingPhase.density())
                        + (1 - fractionalW)*(this->diffProblem.nonWettingPhase.density());
                    rightHandSide_[globalIdxI] += factor*lambda*faceVol*(permeability*gravity);
                }
            }
            // boundary face

            else
            {
                // center of face in global coordinates
                const GlobalPosition& globalPosFace = isIt->intersectionGlobal().global(faceLocal);

                // compute total mobility
                Scalar fractionalW = 1.;
                Scalar lambda = lambdaI;
                if (hasGravity) fractionalW = fractionalWI;

                //get boundary condition for boundary face center
                BoundaryConditions::Flags bctype = this->diffProblem.bctypePress(globalPosFace, *eIt, localPosFace);
                if (bctype == BoundaryConditions::dirichlet)
                {
                    FieldVector<Scalar,dimWorld> distVec(globalPos - globalPosFace);
                    Scalar dist = distVec.two_norm();
                    globalMatrix_[globalIdxI][globalIdxI] -= lambda*faceVol*(kNI*distVec)/(dist*dist);
                    Scalar g = this->diffProblem.dirichletPress(globalPosFace, *eIt, localPosFace);
                    rightHandSide_[globalIdxI] -= lambda*faceVol*g*(kNI*distVec)/(dist*dist);

                    if (hasGravity)
                    {
                        Scalar factor = fractionalW*(this->diffProblem.wettingPhase.density())
                            + (1 - fractionalW)*(this->diffProblem.nonWettingPhase.density());
                        rightHandSide_[globalIdxI] += factor*lambda*faceVol*(kNI*gravity);
                    }
                }
                else
                {
                    Scalar neumannFlux = this->diffProblem.neumannPress(globalPosFace, *eIt, localPosFace);
                    rightHandSide_[globalIdxI] -= faceVol*neumannFlux;
                }
            }
        }
        // end all intersections
    } // end grid traversal
    return;
}

template<class Grid, class Scalar, class VC> void FVDiffusion<Grid, Scalar, VC>::solve()
{
    MatrixAdapter<MatrixType,Vector,Vector> op(globalMatrix_);
    InverseOperatorResult r;

    if (preconditionerName_ == "SeqILU0")
    {
        SeqILU0<MatrixType,Vector,Vector> preconditioner(globalMatrix_, 1.0);
        if (solverName_ == "CG")
        {
            CGSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
            solver.apply((this->diffProblem.variables.pressure), rightHandSide_, r);
        }
        else if (solverName_ == "BiCGSTAB")
        {
            BiCGSTABSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
            solver.apply((this->diffProblem.variables.pressure), rightHandSide_, r);
        }
        else
            DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination "
                       << preconditionerName_<< " and "<< solverName_ << ".");
    }
    else if (preconditionerName_ == "SeqPardiso")
    {
        SeqPardiso<MatrixType,Vector,Vector> preconditioner(globalMatrix_);
        if (solverName_ == "Loop")
        {
            LoopSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
            solver.apply(this->diffProblem.variables.pressure, rightHandSide_, r);
        }
        else
            DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination "
                       << preconditionerName_<< " and "<< solverName_ << ".");
    }
    else
        DUNE_THROW(NotImplemented, "FVDiffusion :: solve : preconditioner "
                   << preconditionerName_ << ".");
    //    printmatrix(std::cout, globalMatrix_, "global stiffness matrix", "row", 11, 3);
    //    printvector(std::cout, rightHandSide_, "right hand side", "row", 200, 1, 3);
    //    printvector(std::cout, (this->diffProblem.variables.pressure), "pressure", "row", 200, 1, 3);
    return;
}

}
#endif
