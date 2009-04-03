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
#include "dumux/diffusion/diffusionproblem.hh"

/**
 * @file
 * @brief  Finite Volume Diffusion Model
 * @author Bernd Flemisch, Jochen Fritz; last changed by Markus Wolff
 */

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

 - Grid         a DUNE grid type
 - Scalar        type used for return values
 */
template<class Grid, class Scalar, class VC, class Problem = DiffusionProblem<
        Grid, Scalar, VC> > class FVDiffusion: public Diffusion<Grid, Scalar,
        VC, Problem>
{

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

    FVDiffusion(Grid& grid,Problem& problem) :
    Diffusion<Grid, Scalar, VC, Problem>(grid, problem),
    A(grid.size(
                    this->level(), 0), grid.size(this->level(), 0), (2*dim+1)
            *grid.size(this->level(), 0), BCRSMatrix<MB>::random),
    f(grid.size(this->level(), 0)), solverName_("BiCGSTAB"),
    preconditionerName_("SeqILU0")
    {
        initializeMatrix();
    }

    FVDiffusion(Grid& grid, Problem& problem, std::string solverName,
            std::string preconditionerName) :
    Diffusion<Grid, Scalar, VC, Problem>(grid, problem),
    A(grid.size(
                    this->level(), 0), grid.size(this->level(), 0), (2*dim+1)
            *grid.size(this->level(), 0), BCRSMatrix<MB>::random),
    f(grid.size(this->level(), 0)), solverName_(solverName),
    preconditionerName_(preconditionerName)
    {
        initializeMatrix();
    }

private:

    MatrixType A;
    BlockVector< FieldVector<Scalar,1> > f;
    std::string solverName_;
    std::string preconditionerName_;
};

template<class Grid, class Scalar, class VC,class Problem> void FVDiffusion<Grid, Scalar, VC, Problem>::initializeMatrix()
{

    const GridView& gridView = this->grid.levelView(this->level());

    // determine matrix row sizes
    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = this->diffProblem.variables.diffMapper.map(*eIt);

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
        A.setrowsize(globalIdxI, rowSize);
    }
    A.endrowsizes();

    // determine position of matrix entries
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = this->diffProblem.variables.diffMapper.map(*eIt);

        // add diagonal index
        A.addindex(globalIdxI, globalIdxI);

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
            int globalIdxJ = this->diffProblem.variables.diffMapper.map(*outside);

            // add off diagonal index
            A.addindex(globalIdxI, globalIdxJ);
        }
    }
    A.endindices();

    return;
}

template<class Grid, class Scalar, class VC,class Problem> void FVDiffusion<Grid, Scalar, VC, Problem>::assemble(const Scalar t=0)
{
    // initialization: set matrix A to zero
    A = 0;

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
        int globalIdxI = this->diffProblem.variables.diffMapper.map(*eIt);

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().integrationElement(localPos)*ReferenceElements<Scalar,dim>::general(gt).volume();

        // set right side to zero
        f[globalIdxI] = volume*this->diffProblem.sourcePress(globalPos, *eIt, localPos);

        // get absolute permeability
        FieldMatrix Ki(this->diffProblem.soil().K(globalPos, *eIt, localPos));

        //compute total mobility
        Scalar lambdaI, fractionalWI;
        Scalar satI = this->diffProblem.variables.saturation[globalIdxI];
        lambdaI = this->diffProblem.materialLaw().mobTotal(satI,globalPos, *eIt, localPos);
        if (hasGravity)
        fractionalWI = this->diffProblem.materialLaw().fractionalW(satI,globalPos, *eIt, localPos);

        IntersectionIterator
        isItEnd = gridView.template iend(*eIt);
        for (IntersectionIterator
                isIt = gridView.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {

            // get geometry type of face
            GeometryType faceGT = isIt->geometryInInside().type();

            // center in face's reference element
            const FieldVector<Scalar,dim-1>&
            faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

            // center of face inside volume reference element
            const LocalPosition& localPosFace = ReferenceElements<Scalar,dim>::general(faceGT).position(isIt->indexInInside(),1);

            // get normal vector
            FieldVector<Scalar,dimWorld> unitOuterNormal
            = isIt->unitOuterNormal(faceLocal);

            // get normal vector scaled with volume
            FieldVector<Scalar,dimWorld> integrationOuterNormal= isIt->integrationOuterNormal(faceLocal);
            integrationOuterNormal*= ReferenceElements<Scalar,dim-1>::general(faceGT).volume();

            // get face volume
            Scalar faceVol = isIt->geometry().volume();

            // compute directed permeability vector Ki.n
            FieldVector<Scalar,dim> Kni(0);
            Ki.umv(unitOuterNormal, Kni);

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = this->diffProblem.variables.diffMapper.map(*neighborPointer);

                // compute factor in neighbor
                GeometryType neighborGT = neighborPointer->geometry().type();
                const LocalPosition& localPosNeighbor = ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

                // neighbor cell center in global coordinates
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

                // distance vector between barycenters
                FieldVector<Scalar,dimWorld>
                distVec = globalPos - globalPosNeighbor;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                // get absolute permeability
                FieldMatrix Kj(this->diffProblem.soil().K(globalPosNeighbor, *neighborPointer, localPosNeighbor));

                // compute vectorized permeabilities
                FieldVector<Scalar,dim> Knj(0);
                Kj.umv(unitOuterNormal, Knj);
                Scalar K_n_i = Kni * unitOuterNormal;
                Scalar K_n_j = Knj * unitOuterNormal;
                Scalar Kn = 2 * K_n_i * K_n_j / (K_n_i + K_n_j);
                // compute permeability tangential to intersection and take arithmetic mean
                FieldVector<Scalar,dim> uON = unitOuterNormal;
                FieldVector<Scalar,dim> K_t_i = Kni - (uON *= K_n_i);
                uON = unitOuterNormal;
                FieldVector<Scalar,dim> K_t_j = Knj - (uON *= K_n_j);
                FieldVector<Scalar,dim> Kt = (K_t_i += K_t_j);
                Kt *= 0.5;
                // Build vectorized averaged permeability
                uON = unitOuterNormal;
                FieldVector<Scalar,dim> K = (Kt += (uON *=Kn));

                //compute total mobility
                Scalar lambdaJ, fractionalWJ;
                Scalar satJ = this->diffProblem.variables.saturation[globalIdxJ];

                lambdaJ = this->diffProblem.materialLaw().mobTotal(satJ,globalPosNeighbor, *neighborPointer, localPosNeighbor);
                if (hasGravity)
                fractionalWJ = this->diffProblem.materialLaw().fractionalW(satJ,globalPosNeighbor, *neighborPointer, localPosNeighbor);

                // compute averaged total mobility
                // CAREFUL: Harmonic weightig can generate zero matrix entries,
                // use arithmetic weighting instead:
                Scalar fractionalW;
                Scalar lambda = 0.5*(lambdaI + lambdaJ);
                if (hasGravity)
                fractionalW = 0.5*(fractionalWI + fractionalWJ);

                // update diagonal entry
                Scalar entry = fabs(lambda*faceVol*(K*distVec)/(dist*dist));
                A[globalIdxI][globalIdxI] += entry;

                // set off-diagonal entry
                A[globalIdxI][globalIdxJ] = -entry;

                if (hasGravity)
                {
                    Scalar factor = fractionalW*(this->diffProblem.wettingPhase.density())
                    + (1 - fractionalW)*(this->diffProblem.nonWettingPhase.density());
                    f[globalIdxI] -= factor*lambda*faceVol*(K*gravity);
                }

                if (this->diffProblem.capillarity)
                {
                    // calculate saturation gradient
                    FieldVector<Scalar,dim> satGradient = distVec;
                    satGradient *= (satJ - satI)/(dist*dist);

                    // arithmetic average of the permeability
                    K = ((Kni + Knj) *= 0.5);

                    // capillary pressure w.r.t. saturation
                    Scalar pCI = this->diffProblem.materialLaw().pC(satI,globalPos, *eIt, localPos);
                    Scalar pCJ = this->diffProblem.materialLaw().pC(satJ,globalPosNeighbor, *neighborPointer, localPosNeighbor);

                    // mobility of the nonwetting phase
                    Scalar lambdaN = 0.5*(this->diffProblem.materialLaw().mobN(1 - satI,globalPos, *eIt, localPos)
                            + this->diffProblem.materialLaw().mobN(1 - satJ,globalPosNeighbor, *neighborPointer, localPosNeighbor));

                    // calculate capillary pressure gradient
                    FieldVector<Scalar,dim> pCGradient = distVec;
                    pCGradient *= -(pCJ - pCI)/(dist*dist);

                    f[globalIdxI] += lambdaN*faceVol*(K*pCGradient);
                }
            }
            // boundary face

            else
            {
                // center of face in global coordinates
                const GlobalPosition& globalPosFace = isIt->geometry().global(faceLocal);

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
                    A[globalIdxI][globalIdxI] -= lambda*faceVol*(Kni*distVec)/(dist*dist);
                    Scalar g = this->diffProblem.dirichletPress(globalPosFace, *eIt, localPosFace);
                    f[globalIdxI] -= lambda*faceVol*g*(Kni*distVec)/(dist*dist);

                    if (hasGravity)
                    {
                        Scalar factor = fractionalW*(this->diffProblem.wettingPhase.density())
                        + (1 - fractionalW)*(this->diffProblem.nonWettingPhase.density());
                        f[globalIdxI] -= factor*lambda*faceVol*(Kni*gravity);
                    }
                    if (this->diffProblem.capillarity)
                    {
                        Scalar satJ = this->diffProblem.dirichletSat(globalPosFace, *eIt, localPosFace);

                        // distance vector between barycenters
                        FieldVector<Scalar,dimWorld>
                        distVec = globalPos - globalPosFace;

                        // compute distance between cell centers
                        Scalar dist = distVec.two_norm();

                        // calculate saturation gradient
                        FieldVector<Scalar,dim> satGradient = distVec;
                        satGradient *= (satJ - satI)/(dist*dist);

                        // capillary pressure w.r.t. saturation
                        Scalar pCI = this->diffProblem.materialLaw().pC(satI,globalPos, *eIt, localPos);
                        Scalar pCJ = this->diffProblem.materialLaw().pC(satJ,globalPosFace, *eIt, localPosFace);

                        // mobility of the nonwetting phase
                        Scalar lambdaN = 0.5*(this->diffProblem.materialLaw().mobN(1 - satI,globalPos, *eIt, localPos)
                                + this->diffProblem.materialLaw().mobN(1 - satJ,globalPosFace, *eIt, localPosFace));

                        // calculate capillary pressure gradient
                        FieldVector<Scalar,dim> pCGradient = distVec;
                        pCGradient *= -(pCJ - pCI)/(dist*dist);

                        f[globalIdxI] += lambdaN*faceVol*(Kni*pCGradient);
                    }
                }
                else
                {
                    Scalar J = this->diffProblem.neumannPress(globalPosFace, *eIt, localPosFace);
                    f[globalIdxI] -= faceVol*J;
                }
            }
        }
        // end all intersections
    } // end grid traversal
    return;
}

template<class Grid, class Scalar, class VC,class Problem> void FVDiffusion<Grid, Scalar, VC, Problem>::solve()
{
    MatrixAdapter<MatrixType,Vector,Vector> op(A);
    InverseOperatorResult r;

    if (preconditionerName_ == "SeqILU0")
    {
        SeqILU0<MatrixType,Vector,Vector> preconditioner(A, 1.0);
        if (solverName_ == "CG")
        {
            CGSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
            solver.apply((this->diffProblem.variables.pressure), f, r);
        }
        else if (solverName_ == "BiCGSTAB")
        {
            BiCGSTABSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 1);
            solver.apply((this->diffProblem.variables.pressure), f, r);
        }
        else
        DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination "
                << preconditionerName_<< " and "<< solverName_ << ".");
    }
    else if (preconditionerName_ == "SeqPardiso")
    {
        SeqPardiso<MatrixType,Vector,Vector> preconditioner(A);
        if (solverName_ == "Loop")
        {
            LoopSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
            solver.apply(this->diffProblem.variables.pressure, f, r);
        }
        else
        DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination "
                << preconditionerName_<< " and "<< solverName_ << ".");
    }
    else
    DUNE_THROW(NotImplemented, "FVDiffusion :: solve : preconditioner "
            << preconditionerName_ << ".");
    //    printmatrix(std::cout, A, "global stiffness matrix", "row", 11, 3);
    //    printvector(std::cout, f, "right hand side", "row", 200, 1, 3);
    //    printvector(std::cout, (this->diffProblem.variables.pressure), "pressure", "row", 200, 1, 3);
    return;
}

}
#endif
