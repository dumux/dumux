// $Id$

#ifndef DUNE_FVDIFFUSION_PW_HH
#define DUNE_FVDIFFUSION_PW_HH

// dune environent:
#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

// dumux environment
#include "dumux/diffusion/diffusionproblem.hh"
#include "dumux/diffusion/diffusion_phasepressure.hh"
#include "dumux/pardiso/pardiso.hh"

/**
 * @file
 * @brief  Finite Volume Diffusion Model
 * @author Bernd Flemisch, Jochen Fritz, Markus Wolff
 */

namespace Dune
{
//! \ingroup diffusion
//! Finite Volume Diffusion Model
/*! Provides a Finite Volume implementation for the evaluation
 * of equations of the form
 * \f$ - \text{div}\, (\lambda K \text{grad}\, p_w + f_n \text{grad}\, p_c + \sum f_\alpha \rho_\alpha g  \text{grad}\, z) = q, \f$,
 * \f$p = p_D\f$ on \f$\Gamma_1\f$, and
 * \f$ \lambda K \text{grad}\, p_w + f_n \text{grad}\, p_c + \sum f_\alpha \rho_\alpha g  \text{grad}\, z \cdot \mathbf{n} = q_N\f$
 * on \f$\Gamma_2\f$. Here,
 * \f$p_w\f$ denotes the wetting phase pressure, \f$K\f$ the absolute permeability,
 * and \f$\lambda\f$ the total mobility, possibly depending on the
 * saturation, \f$q\f$ the source term.
 Template parameters are:

 - Grid         a DUNE grid type
 - Scalar        type used for return values
 */
template<class Grid, class Scalar, class VC, class Problem = DiffusionProblem<Grid, Scalar, VC> > class FVDiffusion: public Diffusion<
        Grid, Scalar, VC, Problem>
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

    void assemble(bool first, const Scalar t);

    void solve();

    void pressure(bool first, const Scalar t=0)
    {
        assemble(first, t);
        solve();
        return;
    }

    void initializeMatrix();

    FVDiffusion(Grid& grid, Problem& problem) :
    Diffusion<Grid, Scalar, VC, Problem>(grid, problem),
    A_(this->diffProblem.variables.gridSizeDiffusion, this->diffProblem.variables.gridSizeDiffusion,
            (2*dim+1) * this->diffProblem.variables.gridSizeDiffusion, BCRSMatrix<MB>::random),
    f_(this->diffProblem.variables.gridSizeDiffusion), solverName_("BiCGSTAB"),
    preconditionerName_("SeqILU0"), gravity_(problem.gravity())
    {
        initializeMatrix();
    }

    FVDiffusion(Grid& grid, Problem& problem, std::string solverName,
            std::string preconditionerName) :
    Diffusion<Grid, Scalar, VC, Problem>(grid, problem),
    A_(this->diffProblem.variables.gridSizeDiffusion, this->diffProblem.variables.gridSizeDiffusion,
            (2*dim+1)*this->diffProblem.variables.gridSizeDiffusion, BCRSMatrix<MB>::random),
    f_(this->diffProblem.variables.gridSizeDiffusion), solverName_(solverName),
    preconditionerName_(preconditionerName), gravity_(problem.gravity())
    {
        initializeMatrix();
    }

private:

    MatrixType A_;
    BlockVector< FieldVector<Scalar,1> > f_;
    std::string solverName_;
    std::string preconditionerName_;

    FieldVector<Scalar,dimWorld> gravity_; // [m/s^2]
};

template<class Grid, class Scalar, class VC, class Problem> void FVDiffusion<Grid, Scalar, VC, Problem>::initializeMatrix()
{

    const GridView& gridView = this->grid.levelView(this->diffProblem.variables.levelDiffusion);

    // determine matrix row sizes
    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = this->diffProblem.variables.indexSetDiffusion.index(*eIt);

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
        A_.setrowsize(globalIdxI, rowSize);
    }
    A_.endrowsizes();

    // determine position of matrix entries
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = this->diffProblem.variables.indexSetDiffusion.index(*eIt);

        // add diagonal index
        A_.addindex(globalIdxI, globalIdxI);

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
            int globalIdxJ = this->diffProblem.variables.indexSetDiffusion.index(*outside);

            // add off diagonal index
            A_.addindex(globalIdxI, globalIdxJ);
        }
    }
    A_.endindices();

    return;
}

template<class Grid, class Scalar, class VC, class Problem> void FVDiffusion<Grid, Scalar, VC, Problem>::assemble(bool first, const Scalar t=0)
{
    // initialization: set matrix A_ to zero
    A_ = 0;
    f_=0;

    const GridView& gridView = this->grid.levelView(this->diffProblem.variables.levelDiffusion);

    Scalar densityW = this->diffProblem.wettingPhase.density();
    Scalar densityNW = this->diffProblem.nonWettingPhase.density();
    Scalar viscosityW = this->diffProblem.wettingPhase.viscosity();
    Scalar viscosityNW = this->diffProblem.nonWettingPhase.viscosity();

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
        int globalIdxI = this->diffProblem.variables.indexSetDiffusion.index(*eIt);

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().integrationElement(localPos)*ReferenceElements<Scalar,dim>::general(gt).volume();

        // set right side to zero
        f_[globalIdxI] = volume*this->diffProblem.sourcePress(globalPos, *eIt, localPos);

        // get absolute permeability
        FieldMatrix permeabilityI(this->diffProblem.soil().K(globalPos, *eIt, localPos));

        //compute total mobility
        Scalar lambdaI, fractionalWI, fractionalNWI;

        // get mobilities and fractional flow factors
        lambdaI = this->diffProblem.variables.mobilityWetting[globalIdxI] + this->diffProblem.variables.mobilityNonWetting[globalIdxI];
        fractionalWI = this->diffProblem.variables.fracFlowFuncWetting[globalIdxI];
        fractionalNWI = this->diffProblem.variables.fracFlowFuncNonWetting[globalIdxI];

        Scalar pcI = this->diffProblem.variables.capillaryPressure[globalIdxI];

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


            // compute directed permeability vector permeabilityI.n
            FieldVector<Scalar,dim> normalPermeabilityI(0);
            permeabilityI.umv(unitOuterNormal, normalPermeabilityI);

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = this->diffProblem.variables.indexSetDiffusion.index(*neighborPointer);
//                std::cout<<"index J = "<<globalIdxJ<<std::endl;

                // compute factor in neighbor
                GeometryType neighborGT = neighborPointer->geometry().type();
                const LocalPosition& localPosNeighbor = ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

                // neighbor cell center in global coordinates
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

                // distance vector between barycenters
                FieldVector<Scalar,dimWorld>
                distVec = globalPosNeighbor - globalPos;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                FieldMatrix permeabilityJ = this->diffProblem.soil().K(globalPosNeighbor,*neighborPointer,localPosNeighbor);

                // compute vectorized permeabilities
                FieldVector<Scalar,dim> normalPermeabilityJ(0);
                permeabilityJ.umv(unitOuterNormal, normalPermeabilityJ);
                // compute permeability normal to intersection and take harmonic mean
                Scalar normalComponentPermeabilityI = normalPermeabilityI * unitOuterNormal;
                Scalar normalComponentPermeabilityJ = normalPermeabilityJ * unitOuterNormal;
                Scalar meanNormalPermeability = 2 * normalComponentPermeabilityI * normalComponentPermeabilityJ / (normalComponentPermeabilityI + normalComponentPermeabilityJ);
                // compute permeability tangential to intersection and take arithmetic mean
                FieldVector<Scalar,dim> normalComponentVector = unitOuterNormal;
                FieldVector<Scalar,dim> tangentialPermeabilityI = normalPermeabilityI - (normalComponentVector *= normalComponentPermeabilityI);
                normalComponentVector = unitOuterNormal;
                FieldVector<Scalar,dim> tangentialPermeabilityJ = normalPermeabilityJ - (normalComponentVector *= normalComponentPermeabilityJ);
                FieldVector<Scalar,dim> meanTangentialPermeability = (tangentialPermeabilityI += tangentialPermeabilityJ);
                meanTangentialPermeability *= 0.5;
                FieldVector<Scalar,dim> meanNormalPermeabilityVector = unitOuterNormal;
                // Build vectorized averaged permeability
                FieldVector<Scalar,dim> permeability = (meanTangentialPermeability += (meanNormalPermeabilityVector *= meanNormalPermeability));

                //compute total mobility
                Scalar lambdaJ, fractionalWJ, fractionalNWJ;

                // get mobilities and fractional flow factors
                lambdaJ = this->diffProblem.variables.mobilityWetting[globalIdxJ] + this->diffProblem.variables.mobilityNonWetting[globalIdxJ];
                fractionalWJ = this->diffProblem.variables.fracFlowFuncWetting[globalIdxJ];
                fractionalNWJ = this->diffProblem.variables.fracFlowFuncNonWetting[globalIdxJ];

                Scalar pcJ = this->diffProblem.variables.capillaryPressure[globalIdxJ];

                // update diagonal entry
                Scalar entry;
                if (first)
                {
                    Scalar lambda = (lambdaI + lambdaJ) / 2;
                    entry = fabs(lambda*faceVol*(permeability*distVec)/(dist*dist));
                    Scalar factor = (fractionalWI + fractionalWJ) * (densityW) / 2 + (fractionalNWI + fractionalNWJ) * (densityNW) / 2;
                    // calculate capillary pressure gradient
                    FieldVector<Scalar,dim> pCGradient = distVec;
                    pCGradient *= (pcI-pcJ)/(dist*dist);
                    Scalar lambdaNW = 0.5*(this->diffProblem.variables.mobilityNonWetting[globalIdxI]+this->diffProblem.variables.mobilityNonWetting[globalIdxJ]);
                    f_[globalIdxI] -= factor * lambda * faceVol * (permeability * gravity_) - lambdaNW*faceVol*(permeability*pCGradient);
                }
                else
                {


                    Scalar potentialW = (unitOuterNormal * distVec) * (this->diffProblem.variables.pressure[globalIdxI] - this->diffProblem.variables.pressure[globalIdxJ]) / (dist * dist);
                    Scalar potentialNW = (unitOuterNormal * distVec) * (this->diffProblem.variables.pressure[globalIdxI]+ pcI - this->diffProblem.variables.pressure[globalIdxJ]-pcJ) / (dist * dist);

                    potentialW += densityW * (unitOuterNormal * gravity_);
                    potentialNW += densityNW * (unitOuterNormal * gravity_);

                    Scalar lambdaW, lambdaNW;

                    if (potentialW >= 0.)
                    {
                        lambdaW = this->diffProblem.variables.mobilityWetting[globalIdxI];
                    }
                    else
                    {
                        lambdaW = this->diffProblem.variables.mobilityWetting[globalIdxJ];
                    }
                    if (potentialNW >= 0.)
                    {
                        lambdaNW = this->diffProblem.variables.mobilityNonWetting[globalIdxI];
                    }
                    else
                    {
                        lambdaNW = this->diffProblem.variables.mobilityNonWetting[globalIdxJ];
                    }

                    entry = (lambdaW + lambdaNW) * fabs(faceVol * (permeability * distVec) / (dist * dist));
                    Scalar rightEntry = (densityW * lambdaW + densityNW * lambdaNW) * faceVol * (permeability * gravity_);

                    // calculate capillary pressure gradient
                    FieldVector<Scalar,dim> pCGradient = distVec;
                    pCGradient *= (pcI-pcJ)/(dist*dist);

                    rightEntry += lambdaNW*faceVol*(permeability*pCGradient);
                    f_[globalIdxI] -= rightEntry;
                }
                // set diagonal entry
                A_[globalIdxI][globalIdxI] += entry;

                // set off-diagonal entry
                A_[globalIdxI][globalIdxJ] = -entry;
            }

            // boundary face
            else
            {
                // center of face in global coordinates
                const GlobalPosition& globalPosFace = isIt->geometry().global(faceLocal);

                //get boundary condition for boundary face center
                BoundaryConditions::Flags bctype = this->diffProblem.bctypePress(globalPosFace, *eIt, localPosFace);

                if (bctype == BoundaryConditions::dirichlet)
                {
                    FieldVector<Scalar,dimWorld> distVec(globalPosFace-globalPos);
                    Scalar dist = distVec.two_norm();

                    Scalar satBound = this->diffProblem.dirichletSat(globalPosFace, *eIt, localPosFace);
                    Scalar pressBound = this->diffProblem.dirichletPress(globalPosFace, *eIt, localPosFace);
                    Scalar pcBound = this->diffProblem.materialLaw().pC(satBound, globalPosFace, *eIt, localPosFace);
                    Scalar pcI = this->diffProblem.variables.capillaryPressure[globalIdxI];

                    if (first)
                    {
                        Scalar lambda = lambdaI;
                        A_[globalIdxI][globalIdxI] += lambda * faceVol * (normalPermeabilityI * distVec) / (dist * dist);
                        f_[globalIdxI] += lambda * faceVol * pressBound * fabs((normalPermeabilityI * distVec) / (dist * dist));
                        f_[globalIdxI] -= (fractionalWI * densityW + fractionalNWI * densityNW) * lambda*faceVol*(normalPermeabilityI*gravity_);
                        // calculate capillary pressure gradient
                        FieldVector<Scalar,dim> pCGradient = distVec;
                        pCGradient *= (pcI-pcBound)/(dist*dist);
                        Scalar lambdaNW = this->diffProblem.variables.mobilityNonWetting[globalIdxI];
                        f_[globalIdxI]-= lambdaNW*faceVol*(normalPermeabilityI*pCGradient);
                    }
                    else
                    {
                        Scalar residualSatW = this->diffProblem.soil().Sr_w(globalPos,*eIt, localPos);
                        Scalar residualSatNW = this->diffProblem.soil().Sr_n(globalPos,*eIt, localPos);
                        Scalar satBoundEff = (satBound-residualSatW)/(1-residualSatW-residualSatNW);

                        // velocities
                        Scalar potentialW = (unitOuterNormal * distVec) * (this->diffProblem.variables.pressure[globalIdxI] - pressBound) / (dist * dist);
                        Scalar potentialNW = (unitOuterNormal * distVec) * (this->diffProblem.variables.pressure[globalIdxI] + pcI - pressBound - pcBound) / (dist * dist);
                        potentialW += densityW * (unitOuterNormal * gravity_);
                        potentialNW += densityNW * (unitOuterNormal * gravity_);

                        Scalar lambdaW, lambdaNW;

                        if (potentialW >= 0.)
                        {
                            lambdaW = this->diffProblem.variables.mobilityWetting[globalIdxI];
                        }
                        else
                        {
                            lambdaW = this->diffProblem.materialLaw().mobW(satBound,globalPosFace, *eIt, localPosFace);
//                            lambdaW = satBoundEff / viscosityW;
                        }
                        if (potentialNW >= 0.)
                        {
                            lambdaNW = this->diffProblem.variables.mobilityNonWetting[globalIdxI];
                        }
                        else
                        {
                            lambdaNW = this->diffProblem.materialLaw().mobN((1-satBound),globalPosFace, *eIt, localPosFace);
                            //lambdaNW = (1 - satBoundEff) / viscosityNW;
                        }

                        Scalar entry = (lambdaW + lambdaNW) * fabs(faceVol * (normalPermeabilityI * distVec) / (dist * dist));
                        Scalar rightEntry = (densityW * lambdaW + densityNW *lambdaNW) * faceVol * (normalPermeabilityI * gravity_);

                        // calculate capillary pressure gradient
                        FieldVector<Scalar,dim> pCGradient = distVec;
                        pCGradient *= (pcI - pcBound)/(dist*dist);
//                        std::cout<<"pcGradBound ="<<pCGradient<<std::endl;

                        rightEntry += lambdaNW*faceVol*(normalPermeabilityI*pCGradient);

                        // set diagonal entry and right hand side entry
                        A_[globalIdxI][globalIdxI] += entry;
                        f_[globalIdxI] += entry * pressBound;
                        f_[globalIdxI] -= rightEntry;
                    }
                }
                else
                {
                    Scalar J = this->diffProblem.neumannPress(globalPosFace, *eIt, localPosFace);
                    f_[globalIdxI] -= faceVol*J;
                }
            }
        }
        // end all intersections
    } // end grid traversal
    return;
}

template<class Grid, class Scalar, class VC, class Problem> void FVDiffusion<Grid, Scalar, VC, Problem>::solve()
{
    MatrixAdapter<MatrixType,Vector,Vector> op(A_);
    InverseOperatorResult r;

    if (preconditionerName_ == "SeqILU0")
    {
        SeqILU0<MatrixType,Vector,Vector> preconditioner(A_, 1.0);
        if (solverName_ == "CG")
        {
            CGSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
            solver.apply((this->diffProblem.variables.pressure), f_, r);
        }
        else if (solverName_ == "BiCGSTAB")
        {
            BiCGSTABSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
            solver.apply((this->diffProblem.variables.pressure), f_, r);
        }
        else
        DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination "
                << preconditionerName_<< " and "<< solverName_ << ".");
    }
    else if (preconditionerName_ == "SeqPardiso")
    {
        SeqPardiso<MatrixType,Vector,Vector> preconditioner(A_);
        if (solverName_ == "Loop")
        {
            LoopSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
            solver.apply(this->diffProblem.variables.pressure, f_, r);
        }
        else
        DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination "
                << preconditionerName_<< " and "<< solverName_ << ".");
    }
    else
    DUNE_THROW(NotImplemented, "FVDiffusion :: solve : preconditioner "
            << preconditionerName_ << ".");
//        printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);
//        printvector(std::cout, f_, "right hand side", "row", 200, 1, 3);
//        printvector(std::cout, (this->diffProblem.variables.pressure), "pressure", "row", 200, 1, 3);
    return;
}

}
#endif
