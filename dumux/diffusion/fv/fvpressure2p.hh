// $Id$

#ifndef DUNE_FVPRESSURE2P_HH
#define DUNE_FVPRESSURE2P_HH

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
#include "dumux/diffusion/diffusion.hh"
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
template<class GridView, class Scalar, class VC,
        class Problem = DiffusionProblem<GridView, Scalar, VC> > class FVPressure2P: public Diffusion<
        GridView, Scalar, VC, Problem>
{
    enum
    {
        dim = GridView::dimension
    };
    enum
    {
        dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = 0, pn = 1, pglobal = 2
    };

typedef    typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
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

    void pressure(bool first = true, const Scalar t=0)
    {
        initializeMaterialLaws();
        assemble(first, t);
        solve();
        if (first = true)
        {
            assemble(false, t);
            solve();
        }
        return;
    }

    void initializeMatrix();

    //constitutive functions are initialized and stored in the variables object
    void initializeMaterialLaws();

    FVPressure2P(GridView& gridView, Problem& problem, std::string pressType) :
    Diffusion<GridView, Scalar, VC, Problem>(gridView, problem),
    A_(problem.variables().gridSizeDiffusion(), problem.variables().gridSizeDiffusion(),
            (2*dim+1) * problem.variables().gridSizeDiffusion(), BCRSMatrix<MB>::random),
    f_(problem.variables().gridSizeDiffusion()), solverName_("BiCGSTAB"),
    preconditionerName_("SeqILU0"), gravity(problem.gravity()),
    pressureType((pressType == "pw") ? 0 : ((pressType == "pn") ? 1 : ((pressType == "pglobal") ? 2 : 999)))
    {
        if (pressureType == 999)
        {
            DUNE_THROW(NotImplemented, "Pressure type not supported!");
        }
        initializeMatrix();
    }

    FVPressure2P(GridView& gridView, Problem& problem, std::string pressType, std::string solverName,
            std::string preconditionerName) :
    Diffusion<GridView, Scalar, VC, Problem>(gridView, problem),
    A_(problem.variables().gridSizeDiffusion(), problem.variables().gridSizeDiffusion(),
            (2*dim+1)*problem.variables().gridSizeDiffusion(), BCRSMatrix<MB>::random),
    f_(problem.variables().gridSizeDiffusion()), solverName_(solverName),
    preconditionerName_(preconditionerName), gravity(problem.gravity()),
    pressureType((pressType == "pw") ? 0 : ((pressType == "pn") ? 1 : ((pressType == "pglobal") ? 2 : 999)))
    {
        if (pressureType == 999)
        {
            DUNE_THROW(NotImplemented, "Pressure type not supported!");
        }
        initializeMatrix();
    }

private:
    MatrixType A_;
    BlockVector< FieldVector<Scalar,1> > f_;
    std::string solverName_;
    std::string preconditionerName_;
protected:
    const FieldVector<Scalar,dimWorld>& gravity; // [m/s^2]
    const int pressureType;
};

template<class GridView, class Scalar, class VC, class Problem> void FVPressure2P<GridView, Scalar, VC, Problem>::initializeMatrix()
{
    // determine matrix row sizes
    ElementIterator eItEnd = this->gridView.template end<0>();
    for (ElementIterator eIt = this->gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = this->diffProblem.variables().indexDiffusion(*eIt);

        // initialize row size
        int rowSize = 1;

        // run through all intersections with neighbors
        IntersectionIterator
        isItEnd = this->gridView.template iend(*eIt);
        for (IntersectionIterator
                isIt = this->gridView.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        if (isIt->neighbor())
        rowSize++;
        A_.setrowsize(globalIdxI, rowSize);
    }
    A_.endrowsizes();

    // determine position of matrix entries
    for (ElementIterator eIt = this->gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = this->diffProblem.variables().indexDiffusion(*eIt);

        // add diagonal index
        A_.addindex(globalIdxI, globalIdxI);

        // run through all intersections with neighbors
        IntersectionIterator
        isItEnd = this->gridView.template iend(*eIt);
        for (IntersectionIterator
                isIt = this->gridView.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        if (isIt->neighbor())
        {
            // access neighbor
            ElementPointer outside = isIt->outside();
            int globalIdxJ = this->diffProblem.variables().indexDiffusion(*outside);

            // add off diagonal index
            A_.addindex(globalIdxI, globalIdxJ);
        }
    }
    A_.endindices();

    return;
}

template<class GridView, class Scalar, class VC, class Problem> void FVPressure2P<GridView, Scalar, VC, Problem>::assemble(bool first, const Scalar t=0)
{
    // initialization: set matrix A_ to zero
    A_ = 0;
    f_=0;

    Scalar densityW = this->diffProblem.wettingPhase().density();
    Scalar densityNW = this->diffProblem.nonWettingPhase().density();

    ElementIterator eItEnd = this->gridView.template end<0>();
    for (ElementIterator eIt = this->gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell geometry type
        GeometryType gt = eIt->geometry().type();

        // cell center in reference element
        const LocalPosition& localPos = ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        const GlobalPosition& globalPos = eIt->geometry().global(localPos);

        // cell index
        int globalIdxI = this->diffProblem.variables().indexDiffusion(*eIt);

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().integrationElement(localPos)*ReferenceElements<Scalar,dim>::general(gt).volume();

        // set right side to zero
        f_[globalIdxI] = volume*this->diffProblem.sourcePress(globalPos, *eIt, localPos);

        // get absolute permeability
        FieldMatrix permeabilityI(this->diffProblem.soil().K(globalPos, *eIt, localPos));

        //compute total mobility
        Scalar lambdaI, fractionalWI, fractionalNWI;

        // get mobilities and fractional flow factors
        lambdaI = this->diffProblem.variables().mobilityWetting()[globalIdxI] + this->diffProblem.variables().mobilityNonWetting()[globalIdxI];
        fractionalWI = this->diffProblem.variables().fracFlowFuncWetting()[globalIdxI];
        fractionalNWI = this->diffProblem.variables().fracFlowFuncNonWetting()[globalIdxI];

        Scalar pcI = this->diffProblem.variables().capillaryPressure()[globalIdxI];

        IntersectionIterator
        isItEnd = this->gridView.template iend(*eIt);
        for (IntersectionIterator
                isIt = this->gridView.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {

            // get geometry type of face
            GeometryType faceGT = isIt->geometryInInside().type();

            // center in face's reference element
            const FieldVector<Scalar,dim-1>&
            faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

            int indexInInside = isIt->indexInInside();

            // center of face inside volume reference element
            const LocalPosition& localPosFace = ReferenceElements<Scalar,dim>::general(faceGT).position(indexInInside,1);

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
                int globalIdxJ = this->diffProblem.variables().indexDiffusion(*neighborPointer);
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
                lambdaJ = this->diffProblem.variables().mobilityWetting()[globalIdxJ] + this->diffProblem.variables().mobilityNonWetting()[globalIdxJ];
                fractionalWJ = this->diffProblem.variables().fracFlowFuncWetting()[globalIdxJ];
                fractionalNWJ = this->diffProblem.variables().fracFlowFuncNonWetting()[globalIdxJ];

                Scalar pcJ = this->diffProblem.variables().capillaryPressure()[globalIdxJ];

                // update diagonal entry
                Scalar entry;
                if (first)
                {
                    Scalar lambda = (lambdaI + lambdaJ) / 2;
                    entry = fabs(lambda*faceVol*(permeability*distVec)/(dist*dist));
                    Scalar factor = (fractionalWI + fractionalWJ) * (densityW) / 2 + (fractionalNWI + fractionalNWJ) * (densityNW) / 2;

                    f_[globalIdxI] -= factor * lambda * faceVol * (permeability * gravity);

                    if (pressureType == pw)
                    {
                        // calculate capillary pressure gradient
                        FieldVector<Scalar,dim> pCGradient = distVec;
                        pCGradient *= (pcI-pcJ)/(dist*dist);
                        Scalar lambdaNW = 0.5*(this->diffProblem.variables().mobilityNonWetting()[globalIdxI]+this->diffProblem.variables().mobilityNonWetting()[globalIdxJ]);

                        f_[globalIdxI] -= lambdaNW*faceVol*(permeability*pCGradient);
                    }
                    if (pressureType == pn)
                    {
                        // calculate capillary pressure gradient
                        FieldVector<Scalar,dim> pCGradient = distVec;
                        pCGradient *= (pcI-pcJ)/(dist*dist);
                        Scalar lambdaW = 0.5*(this->diffProblem.variables().mobilityWetting()[globalIdxI]+this->diffProblem.variables().mobilityWetting()[globalIdxJ]);

                        f_[globalIdxI] += lambdaW*faceVol*(permeability*pCGradient);
                    }
                }
                else
                {
                    Scalar potentialW = 0;
                    Scalar potentialNW = 0;

                    if (pressureType == pw)
                    {
                        potentialW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - this->diffProblem.variables().pressure()[globalIdxJ]) / (dist * dist);
                        potentialNW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI]+ pcI - this->diffProblem.variables().pressure()[globalIdxJ]-pcJ) / (dist * dist);
                    }
                    if (pressureType == pn)
                    {
                        potentialW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - pcI - this->diffProblem.variables().pressure()[globalIdxJ] + pcJ) / (dist * dist);
                        potentialNW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - this->diffProblem.variables().pressure()[globalIdxJ]) / (dist * dist);
                    }
                    if (pressureType == pglobal)
                    {
                        potentialW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - this->diffProblem.variables().pressure()[globalIdxJ] - 0.5 * (fractionalNWI+fractionalNWJ)*(pcI - pcJ)) / (dist * dist);
                        potentialNW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - this->diffProblem.variables().pressure()[globalIdxJ] + 0.5 * (fractionalWI+fractionalWJ)*(pcI - pcJ)) / (dist * dist);
                    }

                    potentialW += densityW * (unitOuterNormal * gravity);
                    potentialNW += densityNW * (unitOuterNormal * gravity);

                    this->diffProblem.variables().potentialWetting()[globalIdxI][indexInInside] = potentialW;
                    this->diffProblem.variables().potentialNonWetting()[globalIdxI][indexInInside] = potentialNW;

                    Scalar lambdaW, lambdaNW;

                    if (potentialW >= 0.)
                    {
                        lambdaW = this->diffProblem.variables().mobilityWetting()[globalIdxI];
                    }
                    else
                    {
                        lambdaW = this->diffProblem.variables().mobilityWetting()[globalIdxJ];
                    }
                    if (potentialNW >= 0.)
                    {
                        lambdaNW = this->diffProblem.variables().mobilityNonWetting()[globalIdxI];
                    }
                    else
                    {
                        lambdaNW = this->diffProblem.variables().mobilityNonWetting()[globalIdxJ];
                    }

                    entry = (lambdaW + lambdaNW) * fabs(faceVol * (permeability * distVec) / (dist * dist));
                    Scalar rightEntry = (densityW * lambdaW + densityNW * lambdaNW) * faceVol * (permeability * gravity);

                    if (pressureType == pw)
                    {
                        // calculate capillary pressure gradient
                        FieldVector<Scalar,dim> pCGradient = distVec;
                        pCGradient *= (pcI-pcJ)/(dist*dist);

                        rightEntry += (lambdaW + lambdaNW) * 0.5 * (fractionalNWI + fractionalNWJ)*faceVol*(permeability*pCGradient);
                    }
                    if (pressureType == pn)
                    {
                        // calculate capillary pressure gradient
                        FieldVector<Scalar,dim> pCGradient = distVec;
                        pCGradient *= (pcI-pcJ)/(dist*dist);

                        rightEntry -= (lambdaW + lambdaNW) * 0.5*(fractionalWI + fractionalWJ) *faceVol*(permeability*pCGradient);
                    }

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
                BoundaryConditions::Flags bcTypeSat = this->diffProblem.bctypeSat(globalPosFace, *eIt, localPosFace);

                if (bctype == BoundaryConditions::dirichlet)
                {
                    FieldVector<Scalar,dimWorld> distVec(globalPosFace-globalPos);
                    Scalar dist = distVec.two_norm();

                    Scalar satBound;
                    if (bcTypeSat == BoundaryConditions::dirichlet)
                    {
                        satBound = this->diffProblem.dirichletSat(globalPosFace, *eIt, localPosFace);
                    }
                    else
                    {
                        satBound = this->diffProblem.variables().saturation()[globalIdxI];
                    }

                    Scalar pressBound = this->diffProblem.dirichletPress(globalPosFace, *eIt, localPosFace);
                    Scalar pcBound = this->diffProblem.materialLaw().pC(satBound, globalPosFace, *eIt, localPosFace);
                    Scalar pcI = this->diffProblem.variables().capillaryPressure()[globalIdxI];

                    if (first)
                    {
                        Scalar lambda = lambdaI;
                        A_[globalIdxI][globalIdxI] += lambda * faceVol * (normalPermeabilityI * distVec) / (dist * dist);
                        f_[globalIdxI] += lambda * faceVol * pressBound * fabs((normalPermeabilityI * distVec) / (dist * dist));
                        f_[globalIdxI] -= (fractionalWI * densityW + fractionalNWI * densityNW) * lambda*faceVol*(normalPermeabilityI*gravity);

                        if (pressureType == pw)
                        {
                            // calculate capillary pressure gradient
                            FieldVector<Scalar,dim> pCGradient = distVec;
                            pCGradient *= (pcI-pcBound)/(dist*dist);
                            Scalar lambdaNW = this->diffProblem.variables().mobilityNonWetting()[globalIdxI];

                            f_[globalIdxI] -= lambdaNW*faceVol*(normalPermeabilityI*pCGradient);
                        }
                        if (pressureType == pn)
                        {
                            // calculate capillary pressure gradient
                            FieldVector<Scalar,dim> pCGradient = distVec;
                            pCGradient *= (pcI-pcBound)/(dist*dist);
                            Scalar lambdaW = this->diffProblem.variables().mobilityWetting()[globalIdxI];

                            f_[globalIdxI] += lambdaW*faceVol*(normalPermeabilityI*pCGradient);
                        }
                    }
                    else
                    {
                        Scalar potentialW = 0;
                        Scalar potentialNW = 0;

                        // potential gradient
                        if (pressureType == pw)
                        {
                            potentialW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - pressBound) / (dist * dist);
                            potentialNW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] + pcI - pressBound - pcBound) / (dist * dist);
                        }
                        if (pressureType == pn)
                        {
                            potentialW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - pcI - pressBound + pcBound) / (dist * dist);
                            potentialNW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - pressBound) / (dist * dist);
                        }
                        if (pressureType == pglobal)
                        {
                            potentialW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - pressBound - fractionalNWI * (pcI - pcBound)) / (dist * dist);
                            potentialNW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - pressBound + fractionalWI * (pcI - pcBound)) / (dist * dist);
                        }

                        potentialW += densityW * (unitOuterNormal * gravity);
                        potentialNW += densityNW * (unitOuterNormal * gravity);

                        this->diffProblem.variables().potentialWetting()[globalIdxI][indexInInside] = potentialW;
                        this->diffProblem.variables().potentialNonWetting()[globalIdxI][indexInInside] = potentialNW;

                        Scalar lambdaW, lambdaNW;

                        if (potentialW >= 0.)
                        {
                            lambdaW = this->diffProblem.variables().mobilityWetting()[globalIdxI];
                        }
                        else
                        {
                            lambdaW = this->diffProblem.materialLaw().mobW(satBound,globalPosFace, *eIt, localPosFace);
                        }
                        if (potentialNW >= 0.)
                        {
                            lambdaNW = this->diffProblem.variables().mobilityNonWetting()[globalIdxI];
                        }
                        else
                        {
                            lambdaNW = this->diffProblem.materialLaw().mobN((1-satBound),globalPosFace, *eIt, localPosFace);
                        }

                        Scalar entry = (lambdaW + lambdaNW) * fabs(faceVol * (normalPermeabilityI * distVec) / (dist * dist));
                        Scalar rightEntry = (densityW * lambdaW + densityNW *lambdaNW) * faceVol * (normalPermeabilityI * gravity);

                        if (pressureType == pw)
                        {
                            // calculate capillary pressure gradient
                            FieldVector<Scalar,dim> pCGradient = distVec;
                            pCGradient *= (pcI - pcBound)/(dist*dist);

                            rightEntry += lambdaNW * faceVol*(normalPermeabilityI*pCGradient);
                        }
                        if (pressureType == pn)
                        {
                            // calculate capillary pressure gradient
                            FieldVector<Scalar,dim> pCGradient = distVec;
                            pCGradient *= (pcI - pcBound)/(dist*dist);

                            rightEntry -= lambdaW *faceVol*(normalPermeabilityI*pCGradient);

                        }

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

template<class GridView, class Scalar, class VC, class Problem> void FVPressure2P<GridView, Scalar, VC, Problem>::solve()
{
    MatrixAdapter<MatrixType,Vector,Vector> op(A_);
    InverseOperatorResult r;

    if (preconditionerName_ == "SeqILU0")
    {
        SeqILU0<MatrixType,Vector,Vector> preconditioner(A_, 1.0);
        if (solverName_ == "CG")
        {
            CGSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
            solver.apply((this->diffProblem.variables().pressure()), f_, r);
        }
        else if (solverName_ == "BiCGSTAB")
        {
            BiCGSTABSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
            solver.apply((this->diffProblem.variables().pressure()), f_, r);
        }
        else
        DUNE_THROW(NotImplemented, "FVPressure2P :: solve : combination "
                << preconditionerName_<< " and "<< solverName_ << ".");
    }
    else if (preconditionerName_ == "SeqPardiso")
    {
        SeqPardiso<MatrixType,Vector,Vector> preconditioner(A_);
        if (solverName_ == "Loop")
        {
            LoopSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
            solver.apply(this->diffProblem.variables().pressure(), f_, r);
        }
        else
        DUNE_THROW(NotImplemented, "FVPressure2P :: solve : combination "
                << preconditionerName_<< " and "<< solverName_ << ".");
    }
    else
    DUNE_THROW(NotImplemented, "FVPressure2P :: solve : preconditioner "
            << preconditionerName_ << ".");
//                printmatrix(std::cout, A_, "global stiffness matrix", "row", 11, 3);
//                printvector(std::cout, f_, "right hand side", "row", 200, 1, 3);
//                printvector(std::cout, (this->diffProblem.variables().pressure()), "pressure", "row", 200, 1, 3);
    return;
}
//constitutive functions are updated once if new saturations are calculated and stored in the variables object
template<class GridView, class Scalar, class VC, class Problem> void FVPressure2P<GridView, Scalar, VC, Problem>::initializeMaterialLaws()
{
    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = this->gridView.template end<0>();
    for (ElementIterator eIt = this->gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // get geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // get cell center in reference element
        const LocalPosition
        &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().global(localPos);

        int globalIdx = this->diffProblem.variables().indexDiffusion(*eIt);

        Scalar sat = this->diffProblem.variables().saturation()[globalIdx];

        std::vector<Scalar> mobilities = this->diffProblem.materialLaw().mob(sat, globalPos, *eIt, localPos);

        // initialize mobilities
        this->diffProblem.variables().mobilityWetting()[globalIdx]= mobilities[0];
        this->diffProblem.variables().mobilityNonWetting()[globalIdx]= mobilities[1];
        this->diffProblem.variables().capillaryPressure()[globalIdx]= this->diffProblem.materialLaw().pC(sat, globalPos, *eIt, localPos);
        this->diffProblem.variables().fracFlowFuncWetting()[globalIdx]= mobilities[0]/(mobilities[0]+mobilities[1]);
        this->diffProblem.variables().fracFlowFuncNonWetting()[globalIdx]= mobilities[1]/(mobilities[0]+mobilities[1]);
    }
    return;
}

}
#endif
