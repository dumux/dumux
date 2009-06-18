// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2007-2009 by Jochen Fritz                                 *
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
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
 * \f[\text{div}\, \boldsymbol{v}_{total} = q.\f]
 * The definition of the total velocity \f$\boldsymbol{v}_total\f$ depends on the kind of pressure chosen. This could be a wetting (w) phase pressure leading to
 * \f[ - \text{div}\,  \left[\lambda \boldsymbol{K} \left(\text{grad}\, p_w + f_n \text{grad}\, p_c + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q, \f]
 * a non-wetting (n) phase pressure yielding
 * \f[ - \text{div}\,  \left[\lambda \boldsymbol{K}  \left(\text{grad}\, p_n - f_w \text{grad}\, p_c + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q, \f]
 * or a global pressure leading to
 * \f[ - \text{div}\, \left[\lambda \boldsymbol{K} \left(\text{grad}\, p_{global} + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right)\right] = q.\f]
 *  Here, \f$p\f$ denotes a pressure, \f$\boldsymbol{K}\f$ the absolute permeability, \f$\lambda\f$ the total mobility, possibly depending on the
 * saturation,\f$f\f$ the fractional flow function of a phase, \f$\rho\f$ a phase density, \f$g\f$ the gravity constant and \f$q\f$ the source term.
 * For all cases, \f$p = p_D\f$ on \f$\Gamma_{Neumann}\f$, and \f$\boldsymbol{v}_{total}  = q_N\f$
 * on \f$\Gamma_{Dirichlet}\f$.
 *
 * Template parameters are:
 *
 - GridView      a DUNE gridview type
 - Scalar        type used for scalar quantities
 - VC            type of a class containing different variables of the model
 - Problem       class defining the physical problem
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
        pw = 0, pn = 1, pglobal = 2, Sw = 0, Sn = 1, other = 999
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
    typedef BCRSMatrix<MB> Matrix;
    typedef BlockVector<FieldVector<Scalar, 1> > Vector;

    //initializes the matrix to store the system of equations
    void initializeMatrix();

    //constitutive functions are initialized and stored in the variables object
    void initializeMaterialLaws();

    //function which assembles the system of equations to be solved
    void assemble(bool first, const Scalar t);

    //solves the system of equations to get the spatial distribution of the pressure
    void solve();

public:
    //! Calculate the pressure.
    /*!
     *  \param t time
     *
     *  Calculates the pressure \f$p\f$ as solution of the boundary value problem
     *  \f[  \text{div}\, \boldsymbol{v} = q, \f]
     *  subject to appropriate boundary conditions.
     */
    void pressure(bool first = true, const Scalar t=0, bool solveTwice = true)
    {
        if (first)
        {
            initializeMaterialLaws();
        }
        assemble(first, t);
        solve();
        if (first && solveTwice)
        {
            assemble(false, t);
            solve();
        }
        return;
    }

    //! Constructs a FVPressure2P object
    /**
     * \param gridView gridView object of type GridView
     * \param problem a problem class object
     * \param pressType a string giving the type of pressure used (could be: pw, pn, pglobal)
     * \param satType a string giving the type of saturation used (could be: Sw, Sn)
     */
    FVPressure2P(GridView& gridView, Problem& problem, std::string pressType, std::string satType) :
    Diffusion<GridView, Scalar, VC, Problem>(gridView, problem),
    A_(problem.variables().gridSizeDiffusion(), problem.variables().gridSizeDiffusion(),
            (2*dim+1) * problem.variables().gridSizeDiffusion(), BCRSMatrix<MB>::random),
    f_(problem.variables().gridSizeDiffusion()), solverName_("BiCGSTAB"),
    preconditionerName_("SeqILU0"), gravity(problem.gravity()),
    pressureType((pressType == "pw") ? pw : ((pressType == "pn") ? pn : ((pressType == "pglobal") ? pglobal : other))),
    saturationType((satType == "Sw") ? Sw : ((satType == "Sn") ? Sn : other))
    {
        if (pressureType == other)
        {
            DUNE_THROW(NotImplemented, "Pressure type not supported!");
        }
        if (saturationType == other)
        {
            DUNE_THROW(NotImplemented, "Saturation type not supported!");
        }
        initializeMatrix();
    }

    //! Constructs a FVPressure2P object
    /**
     * \param gridView gridView object of type GridView
     * \param problem a problem class object
     * \param pressType a string giving the type of pressure used (could be: pw, pn, pglobal)
     * \param satType a string giving the type of saturation used (could be: Sw, Sn)
     * \param solverName a string giving the type of solver used (could be: CG, BiCGSTAB, Loop)
     * \param preconditionerName a string giving the type of the matrix preconditioner used (could be: SeqILU0, SeqPardiso)
     */
    FVPressure2P(GridView& gridView, Problem& problem, std::string pressType, std::string satType, std::string solverName,
            std::string preconditionerName)
    : Diffusion<GridView, Scalar, VC, Problem>(gridView, problem),
    A_(problem.variables().gridSizeDiffusion(), problem.variables().gridSizeDiffusion(),
    		(2*dim+1) * problem.variables().gridSizeDiffusion(), BCRSMatrix<MB>::random),
    f_(problem.variables().gridSizeDiffusion()), solverName_(solverName),
    preconditionerName_(preconditionerName), gravity(problem.gravity()),
    pressureType((pressType == "pw") ? pw : ((pressType == "pn") ? pn : ((pressType == "pglobal") ? pglobal : other))),
    saturationType((satType == "Sw") ? Sw : ((satType == "Sn") ? Sn : other))
    {
    	if (pressureType == other)
    	{
    		DUNE_THROW(NotImplemented, "Pressure type not supported!");
    	}
    	if (saturationType == other)
    	{
    		DUNE_THROW(NotImplemented, "Saturation type not supported!");
    	}
    	initializeMatrix();
    }

private:
    Matrix A_;
    BlockVector< FieldVector<Scalar,1> > f_;
    std::string solverName_;
    std::string preconditionerName_;
protected:
    const FieldVector<Scalar,dimWorld>& gravity; //!< vector including the gravity constant
    const int pressureType; //!< gives kind of pressure used (\f$ 0 = p_w\f$, \f$ 1 = p_n\f$, \f$ 2 = p_{global}\f$)
    const int saturationType; //!< gives kind of saturation used (\f$ 0 = S_w\f$, \f$ 1 = S_n\f$)
};

//initializes the matrix to store the system of equations
template<class GridView, class Scalar, class VC, class Problem>
void FVPressure2P<GridView, Scalar, VC, Problem>::initializeMatrix()
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

//function which assembles the system of equations to be solved
template<class GridView, class Scalar, class VC, class Problem>
void FVPressure2P<GridView, Scalar, VC, Problem>::assemble(bool first, const Scalar t=0)
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

        // get mobilities and fractional flow factors
        Scalar lambdaWI = this->diffProblem.variables().mobilityWetting()[globalIdxI];
        Scalar lambdaNWI = this->diffProblem.variables().mobilityNonWetting()[globalIdxI];
        Scalar lambdaI = lambdaWI+ lambdaNWI;
        Scalar fractionalWI = this->diffProblem.variables().fracFlowFuncWetting()[globalIdxI];
        Scalar fractionalNWI = this->diffProblem.variables().fracFlowFuncNonWetting()[globalIdxI];
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

                // get mobilities and fractional flow factors
                Scalar lambdaWJ = this->diffProblem.variables().mobilityWetting()[globalIdxJ];
                Scalar lambdaNWJ = this->diffProblem.variables().mobilityNonWetting()[globalIdxJ];
                Scalar lambdaJ = lambdaWJ + lambdaNWJ;
                Scalar fractionalWJ = this->diffProblem.variables().fracFlowFuncWetting()[globalIdxJ];
                Scalar fractionalNWJ = this->diffProblem.variables().fracFlowFuncNonWetting()[globalIdxJ];

                Scalar pcJ = this->diffProblem.variables().capillaryPressure()[globalIdxJ];

                // update diagonal entry
                Scalar entry;
                //if we are at the very first iteration a guess of the pressure field has to be calculated to to get phase potentials for the following calculation steps
                if (first)
                {
                    //use central weights at the very first iteration
                    Scalar lambda = (lambdaI + lambdaJ) / 2;

                    //calculate current matrix entry
                    entry = fabs(lambda*faceVol*(permeability*distVec)/(dist*dist));

                    //calculate right hand side
                    Scalar factor = (fractionalWI + fractionalWJ) * (densityW) / 2 + (fractionalNWI + fractionalNWJ) * (densityNW) / 2;
                    f_[globalIdxI] -= factor * lambda * faceVol * (permeability * gravity);

                    if (pressureType == pw)
                    {
                        // calculate capillary pressure gradient
                        FieldVector<Scalar,dim> pCGradient = distVec;
                        pCGradient *= (pcI-pcJ)/(dist*dist);

                        //add capillary pressure term to right hand side
                        f_[globalIdxI] -= 0.5 * (lambdaNWI + lambdaNWJ) *faceVol*(permeability*pCGradient);
                    }
                    if (pressureType == pn)
                    {
                        // calculate capillary pressure gradient
                        FieldVector<Scalar,dim> pCGradient = distVec;
                        pCGradient *= (pcI-pcJ)/(dist*dist);

                        //add capillary pressure term to right hand side
                        f_[globalIdxI] += 0.5* (lambdaWI + lambdaWJ) *faceVol*(permeability*pCGradient);
                    }
                }
                //if a pressure field is knwon from previous calculations phase potentials can be calculated to determine upwind directions
                else
                {
                    //calculate potential gradients
                    Scalar potentialW = 0;
                    Scalar potentialNW = 0;

                    if (pressureType == pw)
                    {
                        potentialW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - this->diffProblem.variables().pressure()[globalIdxJ]) / (dist * dist);
                        potentialNW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - this->diffProblem.variables().pressure()[globalIdxJ]+ pcI - pcJ) / (dist * dist);
                    }
                    if (pressureType == pn)
                    {
                        potentialW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - this->diffProblem.variables().pressure()[globalIdxJ] - pcI + pcJ) / (dist * dist);
                        potentialNW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - this->diffProblem.variables().pressure()[globalIdxJ]) / (dist * dist);
                    }
                    if (pressureType == pglobal)
                    {
                        potentialW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - this->diffProblem.variables().pressure()[globalIdxJ] - 0.5 * (fractionalNWI+fractionalNWJ)*(pcI - pcJ)) / (dist * dist);
                        potentialNW = (unitOuterNormal * distVec) * (this->diffProblem.variables().pressure()[globalIdxI] - this->diffProblem.variables().pressure()[globalIdxJ] + 0.5 * (fractionalWI+fractionalWJ)*(pcI - pcJ)) / (dist * dist);
                    }

                    potentialW += densityW * (unitOuterNormal * gravity);
                    potentialNW += densityNW * (unitOuterNormal * gravity);

                    //store potentials for further calculations (velocity, saturation, ...)
                    this->diffProblem.variables().potentialWetting()[globalIdxI][indexInInside] = potentialW;
                    this->diffProblem.variables().potentialNonWetting()[globalIdxI][indexInInside] = potentialNW;

                    //do the upwinding of the mobility depending on the phase potentials
                    Scalar lambdaW = (potentialW >= 0.) ? lambdaWI : lambdaWJ;
                    Scalar lambdaNW = (potentialNW >= 0.) ? lambdaNWI : lambdaNWJ;

                    //calculate current matrix entry
                    entry = (lambdaW + lambdaNW) * fabs(faceVol * (permeability * distVec) / (dist * dist));

                    //calculate right hand side
                    Scalar rightEntry = (densityW * lambdaW + densityNW * lambdaNW) * faceVol * (permeability * gravity);

                    if (pressureType == pw)
                    {
                        // calculate capillary pressure gradient
                        FieldVector<Scalar,dim> pCGradient = distVec;
                        pCGradient *= (pcI-pcJ)/(dist*dist);

                        //add capillary pressure term to right hand side
                        rightEntry += 0.5 * (lambdaNWI + lambdaNWJ)*faceVol*(permeability*pCGradient);
                    }
                    if (pressureType == pn)
                    {
                        // calculate capillary pressure gradient
                        FieldVector<Scalar,dim> pCGradient = distVec;
                        pCGradient *= (pcI-pcJ)/(dist*dist);

                        //add capillary pressure term to right hand side
                        rightEntry -= 0.5*(lambdaWI + lambdaWJ) *faceVol*(permeability*pCGradient);
                    }

                    //set right hand side
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

                    //determine saturation at the boundary -> if no saturation is known directly at the boundary use the cell saturation
                    Scalar satBound;
                    if (bcTypeSat == BoundaryConditions::dirichlet)
                    {
                        satBound = this->diffProblem.dirichletSat(globalPosFace, *eIt, localPosFace);
                    }
                    else
                    {
                        satBound = this->diffProblem.variables().saturation()[globalIdxI];
                    }

                    //get dirichlet pressure boundary condition
                    Scalar pressBound = this->diffProblem.dirichletPress(globalPosFace, *eIt, localPosFace);
                    Scalar pcBound = 0;
                    Scalar lambdaWBound = 0;
                    Scalar lambdaNWBound = 0;

                    //calculate consitutive relations depending on the kind of saturation used
                    if (saturationType == Sw)
                    {
                        pcBound = this->diffProblem.materialLaw().pC(satBound, globalPosFace, *eIt, localPosFace);
                        std::vector<Scalar> mobilityBound = this->diffProblem.materialLaw().mob(satBound,globalPosFace, *eIt, localPosFace);
                        lambdaWBound = mobilityBound[0];
                        lambdaNWBound = mobilityBound[1];
                    }
                    if (saturationType == Sn)
                    {
                        pcBound = this->diffProblem.materialLaw().pC(1-satBound, globalPosFace, *eIt, localPosFace);
                        std::vector<Scalar> mobilityBound = this->diffProblem.materialLaw().mob(1-satBound,globalPosFace, *eIt, localPosFace);
                        lambdaWBound = mobilityBound[0];
                        lambdaNWBound = mobilityBound[1];
                    }
                    Scalar pcI = this->diffProblem.variables().capillaryPressure()[globalIdxI];

                    //if we are at the very first iteration a guess of the pressure field has to be calculated to to get phase potentials for the following calculation steps
                    if (first)
                    {
                        Scalar lambda = lambdaI;
                        //calculate current matrix entry
                        A_[globalIdxI][globalIdxI] += lambda * faceVol * (normalPermeabilityI * distVec) / (dist * dist);

                        //add dirichlet condition and gravity term to right hand side
                        f_[globalIdxI] += lambda * faceVol * pressBound * fabs((normalPermeabilityI * distVec) / (dist * dist));
                        f_[globalIdxI] -= (fractionalWI * densityW + fractionalNWI * densityNW) * lambda*faceVol*(normalPermeabilityI*gravity);

                        if (pressureType == pw)
                        {
                            // calculate capillary pressure gradient
                            FieldVector<Scalar,dim> pCGradient = distVec;
                            pCGradient *= (pcI-pcBound)/(dist*dist);

                            //add capillary pressure term to right hand side
                            f_[globalIdxI] -= 0.5*(lambdaNWI + lambdaNWBound) *faceVol*(normalPermeabilityI*pCGradient);
                        }
                        if (pressureType == pn)
                        {
                            // calculate capillary pressure gradient
                            FieldVector<Scalar,dim> pCGradient = distVec;
                            pCGradient *= (pcI-pcBound)/(dist*dist);

                            //add capillary pressure term to right hand side
                            f_[globalIdxI] += 0.5*(lambdaWI + lambdaWBound)*faceVol*(normalPermeabilityI*pCGradient);
                        }
                    }
                    //if a pressure field is knwon from previous calculations phase potentials can be calculated to determine upwind directions
                    else
                    {
                        Scalar potentialW = 0;
                        Scalar potentialNW = 0;

                        //calculate potential gradient
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

                        //store potential gradients for further calculations
                        this->diffProblem.variables().potentialWetting()[globalIdxI][indexInInside] = potentialW;
                        this->diffProblem.variables().potentialNonWetting()[globalIdxI][indexInInside] = potentialNW;

                        //do the upwinding of the mobility depending on the phase potentials
                        Scalar lambdaW = (potentialW >= 0.) ? lambdaWI : lambdaWBound;
                        Scalar lambdaNW = (potentialNW >= 0.) ? lambdaNWI : lambdaNWBound;

                        //calculate current matrix entry
                        Scalar entry = (lambdaW + lambdaNW) * fabs(faceVol * (normalPermeabilityI * distVec) / (dist * dist));

                        //calculate right hand side
                        Scalar rightEntry = (densityW * lambdaW + densityNW *lambdaNW) * faceVol * (normalPermeabilityI * gravity);

                        if (pressureType == pw)
                        {
                            // calculate capillary pressure gradient
                            FieldVector<Scalar,dim> pCGradient = distVec;
                            pCGradient *= (pcI - pcBound)/(dist*dist);

                            //add capillary pressure term to right hand side
                            rightEntry += 0.5 * (lambdaNWI + lambdaNWBound) * faceVol*(normalPermeabilityI*pCGradient);
                        }
                        if (pressureType == pn)
                        {
                            // calculate capillary pressure gradient
                            FieldVector<Scalar,dim> pCGradient = distVec;
                            pCGradient *= (pcI - pcBound)/(dist*dist);

                            //add capillary pressure term to right hand side
                            rightEntry -= 0.5 * (lambdaWI + lambdaWBound) *faceVol*(normalPermeabilityI*pCGradient);

                        }

                        // set diagonal entry and right hand side entry
                        A_[globalIdxI][globalIdxI] += entry;
                        f_[globalIdxI] += entry * pressBound;
                        f_[globalIdxI] -= rightEntry;
                    }
                }
                //set neumann boundary condition
                else
                {
                    Scalar J = this->diffProblem.neumannPress(globalPosFace, *eIt, localPosFace);
                    f_[globalIdxI] -= faceVol*J;

                    //Assumes that the phases flow in the same direction at the neumann boundary, which is the direction of the total flux!!!
                    //needed to determine the upwind direction in the saturation equation
                    this->diffProblem.variables().potentialWetting()[globalIdxI][indexInInside] = J;
                    this->diffProblem.variables().potentialNonWetting()[globalIdxI][indexInInside] = J;
                }
            }
        }
        // end all intersections
    } // end grid traversal
    return;
}

//solves the system of equations to get the spatial distribution of the pressure
template<class GridView, class Scalar, class VC, class Problem>
void FVPressure2P<GridView, Scalar, VC, Problem>::solve()
{
	std::cout << "FVPressure2P: solve for pressure" << std::endl;

    MatrixAdapter<Matrix,Vector,Vector> op(A_);
    InverseOperatorResult r;
    double reduction = 1E-12;
    int maxIt = 10000;
    int verboseLevel = 1;
    InverseOperatorResult result;

    if (preconditionerName_ == "SeqILU0")
    {
    	SeqILU0<Matrix,Vector,Vector> preconditioner(A_, 1.0);
    	if (solverName_ == "CG")
    	{
    		CGSolver<Vector> solver(op, preconditioner, reduction, maxIt, verboseLevel);
    		solver.apply(this->diffProblem.variables().pressure(), f_, result);
    	}
    	else if (solverName_ == "BiCGSTAB")
    	{
    		BiCGSTABSolver<Vector> solver(op, preconditioner, reduction, maxIt, verboseLevel);
    		solver.apply(this->diffProblem.variables().pressure(), f_, result);
    	}
    	else
    		DUNE_THROW(NotImplemented, "FVPressure2P :: solve : combination "
    				<< preconditionerName_<< " and "<< solverName_ << ".");
    }
    else if (preconditionerName_ == "SeqPardiso")
    {
    	SeqPardiso<Matrix,Vector,Vector> preconditioner(A_);
    	if (solverName_ == "Loop")
    	{
    		LoopSolver<Vector> solver(op, preconditioner, reduction, maxIt, verboseLevel);
    		solver.apply(this->diffProblem.variables().pressure(), f_, result);
    	}
    	else
    		DUNE_THROW(NotImplemented, "FVPressure2P :: solve : combination "
    				<< preconditionerName_<< " and "<< solverName_ << ".");
    }
    else
    	DUNE_THROW(NotImplemented, "FVPressure2P :: solve : preconditioner "
    			<< preconditionerName_ << ".");

    return;
}
//constitutive functions are updated once if new saturations are calculated and stored in the variables object
template<class GridView, class Scalar, class VC, class Problem>
void FVPressure2P<GridView, Scalar, VC, Problem>::initializeMaterialLaws()
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

        std::vector<Scalar> mobilities(2,0.0);

        if (saturationType == Sw)
        {
            mobilities = this->diffProblem.materialLaw().mob(sat, globalPos, *eIt, localPos);
            this->diffProblem.variables().capillaryPressure()[globalIdx]= this->diffProblem.materialLaw().pC(sat, globalPos, *eIt, localPos);
        }
        else if (saturationType == Sn)
        {
            mobilities = this->diffProblem.materialLaw().mob(1-sat, globalPos, *eIt, localPos);
            this->diffProblem.variables().capillaryPressure()[globalIdx]= this->diffProblem.materialLaw().pC(1-sat, globalPos, *eIt, localPos);
        }
        else
        {
            DUNE_THROW(RangeError, "materialLaws not initialized!");
        }

        // initialize mobilities
        this->diffProblem.variables().mobilityWetting()[globalIdx]= mobilities[0];
        this->diffProblem.variables().mobilityNonWetting()[globalIdx]= mobilities[1];
        this->diffProblem.variables().fracFlowFuncWetting()[globalIdx]= mobilities[0]/(mobilities[0]+mobilities[1]);
        this->diffProblem.variables().fracFlowFuncNonWetting()[globalIdx]= mobilities[1]/(mobilities[0]+mobilities[1]);
    }
    return;
}

}
#endif
