// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
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
#ifndef DUNE_FVWETTINGVELOCITY2P_HH
#define DUNE_FVWETTINGVELOCITY2P_HH

#include "dumux/diffusion/fv/fvpressure2p.hh"
#include "dumux/diffusion/diffusionproblem.hh"

/**
 * @file
 * @brief  Finite Volume Diffusion Model
 * @author Markus Wolff
 */

namespace Dune
{

//! \ingroup diffusion
//! Finite Volume Diffusion Model
/*! Calculates wetting phase velocities from a known pressure field in context of a Finite Volume implementation for the evaluation
 * of equations of the form
 * \f[\text{div}\, \boldsymbol{v}_{total} = q.\f]
 * The wetting or the non-wetting phase pressure has to be given as piecewise constant cell values.
 * The velocity is calculated following  Darcy's law as
 * \f[\boldsymbol{v}_w = \lambda_w \boldsymbol{K} \left(\text{grad}\, p_w + \rho_w g  \text{grad}\, z\right),\f]
 * where, \f$p_w\f$ denotes the wetting phase pressure, \f$\boldsymbol{K}\f$ the absolute permeability, \f$\lambda_w\f$ the wetting phase mobility, \f$\rho_w\f$ the wetting phase density and \f$g\f$ the gravity constant.
 * As in the two-phase pressure equation a total flux depending on a total velocity is considered one has to be careful at neumann flux boundaries. Here, a phase velocity is only uniquely defined, if
 * the saturation is at the maximum (\f$1-S_{rn}\f$, \f$\boldsymbol{v}_{total} = \boldsymbol{v}_w\f$) or at the minimum (\f$ S_{rw} \f$, \f$\boldsymbol{v}_w = 0\f$)
 *
 * Template parameters are:
 *
 - GridView      a DUNE gridview type
 - Scalar        type used for scalar quantities
 - VC            type of a class containing different variables of the model
 - Problem       class defining the physical problem
 */
template<class GridView, class Scalar, class VC, class Problem = DiffusionProblem<GridView, Scalar, VC> >
class FVWettingPhaseVelocity2P: public FVPressure2P<GridView, Scalar, VC, Problem>
{
	typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;

    enum
    {   dim = GridView::dimension};
    enum
    {   dimWorld = GridView::dimensionworld};
    enum
    {
        pw = 0, pn = 1, pglobal = 2
    };

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    //! Constructs a FVWettingPhaseVelocity2P object
    /**
     * \param gridView gridView object of type GridView
     * \param problem a problem class object
     * \param pressType a string giving the type of pressure used (could be: pw, pn, pglobal)
     * \param satType a string giving the type of saturation used (could be: Sw, Sn)
     */
    FVWettingPhaseVelocity2P(GridView& gridView, Problem& problem, std::string pressureType, std::string satType = "Sw")
    : FVPressure2P<GridView,Scalar,VC, Problem>(gridView, problem, pressureType, satType)
    {}
    //! Constructs a FVWettingPhaseVelocity2P object
    /**
     * \param gridView gridView object of type GridView
     * \param problem a problem class object
     * \param pressType a string giving the type of pressure used (could be: pw, pn, pglobal)
     * \param satType a string giving the type of saturation used (could be: Sw, Sn)
     * \param solverName a string giving the type of solver used (could be: CG, BiCGSTAB, Loop)
     * \param preconditionerName a string giving the type of the matrix preconditioner used (could be: SeqILU0, SeqPardiso)
     */
    FVWettingPhaseVelocity2P(GridView& gridView, Problem& problem, std::string pressureType, std::string satType, std::string solverName,
            std::string preconditionerName)
    : FVPressure2P<GridView,Scalar,VC, Problem>(gridView, problem, pressureType, satType, solverName, preconditionerName)
    {}

    //! Calculate the velocity.
    /*!
     *  \param t time
     *
     *
     *  Given the piecewise constant pressure \f$p\f$,
     *  this method calculates the velocity
     *  The method is needed in the IMPES (Implicit Pressure Explicit Saturation) algorithm which is used for a fractional flow formulation
     *  to provide the velocity field required for the solution of the saturation equation.
     */
    void calculateVelocity(const Scalar t=0) const
    {
        // phase densities
        Scalar densityW = this->diffProblem.wettingPhase().density();

        // compute update vector
        ElementIterator eItEnd = this->gridView.template end<0>();
        for (ElementIterator eIt = this->gridView.template begin<0>(); eIt != eItEnd; ++eIt)
        {
            // cell geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // cell center in reference element
            const LocalPosition
            &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

            //
            GlobalPosition globalPos = eIt->geometry().global(localPos);

            // cell index
            int globalIdxI = this->diffProblem.variables().indexDiffusion(*eIt);

            Scalar pressI = this->diffProblem.variables().pressure()[globalIdxI];
            Scalar pcI = this->diffProblem.variables().capillaryPressure()[globalIdxI];
            Scalar lambdaWI = this->diffProblem.variables().mobilityWetting()[globalIdxI];

            // run through all intersections with neighbors and boundary
            IntersectionIterator
            isItEnd = this->gridView.template iend(*eIt);
            for (IntersectionIterator
                    isIt = this->gridView.template ibegin(*eIt); isIt
                    !=isItEnd; ++isIt)
            {
                // local number of facet
                int indexInInside = isIt->indexInInside();

                // get geometry type of face
                Dune::GeometryType faceGT = isIt->geometryInInside().type();

                // center in face's reference element
                const Dune::FieldVector<Scalar,dim-1>&
                faceLocal = Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                // center of face inside volume reference element
                const LocalPosition&
                localPosFace = Dune::ReferenceElements<Scalar,dim>::general(faceGT).position(indexInInside,1);

                Dune::FieldVector<Scalar,dimWorld> unitOuterNormal = isIt->unitOuterNormal(faceLocal);

                // get absolute permeability
                FieldMatrix permeabilityI(this->diffProblem.soil().K(globalPos, *eIt, localPos));
                // compute directed permeability vector permeabilityI.n
                FieldVector<Scalar,dim> normalPermeabilityI(0);
                permeabilityI.umv(unitOuterNormal, normalPermeabilityI);

                // handle interior face
                if (isIt->neighbor())
                {
                    // access neighbor
                    ElementPointer neighborPointer = isIt->outside();
                    int globalIdxJ = this->diffProblem.variables().indexDiffusion(*neighborPointer);

                    // compute factor in neighbor
                    Dune::GeometryType neighborGT = neighborPointer->geometry().type();
                    const LocalPosition&
                    localPosNeighbor = Dune::ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

                    // cell center in global coordinates
                    const GlobalPosition& globalPos = eIt->geometry().global(localPos);

                    // neighbor cell center in global coordinates
                    const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

                    // distance vector between barycenters
                    Dune::FieldVector<Scalar,dimWorld> distVec = globalPosNeighbor - globalPos;

                    // compute distance between cell centers
                    Scalar dist = distVec.two_norm();

                    // get absolute permeability
                    FieldMatrix permeabilityJ(this->diffProblem.soil().K(globalPosNeighbor, *neighborPointer, localPosNeighbor));

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

                    Scalar pressJ = this->diffProblem.variables().pressure()[globalIdxJ];
                    Scalar pcJ = this->diffProblem.variables().capillaryPressure()[globalIdxJ];
                    Scalar lambdaWJ = this->diffProblem.variables().mobilityWetting()[globalIdxJ];

                    Scalar potentialW = this->diffProblem.variables().potentialWetting()[globalIdxI][indexInInside];

                    //do the upwinding of the mobility depending on the phase potentials
                    Scalar lambdaW = (potentialW >= 0.) ? lambdaWI : lambdaWJ;

                    //calculate the gravity term
                    FieldVector<Scalar,dimWorld> velocity(permeability);
                    FieldVector<Scalar,dimWorld> gravityTerm(this->gravity);

                    for (int i = 0;i<dim;i++)
                    {
                        gravityTerm[i] *= permeability[i]*unitOuterNormal[i];
                    }

                    //calculate velocity depending on the pressure used -> use pc = pn - pw
                    if (this->pressureType == pw)
                    {
                        velocity *= lambdaW * (pressI - pressJ) / dist;
                        velocity += (gravityTerm *= (lambdaW * densityW));
                    }
                    if (this->pressureType == pn)
                    {
                        velocity *= lambdaW * (pressI - pressJ)/dist - 0.5*(lambdaWI+lambdaWJ)*(pcI - pcJ) / dist;
                        velocity += (gravityTerm *= (lambdaW * densityW));
                    }
                    if (this->pressureType == pglobal)
                    {
                        DUNE_THROW(NotImplemented, "Pressure type not supported for v_w!");
                    }
                    this->diffProblem.variables().velocity()[globalIdxI][indexInInside] = velocity;

                }

                // handle boundary face
                if (isIt->boundary())
                {
                    // center of face in global coordinates
                    GlobalPosition globalPosFace = isIt->geometry().global(faceLocal);

                    //get boundary type
                    BoundaryConditions::Flags bcTypeSat = this->diffProblem.bctypeSat(globalPosFace, *eIt, localPosFace);
                    BoundaryConditions::Flags bcTypePress = this->diffProblem.bctypePress(globalPosFace, *eIt, localPosFace);

                    // cell center in global coordinates
                    GlobalPosition globalPos = eIt->geometry().global(localPos);

                    // distance vector between barycenters
                    Dune::FieldVector<Scalar,dimWorld> distVec = globalPosFace - globalPos;

                    // compute distance between cell centers
                    Scalar dist = distVec.two_norm();

                    Scalar residualSatW = this->diffProblem.soil().Sr_w(globalPos, *eIt, localPos);
                    Scalar residualSatNW = this->diffProblem.soil().Sr_n(globalPos, *eIt, localPos);

                    Scalar satBound = 0;
                    if (bcTypeSat == BoundaryConditions::dirichlet)
                    {
                        satBound = this->diffProblem.dirichletSat(globalPosFace, *eIt, localPosFace);
                    }
                    else
                    {
                        satBound = this->diffProblem.variables().saturation()[globalIdxI];
                    }

                    if (bcTypePress == BoundaryConditions::dirichlet)
                    {
                        Scalar pressBound = this->diffProblem.dirichletPress(globalPosFace, *eIt, localPosFace);
                        Scalar pcBound = this->diffProblem.materialLaw().pC(satBound, globalPosFace, *eIt, localPosFace);
                        Scalar lambdaWBound = this->diffProblem.materialLaw().mobW(satBound,globalPosFace, *eIt, localPosFace);

                        Scalar potentialW = this->diffProblem.variables().potentialWetting()[globalIdxI][indexInInside];

                        //do the upwinding of the mobility depending on the phase potentials
                        Scalar lambdaW = (potentialW >= 0.) ? lambdaWI : lambdaWBound;

                        //calculate the gravity term
                        FieldVector<Scalar,dimWorld> velocity(normalPermeabilityI);
                        FieldVector<Scalar,dimWorld> gravityTerm(this->gravity);
                        for (int i=0;i<dim;i++)
                        {
                            gravityTerm[i] *= normalPermeabilityI[i]*unitOuterNormal[i];
                        }

                        //calculate velocity depending on the pressure used -> use pc = pn - pw
                        if (this->pressureType == pw)
                        {
                            velocity *= lambdaW * (pressI - pressBound) / dist;
                            velocity += (gravityTerm *= (lambdaW * densityW));
                        }
                        if (this->pressureType == pn)
                        {
                            velocity *= lambdaW * (pressI - pressBound)/dist - 0.5 * (lambdaWI + lambdaWBound) * (pcI - pcBound) / dist;
                            velocity += (gravityTerm *= (lambdaW * densityW));
                        }
                        if (this->pressureType == pglobal)
                        {
                            DUNE_THROW(NotImplemented, "Pressure type not supported for v_w!");
                        }

                        this->diffProblem.variables().velocity()[globalIdxI][indexInInside] = velocity;
                    }
                    else
                    {
                        Scalar J = this->diffProblem.neumannPress(globalPosFace, *eIt, localPosFace);
                        FieldVector<Scalar,dimWorld> velocity(unitOuterNormal);

                        if (J != 0)
                        {
                            if (satBound <= residualSatW)
                            {
                                velocity *= 0;
                            }
                            else if (satBound >= (1-residualSatNW))
                            {
                                velocity *= J;
                            }
                            else
                            {
                                DUNE_THROW(NotImplemented, "v_w can not be determined from v_t at a neumann boundary if both phases are present!");
                            }
                        }
                        else
                        {
                            velocity *= 0;
                        }
                        this->diffProblem.variables().velocity()[globalIdxI][indexInInside] = velocity;
                    }
                }
            }// end all intersections
        }// end grid traversal
//        printvector(std::cout, this->diffProblem.variables().velocity(), "velocity", "row", 4, 1, 3);

        return;
    }
};
}
#endif
