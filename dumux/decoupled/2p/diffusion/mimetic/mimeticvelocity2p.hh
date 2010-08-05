// $Id: fvvelocity2p.hh 3455 2010-04-09 12:11:50Z mwolff $
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
#ifndef DUMUX_FVVELOCITY2P_HH
#define DUMUX_FVVELOCITY2P_HH

/**
 * @file
 * @brief  Finite Volume Diffusion Model
 * @author Markus Wolff
 */

#include <dumux/decoupled/2p/diffusion/fv/fvpressure2p.hh>

namespace Dumux
{
//! \ingroup diffusion
//! Finite Volume Diffusion Model
/*! Calculates non-wetting phase velocities from a known pressure field in context of a Finite Volume implementation for the evaluation
 * of equations of the form
 * \f[\text{div}\, \boldsymbol{v}_{total} = q.\f]
 * The wetting or the non-wetting phase pressure has to be given as piecewise constant cell values.
 * The velocity is calculated following  Darcy's law as
 * \f[\boldsymbol{v}_n = \lambda_n \boldsymbol{K} \left(\text{grad}\, p_n + \rho_n g  \text{grad}\, z\right),\f]
 * where, \f$p_n\f$ denotes the wetting phase pressure, \f$\boldsymbol{K}\f$ the absolute permeability, \f$\lambda_n\f$ the non-wetting phase mobility, \f$\rho_n\f$ the non-wetting phase density and \f$g\f$ the gravity constant.
 * As in the two-phase pressure equation a total flux depending on a total velocity is considered one has to be careful at neumann flux boundaries. Here, a phase velocity is only uniquely defined, if
 * the saturation is at the maximum (\f$1-S_{rw}\f$, \f$\boldsymbol{v}_{total} = \boldsymbol{v}_n\f$) or at the minimum (\f$ S_{rn} \f$, \f$\boldsymbol{v}_n = 0\f$)
 *
 * Template parameters are:
 *
 - GridView a DUNE gridview type
 - Scalar type used for scalar quantities
 - VC type of a class containing different variables of the model
 - Problem class defining the physical problem
 */

template<class TypeTag>
class FVVelocity2P: public FVPressure2P<TypeTag>
{
    typedef FVVelocity2P<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
     typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
     typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
     typedef typename GET_PROP(TypeTag, PTAG(ReferenceElements)) ReferenceElements;
     typedef typename ReferenceElements::Container ReferenceElementContainer;
     typedef typename ReferenceElements::ContainerFaces ReferenceElementFaceContainer;

     typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
     typedef typename SpatialParameters::MaterialLaw MaterialLaw;

     typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

     typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
     typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNW,
        pglobal = Indices::pressureGlobal,
        vw = Indices::velocityW,
        vn = Indices::velocityNW,
        vt = Indices::velocityTotal,
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW,
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    //! Constructs a FVNonWettingPhaseVelocity2P object
    /**
     * \param gridView gridView object of type GridView
     * \param problem a problem class object
     * \param pressureType a string giving the type of pressure used (could be: pw, pn, pglobal)
     * \param satType a string giving the type of saturation used (could be: Sw, Sn)
     */
    FVVelocity2P(Problem& problem)
    : FVPressure2P<TypeTag>(problem)
    {
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableCompressibility)) && velocityType_ == vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Total velocity - global pressure - model cannot be used with compressible fluids!");
        }
        if (velocityType_ != vw && velocityType_ != vn && velocityType_ != vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Velocity type not supported!");
        }
    }
    //! Constructs a FVNonWettingPhaseVelocity2P object
    /**
     * \param gridView gridView object of type GridView
     * \param problem a problem class object
     * \param pressureType a string giving the type of pressure used (could be: pw, pn, pglobal)
     * \param satType a string giving the type of saturation used (could be: Sw, Sn)
     * \param solverName a string giving the type of solver used (could be: CG, BiCGSTAB, Loop)
     * \param preconditionerName a string giving the type of the matrix preconditioner used (could be: Dune::SeqILU0, SeqPardiso)
     */
    FVVelocity2P(Problem& problem, std::string solverName,
            std::string preconditionerName)
    : FVPressure2P<TypeTag>(problem, solverName, preconditionerName)
    {
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableCompressibility)) && velocityType_ == vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Total velocity - global pressure - model cannot be used with compressible fluids!");
        }
        if (velocityType_ != vw && velocityType_ != vn && velocityType_ != vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Velocity type not supported!");
        }
    }

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
    void calculateVelocity();

private:
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, PTAG(VelocityFormulation)); //!< gives kind of velocity used (\f$ 0 = v_w\f$, \f$ 1 = v_n\f$, \f$ 2 = v_t\f$)
};
template<class TypeTag>
void FVVelocity2P<TypeTag>::calculateVelocity()
{
    // compute update vector
    ElementIterator eItEnd = this->problem().gridView().template end<0>();
    for (ElementIterator eIt = this->problem().gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // cell center in reference element
        const LocalPosition
        &localPos = ReferenceElementContainer::general(gt).position(0, 0);

        //
        GlobalPosition globalPos = eIt->geometry().global(localPos);

        // cell index
        int globalIdxI = this->problem().variables().index(*eIt);

        Scalar pressI = this->problem().variables().pressure()[globalIdxI];
        Scalar pcI = this->problem().variables().capillaryPressure(globalIdxI);
        Scalar lambdaWI = this->problem().variables().mobilityWetting(globalIdxI);
        Scalar lambdaNWI = this->problem().variables().mobilityNonwetting(globalIdxI);
        Scalar fractionalWI = this->problem().variables().fracFlowFuncWetting(globalIdxI);
        Scalar fractionalNWI = this->problem().variables().fracFlowFuncNonwetting(globalIdxI);
        Scalar densityWI = this->problem().variables().densityWetting(globalIdxI);
        Scalar densityNWI = this->problem().variables().densityNonwetting(globalIdxI);

        // run through all intersections with neighbors and boundary
        IntersectionIterator
        isItEnd = this->problem().gridView().template iend(*eIt);
        for (IntersectionIterator
                isIt = this->problem().gridView().template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {
            // local number of facet
            int isIndex = isIt->indexInInside();

            // get geometry type of face
            Dune::GeometryType faceGT = isIt->geometryInInside().type();

            // center in face's reference element
            const Dune::FieldVector<Scalar,dim-1>&
            faceLocal = ReferenceElementFaceContainer::general(faceGT).position(0,0);

            // center of face inside volume reference element
            const LocalPosition localPosFace(0);

            Dune::FieldVector<Scalar,dimWorld> unitOuterNormal = isIt->unitOuterNormal(faceLocal);

            // get absolute permeability
            FieldMatrix permeabilityI(this->problem().spatialParameters().intrinsicPermeability(globalPos, *eIt));

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = this->problem().variables().index(*neighborPointer);

                // compute factor in neighbor
                Dune::GeometryType neighborGT = neighborPointer->geometry().type();
                const LocalPosition&
                localPosNeighbor = ReferenceElementContainer::general(neighborGT).position(0,0);

                // cell center in global coordinates
                const GlobalPosition& globalPos = eIt->geometry().global(localPos);

                // neighbor cell center in global coordinates
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

                // distance vector between barycenters
                Dune::FieldVector<Scalar,dimWorld> distVec = globalPosNeighbor - globalPos;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                Dune::FieldVector<Scalar, dimWorld> unitDistVec(distVec);
                unitDistVec /= dist;

                // get absolute permeability
                FieldMatrix permeabilityJ(this->problem().spatialParameters().intrinsicPermeability(globalPosNeighbor, *neighborPointer));

                // compute vectorized permeabilities
                FieldMatrix meanPermeability(0);

                // harmonic mean of permeability
                for (int x = 0;x<dim;x++)
                {
                    for (int y = 0; y < dim;y++)
                    {
                        if (permeabilityI[x][y] && permeabilityJ[x][y])
                        {
                            meanPermeability[x][y]= 2*permeabilityI[x][y]*permeabilityJ[x][y]/(permeabilityI[x][y]+permeabilityJ[x][y]);
                        }
                    }
                }

                Dune::FieldVector<Scalar,dim> permeability(0);
                meanPermeability.mv(unitDistVec,permeability);

                Scalar pressJ = this->problem().variables().pressure()[globalIdxJ];
                Scalar pcJ = this->problem().variables().capillaryPressure(globalIdxJ);
                Scalar lambdaWJ = this->problem().variables().mobilityWetting(globalIdxJ);
                Scalar lambdaNWJ = this->problem().variables().mobilityNonwetting(globalIdxJ);
                Scalar fractionalWJ = this->problem().variables().fracFlowFuncWetting(globalIdxJ);
                Scalar fractionalNWJ = this->problem().variables().fracFlowFuncNonwetting(globalIdxJ);
                Scalar densityWJ = this->problem().variables().densityWetting(globalIdxJ);
                Scalar densityNWJ = this->problem().variables().densityNonwetting(globalIdxJ);

                //calculate potential gradients
                Scalar potentialW = 0;
                Scalar potentialNW = 0;

                potentialW = this->problem().variables().potentialWetting(globalIdxI, isIndex);
                potentialNW = this->problem().variables().potentialNonwetting(globalIdxI, isIndex);

                Scalar densityW = (potentialW> 0.) ? densityWI : densityWJ;
                Scalar densityNW = (potentialNW> 0.) ? densityNWI : densityNWJ;

                densityW = (potentialW == 0.) ? 0.5 * (densityWI + densityWJ) : densityW;
                densityNW = (potentialNW == 0.) ? 0.5 * (densityNWI + densityNWJ) : densityNW;

                switch (this->pressureType)
                {
                case pw:
                {
                    potentialW = (pressI - pressJ) / dist;
                    potentialNW = (pressI - pressJ+ pcI - pcJ) / dist;
                    break;
                }
                case pn:
                {
                    potentialW = (pressI - pressJ - pcI + pcJ) / dist;
                    potentialNW = (pressI - pressJ) / dist;
                    break;
                }
                case pglobal:
                {
                    potentialW = (pressI - pressJ - 0.5 * (fractionalNWI+fractionalNWJ)*(pcI - pcJ)) / dist;
                    potentialNW = (pressI - pressJ + 0.5 * (fractionalWI+fractionalWJ)*(pcI - pcJ)) / dist;
                    break;
                }
                }

                potentialW += densityW * (unitDistVec * this->gravity);//delta z/delta x in unitDistVec[z]
                potentialNW += densityNW * (unitDistVec * this->gravity);

                //store potentials for further calculations (velocity, saturation, ...)
                this->problem().variables().potentialWetting(globalIdxI, isIndex) = potentialW;
                this->problem().variables().potentialNonwetting(globalIdxI, isIndex) = potentialNW;

                //do the upwinding of the mobility depending on the phase potentials
                Scalar lambdaW = (potentialW > 0.) ? lambdaWI : lambdaWJ;
                lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaW;
                Scalar lambdaNW = (potentialNW > 0.) ? lambdaNWI : lambdaNWJ;
                lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWJ) : lambdaNW;
                densityW = (potentialW> 0.) ? densityWI : densityWJ;
                densityW = (potentialW == 0.) ? 0.5 * (densityWI + densityWJ) : densityW;
                densityNW = (potentialNW> 0.) ? densityNWI : densityNWJ;
                densityNW = (potentialNW == 0.) ? 0.5 * (densityNWI + densityNWJ) : densityNW;

                //calculate the gravity term
                Dune::FieldVector<Scalar,dimWorld> velocityW(permeability);
                Dune::FieldVector<Scalar,dimWorld> velocityNW(permeability);
                Dune::FieldVector<Scalar,dimWorld> gravityTermW(unitDistVec);
                Dune::FieldVector<Scalar,dimWorld> gravityTermNW(unitDistVec);

                gravityTermW *= (this->gravity*permeability)*(lambdaW * densityW);
                gravityTermNW *= (this->gravity*permeability)*(lambdaNW * densityNW);

                //calculate velocity depending on the pressure used -> use pc = pn - pw
                switch (this->pressureType)
                {
                case pw:
                {
                    velocityW *= lambdaW * (pressI - pressJ)/dist;
                    velocityNW *= lambdaNW * (pressI - pressJ)/dist + 0.5 * (lambdaNWI + lambdaNWJ) * (pcI - pcJ) / dist;
                    velocityW += gravityTermW;
                    velocityNW += gravityTermNW;
                    break;
                }
                case pn:
                {
                    velocityW *= lambdaW * (pressI - pressJ)/dist - 0.5 * (lambdaWI + lambdaWJ) * (pcI - pcJ) / dist;
                    velocityNW *= lambdaNW * (pressI - pressJ) / dist;
                    velocityW += gravityTermW;
                    velocityNW += gravityTermNW;
                    break;
                }
                case pglobal:
                {
                    this->problem().variables().velocity()[globalIdxI][isIndex] = permeability;
                    this->problem().variables().velocity()[globalIdxI][isIndex]*= (lambdaW + lambdaNW)* (pressI - pressJ )/ dist;
                    this->problem().variables().velocity()[globalIdxI][isIndex] += gravityTermW;
                    this->problem().variables().velocity()[globalIdxI][isIndex] += gravityTermNW;
                    break;
                }
                }

                //store velocities
                switch (velocityType_)
                {
                case vw:
                {
                    this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW;
                    this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityNW;
                    break;
                }
                case vn:
                {
                    this->problem().variables().velocity()[globalIdxI][isIndex] = velocityNW;
                    this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityW;
                    break;
                }
                case vt:
                {
                    switch (this->pressureType)
                    {
                    case pw:
                    {
                        this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW + velocityNW;
                        this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityNW;
                        break;
                    }
                    case pn:
                    {
                        this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW + velocityNW;
                        this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityW;
                        break;
                    }
                    }
                    break;
                }
                }
            }//end intersection with neighbor

            // handle boundary face
            if (isIt->boundary())
            {
                // center of face in global coordinates
                GlobalPosition globalPosFace = isIt->geometry().global(faceLocal);

                //get boundary type
                BoundaryConditions::Flags bcTypeSat = this->problem().bctypeSat(globalPosFace, *isIt);
                BoundaryConditions::Flags bcTypePress = this->problem().bctypePress(globalPosFace, *isIt);

                // cell center in global coordinates
                GlobalPosition globalPos = eIt->geometry().global(localPos);

                // distance vector between barycenters
                Dune::FieldVector<Scalar,dimWorld> distVec = globalPosFace - globalPos;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                Dune::FieldVector<Scalar, dimWorld> unitDistVec(distVec);
                unitDistVec /= dist;

                //multiply with normal vector at the boundary
                Dune::FieldVector<Scalar,dim> permeability(0);
                permeabilityI.mv(unitDistVec, permeability);

                Scalar satBound = 0;
                if (bcTypeSat == BoundaryConditions::dirichlet)
                {
                    satBound = this->problem().dirichletSat(globalPosFace, *isIt);
                }
                else
                {
                    satBound = this->problem().variables().saturation()[globalIdxI];
                }

                if (bcTypePress == BoundaryConditions::dirichlet)
                {
                    //determine phase saturations from primary saturation variable
                    Scalar satW;
                    Scalar satNW;
                    switch (this->saturationType)
                    {
                    case Sw:
                    {
                        satW = satBound;
                        satNW = 1-satBound;
                        break;
                    }
                    case Sn:
                    {
                        satW = 1-satBound;
                        satNW = satBound;
                        break;
                    }
                    default:
                    {
                        DUNE_THROW(Dune::RangeError, "saturation type not implemented");
                    }
                    }
                    Scalar pressBound = this->problem().dirichletPress(globalPosFace, *isIt);
                    Scalar pcBound = MaterialLaw::pC(this->problem().spatialParameters().materialLawParams(globalPos, *eIt), satW);

                    //determine phase pressures from primary pressure variable
                    Scalar pressW=0;
                    Scalar pressNW=0;
                    switch (this->pressureType)
                    {
                    case pw:
                    {
                        pressW = pressBound;
                        pressNW = pressBound + pcBound;
                        break;
                    }
                    case pn:
                    {
                        pressW = pressBound - pcBound;
                        pressNW = pressBound;
                        break;
                    }
                    }

                    //get temperature at current position
                    Scalar temperature = this->problem().temperature(globalPosFace, *eIt);

                    Scalar densityWBound = 0;
                    Scalar densityNWBound = 0;
                    Scalar lambdaWBound = 0;
                    Scalar lambdaNWBound = 0;
                    if (this->compressibility)
                    {
                        FluidState fluidState;
                        fluidState.update(satW, pressW, pressNW, temperature);
                        densityWBound = FluidSystem::phaseDensity(wPhaseIdx, temperature, pressW, fluidState);
                        densityNWBound = FluidSystem::phaseDensity(nPhaseIdx, temperature, pressNW, fluidState);
                        Scalar viscosityWBound = FluidSystem::phaseViscosity(wPhaseIdx, temperature, pressW, fluidState);
                        Scalar viscosityNWBound = FluidSystem::phaseViscosity(nPhaseIdx, temperature, pressNW, fluidState);
                        lambdaWBound = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos, *eIt), satW) / viscosityWBound * densityWBound;
                        lambdaNWBound = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos, *eIt), satW) / viscosityNWBound * densityNWBound;
                    }
                    else
                    {
                        Scalar referencePressure =  this->problem().referencePressure(globalPos, *eIt);
                        FluidState fluidState;
                        fluidState.update(satW, referencePressure, referencePressure, temperature);
                        densityWBound = FluidSystem::phaseDensity(wPhaseIdx, temperature, referencePressure, fluidState);
                        densityNWBound = FluidSystem::phaseDensity(nPhaseIdx, temperature, referencePressure, fluidState);
                        Scalar viscosityWBound = FluidSystem::phaseViscosity(wPhaseIdx, temperature, referencePressure, fluidState);
                        Scalar viscosityNWBound = FluidSystem::phaseViscosity(nPhaseIdx, temperature, referencePressure, fluidState);
                        lambdaWBound = MaterialLaw::krw(this->problem().spatialParameters().materialLawParams(globalPos, *eIt), satW) / viscosityWBound;
                        lambdaNWBound = MaterialLaw::krn(this->problem().spatialParameters().materialLawParams(globalPos, *eIt), satW) / viscosityNWBound;
                    }

                    Scalar potentialW = 0;
                    Scalar potentialNW = 0;

                    potentialW = this->problem().variables().potentialWetting(globalIdxI, isIndex);
                    potentialNW = this->problem().variables().potentialNonwetting(globalIdxI, isIndex);

                    Scalar densityW = (potentialW> 0.) ? densityWI : densityWBound;
                    Scalar densityNW = (potentialNW> 0.) ? densityNWI : densityNWBound;

                    densityW = (potentialW == 0.) ? 0.5 * (densityWI + densityWBound) : densityW;
                    densityNW = (potentialNW == 0.) ? 0.5 * (densityNWI + densityNWBound) : densityNW;

                    //calculate potential gradient
                    switch (this->pressureType)
                    {
                    case pw:
                    {
                        potentialW = (pressI - pressBound) / dist;
                        potentialNW = (pressI + pcI - pressBound - pcBound) / dist;
                        break;
                    }
                    case pn:
                    {
                        potentialW = (pressI - pcI - pressBound + pcBound) / dist;
                        potentialNW = (pressI - pressBound) / dist;
                        break;
                    }
                    case pglobal:
                    {
                        potentialW = (pressI - pressBound - fractionalNWI * (pcI - pcBound)) / dist;
                        potentialNW = (pressI - pressBound + fractionalWI * (pcI - pcBound)) / dist;
                        break;
                    }
                    }

                    potentialW += densityW * (unitDistVec * this->gravity);
                    potentialNW += densityNW * (unitDistVec * this->gravity);

                    //store potential gradients for further calculations
                    this->problem().variables().potentialWetting(globalIdxI, isIndex) = potentialW;
                    this->problem().variables().potentialNonwetting(globalIdxI, isIndex) = potentialNW;

                    //do the upwinding of the mobility depending on the phase potentials
                    Scalar lambdaW = (potentialW > 0.) ? lambdaWI : lambdaWBound;
                    lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWBound) : lambdaW;
                    Scalar lambdaNW = (potentialNW > 0.) ? lambdaNWI : lambdaNWBound;
                    lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWBound) : lambdaNW;
                    densityW = (potentialW> 0.) ? densityWI : densityWBound;
                    densityW = (potentialW == 0.) ? 0.5 * (densityWI + densityWBound) : densityW;
                    densityNW = (potentialNW> 0.) ? densityNWI : densityNWBound;
                    densityNW = (potentialNW == 0.) ? 0.5 * (densityNWI + densityNWBound) : densityNW;

                    //calculate the gravity term
                    Dune::FieldVector<Scalar,dimWorld> velocityW(permeability);
                    Dune::FieldVector<Scalar,dimWorld> velocityNW(permeability);
                    Dune::FieldVector<Scalar,dimWorld> gravityTermW(unitDistVec);
                    Dune::FieldVector<Scalar,dimWorld> gravityTermNW(unitDistVec);

                    gravityTermW *= (this->gravity*permeability)*(lambdaW * densityW);
                    gravityTermNW *= (this->gravity*permeability)*(lambdaNW * densityNW);

                    //calculate velocity depending on the pressure used -> use pc = pn - pw
                    switch (this->pressureType)
                    {
                    case pw:
                    {
                        velocityW *= lambdaW * (pressI - pressBound)/dist;
                        velocityNW *= lambdaNW * (pressI - pressBound)/dist + 0.5 * (lambdaNWI + lambdaNWBound) * (pcI - pcBound) / dist;
                        velocityW += gravityTermW;
                        velocityNW += gravityTermNW;
                        break;
                    }
                    case pn:
                    {
                        velocityW *= lambdaW * (pressI - pressBound)/dist - 0.5 * (lambdaWI + lambdaWBound) * (pcI - pcBound) / dist;
                        velocityNW *= lambdaNW * (pressI - pressBound) / dist;
                        velocityW += gravityTermW;
                        velocityNW += gravityTermNW;
                        break;
                    }
                    case pglobal:
                    {
                        this->problem().variables().velocity()[globalIdxI][isIndex] = permeability;
                        this->problem().variables().velocity()[globalIdxI][isIndex] *= (lambdaW + lambdaNW)* (pressI - pressBound )/ dist;
                        this->problem().variables().velocity()[globalIdxI][isIndex] += gravityTermW;
                        this->problem().variables().velocity()[globalIdxI][isIndex] += gravityTermNW;
                        break;
                    }
                    }

                    //store velocities
                    switch (velocityType_)
                    {
                    case vw:
                    {
                        this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW;
                        this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityNW;
                        break;
                    }
                    case vn:
                    {
                        this->problem().variables().velocity()[globalIdxI][isIndex] = velocityNW;
                        this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityW;
                        break;
                    }
                    case vt:
                    {
                        switch (this->pressureType)
                        {
                        case pw:
                        {
                            this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW + velocityNW;
                            this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityNW;
                            break;
                        }
                        case pn:
                        {
                            this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW + velocityNW;
                            this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityW;
                            break;
                        }
                        }
                        break;
                    }
                    }
                }//end dirichlet boundary

                else
                {
                    std::vector<Scalar> J = this->problem().neumannPress(globalPosFace, *isIt);
                    Dune::FieldVector<Scalar,dimWorld> velocityW(unitDistVec);
                    Dune::FieldVector<Scalar,dimWorld> velocityNW(unitDistVec);

                    velocityW *= J[wPhaseIdx];
                    velocityNW *= J[nPhaseIdx];

                    if (!this->compressibility)
                    {
                        velocityW /= densityWI;
                        velocityNW /= densityNWI;
                    }

                    switch (velocityType_)
                    {
                    case vw:
                    {
                        this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW;
                        this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityNW;
                        break;
                    }
                    case vn:
                    {
                        this->problem().variables().velocity()[globalIdxI][isIndex] = velocityNW;
                        this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityW;
                        break;
                    }
                    case vt:
                    {
                        switch (this->pressureType)
                        {
                        case pw:
                        {
                            this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW + velocityNW;
                            this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityNW;
                            break;
                        }
                        case pn:
                        {
                            this->problem().variables().velocity()[globalIdxI][isIndex] = velocityW + velocityNW;
                            this->problem().variables().velocitySecondPhase()[globalIdxI][isIndex] = velocityW;
                            break;
                        }
                        }
                        break;
                    }
                    }
                }//end neumann boundary
            }//end boundary
        }// end all intersections
    }// end grid traversal
                        printvector(std::cout, this->problem().variables().velocity(), "velocity", "row", 4, 1, 3);
    return;
}
}
#endif
