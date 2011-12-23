// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_FVVELOCITY2P_HH
#define DUMUX_FVVELOCITY2P_HH

/**
 * @file
 * @brief  Velocity Field from a finite volume solution of a pressure equation.
 * @author Markus Wolff
 */

#include <dumux/decoupled/2p/diffusion/fv/fvpressure2p.hh>

namespace Dumux
{
//! \ingroup FV2p
//! \brief Determines the velocity from a finite volume solution of the  pressure equation of the sequential Model (IMPES).
/*! Calculates phase velocities or total velocity from a known pressure field in context of a Finite Volume implementation for the evaluation
 * of equations of the form
 * \f[\text{div}\, \boldsymbol{v}_{total} = q.\f]
 * The wetting or the non-wetting phase pressure, or the global pressure has to be given as piecewise constant cell values.
 * The phase velocities are calculated following  Darcy's law as
 * \f[\boldsymbol{v}_\alpha = \lambda_\alpha \boldsymbol{K} \left(\text{grad}\, p_\alpha + \rho_\alpha g  \text{grad}\, z\right),\f]
 * where \f$p_\alpha\f$ denotes the pressure of phase \f$_\alpha\f$ (wetting or non-wetting), \f$\boldsymbol{K}\f$ the absolute permeability, \f$\lambda_\alpha\f$ the phase mobility, \f$\rho_\alpha\f$ the phase density and \f$g\f$ the gravity constant.
 * The total velocity is either calculated as sum of the phase velocities
 * \f[\boldsymbol{v}_{total} = \boldsymbol{v}_{wetting}+\boldsymbol{v}_{non-wetting},\f]
 * or with a given global pressure
 * \f[\boldsymbol{v}_{total} = \lambda_{total} \boldsymbol{K} \left(\text{grad}\, p_{global} + \sum f_\alpha \rho_\alpha g  \text{grad}\, z\right).\f]
 *
 * \tparam TypeTag The Type Tag
 */

template<class TypeTag>
class FVVelocity2P: public FVPressure2P<TypeTag>
{
    typedef FVVelocity2P<TypeTag> ThisType;
    typedef FVPressure2P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Variables)) Variables;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Indices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(WettingPhase)) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NonwettingPhase)) NonwettingPhase;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(CellData)) CellData;

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

    typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElementContainer;
    typedef Dune::GenericReferenceElements<Scalar, dim - 1> ReferenceElementFaceContainer;
    typedef Dune::GenericReferenceElement<Scalar, dim> ReferenceElement;

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
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        eqIdxPress = Indices::pressEqIdx,
        eqIdxSat = Indices::satEqIdx
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar, dim, dim> FieldMatrix;

public:
    //! Constructs a FVVelocity2P object
    /*!
     * \param problem a problem class object
     */
    FVVelocity2P(Problem& problem) :
            FVPressure2P<TypeTag>(problem), problem_(problem)
    {
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableCompressibility)) && velocityType_ == vt)
        {
            DUNE_THROW(Dune::NotImplemented,
                    "Total velocity - global pressure - model cannot be used with compressible fluids!");
        }
        if (velocityType_ != vw && velocityType_ != vn && velocityType_ != vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Velocity type not supported!");
        }
    }

    //! Calculate the velocity.
    /*!
     *
     *  Given the piecewise constant pressure \f$p\f$,
     *  this method calculates the velocity
     *  The method is needed in the IMPES (Implicit Pressure Explicit Saturation) algorithm which is used for a fractional flow formulation
     *  to provide the velocity field required for the solution of the saturation equation.
     */
    void calculateVelocity();

    void update()
    {
        ParentType::update();

        calculateVelocity();

        return;
    }

    //! \brief Write data files
    /*  \param name file name */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        ParentType::addOutputVtkFields(writer);

        Dune::BlockVector < Dune::FieldVector<Scalar, dim> > &velocity = *(writer.template allocateManagedBuffer<Scalar,
                dim>(problem_.gridView().size(0)));
        Dune::BlockVector < Dune::FieldVector<Scalar, dim> > &velocitySecondPhase =
                *(writer.template allocateManagedBuffer<Scalar, dim>(problem_.gridView().size(0)));

        // compute update vector
        ElementIterator eItEnd = problem_.gridView().template end<0>();
        for (ElementIterator eIt = problem_.gridView().template begin<0>(); eIt != eItEnd; ++eIt)
        {
            // cell index
            int globalIdx = problem_.variables().index(*eIt);

            CellData& cellData = problem_.variables().cellData(globalIdx);

            Dune::FieldVector < Scalar, 2 * dim > fluxW(0);
            Dune::FieldVector < Scalar, 2 * dim > fluxNW(0);
            // run through all intersections with neighbors and boundary
            IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
            for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
            {
                int isIndex = isIt->indexInInside();

                fluxW[isIndex] = isIt->geometry().volume()
                        * (isIt->centerUnitOuterNormal() * cellData.fluxData().velocity(wPhaseIdx, isIndex));
                fluxNW[isIndex] = isIt->geometry().volume()
                        * (isIt->centerUnitOuterNormal() * cellData.fluxData().velocity(nPhaseIdx, isIndex));
            }

            Dune::FieldVector < Scalar, dim > refVelocity(0);
            refVelocity[0] = 0.5 * (fluxW[1] - fluxW[0]);
            refVelocity[1] = 0.5 * (fluxW[3] - fluxW[2]);

            const Dune::FieldVector<Scalar, dim>& localPos =
                    ReferenceElementContainer::general(eIt->geometry().type()).position(0, 0);

            // get the transposed Jacobian of the element mapping
            const FieldMatrix& jacobianInv = eIt->geometry().jacobianInverseTransposed(localPos);
            FieldMatrix jacobianT(jacobianInv);
            jacobianT.invert();

            // calculate the element velocity by the Piola transformation
            Dune::FieldVector < Scalar, dim > elementVelocity(0);
            jacobianT.umtv(refVelocity, elementVelocity);
            elementVelocity /= eIt->geometry().integrationElement(localPos);

            velocity[globalIdx] = elementVelocity;

            refVelocity = 0;
            refVelocity[0] = 0.5 * (fluxNW[1] - fluxNW[0]);
            refVelocity[1] = 0.5 * (fluxNW[3] - fluxNW[2]);

            // calculate the element velocity by the Piola transformation
            elementVelocity = 0;
            jacobianT.umtv(refVelocity, elementVelocity);
            elementVelocity /= eIt->geometry().integrationElement(localPos);

            velocitySecondPhase[globalIdx] = elementVelocity;
        }

        //switch velocities
        switch (velocityType_)
        {
        case vw:
        {
            writer.attachCellData(velocity, "wetting-velocity", dim);
            writer.attachCellData(velocitySecondPhase, "non-wetting-velocity", dim);
            break;
        }
        case vn:
        {
            writer.attachCellData(velocity, "non-wetting-velocity", dim);
            writer.attachCellData(velocitySecondPhase, "wetting-velocity", dim);
            break;
        }
        case vt:
        {
            writer.attachCellData(velocity, "total velocity", dim);
            break;
        }
        }

        return;
    }

private:
    Problem& problem_;
    static const int velocityType_ = GET_PROP_VALUE(TypeTag, PTAG(VelocityFormulation)); //!< gives kind of velocity used (\f$ 0 = v_w\f$, \f$ 1 = v_n\f$, \f$ 2 = v_t\f$)
    static const bool compressibility_ = GET_PROP_VALUE(TypeTag, PTAG(EnableCompressibility));
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PTAG(PressureFormulation)); //!< gives kind of pressure used (\f$p_w\f$, \f$p_n\f$, \f$p_{global}\f$)
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, PTAG(SaturationFormulation)); //!< gives kind of saturation used (\f$S_w\f$, \f$S_n\f$)
};

template<class TypeTag>
void FVVelocity2P<TypeTag>::calculateVelocity()
{
    BoundaryTypes bcType;

    ElementIterator eItBegin = problem_.gridView().template begin<0>();
    Scalar densityW = 0;
    Scalar densityNW = 0;
    Scalar viscosityW = 0;
    Scalar viscosityNW = 0;

    if (!compressibility_)
    {
        Scalar temp = problem_.temperature(*eItBegin);
        Scalar pRef = problem_.referencePressure(*eItBegin);
        densityW = WettingPhase::density(temp, pRef);
        densityNW = NonwettingPhase::density(temp, pRef);
        viscosityW = WettingPhase::viscosity(temp, pRef);
        viscosityNW = NonwettingPhase::viscosity(temp, pRef);
    }

    // compute update vector
    ElementIterator eItEnd = problem_.gridView().template end<0>();
    for (ElementIterator eIt = eItBegin; eIt != eItEnd; ++eIt)
    {
#if HAVE_MPI
        if (eIt->partitionType() == Dune::GhostEntity || eIt->partitionType() == Dune::OverlapEntity)
        {
            continue;
        }
#endif

        // cell index
        int globalIdxI = problem_.variables().index(*eIt);

        CellData& cellDataI = problem_.variables().cellData(globalIdxI);

        GlobalPosition globalPos = eIt->geometry().center();

        this->storePressureSolution(globalIdxI, cellDataI);
        Scalar pcI = cellDataI.capillaryPressure();
        Scalar lambdaWI = cellDataI.mobility(wPhaseIdx);
        Scalar lambdaNWI = cellDataI.mobility(nPhaseIdx);
        Scalar fractionalWI = cellDataI.fracFlowFunc(wPhaseIdx);
        Scalar fractionalNWI = cellDataI.fracFlowFunc(nPhaseIdx);

        // run through all intersections with neighbors and boundary
        IntersectionIterator isItEnd = problem_.gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem_.gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            // local number of facet
            int isIndex = isIt->indexInInside();

            Dune::FieldVector < Scalar, dimWorld > unitOuterNormal = isIt->centerUnitOuterNormal();

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = problem_.variables().index(*neighborPointer);

                CellData& cellDataJ = problem_.variables().cellData(globalIdxJ);

                // neighbor cell center in global coordinates
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

                // distance vector between barycenters
                Dune::FieldVector < Scalar, dimWorld > distVec = globalPosNeighbor - globalPos;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

                // compute vectorized permeabilities
                FieldMatrix meanPermeability(0);

                problem_.spatialParameters().meanK(meanPermeability,
                        problem_.spatialParameters().intrinsicPermeability(*eIt),
                        problem_.spatialParameters().intrinsicPermeability(*neighborPointer));

                Dune::FieldVector < Scalar, dim > permeability(0);
                meanPermeability.mv(unitOuterNormal, permeability);

                this->storePressureSolution(globalIdxJ, cellDataJ);
                Scalar pcJ = cellDataJ.capillaryPressure();
                Scalar lambdaWJ = cellDataJ.mobility(wPhaseIdx);
                Scalar lambdaNWJ = cellDataJ.mobility(nPhaseIdx);
                Scalar fractionalWJ = cellDataJ.fracFlowFunc(wPhaseIdx);
                Scalar fractionalNWJ = cellDataJ.fracFlowFunc(nPhaseIdx);

                //calculate potential gradients
                Scalar potentialW = 0;
                Scalar potentialNW = 0;

                potentialW = cellDataI.fluxData().potential(wPhaseIdx, isIndex);
                potentialNW = cellDataI.fluxData().potential(nPhaseIdx, isIndex);

                if (compressibility_)
                {
                    densityW = (potentialW > 0.) ? cellDataI.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
                    densityNW = (potentialNW > 0.) ? cellDataI.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);

                    densityW =
                            (potentialW == 0) ? 0.5 * (cellDataI.density(wPhaseIdx) + cellDataJ.density(wPhaseIdx)) :
                                    densityW;
                    densityNW =
                            (potentialNW == 0) ? 0.5 * (cellDataI.density(nPhaseIdx) + cellDataJ.density(nPhaseIdx)) :
                                    densityNW;
                }

                if (pressureType_ == pglobal)
                {
                    potentialW = (cellDataI.globalPressure() - cellDataJ.globalPressure()
                            - 0.5 * (fractionalNWI + fractionalNWJ) * (pcI - pcJ));
                    potentialNW = (cellDataI.globalPressure() - cellDataJ.globalPressure()
                            + 0.5 * (fractionalWI + fractionalWJ) * (pcI - pcJ));
                }
                else
                {
                    potentialW = (cellDataI.pressure(wPhaseIdx) - cellDataJ.pressure(wPhaseIdx));
                    potentialNW = (cellDataI.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx));
                }

                potentialW += densityW * (distVec * this->gravity); //delta z/delta x in unitOuterNormal[z]
                potentialNW += densityNW * (distVec * this->gravity);

                //store potentials for further calculations (velocity, saturation, ...)
                cellDataI.fluxData().setPotential(wPhaseIdx, isIndex, potentialW);
                cellDataI.fluxData().setPotential(nPhaseIdx, isIndex, potentialNW);

                //do the upwinding of the mobility depending on the phase potentials
                Scalar lambdaW = (potentialW > 0.) ? lambdaWI : lambdaWJ;
                lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaW;
                Scalar lambdaNW = (potentialNW > 0.) ? lambdaNWI : lambdaNWJ;
                lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWJ) : lambdaNW;

                if (compressibility_)
                {
                    densityW = (potentialW > 0.) ? cellDataI.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
                    densityNW = (potentialNW > 0.) ? cellDataI.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);

                    densityW =
                            (potentialW == 0) ? 0.5 * (cellDataI.density(wPhaseIdx) + cellDataJ.density(wPhaseIdx)) :
                                    densityW;
                    densityNW =
                            (potentialNW == 0) ? 0.5 * (cellDataI.density(nPhaseIdx) + cellDataJ.density(nPhaseIdx)) :
                                    densityNW;
                }

                //calculate the gravity term
                Dune::FieldVector < Scalar, dimWorld > velocityW(permeability);
                Dune::FieldVector < Scalar, dimWorld > velocityNW(permeability);
                Dune::FieldVector < Scalar, dimWorld > gravityTermW(unitOuterNormal);
                Dune::FieldVector < Scalar, dimWorld > gravityTermNW(unitOuterNormal);

                gravityTermW *= (this->gravity * permeability) * (lambdaW * densityW);
                gravityTermNW *= (this->gravity * permeability) * (lambdaNW * densityNW);

                //calculate velocity depending on the pressure used -> use pc = pn - pw
                switch (pressureType_)
                {
                case pw:
                {
                    velocityW *= lambdaW * (cellDataI.pressure(wPhaseIdx) - cellDataJ.pressure(wPhaseIdx)) / dist;
                    velocityNW *= lambdaNW * (cellDataI.pressure(wPhaseIdx) - cellDataJ.pressure(wPhaseIdx)) / dist
                            + 0.5 * (lambdaNWI + lambdaNWJ) * (pcI - pcJ) / dist;
                    velocityW += gravityTermW;
                    velocityNW += gravityTermNW;
                    break;
                }
                case pn:
                {
                    velocityW *= lambdaW * (cellDataI.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx)) / dist
                            - 0.5 * (lambdaWI + lambdaWJ) * (pcI - pcJ) / dist;
                    velocityNW *= lambdaNW * (cellDataI.pressure(nPhaseIdx) - cellDataJ.pressure(nPhaseIdx)) / dist;
                    velocityW += gravityTermW;
                    velocityNW += gravityTermNW;
                    break;
                }
                case pglobal:
                {
                    velocityW *= (lambdaW + lambdaNW) * (cellDataI.globalPressure() - cellDataJ.globalPressure())
                            / dist;
                    velocityW += gravityTermW;
                    velocityW += gravityTermNW;
                    velocityNW = 0;
                    break;
                }
                }

                //store velocities
                cellDataI.fluxData().velocity(wPhaseIdx, isIndex) = velocityW;
                cellDataI.fluxData().velocity(nPhaseIdx, isIndex) = velocityNW;
            } //end intersection with neighbor

            // handle boundary face
            if (isIt->boundary())
            {
                //get boundary type
                problem_.boundaryTypes(bcType, *isIt);
                PrimaryVariables boundValues(0.0);

                if (bcType.isDirichlet(eqIdxPress))
                {
                    problem_.dirichlet(boundValues, *isIt);

                    // center of face in global coordinates
                    GlobalPosition globalPosFace = isIt->geometry().center();

                    // distance vector between barycenters
                    Dune::FieldVector < Scalar, dimWorld > distVec = globalPosFace - globalPos;

                    // compute distance between cell centers
                    Scalar dist = distVec.two_norm();

                    //permeability vector at boundary
                    // compute vectorized permeabilities
                    FieldMatrix meanPermeability(0);

                    problem_.spatialParameters().meanK(meanPermeability,
                            problem_.spatialParameters().intrinsicPermeability(*eIt));

                    Dune::FieldVector < Scalar, dim > permeability(0);
                    meanPermeability.mv(unitOuterNormal, permeability);

                    Scalar satBound = 0;
                    if (bcType.isDirichlet(eqIdxSat))
                    {
                        satBound = boundValues[saturationIdx];
                    }
                    else
                    {
                        satBound = problem_.variables().primaryVariablesGlobal(eqIdxSat)[globalIdxI];
                    }

                    //determine phase saturations from primary saturation variable
                    Scalar satW;
                    Scalar satNW;
                    switch (saturationType_)
                    {
                    case Sw:
                    {
                        satW = satBound;
                        satNW = 1 - satBound;
                        break;
                    }
                    case Sn:
                    {
                        satW = 1 - satBound;
                        satNW = satBound;
                        break;
                    }
                    default:
                    {
                        DUNE_THROW(Dune::RangeError, "saturation type not implemented");
                        break;
                    }
                    }
                    Scalar pressBound = boundValues[pressureIdx];
                    Scalar pcBound = MaterialLaw::pC(problem_.spatialParameters().materialLawParams(*eIt), satW);

                    //determine phase pressures from primary pressure variable
                    Scalar pressWBound = 0;
                    Scalar pressNWBound = 0;
                    switch (pressureType_)
                    {
                    case pw:
                    {
                        pressWBound = pressBound;
                        pressNWBound = pressBound + pcBound;
                        break;
                    }
                    case pn:
                    {
                        pressWBound = pressBound - pcBound;
                        pressNWBound = pressBound;
                        break;
                    }
                    }

                    //get temperature at current position
                    Scalar temperature = problem_.temperature(*eIt);

                    Scalar densityWBound = densityW;
                    Scalar densityNWBound = densityNW;
                    Scalar viscosityWBound = viscosityW;
                    Scalar viscosityNWBound =  viscosityNW;

                    if (compressibility_)
                    {
                        FluidState fluidState;
                        fluidState.setSaturation(wPhaseIdx, satW);
                        fluidState.setSaturation(nPhaseIdx, satNW);
                        fluidState.setTemperature(temperature);
                        fluidState.setPressure(wPhaseIdx, pressWBound);
                        fluidState.setPressure(nPhaseIdx, pressNWBound);
                        densityWBound = FluidSystem::density(fluidState, wPhaseIdx);
                        densityNWBound = FluidSystem::density(fluidState, nPhaseIdx);
                        viscosityWBound = FluidSystem::viscosity(fluidState, wPhaseIdx)/densityWBound;
                        viscosityNWBound = FluidSystem::viscosity(fluidState, nPhaseIdx)/densityNWBound;
                    }



                    Scalar lambdaWBound = MaterialLaw::krw(problem_.spatialParameters().materialLawParams(*eIt), satW)
                            / viscosityWBound;
                    Scalar lambdaNWBound = MaterialLaw::krn(problem_.spatialParameters().materialLawParams(*eIt), satW)
                            / viscosityNWBound;

                    Scalar potentialW = 0;
                    Scalar potentialNW = 0;

                    potentialW = cellDataI.fluxData().potential(wPhaseIdx, isIndex);
                    potentialNW = cellDataI.fluxData().potential(nPhaseIdx, isIndex);

                    if (compressibility_)
                    {
                        densityW = (potentialW > 0.) ? cellDataI.density(wPhaseIdx) : densityWBound;
                        densityW = (potentialW == 0) ? 0.5 * (cellDataI.density(wPhaseIdx) + densityWBound) : densityW;
                        densityNW = (potentialNW > 0.) ? cellDataI.density(nPhaseIdx) : densityNWBound;
                        densityNW = (potentialNW == 0) ? 0.5 * (cellDataI.density(nPhaseIdx) + densityNWBound) : densityNW;
                    }

                    //calculate potential gradient
                    if (pressureType_ == pglobal)
                    {
                        potentialW = (cellDataI.globalPressure() - pressBound - fractionalNWI * (pcI - pcBound));
                        potentialNW = (cellDataI.globalPressure() - pressBound + fractionalWI * (pcI - pcBound));
                    }
                    else
                    {
                        potentialW = (cellDataI.pressure(wPhaseIdx) - pressWBound);
                        potentialNW = (cellDataI.pressure(nPhaseIdx) - pressNWBound);
                    }

                    potentialW += densityW * (distVec * this->gravity);
                    potentialNW += densityNW * (distVec * this->gravity);

                    //store potential gradients for further calculations
                    cellDataI.fluxData().setPotential(wPhaseIdx, isIndex, potentialW);
                    cellDataI.fluxData().setPotential(nPhaseIdx, isIndex, potentialNW);

                    //do the upwinding of the mobility depending on the phase potentials
                    Scalar lambdaW = (potentialW > 0.) ? lambdaWI : lambdaWBound;
                    lambdaW = (potentialW == 0) ? 0.5 * (lambdaWI + lambdaWBound) : lambdaW;
                    Scalar lambdaNW = (potentialNW > 0.) ? lambdaNWI : lambdaNWBound;
                    lambdaNW = (potentialNW == 0) ? 0.5 * (lambdaNWI + lambdaNWBound) : lambdaNW;

                    if (compressibility_)
                    {
                        densityW = (potentialW > 0.) ? cellDataI.density(wPhaseIdx) : densityWBound;
                        densityW = (potentialW == 0) ? 0.5 * (cellDataI.density(wPhaseIdx) + densityWBound) : densityW;
                        densityNW = (potentialNW > 0.) ? cellDataI.density(nPhaseIdx) : densityNWBound;
                        densityNW = (potentialNW == 0) ? 0.5 * (cellDataI.density(nPhaseIdx) + densityNWBound) : densityNW;
                    }

                    //calculate the gravity term
                    Dune::FieldVector < Scalar, dimWorld > velocityW(permeability);
                    Dune::FieldVector < Scalar, dimWorld > velocityNW(permeability);
                    Dune::FieldVector < Scalar, dimWorld > gravityTermW(unitOuterNormal);
                    Dune::FieldVector < Scalar, dimWorld > gravityTermNW(unitOuterNormal);

                    gravityTermW *= (this->gravity * permeability) * (lambdaW * densityW);
                    gravityTermNW *= (this->gravity * permeability) * (lambdaNW * densityNW);

                    //calculate velocity depending on the pressure used -> use pc = pn - pw
                    switch (pressureType_)
                    {
                    case pw:
                    {
                        velocityW *= lambdaW * (cellDataI.pressure(wPhaseIdx) - pressBound) / dist;
                        velocityNW *= lambdaNW * (cellDataI.pressure(wPhaseIdx) - pressBound) / dist
                                + 0.5 * (lambdaNWI + lambdaNWBound) * (pcI - pcBound) / dist;
                        velocityW += gravityTermW;
                        velocityNW += gravityTermNW;
                        break;
                    }
                    case pn:
                    {
                        velocityW *= lambdaW * (cellDataI.pressure(nPhaseIdx) - pressBound) / dist
                                - 0.5 * (lambdaWI + lambdaWBound) * (pcI - pcBound) / dist;
                        velocityNW *= lambdaNW * (cellDataI.pressure(nPhaseIdx) - pressBound) / dist;
                        velocityW += gravityTermW;
                        velocityNW += gravityTermNW;
                        break;
                    }
                    case pglobal:
                    {
                        velocityW *= (lambdaW + lambdaNW) * (cellDataI.globalPressure() - pressBound) / dist;
                        velocityW += gravityTermW;
                        velocityW += gravityTermNW;
                        velocityNW = 0;
                        break;
                    }
                    }

                    //store velocities
                    cellDataI.fluxData().velocity(wPhaseIdx, isIndex) = velocityW;
                    cellDataI.fluxData().velocity(nPhaseIdx, isIndex) = velocityNW;

                } //end dirichlet boundary

                else if (bcType.isNeumann(eqIdxPress))
                {
                    problem_.neumann(boundValues, *isIt);

                    Dune::FieldVector < Scalar, dimWorld > velocityW(unitOuterNormal);
                    Dune::FieldVector < Scalar, dimWorld > velocityNW(unitOuterNormal);

                    velocityW *= boundValues[wPhaseIdx];
                    velocityNW *= boundValues[nPhaseIdx];

                    if (!compressibility_)
                    {
                        velocityW /= densityW;
                        velocityNW /= densityNW;
                    }

                    cellDataI.fluxData().velocity(wPhaseIdx, isIndex) = velocityW;
                    cellDataI.fluxData().velocity(nPhaseIdx, isIndex) = velocityNW;
                } //end neumann boundary
                else
                {
                    DUNE_THROW(Dune::NotImplemented, "No valid boundary condition type defined for pressure equation!");
                }
            } //end boundary
        } // end all intersections
    } // end grid traversal
//                        printvector(std::cout, cellDataI.fluxData().velocity(), "velocity", "row", 4, 1, 3);
    return;
}
}
#endif
