// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
#ifndef DUMUX_FVVELOCITY2P_ADAPTIVE_HH
#define DUMUX_FVVELOCITY2P_ADAPTIVE_HH

/**
 * @file
 * @brief  Velocity Field from a finite volume solution of a pressure equation.
 * @author Markus Wolff
 */

#include <dumux/decoupled/2p/diffusion/fv/fvvelocity2p.hh>

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
class FVVelocity2PAdaptive: public FVVelocity2P<TypeTag>
{
    typedef FVVelocity2P<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
     typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
     typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;


     typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
     typedef typename SpatialParameters::MaterialLaw MaterialLaw;

     typedef typename GET_PROP_TYPE(TypeTag, TwoPIndices) Indices;

     typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
     typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

     typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
     typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVariables PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef Dune::GenericReferenceElements<Scalar, dim> ReferenceElementContainer;

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
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx, numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    //! Constructs a FVVelocity2PAdaptive object
    /*!
     * \param problem a problem class object
     */
    FVVelocity2PAdaptive(Problem& problem)
    : ParentType(problem), problem_(problem), gravity_(problem.gravity())
    {
    	// todo: kompatibilität prüfen
        if (GET_PROP_VALUE(TypeTag, EnableCompressibility) && velocityType_ == vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Total velocity - global pressure - model cannot be used with compressible fluids!");
        }
        if (velocityType_ != vw && velocityType_ != vn && velocityType_ != vt)
        {
            DUNE_THROW(Dune::NotImplemented, "Velocity type not supported!");
        }

        if (!compressibility_)
        {
            const Element& element = *(problem_.gridView().template begin<0> ());
            FluidState fluidState;
            fluidState.setPressure(wPhaseIdx, problem_.referencePressure(element));
            fluidState.setPressure(nPhaseIdx, problem_.referencePressure(element));
            fluidState.setTemperature(problem_.temperature(element));
            fluidState.setSaturation(wPhaseIdx, 1.);
            fluidState.setSaturation(nPhaseIdx, 0.);
            density_[wPhaseIdx] = FluidSystem::density(fluidState, wPhaseIdx);
            density_[nPhaseIdx] = FluidSystem::density(fluidState, nPhaseIdx);
            viscosity_[wPhaseIdx] = FluidSystem::viscosity(fluidState, wPhaseIdx);
            viscosity_[nPhaseIdx] = FluidSystem::viscosity(fluidState, nPhaseIdx);
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
    void calculateVelocity(const Intersection& intersection, CellData& cellDataI);

    bool calculateVelocityInTransport()
    {
        return false;
    }

private:
    Problem& problem_;
    const GlobalPosition& gravity_; //!< vector including the gravity constant
    Scalar density_[numPhases];
    Scalar viscosity_[numPhases];

    static const int velocityType_ = GET_PROP_VALUE(TypeTag, VelocityFormulation); //!< gives kind of velocity used (\f$ 0 = v_w\f$, \f$ 1 = v_n\f$, \f$ 2 = v_t\f$)
    static const bool compressibility_ = GET_PROP_VALUE(TypeTag, EnableCompressibility);
    static const int pressureType_ = GET_PROP_VALUE(TypeTag, PressureFormulation); //!< gives kind of pressure used (\f$p_w\f$, \f$p_n\f$, \f$p_{global}\f$)
    static const int saturationType_ = GET_PROP_VALUE(TypeTag, SaturationFormulation); //!< gives kind of saturation used (\f$S_w\f$, \f$S_n\f$)
};

template<class TypeTag>
void FVVelocity2PAdaptive<TypeTag>::calculateVelocity(const Intersection& intersection, CellData& cellDataI)
{
    ElementPointer elementI = intersection.inside();
    ElementPointer elementJ = intersection.outside();

    if (elementI->level() == elementJ->level())
    {
        ParentType::calculateVelocity(intersection, cellDataI);
    }
    else if (elementJ->level() == elementI->level() + 1)
    {
        CellData& cellDataJ = problem_.variables().cellData(problem_.variables().index(*elementJ));

        // get global coordinates of cell centers
        const GlobalPosition& globalPosI = elementI->geometry().center();
        const GlobalPosition& globalPosJ = elementJ->geometry().center();

        // get mobilities and fractional flow factors
        Scalar lambdaWI = cellDataI.mobility(wPhaseIdx);
        Scalar lambdaNWI = cellDataI.mobility(nPhaseIdx);
        Scalar fractionalWI = cellDataI.fracFlowFunc(wPhaseIdx);
        Scalar fractionalNWI = cellDataI.fracFlowFunc(nPhaseIdx);
        Scalar lambdaWJ = cellDataJ.mobility(wPhaseIdx);
        Scalar lambdaNWJ = cellDataJ.mobility(nPhaseIdx);
        Scalar fractionalWJ = cellDataJ.fracFlowFunc(wPhaseIdx);
        Scalar fractionalNWJ = cellDataJ.fracFlowFunc(nPhaseIdx);

        // get capillary pressure
        Scalar pcI = cellDataI.capillaryPressure();
        Scalar pcJ = cellDataJ.capillaryPressure();

        //get face index
        int isIndexI = intersection.indexInInside();

        //get face normal
        const Dune::FieldVector<Scalar, dim>& unitOuterNormal = intersection.centerUnitOuterNormal();

        // Count number of hanging nodes
        // not really necessary
        int isIndexJ = intersection.indexInOutside();

        int globalIdxK = 0;
        ElementPointer elementK = intersection.outside();
        // We are looking for two things:
        // IsIndexJ, the index of the interface from the neighbor-cell point of view
        // GlobalIdxK, the index of the third cell
        // for efficienty this is done in one IntersectionIterator-Loop

        // Intersectioniterator around cell I
        IntersectionIterator isItEndI = problem_.gridView().iend(*elementI);
        for (IntersectionIterator isItI = problem_.gridView().ibegin(*elementI); isItI != isItEndI; ++isItI)
        {
            if (isItI->neighbor())
            {
                ElementPointer neighborPointer2 = isItI->outside();

                // make sure we do not chose elemntI as third element
                // -> faces with hanging node have more than one intersection but only one face index!
                if (neighborPointer2 != elementJ && isItI->indexInInside() == isIndexI)
                {
                    globalIdxK = problem_.variables().index(*neighborPointer2);
                    elementK = neighborPointer2;
                    break;
                }
            }
        }

        CellData& cellDataK = problem_.variables().cellData(globalIdxK);

        // face global coordinates
        const GlobalPosition& globalPosInterface = intersection.geometry().center();

        Dune::FieldVector < Scalar, dimWorld > distVec = globalPosInterface - globalPosI;
        Scalar lI = distVec * unitOuterNormal;
        distVec = globalPosJ - globalPosInterface;
        Scalar lJ = distVec * unitOuterNormal;
        Scalar l = lI + lJ;

        FieldMatrix permeabilityI(0);
        FieldMatrix permeabilityJ(0);
        FieldMatrix permeabilityK(0);

        problem_.spatialParameters().meanK(permeabilityI,
                problem_.spatialParameters().intrinsicPermeability(*elementI));
        problem_.spatialParameters().meanK(permeabilityJ,
                problem_.spatialParameters().intrinsicPermeability(*elementJ));
        problem_.spatialParameters().meanK(permeabilityK,
                problem_.spatialParameters().intrinsicPermeability(*elementK));

        // Calculate permeablity component normal to interface
        Scalar kI, kJ, kK, ng, kMean; //, gI, gJ, gK;
        Dune::FieldVector < Scalar, dim > permI(0);
        Dune::FieldVector < Scalar, dim > permJ(0);
        Dune::FieldVector < Scalar, dim > permK(0);

        permeabilityI.mv(unitOuterNormal, permI);
        permeabilityJ.mv(unitOuterNormal, permJ);
        permeabilityK.mv(unitOuterNormal, permK);

        // kI,kJ,kK = (n^T)Kn
        kI = unitOuterNormal * permI;
        kJ = unitOuterNormal * permJ;
        kK = unitOuterNormal * permK;
        kMean = kI * kJ * kK * l / (kJ * kK * lI + kI * (kJ + kK) / 2 * lJ);

        ng = gravity_ * unitOuterNormal;

        Scalar fractionalWK = cellDataK.fracFlowFunc(wPhaseIdx);
        Scalar fractionalNWK = cellDataK.fracFlowFunc(nPhaseIdx);

        Scalar pcK = cellDataK.capillaryPressure();
        Scalar pcJK = (pcJ + pcK) / 2;

        // calculate potential gradients
        // reuse potentials from fvpressure2padaptive
        Scalar potentialWIJ = cellDataI.fluxData().potential(wPhaseIdx, isIndexI);
        Scalar potentialNWIJ = cellDataI.fluxData().potential(nPhaseIdx, isIndexI);
        Scalar potentialWIK = potentialWIJ;
        Scalar potentialNWIK = potentialNWIJ;
        // preliminary potential. The "real" ones are found below

        // Comment: reversed weighting is plausible, too (swap lJ and lI)
        Scalar rhoMeanWIJ = density_[wPhaseIdx];
        Scalar rhoMeanNWIJ = density_[nPhaseIdx];
        Scalar rhoMeanWIK = density_[wPhaseIdx];
        Scalar rhoMeanNWIK = density_[nPhaseIdx];

        if (compressibility_)
        {
            rhoMeanWIJ = (lJ * cellDataI.density(wPhaseIdx) + lI * cellDataJ.density(wPhaseIdx)) / l;
            rhoMeanNWIJ = (lJ * cellDataI.density(nPhaseIdx) + lI * cellDataJ.density(nPhaseIdx)) / l;
            rhoMeanWIK = (lJ * cellDataI.density(wPhaseIdx) + lI * cellDataK.density(wPhaseIdx)) / l;
            rhoMeanNWIK = (lJ * cellDataI.density(nPhaseIdx) + lI * cellDataK.density(nPhaseIdx)) / l;
        }

        Scalar fMeanWIJ = (lJ * fractionalWI + lI * fractionalWJ) / l;
        Scalar fMeanNWIJ = (lJ * fractionalNWI + lI * fractionalNWJ) / l;
        Scalar fMeanWIK = (lJ * fractionalWI + lI * fractionalWK) / l;
        Scalar fMeanNWIK = (lJ * fractionalNWI + lI * fractionalNWK) / l;

        Scalar densityWIJ = density_[wPhaseIdx];
        Scalar densityNWIJ = density_[nPhaseIdx];
        Scalar densityWIK = density_[wPhaseIdx];
        Scalar densityNWIK = density_[nPhaseIdx];

        if (compressibility_)
        {
        // Upwinding for finding the upwinding direction
            densityWIJ = (potentialWIJ > 0.) ? cellDataI.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
            densityNWIJ = (potentialNWIJ > 0.) ? cellDataI.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);
            densityWIK = (potentialWIK > 0.) ? cellDataI.density(wPhaseIdx) : cellDataK.density(wPhaseIdx);
            densityNWIK = (potentialNWIK > 0.) ? cellDataI.density(nPhaseIdx) : cellDataK.density(nPhaseIdx);

            densityWIJ = (potentialWIJ == 0.) ? rhoMeanWIJ : densityWIJ;
            densityNWIJ = (potentialNWIJ == 0.) ? rhoMeanNWIJ : densityNWIJ;
            densityWIK = (potentialWIK == 0.) ? rhoMeanWIK : densityWIK;
            densityNWIK = (potentialNWIK == 0.) ? rhoMeanNWIK : densityNWIK;
        }

        Scalar fractionalWIJ = (potentialWIJ > 0.) ? fractionalWI : fractionalWJ;
        Scalar fractionalNWIJ = (potentialNWIJ > 0.) ? fractionalNWI : fractionalNWJ;
        Scalar fractionalWIK = (potentialWIK > 0.) ? fractionalWI : fractionalWK;
        Scalar fractionalNWIK = (potentialNWIK > 0.) ? fractionalNWI : fractionalNWK;

        fractionalWIJ = (potentialWIJ == 0.) ? fMeanWIJ : fractionalWIJ;
        fractionalNWIJ = (potentialNWIJ == 0.) ? fMeanNWIJ : fractionalNWIJ;
        fractionalWIK = (potentialWIK == 0.) ? fMeanWIK : fractionalWIK;
        fractionalNWIK = (potentialNWIK == 0.) ? fMeanNWIK : fractionalNWIK;

        switch (pressureType_)
        {
        case pglobal:
        {
            Scalar pressJK = (cellDataJ.globalPressure() + cellDataK.globalPressure()) / 2;

            potentialWIJ = (cellDataI.globalPressure() - fMeanNWIJ * pcI - (pressJK - fMeanNWIJ * pcJK)) / l
                    + (densityWIJ - lJ / l * (kI + kK) / kI * (densityWIK - densityWIJ) / 2) * ng;
            potentialNWIJ = (cellDataI.globalPressure() + fMeanWIJ * pcI - (pressJK + fMeanWIJ * pcJK)) / l
                    + (densityNWIJ - lJ / l * (kI + kK) / kI * (densityNWIK - densityNWIJ) / 2) * ng;
            potentialWIK = (cellDataI.globalPressure() - fMeanNWIK * pcI - (pressJK - fMeanNWIK * pcJK)) / l
                    + (densityWIK - lJ / l * (kI + kK) / kI * (densityWIJ - densityWIK) / 2) * ng;
            potentialNWIK = (cellDataI.globalPressure() + fMeanWIK * pcI - (pressJK + fMeanWIK * pcJK)) / l
                    + (densityNWIK - lJ / l * (kI + kK) / kI * (densityNWIJ - densityNWIK) / 2) * ng;
            break;
        }
        default:
        {
            potentialWIJ = (cellDataI.pressure(wPhaseIdx)
                    - 0.5 * (cellDataJ.pressure(wPhaseIdx) + cellDataK.pressure(wPhaseIdx))) / l
                    + (densityWIJ - lJ / l * (kI + kK) / kI * (densityWIK - densityWIJ) / 2) * ng;
            potentialNWIJ = (cellDataI.pressure(nPhaseIdx)
                    - (0.5 * (cellDataJ.pressure(nPhaseIdx) + cellDataK.pressure(nPhaseIdx)))) / l
                    + (densityNWIJ - lJ / l * (kI + kK) / kI * (densityNWIK - densityNWIJ) / 2) * ng;
            potentialWIK = (cellDataI.pressure(wPhaseIdx)
                    - 0.5 * (cellDataJ.pressure(wPhaseIdx) + cellDataK.pressure(wPhaseIdx))) / l
                    + (densityWIK - lJ / l * (kI + kK) / kI * (densityWIJ - densityWIK) / 2) * ng;
            potentialNWIK = (cellDataI.pressure(nPhaseIdx)
                    - (0.5 * (cellDataJ.pressure(nPhaseIdx) + cellDataK.pressure(nPhaseIdx)))) / l
                    + (densityNWIK - lJ / l * (kI + kK) / kI * (densityNWIJ - densityNWIK) / 2) * ng;
            break;
        }
        }

        //store potentials for further calculations (velocity, saturation, ...)
        // these quantities only have correct sign (needed for upwinding)
        // potentials are defined slightly different for adaptive scheme
        cellDataI.fluxData().addPotential(wPhaseIdx, isIndexI, potentialWIJ);
        cellDataI.fluxData().addPotential(nPhaseIdx, isIndexI, potentialNWIJ);
        cellDataJ.fluxData().setPotential(wPhaseIdx, isIndexJ, -potentialWIJ);
        cellDataJ.fluxData().setPotential(nPhaseIdx, isIndexJ, -potentialNWIJ);

        //do the upwinding of the mobility depending on the phase potentials
        Scalar lambdaWIJ = (potentialWIJ > 0.) ? lambdaWI : lambdaWJ;
        lambdaWIJ = (potentialWIJ == 0) ? 0.5 * (lambdaWI + lambdaWJ) : lambdaWIJ;
        Scalar lambdaNWIJ = (potentialNWIJ > 0.) ? lambdaNWI : lambdaNWJ;
        lambdaNWIJ = (potentialNWIJ == 0) ? 0.5 * (lambdaNWI + lambdaNWJ) : lambdaNWIJ;

        if (compressibility_)
        {
            densityWIJ = (potentialWIJ > 0.) ? cellDataI.density(wPhaseIdx) : cellDataJ.density(wPhaseIdx);
            densityNWIJ = (potentialNWIJ > 0.) ? cellDataI.density(nPhaseIdx) : cellDataJ.density(nPhaseIdx);
            densityWIJ = (potentialWIJ == 0) ? rhoMeanWIJ : densityWIJ;
            densityNWIJ = (potentialNWIJ == 0) ? rhoMeanNWIJ : densityNWIJ;
            densityWIK = (potentialWIK > 0.) ? cellDataI.density(wPhaseIdx) : cellDataK.density(wPhaseIdx);
            densityNWIK = (potentialNWIK > 0.) ? cellDataI.density(nPhaseIdx) : cellDataK.density(nPhaseIdx);
            densityWIK = (potentialWIK == 0) ? rhoMeanWIK : densityWIK;
            densityNWIK = (potentialNWIK == 0) ? rhoMeanNWIK : densityNWIK;
        }

        //calculate velocities and the gravity term
        Dune::FieldVector < Scalar, dimWorld > velocityW(unitOuterNormal);
        Dune::FieldVector < Scalar, dimWorld > velocityNW(unitOuterNormal);
        Dune::FieldVector < Scalar, dimWorld > gravityTermW(unitOuterNormal);
        Dune::FieldVector < Scalar, dimWorld > gravityTermNW(unitOuterNormal);

        gravityTermW *= lambdaWIJ * kMean * ng;
        gravityTermW *= densityWIJ - (lJ / l) * (kI + kK) / kI * (densityWIK - densityWIJ) / 2;
        gravityTermNW *= lambdaNWIJ * kMean * ng;
        gravityTermNW *= densityNWIJ - (lJ / l) * (kI + kK) / kI * (densityNWIK - densityNWIJ) / 2;

        switch (pressureType_)
        {
        case pw:
        {
            velocityW *= lambdaWIJ * kMean * (cellDataI.pressure(wPhaseIdx) - (cellDataJ.pressure(wPhaseIdx) + cellDataK.pressure(wPhaseIdx)) / 2.0) / l;
            velocityNW *= lambdaNWIJ * kMean * (cellDataI.pressure(nPhaseIdx) - (cellDataJ.pressure(nPhaseIdx) + cellDataK.pressure(nPhaseIdx)) / 2.0) / l;

            velocityW += gravityTermW;
            velocityNW += gravityTermNW;
            break;
        }
        case pn:
        {
            velocityNW *= lambdaNWIJ * kMean * (cellDataI.pressure(nPhaseIdx) - (cellDataJ.pressure(nPhaseIdx) + cellDataK.pressure(nPhaseIdx)) / 2.0) / l;
            velocityW *= lambdaWIJ * kMean * (cellDataI.pressure(wPhaseIdx) - (cellDataJ.pressure(wPhaseIdx) + cellDataK.pressure(wPhaseIdx)) / 2.0) / l;

            velocityW += gravityTermW;
            velocityNW += gravityTermNW;
            break;
        }
        case pglobal:
        {
            Scalar pressJK = (cellDataJ.globalPressure() + cellDataK.globalPressure()) / 2;

            velocityW *= lambdaWIJ * kMean * (cellDataI.globalPressure() - pressJK) / l;
            velocityW += gravityTermW;
            velocityW += gravityTermNW;
            velocityNW = 0;
            break;
        }
        }

        cellDataJ.fluxData().setVelocity(wPhaseIdx, isIndexJ, velocityW);
        cellDataJ.fluxData().setVelocity(nPhaseIdx, isIndexJ, velocityNW);
        cellDataJ.fluxData().setVelocityMarker(isIndexJ);

        //times 0.5 because cell face with hanging node is called twice! Do not set marker because it should be called twice!
        velocityW *= 0.5;
        velocityNW *= 0.5;
        cellDataI.fluxData().addVelocity(wPhaseIdx, isIndexI, velocityW);
        cellDataI.fluxData().addVelocity(nPhaseIdx, isIndexI, velocityNW);
    }

    return;
}
}
#endif
