// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadFlux
 * \copydoc Dumux::BedloadFlux
 */
#ifndef DUMUX_FLUX_BEDLOAD_HH
#define DUMUX_FLUX_BEDLOAD_HH

#include "extendvector.hh"
#include "locallaxfriedrichs.hh"
#include "rotatevector.hh"
#include "slopeeffect.hh"

namespace Dumux {

/*!
 * \ingroup BedloadFlux
 * \brief Contains all quantities needed to define and solve the Riemann state
 */
template<typename Scalar>
struct RiemannState {
    Scalar bedloadDischargeLeft;
    Scalar bedloadDischargeRight;
    Scalar bedSurfaceLeft;
    Scalar bedSurfaceRight;
    Scalar bottomActiveLayerLeft;
    Scalar bottomActiveLayerRight;
    Scalar cellAreaLeft;
    Scalar cellAreaRight;
    Scalar grainDensity;
    Scalar gravity;
    Scalar lowerMobilityLimitLeft;
    Scalar lowerMobilityLimitRight;
    Scalar massFractionLeft;
    Scalar massFractionRight;
    Scalar porosity;
    Scalar sedimentMassLeft;
    Scalar sedimentMassRight;
    Scalar upperMobilityLimitLeft;
    Scalar upperMobilityLimitRight;
    Scalar velocityLeft;
    Scalar velocityRight;
    Scalar velocityScalarLeft;
    Scalar velocityScalarRight;
    Scalar waterDepthLeft;
    Scalar waterDepthRight;
};

/*!
 * \ingroup BedloadFlux
 * \brief Compute the bedload flux at a cell face by solving a Riemann problem.
 */
template<class NumEqVector>
class BedloadFlux
{

public:

    using Cache = FluxVariablesCaching::EmptyAdvectionCache;
    using CacheFiller = FluxVariablesCaching::EmptyCacheFiller;
    using Scalar = typename NumEqVector::field_type;

    /*!
     * \brief Prepare the Riemann problem and solve it by the local Lax-Friedrichs solver.
     *
     * \param problem The object specifying the problem which ought to be simulate
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables associated with the element stencil
     * \param scvf The sub control volume face
     *
     * \return The bedload flux \f$[kg/s]\f$. The first three components are zero, since they are related to the shallow water flux.
     */
    template<class Problem, class FVElementGeometry, class ElementVolumeVariables>
    static NumEqVector flux(const Problem& problem,
                            const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const typename FVElementGeometry::SubControlVolumeFace& scvf)
    {
        using std::max;
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& nxy = scvf.unitOuterNormal();

        // prepare the quantities needed for the slope effect
        Dune::FieldVector<Scalar, 2> bedShearStressLeft = Sediment::getExtendedVector(insideVolVars.bottomShearStress(), insideVolVars.tanAngleSecondaryCurrents());
        Dune::FieldVector<Scalar, 2> bedShearStressRight = Sediment::getExtendedVector(outsideVolVars.bottomShearStress(), outsideVolVars.tanAngleSecondaryCurrents());

        // compute the bed surface gradient
        const auto cellCenterToCellCenter = outsideScv.center() - insideScv.center();
        const Scalar distance = cellCenterToCellCenter.two_norm();
        const auto cellCenterToCellCenterNormed = cellCenterToCellCenter / distance;
        const Dune::FieldVector<Scalar, 2> bedSlope = (outsideVolVars.bedSurface() - insideVolVars.bedSurface()) / distance * cellCenterToCellCenterNormed;

        // get the flow velocity, modified by the secondary currents angle
        Dune::FieldVector<Scalar, 2> velocityLeft = Sediment::getExtendedVector(insideVolVars.velocity(), insideVolVars.tanAngleSecondaryCurrents());
        Dune::FieldVector<Scalar, 2> velocityRight = Sediment::getExtendedVector(outsideVolVars.velocity(), outsideVolVars.tanAngleSecondaryCurrents());

        // fill the riemann state
        RiemannState<Scalar> riemannState;
        riemannState.waterDepthLeft = insideVolVars.waterDepth();
        riemannState.waterDepthRight = outsideVolVars.waterDepth();
        riemannState.bedSurfaceLeft = insideVolVars.bedSurface();
        riemannState.bedSurfaceRight = outsideVolVars.bedSurface();
        riemannState.gravity = problem.spatialParams().gravity();
        riemannState.porosity = problem.spatialParams().porosity();
        riemannState.cellAreaLeft = insideScv.volume();
        riemannState.cellAreaRight = outsideScv.volume();
        riemannState.bottomActiveLayerLeft = insideVolVars.bottomActiveLayer();
        riemannState.bottomActiveLayerRight = outsideVolVars.bottomActiveLayer();
        Scalar fluxLimiterUpperSedimentThickness = problem.spatialParams().fluxLimiterUpperSedimentThickness();
        Scalar fluxLimiterLowerSedimentThickness = problem.spatialParams().fluxLimiterLowerSedimentThickness();

        NumEqVector localFlux(0.0);
        static const int nGrainClasses = getParam<int>("Sediment.NumberGrainClasses");
        for (int i = 0; i < nGrainClasses; i++)
        {
            // calculate the deviation of the bedload discharge due to the bed slope
            Scalar angleSlopeEffectLeft = calculateDirectionCorrection(bedSlope, bedShearStressLeft, velocityLeft, problem.spatialParams().representativeGrainDiameter(i));
            Scalar angleSlopeEffectRight = calculateDirectionCorrection(bedSlope, bedShearStressRight, velocityRight, problem.spatialParams().representativeGrainDiameter(i));

            // rotate the bedload discharge by the angle caused by slope effect.
            auto bedloadDischargeLeft = Sediment::getRotatedVector(insideVolVars.bedloadDischarge(i), angleSlopeEffectLeft);
            auto bedloadDischargeRight = Sediment::getRotatedVector(outsideVolVars.bedloadDischarge(i), angleSlopeEffectRight);

            // correct the magnitude of the bedload discharge due to the bed slope
            correctBedloadMagnitude(bedloadDischargeLeft, bedSlope);
            correctBedloadMagnitude(bedloadDischargeRight, bedSlope);

            // rotate the velocity by the angle caused by slope effect
            Dune::FieldVector<Scalar, 2> velocityLeftRotated = Sediment::getRotatedVector(velocityLeft, angleSlopeEffectLeft);
            Dune::FieldVector<Scalar, 2> velocityRightRotated = Sediment::getRotatedVector(velocityRight, angleSlopeEffectRight);

            riemannState.grainDensity = problem.spatialParams().grainDensity(i);
            riemannState.velocityScalarLeft = velocityLeftRotated.two_norm();
            riemannState.velocityScalarRight = velocityRightRotated.two_norm();
            riemannState.velocityLeft = velocityLeftRotated[0] * nxy[0] + velocityLeftRotated[1] * nxy[1];
            riemannState.velocityRight = velocityRightRotated[0] * nxy[0] + velocityRightRotated[1] * nxy[1];
            riemannState.sedimentMassLeft = insideVolVars.sedimentMass(i);
            riemannState.sedimentMassRight = outsideVolVars.sedimentMass(i);
            riemannState.massFractionLeft = insideVolVars.massFraction(i);
            riemannState.massFractionRight = outsideVolVars.massFraction(i);
            riemannState.bedloadDischargeLeft = (bedloadDischargeLeft[0] * nxy[0] + bedloadDischargeLeft[1] * nxy[1]) * insideVolVars.massFraction(i);
            riemannState.bedloadDischargeRight = (bedloadDischargeRight[0] * nxy[0] + bedloadDischargeRight[1] * nxy[1]) * outsideVolVars.massFraction(i);
            // mobilityLimit has the unit [kg/m^2]
            const Scalar upperMobilityLimit = fluxLimiterUpperSedimentThickness * (1 - problem.spatialParams().porosity()) * problem.spatialParams().grainDensity(i);
            const Scalar lowerMobilityLimit = fluxLimiterLowerSedimentThickness * (1 - problem.spatialParams().porosity()) * problem.spatialParams().grainDensity(i);
            // Scale the flux limiter values for the limiting with the mass mass fraction. To avoid very small values the scaling is limited by 1/nGrainClasses
            riemannState.upperMobilityLimitLeft = upperMobilityLimit * max(insideVolVars.massFraction(i), 1.0/nGrainClasses);
            riemannState.lowerMobilityLimitLeft = lowerMobilityLimit * max(insideVolVars.massFraction(i), 1.0/nGrainClasses);
            riemannState.upperMobilityLimitRight = upperMobilityLimit * max(outsideVolVars.massFraction(i), 1.0/nGrainClasses);
            riemannState.lowerMobilityLimitRight = lowerMobilityLimit * max(outsideVolVars.massFraction(i), 1.0/nGrainClasses);
            auto riemannFlux = Sediment::localLaxFriedrichs<Scalar>(riemannState);
            localFlux[i] = riemannFlux * scvf.area() * problem.spatialParams().grainDensity(i);
        }

        return localFlux;
    }
};

} // end namespace Dumux

#endif
