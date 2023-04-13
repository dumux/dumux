// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \ingroup PoreNetworkModels
 * \brief Implementation of the single-phase transmissibility laws for throats
 */
#ifndef DUMUX_PNM_THROAT_TRANSMISSIBILITY_1P_HH
#define DUMUX_PNM_THROAT_TRANSMISSIBILITY_1P_HH

#include <dumux/common/parameters.hh>
#include <dumux/porenetwork/common/throatproperties.hh>
#include "emptycache.hh"

namespace Dumux::PoreNetwork {

/*!
 * \ingroup Fluidmatrixinteractions
 * \ingroup PoreNetworkModels
 * \brief Collection of single-phase flow throat transmissibilities based on
 *        Bruus, H. (2011). Acoustofluidics 1: Governing equations in microfluidics. Lab on a Chip, 11(22), 3742-3751.
 *        https://backend.orbit.dtu.dk/ws/portalfiles/portal/5900070/rsc%5B1%5D.pdf
 */
template<class Scalar>
class TransmissibilityBruus
{
public:

    using SinglePhaseCache = EmptyCache;

    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static Scalar singlePhaseTransmissibility(const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                              const ElementVolumeVariables& elemVolVars,
                                              const FluxVariablesCache& fluxVarsCache,
                                              const int phaseIdx)
    {
        const auto shape = fluxVarsCache.throatCrossSectionShape();
        const Scalar throatInscribedRadius = fluxVarsCache.throatInscribedRadius();
        const Scalar throatLength = fluxVarsCache.throatLength();
        const Scalar area = fluxVarsCache.throatCrossSectionalArea();
        return singlePhaseTransmissibility(shape, throatInscribedRadius, throatLength, area);
    }

    //! Returns the conductivity of a throat when only one phase is present.
    static Scalar singlePhaseTransmissibility(const Throat::Shape shape,
                                              const Scalar inscribedRadius,
                                              const Scalar throatLength,
                                              const Scalar area)
    { return 1.0 / rHydThroat_(shape, inscribedRadius, throatLength, area); }

protected:

    static Scalar rHydThroat_(const Throat::Shape shape,
                              const Scalar radius,
                              const Scalar length,
                              const Scalar area)
    {
        switch(shape)
        {
            case Throat::Shape::square:
            {
                const Scalar sideLength = 2.0*radius;
                return 28.4*length * 1.0/(sideLength*sideLength*sideLength*sideLength);
            }
            case Throat::Shape::circle:
            {
                return 8.0/M_PI * length * 1.0/(radius*radius*radius*radius);
            }
            case Throat::Shape::equilateralTriangle:
            {
                using std::sqrt;
                static Scalar sqrt3 = sqrt(3.0);
                const Scalar sideLength = 6.0/sqrt3 * radius;
                return 320.0/sqrt3 * length * 1.0/(sideLength*sideLength*sideLength*sideLength);
            }
            case Throat::Shape::twoPlates:
            {
                // the distance between the two parallel plates
                const Scalar width = 2*radius;
                return 12.0/(width*width*width) * length;
            }
            case Throat::Shape::rectangle:
            {
                // requires width >> height for good accuracy
                Scalar height = 2.0*radius;
                Scalar width = area/height;

                using std::swap;
                if (width < height)
                    swap(width, height);

                return 12.0*length / (1.0 - 0.63*(height/width)) * 1.0/(height*height*height*width);
            }
            default: DUNE_THROW(Dune::InvalidStateException, "Throat geometry not supported");
        }
    }
};

/*!
 * \ingroup Fluidmatrixinteractions
 * \ingroup PoreNetworkModels
 * \brief Single-phase flow throat transmissibility based on Patzek & Silin (2001) https://doi.org/10.1006/jcis.2000.7413
 */
template<class Scalar, bool considerPoreResistance = false, bool interpolateK = false>
class TransmissibilityPatzekSilin
{
    static_assert(!interpolateK,  "Interpolation of k not implemented");
public:

    using SinglePhaseCache = EmptyCache;

    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static Scalar singlePhaseTransmissibility(const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                              const ElementVolumeVariables& elemVolVars,
                                              const FluxVariablesCache& fluxVarsCache,
                                              const int phaseIdx)
    {
        assert(fluxVarsCache.throatCrossSectionShape() != Throat::Shape::twoPlates && "TwoPlates not supported. Use TransmissibilityBruus instead!");

        const Scalar shapeFactor = fluxVarsCache.throatShapeFactor();
        const Scalar area = fluxVarsCache.throatCrossSectionalArea();
        const Scalar throatLength = fluxVarsCache.throatLength();
        const Scalar throatTransmissibility = singlePhaseTransmissibility(shapeFactor, throatLength, area);

        if constexpr (!considerPoreResistance)
            return throatTransmissibility;
        else
        {
            static const bool considerPoreResistanceOnRuntime = getParamFromGroup<bool>(problem.paramGroup(), "Transmissibility.ConsiderPoreResistance", true);
            if (!considerPoreResistanceOnRuntime)
                return throatTransmissibility;

            const auto& scv0 = fvGeometry.scv(scvf.insideScvIdx());
            const auto& scv1 = fvGeometry.scv(scvf.outsideScvIdx());

            const auto& spatialParams = problem.spatialParams();
            const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);

            // TODO maybe include this in fluxVarsCache if this is general enough
            const Scalar poreLength0 = spatialParams.poreLength(element, scv0, elemSol);
            const Scalar poreLength1 = spatialParams.poreLength(element, scv1, elemSol);

            const Scalar poreShapeFactor0 = spatialParams.poreShapeFactor(element, scv0, elemSol);
            const Scalar poreShapeFactor1 = spatialParams.poreShapeFactor(element, scv1, elemSol);

            const Scalar poreCrossSectionalArea0 = spatialParams.poreCrossSectionalArea(element, scv0, elemSol);
            const Scalar poreCrossSectionalArea1 = spatialParams.poreCrossSectionalArea(element, scv1, elemSol);

            const Scalar poreTransmissibility0 = singlePhaseTransmissibility(poreShapeFactor0, poreLength0, poreCrossSectionalArea0);
            const Scalar poreTransmissibility1 = singlePhaseTransmissibility(poreShapeFactor1, poreLength1, poreCrossSectionalArea1);

            return 1 / (1.0/throatTransmissibility + 1.0/poreTransmissibility0 + 1.0/poreTransmissibility1);
        }
    }

    //! Returns the conductivity of a throat when only one phase is present. See Patzek & Silin (2001)
    static Scalar singlePhaseTransmissibility(const Scalar shapeFactor,
                                              const Scalar length,
                                              const Scalar area)
    {
        const Scalar k = k_(shapeFactor);
        return k * area*area * shapeFactor / length;
    }

private:
    static Scalar k_(const Scalar shapeFactor)
    {
        if (shapeFactor <= Throat::shapeFactorEquilateralTriangle<Scalar>())
            return 0.6; // == 3/5
        else if (shapeFactor <= Throat::shapeFactorSquare<Scalar>())
            return 0.5623;
        else // circle
            return 0.5;
    }

    // TODO interpolation
};


//! Used by Joeakar-Niasar
template<class Scalar>
class TransmissibilityAzzamDullien
{
public:

    using SinglePhaseCache = EmptyCache;

    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static Scalar singlePhaseTransmissibility(const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                              const ElementVolumeVariables& elemVolVars,
                                              const FluxVariablesCache& fluxVarsCache,
                                              const int phaseIdx)
    {
        const Scalar throatInscribedRadius = fluxVarsCache.throatInscribedRadius();
        const Scalar throatLength = fluxVarsCache.throatLength();
        return singlePhaseTransmissibility(throatInscribedRadius, throatLength);
    }

    //! Returns the conductivity of a throat when only one phase is present.
    static Scalar singlePhaseTransmissibility(const Scalar radius,
                                              const Scalar length)
    {
        const Scalar rEff= std::sqrt(4.0/M_PI)*radius ;
        return M_PI/(8.0*length) *rEff*rEff*rEff*rEff ;
    }
};

} // end namespace Dumux::Porenetwork

#endif
