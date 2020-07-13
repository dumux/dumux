// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Implementation of the single-phase transmissibility laws for throats
 */
#ifndef DUMUX_PNM_THROAT_TRANSMISSIBILITY_1P_HH
#define DUMUX_PNM_THROAT_TRANSMISSIBILITY_1P_HH

#include <dumux/porenetworkflow/common/throatproperties.hh>
#include "emptycache.hh"

namespace Dumux {

/*!
 * \ingroup PoreNetworkModel
 * \brief Collection of single-phase flow throat transmissibilities based on
 *        Bruus, H. (2011). Acoustofluidics 1: Governing equations in microfluidics. Lab on a Chip, 11(22), 3742-3751.
 *        https://backend.orbit.dtu.dk/ws/portalfiles/portal/5900070/rsc%5B1%5D.pdf
 */
template<class Scalar>
class TransmissibilityBruus
{
public:

    using Cache = EmptyCache;

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
        const Scalar throatRadius = fluxVarsCache.throatRadius();
        const Scalar throatLength = fluxVarsCache.throatLength();
        const Scalar area = fluxVarsCache.throatCrossSectionalArea();
        return singlePhaseTransmissibility(shape, throatRadius, throatLength, area);
    }

     //! Returns the conductivity of a throat when only one phase is present.
    template<class Element, class FVElementGeometry, class FluxVariablesCache>
    [[deprecated("Use singlePhaseTransmissibility(shape, inscribedRadius, throatLength, area) instead. Will be removed soon.")]]
    static Scalar singlePhaseTransmissibility(const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                              const FluxVariablesCache& fluxVarsCache)
    {
        const auto shape = fluxVarsCache.shape();
        const Scalar inscribedRadius = fluxVarsCache.throatRadius();
        const Scalar throatLength = fluxVarsCache.throatLength();
        const Scalar area = fluxVarsCache.throatCrossSectionalArea();
        return singlePhaseTransmissibility(shape, inscribedRadius, throatLength, area);
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
                static constexpr Scalar sqrt3 = sqrt(3.0);
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

template<class Scalar, bool interpolateK = false>
class TransmissibilityPatzekSilin
{
    static_assert(!interpolateK,  "Interpolation of k not implemented");
public:

    using Cache = EmptyCache;

    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static Scalar singlePhaseTransmissibility(const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                              const ElementVolumeVariables& elemVolVars,
                                              const FluxVariablesCache& fluxVarsCache,
                                              const int phaseIdx)
    {
        const auto shapeFactor = fluxVarsCache.throatShapeFactor();
        const Scalar area = fluxVarsCache.throatCrossSectionalArea();
        const Scalar throatLength = fluxVarsCache.throatLength();
        return singlePhaseTransmissibility(shapeFactor, throatLength, area);
    }

    //! Returns the conductivity of a throat when only one phase is present. See Patzek & Silin (2001)
   template<class Element, class FVElementGeometry, class FluxVariablesCache>
   [[deprecated("Use singlePhaseTransmissibility(shapeFactor, throatLength, area) instead. Will be removed soon.")]]
   static Scalar singlePhaseTransmissibility(const Element& element,
                                             const FVElementGeometry& fvGeometry,
                                             const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                             const FluxVariablesCache& fluxVarsCache)
   {
       const auto shapeFactor = fluxVarsCache.shapeFactor();
       const Scalar area = fluxVarsCache.throatCrossSectionalArea();
       const Scalar throatLength = fluxVarsCache.throatLength();
       return singlePhaseTransmissibility(shapeFactor, throatLength, area);
   }

     //! Returns the conductivity of a throat when only one phase is present. See Patzek & Silin (2001)
    static Scalar singlePhaseTransmissibility(const Scalar shapeFactor,
                                              const Scalar throatLength,
                                              const Scalar area)
    {
        const Scalar k = k_(shapeFactor);
        return k * area*area * shapeFactor / throatLength;
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


//! Used by Joeakar-Niasar, probably wrong, TODO: check and maybe remove
template<class Scalar>
class TransmissibilityAzzamDullien
{
public:

    using Cache = EmptyCache;

    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class FluxVariablesCache>
    static Scalar singlePhaseTransmissibility(const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                              const ElementVolumeVariables& elemVolVars,
                                              const FluxVariablesCache& fluxVarsCache,
                                              const int phaseIdx)
    {
        const Scalar throatRadius = fluxVarsCache.throatRadius();
        const Scalar throatLength = fluxVarsCache.throatLength();
        return singlePhaseTransmissibility(throatRadius, throatLength);
    }

    //! Returns the conductivity of a throat when only one phase is present.
    template<class Element, class FVElementGeometry,  class FluxVariablesCache>
    [[deprecated("Use singlePhaseTransmissibility(shapeFactor, inscribedRadius, throatLength, area) instead. Will be removed soon.")]]
    static Scalar singlePhaseTransmissibility(const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                              const FluxVariablesCache& fluxVarsCache)
    {
        const Scalar throatRadius = fluxVarsCache.throatRadius();
        const Scalar throatLength = fluxVarsCache.throatLength();
        return singlePhaseTransmissibility(throatRadius, throatLength);
    }

    //! Returns the conductivity of a throat when only one phase is present.
    static Scalar singlePhaseTransmissibility(const Scalar radius,
                                              const Scalar length)
    {
        const Scalar rEff= std::sqrt(4.0/M_PI)*radius ;
        return M_PI/(8.0*length) *rEff*rEff*rEff*rEff ;
    }
};

}

#endif // DUMUX_PNM_THROAT_TRANSMISSIBILITY_1P_HH
