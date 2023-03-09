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
 * \brief This file contains the data which is required to calculate
 *        diffusive heat fluxes with Fourier's law.
 */
#ifndef DUMUX_FLUX_PNM_GRAIN_FOURIERS_LAW_HH
#define DUMUX_FLUX_PNM_GRAIN_FOURIERS_LAW_HH

#include <dumux/common/math.hh>

namespace Dumux::PoreNetwork {

/*!
 * \brief Specialization of Fourier's Law for the pore-network SOLID model.
 * \note See Koch et al (2021) https://doi.org/10.1007/s11242-021-01602-5
 *       and Khan et al (2019) https://doi.org/10.1016/j.compchemeng.2018.12.025 
 */
template <class Scalar>
struct TruncatedPyramidGrainFouriersLaw
{
    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class ElementFluxVariablesCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const typename FVElementGeometry::SubControlVolumeFace& scvf,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        const Scalar topSideLength = 2.0*fluxVarsCache.throatRadius(); // maybe contact radius?

        // We assume that the distance between pore centroid and throat
        // centroid (i.e., the height of the pyramid) equals the inscribed pore radius.
        const Scalar insideHeight = insideVolVars.poreRadius();
        const Scalar outsideHeight = outsideVolVars.poreRadius();

        auto getPyramidBaseLengthFromVolume = [&](const Scalar v, const Scalar h)
        {
            const Scalar b = topSideLength;
            using std::sqrt;
            return 0.5*sqrt(3.0) * sqrt(-(b*b*h-4.0*v)/h) -0.5*b;
        };

        // the pyramid base length of the inside pore
        const Scalar insideBaseSideLength = [&]()
        {
            static const bool useAdaptedVolume = getParamFromGroup<bool>(problem.paramGroup(), "GrainFouriersLaw.UseAdaptedVolumeForPyramid", false);

            if (useAdaptedVolume)
                return getPyramidBaseLengthFromVolume(0.5*insideVolVars.poreVolume(), insideHeight);
            else
                return 2.0 * insideVolVars.poreRadius();
        }();

        // the pyramid base length of the outside pore
        const Scalar outsideBaseSideLength = [&]()
        {
            static const bool useAdaptedVolume = getParamFromGroup<bool>(problem.paramGroup(), "GrainFouriersLaw.UseAdaptedVolumeForPyramid", false);

            if (useAdaptedVolume)
                return getPyramidBaseLengthFromVolume(0.5*outsideVolVars.poreVolume(), outsideHeight);
            else
                return 2.0 * outsideVolVars.poreRadius();
        }();

        auto insideThermalConducitivity = insideVolVars.solidThermalConductivity();
        auto outsideThermalConducitivity = outsideVolVars.solidThermalConductivity();

        const Scalar gInside = 4.0*insideThermalConducitivity *  0.5*topSideLength * 0.5*insideBaseSideLength / insideHeight;
        const Scalar gOutside = 4.0*outsideThermalConducitivity *  0.5*topSideLength * 0.5*outsideBaseSideLength / outsideHeight;
        const Scalar gThroat = Dumux::harmonicMean(insideThermalConducitivity, outsideThermalConducitivity)
                             * fluxVarsCache.grainContactArea() / fluxVarsCache.throatLength();

        const Scalar g = 1.0 / (1.0/gInside + 1.0/gOutside + 1.0/gThroat);

        // calculate the temperature gradient
        const Scalar deltaT = insideVolVars.temperature() - outsideVolVars.temperature();

        return deltaT * g;
    }
};

/*!
 * \brief Specialization of Fourier's Law for the pore-network SOLID model.
 * \note See Koch et al (2021) https://doi.org/10.1007/s11242-021-01602-5
 */
template <class Scalar>
struct SphereCapGrainFouriersLaw
{
    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class ElementFluxVariablesCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const typename FVElementGeometry::SubControlVolumeFace& scvf,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        auto gSphereCap = [&](const auto& scv, const auto& volVars)
        {
            const Scalar R = problem.spatialParams.extendedPoreRadius(scv.dofIndex());
            const Scalar lambda = volVars.solidThermalConductivity();
            const Scalar rC = volVars.poreRadius();

            return (lambda*M_PI*R) / std::atanh(rC/R);
        };

        auto insideThermalConducitivity = insideVolVars.solidThermalConductivity();
        auto outsideThermalConducitivity = outsideVolVars.solidThermalConductivity();

        const Scalar gInside = gSphereCap(insideScv, insideVolVars);
        const Scalar gOutside = gSphereCap(outsideScv, outsideVolVars);
        const Scalar gThroat = Dumux::harmonicMean(insideThermalConducitivity, outsideThermalConducitivity)
                             * fluxVarsCache.grainContactArea() / fluxVarsCache.throatLength();

        const Scalar g = 1.0 / (1.0/gInside + 1.0/gOutside + 1.0/gThroat);

        // calculate the temperature gradient
        const Scalar deltaT = insideVolVars.temperature() - outsideVolVars.temperature();

        return deltaT * g;
    }

};

} // end namespace Dumux::PoreNetwork

#endif
