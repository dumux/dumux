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
 * \ingroup StaggeredFlux
 * \brief This file contains the data which is required to calculate
 *        diffusive molar fluxes due to molecular diffusion with Fick's law.
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FICKS_LAW_HH

#include <numeric>
#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/fluxvariablescaching.hh>
#include <dumux/flux/fickiandiffusioncoefficients.hh>
#include <dumux/flux/referencesystemformulation.hh>


namespace Dumux {

// forward declaration
template<class TypeTag, DiscretizationMethod discMethod, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation;

/*!
 * \ingroup StaggeredFlux
 * \brief Specialization of Fick's Law for the staggered free flow method.
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem>
class FicksLawImplementation<TypeTag, DiscretizationMethod::staggered, referenceSystem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    static constexpr int numComponents = ModelTraits::numFluidComponents();
    static constexpr int numPhases = ModelTraits::numFluidPhases();

    using NumEqVector = Dune::FieldVector<Scalar, numComponents>;

    static_assert(ModelTraits::numFluidPhases() == 1, "Only one phase supported!");

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::staggered;
    //return the reference system
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    //! state the type for the corresponding cache
    //! We don't cache anything for this law
    using Cache = FluxVariablesCaching::EmptyDiffusionCache;

    using DiffusionCoefficientsContainer = FickianDiffusionCoefficients<Scalar, numPhases, numComponents>;

    /*!
     * \brief Returns the diffusive fluxes of all components within
     *        a fluid phase across the given sub-control volume face.
     *        The computed fluxes are given in mole/s or kg/s, depending
     *        on the template parameter ReferenceSystemFormulation.
     */
    template<class Problem, class ElementVolumeVariables>
    static NumEqVector flux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace &scvf)
    {
        NumEqVector flux(0.0);

        // There is no diffusion over outflow boundaries (grad x == 0).
        // We assume that if an outflow BC is set for the first transported component, this
        // also holds for all other components.
        if (scvf.boundary() && problem.boundaryTypes(element, scvf).isOutflow(Indices::conti0EqIdx + 1))
            return flux;

        const int phaseIdx = 0;

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        const Scalar insideDistance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
        const Scalar insideDensity = massOrMolarDensity(insideVolVars, referenceSystem, phaseIdx);

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            if (compIdx == FluidSystem::getMainComponent(phaseIdx))
                continue;

            const Scalar massOrMoleFractionInside = massOrMoleFraction(insideVolVars, referenceSystem, phaseIdx, compIdx);
            const Scalar massOrMoleFractionOutside =  massOrMoleFraction(outsideVolVars, referenceSystem, phaseIdx, compIdx);
            const Scalar insideD = getEffectiveDiffusionCoefficient_(insideVolVars, phaseIdx, compIdx) * insideVolVars.extrusionFactor();

            if (scvf.boundary())
            {
                flux[compIdx] = insideDensity * insideD
                                * (massOrMoleFractionInside - massOrMoleFractionOutside) / insideDistance;
            }
            else
            {
                const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                const Scalar outsideD = getEffectiveDiffusionCoefficient_(outsideVolVars, phaseIdx, compIdx)
                                      * outsideVolVars.extrusionFactor();
                const Scalar outsideDistance = (outsideScv.dofPosition() - scvf.ipGlobal()).two_norm();
                const Scalar outsideDensity = massOrMolarDensity(outsideVolVars, referenceSystem, phaseIdx);

                const Scalar avgDensity = 0.5*(insideDensity + outsideDensity);
                const Scalar avgD = harmonicMean(insideD, outsideD, insideDistance, outsideDistance);

                flux[compIdx] = avgDensity * avgD
                                * (massOrMoleFractionInside - massOrMoleFractionOutside) / (insideDistance + outsideDistance);
            }
        }

        // Fick's law (for binary systems) states that the net flux of mass within the bulk phase has to be zero:
        const Scalar cumulativeFlux = std::accumulate(flux.begin(), flux.end(), 0.0);
        flux[FluidSystem::getMainComponent(0)] = -cumulativeFlux;

        flux *= Extrusion::area(scvf);

        return flux;
    }

private:
    static Scalar getEffectiveDiffusionCoefficient_(const VolumeVariables& volVars, const int phaseIdx, const int compIdx)
    {
        return volVars.effectiveDiffusionCoefficient(phaseIdx, FluidSystem::getMainComponent(phaseIdx), compIdx);
    }
};
} // end namespace Dumux

#endif
