// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *        diffusive molar fluxes due to molecular diffusion with Fick's law.
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FICKS_LAW_HH

#include <numeric>
#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fluxvariablescaching.hh>

namespace Dumux
{
// forward declaration
template<class TypeTag, DiscretizationMethod discMethod>
class FicksLawImplementation;

/*!
 * \ingroup StaggeredFicksLaw
 * \brief Specialization of Fick's Law for the staggered free flow method.
 */
template <class TypeTag>
class FicksLawImplementation<TypeTag, DiscretizationMethod::staggered >
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using Indices = typename ModelTraits::Indices;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    static constexpr int numComponents = ModelTraits::numComponents();
    using NumEqVector = Dune::FieldVector<Scalar, numComponents>;

    static_assert(ModelTraits::numPhases() == 1, "Only one phase supported!");

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::staggered;

    //! state the type for the corresponding cache
    //! We don't cache anything for this law
    using Cache = FluxVariablesCaching::EmptyDiffusionCache;

    template<class Problem>
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

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        const Scalar insideDistance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
        const Scalar insideMolarDensity = insideVolVars.molarDensity();

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            if (compIdx == FluidSystem::getMainComponent(0))
                continue;

            const Scalar insideMoleFraction = insideVolVars.moleFraction(compIdx);
            const Scalar outsideMoleFraction = outsideVolVars.moleFraction(compIdx);

            const Scalar insideD = insideVolVars.effectiveDiffusivity(0, compIdx) * insideVolVars.extrusionFactor();

            if (scvf.boundary())
            {
                flux[compIdx] = insideMolarDensity * insideD
                                * (insideMoleFraction - outsideMoleFraction) / insideDistance;
            }
            else
            {
                const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                const Scalar outsideD = outsideVolVars.effectiveDiffusivity(0, compIdx) * outsideVolVars.extrusionFactor();
                const Scalar outsideDistance = (outsideScv.dofPosition() - scvf.ipGlobal()).two_norm();
                const Scalar outsideMolarDensity = outsideVolVars.molarDensity();

                const Scalar avgDensity = 0.5*(insideMolarDensity + outsideMolarDensity);
                const Scalar avgD = harmonicMean(insideD, outsideD, insideDistance, outsideDistance);

                flux[compIdx] = avgDensity * avgD
                                * (insideMoleFraction - outsideMoleFraction) / (insideDistance + outsideDistance);
            }
        }

        // Fick's law (for binary systems) states that the net flux of moles within the bulk phase has to be zero:
        // If a given amount of molecules A travel into one direction, the same amount of molecules B have to
        // go into the opposite direction.
        const Scalar cumulativeFlux = std::accumulate(flux.begin(), flux.end(), 0.0);
        flux[FluidSystem::getMainComponent(0)] = -cumulativeFlux;

        flux *= scvf.area();

        return flux;
    }
};
} // end namespace

#endif
