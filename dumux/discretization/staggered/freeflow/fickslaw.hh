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
 *        diffusive mass fluxes due to molecular diffusion with Fick's law.
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FICKS_LAW_HH

#include <numeric>
#include <dune/common/float_cmp.hh>

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
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    static constexpr int numComponents = ModelTraits::numComponents();
    static constexpr int phaseIdx = Indices::phaseIdx;
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    static_assert(ModelTraits::numPhases() == 1, "Only one phase allowed supported!");

    enum {
        pressureIdx = Indices::pressureIdx,
        conti0EqIdx = Indices::conti0EqIdx,
        mainCompIdx = Indices::mainCompIdx,
        replaceCompEqIdx = Indices::replaceCompEqIdx,
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::staggered;

    //! state the type for the corresponding cache
    //! We don't cache anything for this law
    using Cache = FluxVariablesCaching::EmptyDiffusionCache;

    static CellCenterPrimaryVariables diffusiveFluxForCellCenter(const Problem& problem,
                                                           const FVElementGeometry& fvGeometry,
                                                           const ElementVolumeVariables& elemVolVars,
                                                           const SubControlVolumeFace &scvf)
    {
        CellCenterPrimaryVariables flux(0.0);

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        const Scalar insideMolarDensity = insideVolVars.molarDensity();

        for(int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            if(compIdx == mainCompIdx)
                continue;

            // get equation index
            const auto eqIdx = conti0EqIdx + compIdx;

            if(scvf.boundary())
            {
                const auto bcTypes = problem.boundaryTypesAtPos(scvf.center());
                if(bcTypes.isOutflow(eqIdx) && eqIdx != pressureIdx)
                    return flux;
            }

            const Scalar tij = transmissibility_(problem, fvGeometry, elemVolVars, scvf, compIdx);
            const Scalar insideMoleFraction = insideVolVars.moleFraction(phaseIdx, compIdx);

            const Scalar outsideMolarDensity = scvf.boundary() ? insideVolVars.molarDensity() : outsideVolVars.molarDensity();
            const Scalar avgDensity = 0.5*(insideMolarDensity + outsideMolarDensity);
            const Scalar outsideMoleFraction = outsideVolVars.moleFraction(phaseIdx, compIdx);
            flux[compIdx] = avgDensity * tij * (insideMoleFraction - outsideMoleFraction);
        }

        if(!(useMoles && replaceCompEqIdx == mainCompIdx))
        {
            const Scalar cumulativeFlux = std::accumulate(flux.begin(), flux.end(), 0.0);
            flux[mainCompIdx] = - cumulativeFlux;
        }

        if(useMoles && replaceCompEqIdx <= numComponents)
            flux[replaceCompEqIdx] = 0.0;

        // Fick's law (for binary systems) states that the net flux of moles within the bulk phase has to be zero:
        // If a given amount of molecules A travel into one direction, the same amount of molecules B have to
        // go into the opposite direction. This, however, means that the net mass flux whithin the bulk phase
        // (given different molecular weigths of A and B) does not have to be zero. We account for this:
        if(!useMoles)
        {
            //convert everything to a mass flux
            for(int compIdx = 0; compIdx < numComponents; ++compIdx)
                flux[compIdx] *= FluidSystem::molarMass(compIdx);


            if(replaceCompEqIdx < numComponents)
            {
                for(int compIdx = 0; compIdx < numComponents; ++compIdx)
                    flux[replaceCompEqIdx] += (compIdx != replaceCompEqIdx) ? flux[compIdx] : 0.0;
            }
        }

        return flux;
    }

    static Scalar transmissibility_(const Problem& problem,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolumeFace& scvf,
                                    const int compIdx)
    {
        Scalar tij = 0.0;
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        const Scalar insideDistance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
        const Scalar insideD = insideVolVars.diffusionCoefficient(phaseIdx, compIdx);
        const Scalar ti = calculateOmega_(insideDistance, insideD, 1.0);

        if(scvf.boundary())
            tij = scvf.area() * ti;
        else
        {
            const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
            const Scalar outsideDistance = (outsideScv.dofPosition() - scvf.ipGlobal()).two_norm();
            const Scalar outsideD = outsideVolVars.diffusionCoefficient(phaseIdx, compIdx);
            const Scalar tj = calculateOmega_(outsideDistance, outsideD, 1.0);

            tij = scvf.area()*(ti * tj)/(ti + tj);
        }
        return tij;
    }

    static Scalar calculateOmega_(const Scalar distance,
                                  const Scalar D,
                                  const Scalar extrusionFactor)
    {
        Scalar omega = D / distance;
        omega *= extrusionFactor;

        return omega;
    }

};
} // end namespace

#endif
