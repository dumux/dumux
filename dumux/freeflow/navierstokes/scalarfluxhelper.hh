// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \brief Navier Stokes scalar boundary flux helper
 */
#ifndef DUMUX_NAVIERSTOKES_SCALAR_BOUNDARY_FLUXHELPER_HH
#define DUMUX_NAVIERSTOKES_SCALAR_BOUNDARY_FLUXHELPER_HH

#include <dune/common/float_cmp.hh>
#include <dune/common/std/type_traits.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {
// helper structs and functions detecting if the VolumeVariables belong to a non-isothermal model
template <class Indices>
using NonisothermalDetector = decltype(std::declval<Indices>().energyEqIdx);

template<class Indices>
static constexpr bool isNonIsothermal()
{ return Dune::Std::is_detected<NonisothermalDetector, Indices>::value; }

} // end namespace Detail
#endif

/*!
 * \ingroup NavierStokesModel
 * \brief Navier Stokes scalar boundary flux helper
 */
template<class AdvectiveFlux>
struct NavierStokesScalarBoundaryFluxHelper
{
    /*!
     * \brief Return the area-specific, weighted advective flux of a scalar quantity.
     */
    template<class VolumeVariables, class SubControlVolumeFace, class Scalar, class UpwindTerm>
    static Scalar advectiveScalarUpwindFlux(const VolumeVariables& insideVolVars,
                                            const VolumeVariables& outsideVolVars,
                                            const SubControlVolumeFace& scvf,
                                            const Scalar volumeFlux,
                                            const Scalar upwindWeight,
                                            UpwindTerm upwindTerm)
    {
        using std::signbit;
        const bool insideIsUpstream = !signbit(volumeFlux);

        const auto& upstreamVolVars = insideIsUpstream ? insideVolVars : outsideVolVars;
        const auto& downstreamVolVars = insideIsUpstream ? outsideVolVars : insideVolVars;

        return (upwindWeight * upwindTerm(upstreamVolVars) +
               (1.0 - upwindWeight) * upwindTerm(downstreamVolVars))
               * volumeFlux;
    }

    template<class Indices, class NumEqVector, class UpwindFunction>
    static void addModelSpecificAdvectiveFlux(NumEqVector& flux,
                                              const UpwindFunction& upwind)
    {
        // add advective fluxes based on physical type of model
        AdvectiveFlux::addAdvectiveFlux(flux, upwind);

        // for non-isothermal models, add the energy flux
        if constexpr (Detail::isNonIsothermal<Indices>())
        {
            auto upwindTerm = [](const auto& volVars) { return volVars.density()*volVars.enthalpy(); };
            flux[Indices::energyEqIdx] = upwind(upwindTerm);
        }
    }

    /*!
     * \brief Return the area-specific outflow fluxes for all scalar balance equations.
     *        The values specified in outsideBoundaryPriVars are used in case of flow reversal.
     */
    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables>
    static auto scalarOutflowFlux(const Problem& problem,
                                  const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                  const ElementVolumeVariables& elemVolVars,
                                  typename ElementVolumeVariables::VolumeVariables::PrimaryVariables&& outsideBoundaryPriVars,
                                  const typename ElementVolumeVariables::VolumeVariables::PrimaryVariables::value_type upwindWeight = 1.0)
    {
        using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
        using PrimaryVariables = typename VolumeVariables::PrimaryVariables;
        using NumEqVector = PrimaryVariables;
        NumEqVector flux;
        const auto velocity = problem.faceVelocity(element,fvGeometry, scvf);
        const auto volumeFlux = velocity * scvf.unitOuterNormal();
        using std::signbit;
        const bool insideIsUpstream = !signbit(volumeFlux);
        const VolumeVariables& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const VolumeVariables& outsideVolVars = [&]()
        {
            // only use the inside volVars for "true" outflow conditions (avoid constructing the outside volVars)
            if (insideIsUpstream && Dune::FloatCmp::eq(upwindWeight, 1.0, 1e-6))
                return insideVolVars;
            else
            {
                // construct outside volVars from the given priVars for situations of flow reversal
                VolumeVariables boundaryVolVars;
                boundaryVolVars.update(elementSolution<FVElementGeometry>(std::forward<PrimaryVariables>(outsideBoundaryPriVars)),
                                       problem,
                                       element,
                                       fvGeometry.scv(scvf.insideScvIdx()));
                return boundaryVolVars;
            }
        }();

        auto upwindFuntion = [&](const auto& upwindTerm)
        {
            return advectiveScalarUpwindFlux(insideVolVars, outsideVolVars, scvf, volumeFlux, upwindWeight, upwindTerm);
        };

        addModelSpecificAdvectiveFlux<typename VolumeVariables::Indices>(flux, upwindFuntion);

        return flux;
    }

    /*!
     * \brief Return the area-specific outflow fluxes for all scalar balance equations.
     *        This should only be used of flow reversal does never occur.
     *        A (deactivable) warning is emitted otherwise.
     */
    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables>
    static auto scalarOutflowFlux(const Problem& problem,
                                  const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                  const ElementVolumeVariables& elemVolVars)
    {
        using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
        using NumEqVector = typename VolumeVariables::PrimaryVariables;
        NumEqVector flux;
        const auto velocity = problem.faceVelocity(element,fvGeometry, scvf);
        const auto volumeFlux = velocity * scvf.unitOuterNormal();
        using std::signbit;
        const bool insideIsUpstream = !signbit(volumeFlux);
        const VolumeVariables& insideVolVars = elemVolVars[scvf.insideScvIdx()];

        if constexpr (VolumeVariables::FluidSystem::isCompressible(0/*phaseIdx*/) /*TODO viscosityIsConstant*/ || NumEqVector::size() > 1)
        {
            static const bool verbose = getParamFromGroup<bool>(problem.paramGroup(), "Flux.EnableOutflowReversalWarning", true);
            using std::abs;
            if (verbose && !insideIsUpstream && abs(volumeFlux) > 1e-10)
            {
                std::cout << "velo " << velocity << ", flux " << volumeFlux << std::endl;
                std::cout << "\n ********** WARNING ********** \n\n"
                "Outflow condition set at " << scvf.center() << " might be invalid due to flow reversal. "
                "Consider using \n"
                "outflowFlux(problem, element, fvGeometry, scvf, elemVolVars, outsideBoundaryPriVars, upwindWeight) \n"
                "instead where you can specify primary variables for inflow situations.\n"
                "\n ***************************** \n" << std::endl;
            }
        }

        auto upwindFuntion = [&](const auto& upwindTerm)
        {
            return advectiveScalarUpwindFlux(insideVolVars, insideVolVars, scvf, volumeFlux, 1.0 /*upwindWeight*/, upwindTerm);
        };

        addModelSpecificAdvectiveFlux<typename VolumeVariables::Indices>(flux, upwindFuntion);

        return flux;
    }
};

} // end namespace Dumux

#endif
