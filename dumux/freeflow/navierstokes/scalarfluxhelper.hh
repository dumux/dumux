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

#ifndef DUMUX_NAVIERSTOKES_SCALAR_BOUNDARY_FLUXHELPER_HH
#define DUMUX_NAVIERSTOKES_SCALAR_BOUNDARY_FLUXHELPER_HH


#include <type_traits>
#include <dune/common/float_cmp.hh>
#include <dune/common/std/type_traits.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/elementsolution.hh>

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

struct NavierStokesScalarFluxHelper
{

    struct UpwindTerms
    {
        static auto totalMassFlux()
        { return [](const auto& volVars) { return volVars.density(); }; }

        static auto totalMoleFlux()
        { return [](const auto& volVars) { return volVars.molarDensity(); }; }

        static auto componentMoleFlux(const int compIdx)
        { return [compIdx](const auto& volVars) { return volVars.molarDensity()*volVars.moleFraction(compIdx); }; }

        static auto componentMassFlux(const int compIdx)
        { return [compIdx](const auto& volVars) { return volVars.density()*volVars.massFraction(compIdx); }; }

        static auto energyFlux()
        { return [](const auto& volVars) { return volVars.density()*volVars.enthalpy(); }; }
    };

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
};



template<class ModelTraits>
struct OnePFlux
{
    template<class NumEqVector, class UpwindFunction>
    static void addAdvectiveFlux(NumEqVector& flux,
                                 const UpwindFunction& upwind)
    {
        auto upwindTerm = NavierStokesScalarFluxHelper::UpwindTerms::totalMassFlux();
        flux[ModelTraits::Indices::conti0EqIdx] = upwind(upwindTerm);
    }
};

// template<class ModelTraits, class VolumeVariables>
// struct ModelSpecificFlux<ModelTraits, VolumeVariables, true>
template<class ModelTraits>
struct OnePNCFlux
{
    template<class NumEqVector, class UpwindFunction>
    static void addAdvectiveFlux(NumEqVector& flux,
                                 const UpwindFunction& upwind)
    {
        static constexpr bool useMoles = ModelTraits::useMoles();
        static constexpr auto numComponents = ModelTraits::numFluidComponents();
        static constexpr auto replaceCompEqIdx = ModelTraits::replaceCompEqIdx();
        static constexpr bool useTotalMoleOrMassBalance = replaceCompEqIdx < numComponents;

        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            // get equation index
            const auto eqIdx = ModelTraits::Indices::conti0EqIdx + compIdx;

            if (eqIdx == replaceCompEqIdx)
                continue;

            auto upwindTerm = [&]()
            {
                if constexpr (useMoles)
                    return NavierStokesScalarFluxHelper::UpwindTerms::componentMoleFlux(compIdx);
                else
                    return NavierStokesScalarFluxHelper::UpwindTerms::componentMassFlux(compIdx);
            }();

            flux[eqIdx] = upwind(upwindTerm);
        }

        // in case one balance is substituted by the total mole balance
        if constexpr(useTotalMoleOrMassBalance)
        {

            auto upwindTerm = [&]()
            {
                if constexpr (useMoles)
                    return NavierStokesScalarFluxHelper::UpwindTerms::totalMoleFlux();
                else
                    return NavierStokesScalarFluxHelper::UpwindTerms::totalMassFlux();
            }();

            flux[replaceCompEqIdx] = upwind(upwindTerm);
        }
    }
};

template<class ModelTraits, template<class> class ModelSpecificFlux>
struct NavierStokesScalarBoundaryFluxHelper
{



    template<class NumEqVector, class UpwindFunction>
    static void addModelSpecificAdvectiveFlux(NumEqVector& flux,
                                              const UpwindFunction& upwind)
    {
        // add advective fluxes based on physical type of model
        ModelSpecificFlux<ModelTraits>::addAdvectiveFlux(flux, upwind);

        // for non-isothermal models, add the energy flux
        if constexpr (Detail::isNonIsothermal<typename ModelTraits::Indices>())
        {
            auto upwindTerm = [](const auto& volVars) { return volVars.density()*volVars.enthalpy(); };
            flux[ModelTraits::Indices::energyEqIdx] = upwind(upwindTerm);
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
            return NavierStokesScalarFluxHelper::advectiveScalarUpwindFlux(insideVolVars, outsideVolVars, scvf, volumeFlux, upwindWeight, upwindTerm);
        };

        addModelSpecificAdvectiveFlux(flux, upwindFuntion);

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
            return NavierStokesScalarFluxHelper::advectiveScalarUpwindFlux(insideVolVars, insideVolVars, scvf, volumeFlux, 1.0 /*upwindWeight*/, upwindTerm);
        };

        addModelSpecificAdvectiveFlux(flux, upwindFuntion);

        return flux;
    }


    // template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class ElementFluxVariablesCache, class Scalar>
    // static auto dirichletWithZeroPressureGradient(const Problem& problem,
    //                                               const Element& element,
    //                                               const FVElementGeometry& fvGeometry,
    //                                               const typename FVElementGeometry::SubControlVolumeFace& scvf,
    //                                               const ElementVolumeVariables& elemVolVars,
    //                                               const ElementFluxVariablesCache& elemFluxVarsCache,
    //                                               PrimaryVariables&& outsideBoundaryPriVars,)
    // {
    //
    // }
};

} // end namespace Dumux

#endif
