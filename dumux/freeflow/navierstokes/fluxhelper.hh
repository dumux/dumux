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

#ifndef DUMUX_NAVIERSTOKES_BOUNDARY_FLUXHELPER_HH
#define DUMUX_NAVIERSTOKES_BOUNDARY_FLUXHELPER_HH

#include <type_traits>
#include <dune/common/float_cmp.hh>
#include <dune/common/std/type_traits.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/elementsolution.hh>
#include <dumux/freeflow/navierstokes/momentum/velocitygradients.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {
// helper structs and functions detecting if the VolumeVariables belong to a non-isothermal model
template <typename T>
using NonisothermalDetector = decltype(std::declval<T>().energyEqIdx);

template<class T>
static constexpr bool isNonIsothermal()
{ return Dune::Std::is_detected<NonisothermalDetector, T>::value; }


} // end namespace Detail
#endif


struct NavierStokesUpwindTerms
{
    static auto totalMassFlux()
    { return [](const auto& volVars) { return volVars.density(); }; }

    static auto comoponentMoleFlux(const int compIdx)
    { return [compIdx](const auto& volVars) { return volVars.molarDensity()*volVars.moleFraction(compIdx); }; }

    static auto comoponentMassFlux(const int compIdx)
    { return [compIdx](const auto& volVars) { return volVars.density()*volVars.massFraction(compIdx); }; }

    static auto energyFlux()
    { return [](const auto& volVars) { return volVars.density()*volVars.enthalpy(); }; }
};

struct NavierStokesMassOnePModelTraits;




template<class ModelTraits>
class NavierStokesBoundaryFluxHelper
{


public:

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

    template<class Traits, class NumEqVector, class UpwindFunction>
    static void addModelSpecificFlux(NumEqVector& flux,
                                     const UpwindFunction& upwind)
    {
        // for non-isothermal models, first add the fluxes of the underlying
        // isothermal model, then add the energy flux
        if constexpr (Detail::isNonIsothermal<typename Traits::Indices>())
        {
            addModelSpecificFlux<typename Traits::IsothermalTraits>(flux, upwind);

            auto upwindTerm = NavierStokesUpwindTerms::energyFlux();
            flux[Traits::Indices::energyEqIdx] = upwind(upwindTerm);
        }

        // the 1p model: add the total mass flux
        if constexpr (std::is_same_v<Traits, NavierStokesMassOnePModelTraits>)
        {
            auto upwindTerm = NavierStokesUpwindTerms::totalMassFlux();
            flux[Traits::Indices::conti0EqIdx] = upwind(upwindTerm);
        }

        // TODO 1pnc

        // TODO RANS

        // TODO maybe use tag dispatch istead of constexpr if?
    }

    /*!
     * \brief Return the area-specific outflow fluxes for all scalar balance equations.
     *        The values specified in outsideBoundaryPriVars are used in case of flow reversal.
     */
    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class PrimaryVariables, class Scalar>
    static auto scalarOutflowFlux(const Problem& problem,
                                  const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                  const ElementVolumeVariables& elemVolVars,
                                  PrimaryVariables&& outsideBoundaryPriVars,
                                  const Scalar upwindWeight = 1.0)
    {
        using NumEqVector = PrimaryVariables;
        using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
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

        addModelSpecificFlux<ModelTraits>(flux, upwindFuntion);

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
            if (verbose && !insideIsUpstream)
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

        addModelSpecificFlux<ModelTraits>(flux, upwindFuntion);

        return flux;
    }

    template<class Problem, class Element, class FVElementGeometry, class ElementVolumeVariables, class ElementFluxVariablesCache, class Scalar>
    static auto fixedPressureMomentumFlux(const Problem& problem,
                                          const Element& element,
                                          const FVElementGeometry& fvGeometry,
                                          const typename FVElementGeometry::SubControlVolumeFace& scvf,
                                          const ElementVolumeVariables& elemVolVars,
                                          const ElementFluxVariablesCache& elemFluxVarsCache,
                                          const Scalar pressure,
                                          const bool zeroNormalVelocityGradient = true)
    {
        // TODO density upwinding?
        static_assert(FVElementGeometry::GridGeometry::discMethod == DiscretizationMethod::fcstaggered);
        using NumEqVector = decltype(problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf));
        NumEqVector flux;

        if (scvf.isFrontal())
        {
            // pressure contribution
            flux[scvf.directionIndex()] = (pressure - problem.referencePressure(element, fvGeometry, scvf)) * scvf.directionSign();

            if (problem.enableInertiaTerms())
            {
                const auto v = elemVolVars[scvf.insideScvIdx()].velocity();
                flux[scvf.directionIndex()] += v*v * problem.density(element, fvGeometry, scvf) * scvf.directionSign();
            }
        }

        // if no zero velocity gradient is desired we are done here, ....
        if (!zeroNormalVelocityGradient)
            return flux;

        // ..., otherwise, make sure the flow does not diverge by accounting for the off-diagonal entries of the stress tensor
        // TODO normal viscous term?

        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

        if (scvf.isLateral())
        {
            // viscous terms
            const Scalar mu = problem.effectiveViscosity(element, fvGeometry, scvf);

            // lateral face normal to boundary (integration point touches boundary)
            if (scv.boundary())
                flux[scvf.directionIndex()] -= mu * StaggeredVelocityGradients::velocityGradIJ(fvGeometry, scvf, elemVolVars)
                                               * scvf.directionSign();
            // lateral face coinciding with boundary
            else if (scvf.boundary())
                flux[scv.directionIndex()] -= mu * StaggeredVelocityGradients::velocityGradJI(fvGeometry, scvf, elemVolVars)
                                              * scvf.directionSign();

            // advective terms
            if (problem.enableInertiaTerms())
            {
                const auto transportingVelocity = [&]()
                {
                    const auto& orthogonalScvf = fvGeometry.lateralOrthogonalScvf(scvf);
                    const auto innerTransportingVelocity = elemVolVars[orthogonalScvf.insideScvIdx()].velocity();

                    if (scvf.boundary())
                    {
                        if (const auto bcTypes = problem.boundaryTypes(element, scvf); bcTypes.isDirichlet(scvf.directionIndex()))
                            return problem.dirichlet(element, scvf)[scvf.directionIndex()];
                        else
                            return
                                innerTransportingVelocity; // fallback
                    }
                    else
                    {
                        static const bool useOldScheme = getParam<bool>("FreeFlow.UseOldTransportingVelocity", true); // TODO how to deprecate?
                        if (useOldScheme)
                            return innerTransportingVelocity;
                        else
                        {
                            // average the transporting velocity by weighting with the scv volumes
                            const auto insideVolume = fvGeometry.scv(orthogonalScvf.insideScvIdx()).volume();
                            const auto outsideVolume = fvGeometry.scv(orthogonalScvf.outsideScvIdx()).volume();
                            const auto outerTransportingVelocity = elemVolVars[orthogonalScvf.outsideScvIdx()].velocity();
                            return (insideVolume*innerTransportingVelocity + outsideVolume*outerTransportingVelocity) / (insideVolume + outsideVolume);
                        }
                    }
                }();

                // lateral face normal to boundary (integration point touches boundary)
                if (scv.boundary())
                {
                    const auto innerVelocity = elemVolVars[scvf.insideScvIdx()].velocity();
                    const auto outerVelocity = elemVolVars[scvf.outsideScvIdx()].velocity();
                    const auto rho = problem.getInsideAndOutsideDensity(element, fvGeometry, scvf);

                    const bool selfIsUpstream = scvf.directionSign() == sign(transportingVelocity);

                    const auto insideMomentum = innerVelocity * rho.first;
                    const auto outsideMomentum = outerVelocity * rho.second;

                    static const auto upwindWeight = getParamFromGroup<Scalar>(problem.paramGroup(), "Flux.UpwindWeight");

                    const auto transportedMomentum =  selfIsUpstream ? (upwindWeight * insideMomentum + (1.0 - upwindWeight) * outsideMomentum)
                                                                     : (upwindWeight * outsideMomentum + (1.0 - upwindWeight) * insideMomentum);

                    flux[scvf.directionIndex()] += transportingVelocity * transportedMomentum * scvf.directionSign();
                }

                // lateral face coinciding with boundary
                else if (scvf.boundary())
                {
                    const auto insideDensity = problem.density(element, fvGeometry.scv(scvf.insideScvIdx()));
                    const auto innerVelocity = elemVolVars[scvf.insideScvIdx()].velocity();
                    flux[scv.directionIndex()] += innerVelocity * transportingVelocity * insideDensity * scvf.directionSign();
                }
            }
        }


        return flux;
    }

};

} // end namespace Dumux

#endif
