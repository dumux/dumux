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
 * \ingroup CCTpfaDiscretization
 * \brief Darcy's law for cell-centered finite volume schemes with two-point flux approximation
 */
#ifndef MY_DUMUX_DISCRETIZATION_CC_TPFA_DARCYS_LAW_HH
#define MY_DUMUX_DISCRETIZATION_CC_TPFA_DARCYS_LAW_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/discretization/cellcentered/tpfa/darcyslaw.hh>

namespace Dumux {

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Darcy's law for cell-centered finite volume schemes with two-point flux approximation
 * \note Darcy's law is speialized for network and surface grids (i.e. if grid dim < dimWorld)
 * \tparam Scalar the scalar type for scalar physical quantities
 * \tparam FVGridGeometry the grid geometry
 * \tparam isNetwork whether we are computing on a network grid embedded in a higher world dimension
 */
template<class Scalar, class FVGridGeometry, bool isNetwork = false>
class MyCCTpfaDarcysLaw;

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Specialization of the CCTpfaDarcysLaw grids where dim=dimWorld
 */
template<class ScalarType, class FVGridGeometry>
class MyCCTpfaDarcysLaw<ScalarType, FVGridGeometry, /*isNetwork*/ false>
{
    using ThisType = CCTpfaDarcysLaw<ScalarType, FVGridGeometry, /*isNetwork*/ false>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

  public:
    //! state the scalar type of the law
    using Scalar = ScalarType;

    //! state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;

    //! state the type for the corresponding cache
    using Cache = TpfaDarcysLawCache<ThisType, FVGridGeometry>;

    //! Compute the advective flux
    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");

        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        if (enableGravity)
        {
            // do averaging for the density over all neighboring elements
            const auto rho = scvf.boundary() ? outsideVolVars.density(phaseIdx)
                                             : (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;

            // Obtain inside and outside pressures
            const auto pInside = insideVolVars.pressure(phaseIdx);
            const auto pOutside = outsideVolVars.pressure(phaseIdx);

            const auto& tij = fluxVarsCache.advectionTij();
            const auto& g = problem.gravityAtPos(scvf.ipGlobal());

            //! compute alpha := n^T*K*g
            const auto alpha_inside = vtmv(scvf.unitOuterNormal(), insideVolVars.permeability(), g)*insideVolVars.extrusionFactor();

            Scalar flux = tij*(pInside - pOutside) + rho*scvf.area()*alpha_inside;

            //! On interior faces we have to add K-weighted gravitational contributions
            if (!scvf.boundary())
            {
                const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                const auto outsideK = outsideVolVars.permeability();
                const auto outsideTi = computeTpfaTransmissibility(scvf, outsideScv, outsideK, outsideVolVars.extrusionFactor());
                const auto alpha_outside = vtmv(scvf.unitOuterNormal(), outsideK, g)*outsideVolVars.extrusionFactor();

                flux += rho*tij/outsideTi*(alpha_inside - alpha_outside);
            }

            return flux;
        }
        else
        {
            // Obtain inside and outside pressures
            const auto pInside = insideVolVars.pressure(phaseIdx);
            const auto pOutside = outsideVolVars.pressure(phaseIdx);

            // return flux
            auto flux = fluxVarsCache.advectionTij()*(pInside - pOutside);

            if (std::abs(scvf.unitOuterNormal()*GlobalPosition({0.0, 0.0, 1.0})) < 1e-10)
            {
                static const Scalar k = getParam<Scalar>("SpatialParams.PermeabilityTissue");


                // additional fluxes due to point sources
                auto key = std::make_pair(problem.fvGridGeometry().elementMapper().index(element), 0);
                auto keyOutside = std::make_pair(problem.fvGridGeometry().elementMapper().index(problem.fvGridGeometry().element(scvf.outsideScvIdx())), 0);
                if (problem.pointSourceMap().count(key) || problem.pointSourceMap().count(keyOutside))
                {
                    // call the solDependent function. Herein the user might fill/add values to the point sources
                    // we make a copy of the local point sources here
                    auto pointSources = problem.pointSourceMap().at(key);

                    for (auto&& pointSource : pointSources)
                    {
                        auto dir = ((scvf.center())*scvf.unitOuterNormal())/scvf.center().two_norm();
                        if (std::abs(dir) > 1e-15)
                            flux += -k*problem.exactGradient(scvf.center())*scvf.area()*insideVolVars.extrusionFactor()*dir;
                    }
                }
            }

            return flux;
        }
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibility will be computed and stored using the method below.
    template<class Problem, class ElementVolumeVariables>
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const SubControlVolumeFace& scvf)
    {
        Scalar tij;

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];

        const Scalar ti = computeTpfaTransmissibility(scvf, insideScv,
                                                      getPermeability_(problem, insideVolVars, scvf.ipGlobal()),
                                                      insideVolVars.extrusionFactor());

        // on the boundary (dirichlet) we only need ti
        if (scvf.boundary())
            tij = scvf.area()*ti;

        // otherwise we compute a tpfa harmonic mean
        else
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            // as we assemble fluxes from the neighbor to our element the outside index
            // refers to the scv of our element, so we use the scv method
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const Scalar tj = -1.0*computeTpfaTransmissibility(scvf, outsideScv,
                                                               getPermeability_(problem, outsideVolVars, scvf.ipGlobal()),
                                                               outsideVolVars.extrusionFactor());

            // harmonic mean (check for division by zero!)
            // TODO: This could lead to problems!? Is there a better way to do this?
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = scvf.area()*(ti * tj)/(ti + tj);
        }

        return tij;
    }

private:
    template<class Problem, class VolumeVariables,
             std::enable_if_t<!Problem::SpatialParams::evaluatePermeabilityAtScvfIP(), int> = 0>
    static decltype(auto) getPermeability_(const Problem& problem,
                                           const VolumeVariables& volVars,
                                           const GlobalPosition& scvfIpGlobal)
    { return volVars.permeability(); }

    template<class Problem, class VolumeVariables,
             std::enable_if_t<Problem::SpatialParams::evaluatePermeabilityAtScvfIP(), int> = 0>
    static decltype(auto) getPermeability_(const Problem& problem,
                                           const VolumeVariables& volVars,
                                           const GlobalPosition& scvfIpGlobal)
    { return problem.spatialParams().permeabilityAtPos(scvfIpGlobal); }
};

} // end namespace Dumux

#endif
