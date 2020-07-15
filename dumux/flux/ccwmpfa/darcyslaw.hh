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
 * \ingroup CCWMpfaFlux
 * \brief Darcy's law for cell-centered finite volume schemes with two-point flux approximation
 */
#ifndef DUMUX_DISCRETIZATION_CC_WMPFA_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_WMPFA_DARCYS_LAW_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

namespace Dumux {

//! forward declaration of the method-specific implementation
template<class TypeTag, DiscretizationMethod discMethod>
class DarcysLawImplementation;

/*!
 * \ingroup CCWMpfaFlux
 * \brief Class that fills the cache corresponding to wmpfa Darcy's Law
 */
template<class GridGeometry>
class WMpfaDarcysLawCacheFiller
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

public:
    //! Function to fill a WMpfaDarcysLawCache of a given scvf
    //! This interface has to be met by any advection-related cache filler class
    //! TODO: Probably get cache type out of the filler
    template<class FluxVariablesCache, class Problem, class ElementVolumeVariables, class FluxVariablesCacheFiller>
    static void fill(FluxVariablesCache& scvfFluxVarsCache,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolumeFace& scvf,
                     const FluxVariablesCacheFiller& fluxVarsCacheFiller)
    {
        scvfFluxVarsCache.updateAdvection(fluxVarsCacheFiller.faceDataHandle());
    }
};

/*!
 * \ingroup CCWMpfaFlux
 * \brief The cache corresponding to tpfa Darcy's Law
 */
template<class GridGeometry, class ElementFluxVariablesCache>
class WMpfaDarcysLawCache
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

    using DataHandle = typename ElementFluxVariablesCache::FaceDataHandle;
    using AdvectionDataHandle = typename DataHandle::AdvectionHandle;

public:
    using Filler = WMpfaDarcysLawCacheFiller<GridGeometry>;

    void updateAdvection(const DataHandle& dataHandle)
    {
        handle_ = &dataHandle.advectionHandle();
    }

    const AdvectionDataHandle& dataHandle() const
    { return *handle_; }

private:
    const AdvectionDataHandle* handle_;
};

/*!
 * \ingroup CCWMpfaFlux
 * \brief Specialization of the CCWMpfaDarcysLaw grids where dim=dimWorld
 */
template<class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethod::ccwmpfa>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;

    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

  public:
    //! state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::ccwmpfa;

    //! state the type for the corresponding cache
    using Cache = WMpfaDarcysLawCache<GridGeometry, ElementFluxVariablesCache>;

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
        const auto& dataHandle = fluxVarsCache.dataHandle();

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        const auto& subFluxDataIJ = dataHandle.subFluxData();

        if (enableGravity)
        {
            // do averaging for the density over all neighboring elements
            const auto rho = scvf.boundary() ? outsideVolVars.density(phaseIdx)
                                             : (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;

            const auto& g = problem.spatialParams().gravity(scvf.ipGlobal());

            auto pseudoPotential = [&elemVolVars, phaseIdx, &g, &rho](const auto& volVars, const auto& x)
                                    { return volVars.pressure(phaseIdx) - rho*(g*x);};

            Scalar fluxIJ = 0.0;
            Scalar wIJ = 0.5;
            std::for_each(subFluxDataIJ.cbegin(), subFluxDataIJ.cend(),
                            [&fluxIJ, &elemVolVars, &pseudoPotential](const auto& e)
                             { fluxIJ += e.coefficient * pseudoPotential(elemVolVars[e.index], e.position); } );

            Scalar fluxJI = 0.0;
            Scalar wJI = 0.5;
            if(!scvf.boundary())
            {
                const auto& flippedScvf = fvGeometry.flipScvf(scvf.index());
                const auto& dataHandleJ = elemFluxVarsCache[flippedScvf].dataHandle();
                const auto& subFluxDataJI = dataHandleJ.subFluxData();
                std::for_each(subFluxDataJI.cbegin(), subFluxDataJI.cend(),
                                [&fluxJI, &elemVolVars, &pseudoPotential](const auto& e)
                                { fluxJI += e.coefficient * pseudoPotential(elemVolVars[e.index], e.position); } );
            }
            else
            {
                wIJ = 1.0;
                wJI = 0.0;
            }

            return scvf.area() * (wIJ*fluxIJ - wJI*fluxJI);
        }
        else
        {
            auto pressure = [&elemVolVars, phaseIdx](const auto& volVars)
                                    { return volVars.pressure(phaseIdx);};

            Scalar fluxIJ = 0.0;
            Scalar wIJ = 0.5;
            std::for_each(subFluxDataIJ.cbegin(), subFluxDataIJ.cend(),
                            [&fluxIJ, &elemVolVars, &pressure](const auto& e)
                             { fluxIJ += e.coefficient * pressure(elemVolVars[e.index]); } );

            Scalar fluxJI = 0.0;
            Scalar wJI = 0.5;
            if(!scvf.boundary())
            {
                const auto& flippedScvf = fvGeometry.flipScvf(scvf.index());
                const auto& dataHandleJ = elemFluxVarsCache[flippedScvf].dataHandle();
                const auto& subFluxDataJI = dataHandleJ.subFluxData();
                std::for_each(subFluxDataJI.cbegin(), subFluxDataJI.cend(),
                                [&fluxJI, &elemVolVars, &pressure](const auto& e)
                                { fluxJI += e.coefficient * pressure(elemVolVars[e.index]); } );
            }
            else
            {
                wIJ = 1.0;
                wJI = 0.0;
            }

            return scvf.area() * (wIJ*fluxIJ - wJI*fluxJI);
        }

        return 0.0;
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
