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
#include <dumux/discretization/cellcentered/wmpfa/methods.hh>
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

template<WMpfaMethod method>
class WMpfaDarcysLawsFluxCalculator;

template<>
class WMpfaDarcysLawsFluxCalculator<WMpfaMethod::avgmpfa>
{
public:

    template<class Scalar, class SubControlVolumeFace, class ElementVolumeVariables, class FluxData, class SolValues>
    static void flux(Scalar& flux,
                     const SubControlVolumeFace& scvf,
                     const ElementVolumeVariables& elemVolVars,
                     const FluxData& dIJ,
                     const FluxData& dJI,
                     const SolValues& solValues)
    {
        Scalar fluxIJ = 0.0;
        Scalar fluxJI = 0.0;

        fluxIJ += dIJ[0].coefficient * solValues(elemVolVars[dIJ[0].index], dIJ[0].position);
        std::for_each(dIJ.cbegin()+1, dIJ.cend(),
                        [&fluxIJ, &elemVolVars, &solValues](const auto& e)
                        { fluxIJ -= e.coefficient * solValues(elemVolVars[e.index], e.position); } );
        fluxJI += dJI[0].coefficient * solValues(elemVolVars[dJI[0].index], dJI[0].position);
        std::for_each(dJI.cbegin()+1, dJI.cend(),
                        [&fluxJI, &elemVolVars, &solValues](const auto& e)
                        { fluxJI -= e.coefficient * solValues(elemVolVars[e.index], e.position); } );

        flux = scvf.area()*(0.5*fluxIJ - 0.5*fluxJI);
    }

    template <class Scalar, class SubControlVolumeFace, class ElementVolumeVariables, class FluxData, class SolValues>
    static void boundaryFlux(Scalar& flux,
                             const SubControlVolumeFace& scvf,
                             const ElementVolumeVariables& elemVolVars,
                             const FluxData& dIJ,
                             const SolValues& solValues)
    {
        Scalar fluxIJ = 0.0;
        fluxIJ += dIJ[0].coefficient * solValues(elemVolVars[dIJ[0].index], dIJ[0].position);
        std::for_each(dIJ.cbegin()+1, dIJ.cend(),
                        [&fluxIJ, &elemVolVars, &solValues](const auto& e)
                        { fluxIJ -= e.coefficient * solValues(elemVolVars[e.index], e.position); } );

        flux = scvf.area()*fluxIJ;
    }
};

template<>
class WMpfaDarcysLawsFluxCalculator<WMpfaMethod::nltpfa>
{
public:

    template<class Scalar, class SubControlVolumeFace, class ElementVolumeVariables, class FluxData, class SolValues>
    static void flux(Scalar& flux,
                     const SubControlVolumeFace& scvf,
                     const ElementVolumeVariables& elemVolVars,
                     const FluxData& dIJ,
                     const FluxData& dJI,
                     const SolValues& solValues)
    {
        Scalar lambdaIJ = 0.0;
        Scalar lambdaJI = 0.0;
        Scalar tIJ = 0.0;
        Scalar tJI = 0.0;
        std::for_each(dIJ.cbegin()+1, dIJ.cend(),
                        [&lambdaIJ, &tJI, &scvf, &elemVolVars, &solValues](const auto& e)
                        {
                            (e.index == scvf.outsideScvIdx())
                            ? tJI +=  e.coefficient
                            : lambdaIJ += e.coefficient * solValues(elemVolVars[e.index], e.position);
                        } );

        std::for_each(dJI.cbegin()+1, dJI.cend(),
                        [&lambdaJI, &tIJ, &scvf, &elemVolVars, &solValues](const auto& e)
                        {
                            (e.index == scvf.insideScvIdx())
                            ? tIJ +=  e.coefficient
                            : lambdaJI += e.coefficient * solValues(elemVolVars[e.index], e.position);
                        } );

        auto calcAbs = [](Scalar a){ return std::abs(a) + 1e-8; };

        Scalar wIJ = (calcAbs(lambdaJI))/(calcAbs(lambdaIJ) + calcAbs(lambdaJI));
        Scalar wJI = (calcAbs(lambdaIJ))/(calcAbs(lambdaIJ) + calcAbs(lambdaJI));

        tIJ *= wJI;
        tIJ += wIJ*dIJ[0].coefficient;
        tJI *= wIJ;
        tJI += wJI*dJI[0].coefficient;

        flux = scvf.area()*( tIJ * solValues(elemVolVars[dIJ[0].index], dIJ[0].position)
                            -tJI * solValues(elemVolVars[dJI[0].index], dJI[0].position));
    }

    template<class Scalar, class SubControlVolumeFace, class ElementVolumeVariables, class FluxData, class SolValues>
    static void boundaryFlux(Scalar& flux,
                             const SubControlVolumeFace& scvf,
                             const ElementVolumeVariables& elemVolVars,
                             const FluxData& dIJ,
                             const SolValues& solValues)
    {
        Scalar fluxIJ = 0.0;
        fluxIJ += dIJ[0].coefficient * solValues(elemVolVars[dIJ[0].index], dIJ[0].position);
        std::for_each(dIJ.cbegin()+1, dIJ.cend(),
                        [&fluxIJ, &elemVolVars, &solValues](const auto& e)
                        { fluxIJ -= e.coefficient * solValues(elemVolVars[e.index], e.position); } );

        flux = scvf.area()*fluxIJ;
    }
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

    static constexpr WMpfaMethod method = getPropValue<TypeTag, Properties::DiscretizationSubmethod>();

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
        using FluxCalculator = WMpfaDarcysLawsFluxCalculator<method>;

        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");

        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& dataHandle = fluxVarsCache.dataHandle();

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        const auto& subFluxDataIJ = dataHandle.subFluxData();

        Scalar flux = 0.0;

        if (enableGravity)
        {
            // do averaging for the density over all neighboring elements
            const auto rho = scvf.boundary() ? outsideVolVars.density(phaseIdx)
                                             : (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;

            const auto& g = problem.spatialParams().gravity(scvf.ipGlobal());

            auto pseudoPotential = [&elemVolVars, phaseIdx, &g, &rho](const auto& volVars, const auto& x)
                                        { return volVars.pressure(phaseIdx) - rho*(g*x);};

            if(!scvf.boundary())
            {
                const auto& flippedScvf = fvGeometry.flipScvf(scvf.index());
                const auto& dataHandleJ = elemFluxVarsCache[flippedScvf].dataHandle();
                const auto& subFluxDataJI = dataHandleJ.subFluxData();
                FluxCalculator::flux(flux, scvf, elemVolVars, subFluxDataIJ, subFluxDataJI, pseudoPotential);
            }
            else
                FluxCalculator::boundaryFlux(flux, scvf, elemVolVars, subFluxDataIJ, pseudoPotential);
        }
        else
        {
            auto pressure = [&elemVolVars, phaseIdx](const auto& volVars, const auto& x)
                                { return volVars.pressure(phaseIdx);};

            if(!scvf.boundary())
            {
                const auto& flippedScvf = fvGeometry.flipScvf(scvf.index());
                const auto& dataHandleJ = elemFluxVarsCache[flippedScvf].dataHandle();
                const auto& subFluxDataJI = dataHandleJ.subFluxData();
                FluxCalculator::flux(flux, scvf, elemVolVars, subFluxDataIJ, subFluxDataJI, pressure);
            }
            else
                FluxCalculator::boundaryFlux(flux, scvf, elemVolVars, subFluxDataIJ, pressure);
        }

        return flux;
    }
};

} // end namespace Dumux

#endif
