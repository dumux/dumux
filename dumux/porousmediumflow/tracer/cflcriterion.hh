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
 *
 * \brief Velocity output for porous media models
 */
#ifndef DUMUX_CFL_TRACER_HH
#define DUMUX_CFL_TRACER_HH

#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux {

/*!
 * \brief CFL criterion for tracer model
 */
template<class GridVariables, class FluxVariables, class ModelTraits>
class CFLTracer
{
    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using VolumeVariables = typename GridVariables::VolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using FluidSystem = typename VolumeVariables::FluidSystem;
    using Scalar = typename GridVariables::Scalar;

    // TODO should be possible to get this
    using Problem = typename std::decay_t<decltype(std::declval<GridVolumeVariables>().problem())>;
    using BoundaryTypes = typename std::decay_t<decltype(std::declval<GridVolumeVariables>().problem()
                                                         .boundaryTypes(std::declval<Element>(), std::declval<SubControlVolumeFace>()))>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr bool isBox = FVGridGeometry::discMethod == DiscretizationMethod::box;
    static constexpr int dofCodim = isBox ? dim : 0;

    static constexpr auto numEq = ModelTraits::numEq();
    static constexpr bool useMoles = ModelTraits::useMoles();

    using NumEqVector = Dune::FieldVector<Scalar, numEq>;

public:

    /*!
     * \brief Constructor initializes the static data with the initial solution.
     *
     * \param gridVariables The grid variables
     */
    CFLTracer(const GridVariables& gridVariables)
    : problem_(gridVariables.curGridVolVars().problem())
    , fvGridGeometry_(gridVariables.fvGridGeometry())
    , gridVariables_(gridVariables)
    {}

    template <class SolutionVector>
    Scalar suggestTimeStepSize(SolutionVector& sol)
    {
        Scalar deltaT = std::numeric_limits<Scalar>::infinity();

        for (const auto& element : elements(fvGridGeometry_.gridView(), Dune::Partitions::interior))
        {
            const auto eIdxGlobal = fvGridGeometry_.elementMapper().index(element);

            auto fvGeometry = localView(fvGridGeometry_);
            fvGeometry.bind(element);

            auto elemVolVars = localView(gridVariables_.curGridVolVars());
            elemVolVars.bind(element, fvGeometry, sol);

            deltaT = std::min(deltaT,suggestLocalTimeStepSize(elemVolVars, fvGeometry, element));
        }

        return deltaT;
    }

private:
    //Todo Implement also boundary handling
    Scalar suggestLocalTimeStepSize(const ElementVolumeVariables& elemVolVars,
                                    const FVElementGeometry& fvGeometry,
                                    const Element& element)
    {
        const auto geometry = element.geometry();

        // bind the element flux variables cache
        auto elemFluxVarsCache = localView(gridVariables_.gridFluxVarsCache());
        elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

        const auto& scv = fvGeometry.scv(fvGridGeometry_.elementMapper().index(element));
        const auto& volVar = elemVolVars[scv];

        NumEqVector flux(0.0);
        Scalar deltaT = std::numeric_limits<Scalar>::infinity();

        using std::abs;
        for (auto&& scvf : scvfs(fvGeometry))
        {
            FluxVariables fluxVars;
            fluxVars.init(problem_, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

            // formulation with mole balances
            if (useMoles)
            {
                for (int compIdx = 0; compIdx < numEq; ++compIdx)
                {
                    // the physical quantities for which we perform upwinding
                    auto upwindTerm = [compIdx](const VolumeVariables& volVars)
                    { return volVars.molarDensity(); };

                    // advective fluxes
                    flux[compIdx] += abs(fluxVars.advectiveFlux(0, upwindTerm));
                }
            }
            // formulation with mass balances
            else
            {
                for (int compIdx = 0; compIdx < numEq; ++compIdx)
                {
                    // the physical quantities for which we perform upwinding
                    auto upwindTerm = [compIdx](const VolumeVariables& volVars)
                    { return volVars.density(); };

                    // advective fluxes
                    flux[compIdx] += abs(fluxVars.advectiveFlux(0, upwindTerm));
                }
            }
        }

        using std::min;
        for (int compIdx = 0; compIdx < numEq; ++compIdx)
        {
            Scalar density = volVar.density(0);
            if(useMoles)
                density = volVar.molarDensity(0);

            // The factor 2 is needed such the the CFL criterion ist exact for constant velocities on cubic grids
            deltaT = min(deltaT,2.0*(element.geometry().volume()*volVar.porosity()*density)/flux[compIdx]);
        }

        return deltaT;
    }

    // The following SFINAE enable_if usage allows compilation, even if only a
    //
    // boundaryTypes(const Element&, const scv&)
    //
    // is provided in the problem file. In that case, the compiler cannot detect
    // (without additional measures like "using...") the signature
    //
    // boundaryTypes(const Element&, const scvf&)
    //
    // in the problem base class. Therefore, calls to this method trigger a
    // compiler error. However, that call is needed for calculating velocities
    // if the cell-centered discretization is used. By proceeding as in the
    // following lines, that call will only be compiled if cell-centered
    // actually is used.
    template <bool enable = isBox, typename std::enable_if_t<!enable, int> = 0>
    BoundaryTypes problemBoundaryTypes_(const Element& element, const SubControlVolumeFace& scvf) const
    { return problem_.boundaryTypes(element, scvf); }

    //! we should never call this method for box models
    template <bool enable = isBox, typename std::enable_if_t<enable, int> = 0>
    BoundaryTypes problemBoundaryTypes_(const Element& element, const SubControlVolumeFace& scvf) const
    { return BoundaryTypes(); }

    const Problem& problem_;
    const FVGridGeometry& fvGridGeometry_;
    const GridVariables& gridVariables_;
};

} // end namespace Dumux

#endif
