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
 * \ingroup CCWMpfaDiscretization
 * \brief Flux variable caches on a gridview
 */
#ifndef DUMUX_DISCRETIZATION_CCWMPFA_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_CCWMPFA_GRID_FLUXVARSCACHE_HH

// make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>
#include "elementfluxvariablescache.hh"

namespace Dumux {

template<class ModelTraits>
struct DataHandlePhysicsTraits
{
    static constexpr bool enableAdvection = ModelTraits::enableAdvection();
    static constexpr bool enableMolecularDiffusion = ModelTraits::enableMolecularDiffusion();
    static constexpr bool enableHeatConduction = ModelTraits::enableEnergyBalance();

    static constexpr int numPhases = ModelTraits::numFluidPhases();
    static constexpr int numComponents = ModelTraits::numFluidComponents();
};

template<class P,
         class FVC,
         class FVCF,
         class FDH,
         class IOP>
struct CCWMpfaDefaultGridFluxVariablesCacheTraits
{
    using Problem = P;
    using FluxVariablesCache = FVC;
    using FluxVariablesCacheFiller = FVCF;

    using InterpolationOperator = IOP;
    using FaceDataHandle = FDH;

    template<class GridFluxVariablesCache, bool cachingEnabled>
    using LocalView = CCWMpfaElementFluxVariablesCache<GridFluxVariablesCache, cachingEnabled>;
};

/*!
 * \ingroup CCWMpfaDiscretization
 * \brief Flux variable caches on a gridview
 * \note The class is specialized for a version with and without grid caching
 */
template<class Traits, bool cachingEnabled>
class CCWMpfaGridFluxVariablesCache;

/*!
 * \ingroup CCWMpfaDiscretization
 * \brief Flux variable caches on a gridview with grid caching enabled
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 */
template<class TheTraits>
class CCWMpfaGridFluxVariablesCache<TheTraits, true>
{
    using Problem = typename TheTraits::Problem;
    using ThisType = CCWMpfaGridFluxVariablesCache<TheTraits, true>;

    //! the flux variable cache filler type
    using FluxVariablesCacheFiller = typename TheTraits::FluxVariablesCacheFiller;
public:
    //! export the Traits
    using Traits = TheTraits;

    //! export the data handle types used
    using InterpolationOperator = typename Traits::InterpolationOperator;
    using FaceDataHandle = typename Traits::FaceDataHandle;

    //! export the flux variable cache type
    using FluxVariablesCache = typename Traits::FluxVariablesCache;

    //! make it possible to query if caching is enabled
    static constexpr bool cachingEnabled = true;

    //! export the type of the local view
    using LocalView = typename Traits::template LocalView<ThisType, cachingEnabled>;

    //! The constructor
    CCWMpfaGridFluxVariablesCache(const Problem& problem)
    : problemPtr_(&problem)
    {}

    //! When global caching is enabled, precompute transmissibilities for all scv faces
    template<class GridGeometry, class GridVolumeVariables, class SolutionVector>
    void update(const GridGeometry& gridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false)
    {
        // Update only if the filler puts solution-dependent
        // stuff into the caches or if update is enforced
        if (FluxVariablesCacheFiller::isSolDependent || forceUpdate)
        {
            // clear previous data if forced update is desired
            if (forceUpdate)
            {
                clear_();

                dataStorage_.interpolationOperators.resize(gridGeometry.gridView().size(0));
                dataStorage_.facesDataHandles.resize(gridGeometry.numScvf());

                // reserve memory estimate for caches, interaction volumes and corresponding data
                fluxVarsCache_.resize(gridGeometry.numScvf());
            }

            // instantiate helper class to fill the caches
            FluxVariablesCacheFiller filler(problem());

            for (const auto& element : elements(gridGeometry.gridView()))
            {
                auto fvGeometry = localView(gridGeometry);
                fvGeometry.bind(element);

                auto elemVolVars = localView(gridVolVars);
                elemVolVars.bind(element, fvGeometry, sol);

                auto& intOperator = dataStorage_.interpolationOperators[gridGeometry.elementMapper().index(element)];
                intOperator.bind(element, fvGeometry, elemVolVars);

                // Prepare all caches of the scvfs inside the corresponding interaction volume. Skip
                // those ivs that are touching a boundary, we only store the data on interior ivs here.
                for (const auto& scvf : scvfs(fvGeometry))
                    //if (!isEmbeddedInBoundaryIV_(scvf, gridGeometry) && !fluxVarsCache_[scvf.index()].isUpdated())
                        filler.fill(*this, fluxVarsCache_[scvf.index()], dataStorage_, element, fvGeometry, elemVolVars, scvf, forceUpdate);
            }
        }
    }

    template<class FVElementGeometry, class ElementVolumeVariables>
    void updateElement(const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars)
    {
        // // Update only if the filler puts
        // // solution-dependent stuff into the caches
        // if (FluxVariablesCacheFiller::isSolDependent)
        // {
        //     const auto& gridGeometry = fvGeometry.gridGeometry();
        //     const auto& assemblyMapI = gridGeometry.connectivityMap()[gridGeometry.elementMapper().index(element)];

        //     // helper class to fill flux variables caches
        //     FluxVariablesCacheFiller filler(problem());

        //     // first, set all the caches to "outdated"
        //     for (const auto& scvf : scvfs(fvGeometry))
        //         fluxVarsCache_[scvf.index()].setUpdateStatus(false);
        //     for (const auto& dataJ : assemblyMapI)
        //         for (const auto scvfIdx : dataJ.scvfsJ)
        //             fluxVarsCache_[scvfIdx].setUpdateStatus(false);

        //     // go through the caches maybe update them
        //     for (const auto& scvf : scvfs(fvGeometry))
        //     {
        //         auto& scvfCache = fluxVarsCache_[scvf.index()];
        //         if (!isEmbeddedInBoundaryIV_(scvf, gridGeometry) && !scvfCache.isUpdated())
        //             filler.fill(*this, scvfCache, ivDataStorage_, element, fvGeometry, elemVolVars, scvf);
        //     }

        //     for (const auto& dataJ : assemblyMapI)
        //     {
        //         const auto elementJ = gridGeometry.element(dataJ.globalJ);
        //         for (const auto scvfIdx : dataJ.scvfsJ)
        //         {
        //             auto& scvfCache = fluxVarsCache_[scvfIdx];
        //             const auto& scvf = fvGeometry.scvf(scvfIdx);
        //             if (!isEmbeddedInBoundaryIV_(scvf, gridGeometry) && !scvfCache.isUpdated())
        //                 filler.fill(*this, scvfCache, ivDataStorage_, elementJ, fvGeometry, elemVolVars, scvf);
        //         }
        //     }
        // }
    }

    //! access operators in the case of caching
    template<class SubControlVolumeFace>
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[scvf.index()]; }

    //! access operators in the case of caching
    template<class SubControlVolumeFace>
    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[scvf.index()]; }

    template<class SubControlVolume>
    const InterpolationOperator& interpolationOperator(const SubControlVolume& scv) const
    { return dataStorage_.interpolationOperators[scv.elementIndex()]; }

    template<class SubControlVolumeFace>
    const FaceDataHandle& faceDataHandle(const SubControlVolumeFace& scvf) const
    { return dataStorage_.facesDataHandles[scvf.index()]; }

    const Problem& problem() const
    { return *problemPtr_; }

private:

    //! clear all containers
    void clear_()
    {
        fluxVarsCache_.clear();
        dataStorage_.interpolationOperators.clear();
        dataStorage_.facesDataHandles.clear();
    }

    const Problem* problemPtr_;
    std::vector<FluxVariablesCache> fluxVarsCache_;

    // stored interaction volumes and handles
    using DS = DataStorage<InterpolationOperator,
                                    FaceDataHandle>;
    DS dataStorage_;
};

} // end namespace Dumux

#endif
