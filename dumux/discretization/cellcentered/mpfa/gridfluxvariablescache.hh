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
 * \ingroup CCMpfaDiscretization
 * \brief Flux variable caches on a gridview
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_GRID_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_CCMPFA_GRID_FLUXVARSCACHE_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/interactionvolumedatahandle.hh>
#include <dumux/discretization/cellcentered/mpfa/fluxvariablescachefiller.hh>

//! make the local view function available whenever we use this class
#include <dumux/discretization/localview.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Flux variable caches on a gridview
 * \note The class is specialized for a version with and without grid caching
 */
template<class TypeTag, bool EnableGridFluxVariablesCache>
class CCMpfaGridFluxVariablesCache;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Flux variable caches on a gridview with grid caching enabled
 * \note The flux caches of the gridview are stored which is memory intensive but faster
 */
template<class TypeTag>
class CCMpfaGridFluxVariablesCache<TypeTag, true>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using FluxVariablesCacheFiller = CCMpfaFluxVariablesCacheFiller<TypeTag>;
    using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
    using PrimaryMatVecTraits = typename PrimaryInteractionVolume::Traits::MatVecTraits;
    using SecondaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, SecondaryInteractionVolume);
    using SecondaryMatVecTraits = typename SecondaryInteractionVolume::Traits::MatVecTraits;

    //! physics traits class to define the data handles
    struct PhysicsTraits
    {
        static constexpr bool enableAdvection = GET_PROP_VALUE(TypeTag, EnableAdvection);
        static constexpr bool enableMolecularDiffusion = GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion);
        static constexpr bool enableHeatConduction = GET_PROP_VALUE(TypeTag, EnableEnergyBalance);

        static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
        static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    };

public:
    // export the data handle types used
    using PrimaryIvDataHandle = InteractionVolumeDataHandle< PrimaryMatVecTraits, PhysicsTraits >;
    using SecondaryIvDataHandle = InteractionVolumeDataHandle< SecondaryMatVecTraits, PhysicsTraits >;

    //! export the type of the local view
    using LocalView = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

    //! The constructor
    CCMpfaGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    //! When global caching is enabled, precompute transmissibilities for all scv faces
    void update(const FVGridGeometry& fvGridGeometry,
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

                const auto& gridIvIndexSets = fvGridGeometry.gridInteractionVolumeIndexSets();
                const auto numPrimaryIvs = gridIvIndexSets.numPrimaryInteractionVolumes();
                const auto numSecondaryIVs = gridIvIndexSets.numSecondaryInteractionVolumes();
                primaryInteractionVolumes_.reserve(numPrimaryIvs);
                secondaryInteractionVolumes_.reserve(numSecondaryIVs);
                primaryIvDataHandles_.reserve(numPrimaryIvs);
                secondaryIvDataHandles_.reserve(numSecondaryIVs);

                // reserve memory estimate for caches, interaction volumes and corresponding data
                fluxVarsCache_.resize(fvGridGeometry.numScvf());
            }

            // instantiate helper class to fill the caches
            FluxVariablesCacheFiller filler(problem());

            // set all the caches to "outdated"
            for (auto& cache : fluxVarsCache_)
                cache.setUpdateStatus(false);

            for (const auto& element : elements(fvGridGeometry.gridView()))
            {
                auto fvGeometry = localView(fvGridGeometry);
                fvGeometry.bind(element);

                auto elemVolVars = localView(gridVolVars);
                elemVolVars.bind(element, fvGeometry, sol);

                // prepare all the caches of the scvfs inside the corresponding interaction volume
                for (const auto& scvf : scvfs(fvGeometry))
                    if (!fluxVarsCache_[scvf.index()].isUpdated())
                        filler.fill(*this, fluxVarsCache_[scvf.index()], element, fvGeometry, elemVolVars, scvf, forceUpdate);
            }
        }
    }

    void updateElement(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars)
    {
        // Update only if the filler puts
        // solution-dependent stuff into the caches
        if (FluxVariablesCacheFiller::isSolDependent)
        {
            const auto& fvGridGeometry = fvGeometry.fvGridGeometry();
            const auto& assemblyMapI = fvGridGeometry.connectivityMap()[fvGridGeometry.elementMapper().index(element)];

            // helper class to fill flux variables caches
            FluxVariablesCacheFiller filler(problem());

            // first, set all the caches to "outdated"
            for (const auto& scvf : scvfs(fvGeometry))
                fluxVarsCache_[scvf.index()].setUpdateStatus(false);
            for (const auto& dataJ : assemblyMapI)
                for (const auto scvfIdx : dataJ.scvfsJ)
                    fluxVarsCache_[scvfIdx].setUpdateStatus(false);

            // go through the caches maybe update them
            for (const auto& scvf : scvfs(fvGeometry))
            {
                auto& scvfCache = fluxVarsCache_[scvf.index()];
                if (!scvfCache.isUpdated())
                    filler.fill(*this, scvfCache, element, fvGeometry, elemVolVars, scvf);
            }

            for (const auto& dataJ : assemblyMapI)
            {
                const auto elementJ = fvGridGeometry.element(dataJ.globalJ);
                for (const auto scvfIdx : dataJ.scvfsJ)
                {
                    auto& scvfCache = fluxVarsCache_[scvfIdx];
                    if (!scvfCache.isUpdated())
                    {
                        const auto& scvf = fvGeometry.scvf(scvfIdx);
                        filler.fill(*this, scvfCache, elementJ, fvGeometry, elemVolVars, scvf);
                    }
                }
            }
        }
    }

    //! access operators in the case of caching
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[scvf.index()]; }

    //! access operators in the case of caching
    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[scvf.index()]; }

    //! access to the stored interaction volumes
    const std::vector<PrimaryInteractionVolume>& primaryInteractionVolumes() const
    { return primaryInteractionVolumes_; }

    //! access to the stored interaction volumes
    std::vector<PrimaryInteractionVolume>& primaryInteractionVolumes()
    { return primaryInteractionVolumes_; }

    //! access to the stored data handles
    const std::vector<PrimaryIvDataHandle>& primaryDataHandles() const
    { return primaryIvDataHandles_; }

    //! access to the stored data handles
    std::vector<PrimaryIvDataHandle>& primaryDataHandles()
    { return primaryIvDataHandles_; }

    //! access to the stored interaction volumes
    const std::vector<SecondaryInteractionVolume>& secondaryInteractionVolumes() const
    { return secondaryInteractionVolumes_; }

    //! access to the stored interaction volumes
    std::vector<SecondaryInteractionVolume>& secondaryInteractionVolumes()
    { return secondaryInteractionVolumes_; }

    //! access to the stored data handles
    const std::vector<SecondaryIvDataHandle>& secondaryDataHandles() const
    { return secondaryIvDataHandles_; }

    //! access to the stored data handles
    std::vector<SecondaryIvDataHandle>& secondaryDataHandles()
    { return secondaryIvDataHandles_; }

    const Problem& problem() const
    { return *problemPtr_; }

private:

    //! clear all containers
    void clear_()
    {
        fluxVarsCache_.clear();
        primaryInteractionVolumes_.clear();
        secondaryInteractionVolumes_.clear();
        primaryIvDataHandles_.clear();
        secondaryIvDataHandles_.clear();
    }

    const Problem* problemPtr_;
    std::vector<FluxVariablesCache> fluxVarsCache_;

    // store the interaction volumes and handles
    std::vector<PrimaryInteractionVolume> primaryInteractionVolumes_;
    std::vector<SecondaryInteractionVolume> secondaryInteractionVolumes_;
    std::vector<PrimaryIvDataHandle> primaryIvDataHandles_;
    std::vector<SecondaryIvDataHandle> secondaryIvDataHandles_;
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Flux variable caches on a gridview with grid caching disabled
 */
template<class TypeTag>
class CCMpfaGridFluxVariablesCache<TypeTag, false>
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

public:
    //! export the type of the local view
    using LocalView = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

    //! export the data handle types used by the local view
    using PrimaryIvDataHandle = typename LocalView::PrimaryIvDataHandle;
    using SecondaryIvDataHandle = typename LocalView::SecondaryIvDataHandle;

    //! The constructor
    CCMpfaGridFluxVariablesCache(const Problem& problem) : problemPtr_(&problem) {}

    //! When global flux variables caching is disabled, we don't need to update the cache
    void update(const FVGridGeometry& fvGridGeometry,
                const GridVolumeVariables& gridVolVars,
                const SolutionVector& sol,
                bool forceUpdate = false) {}

    //! When global flux variables caching is disabled, we don't need to update the cache
    void updateElement(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars) {}

    const Problem& problem() const
    { return *problemPtr_; }

private:
    const Problem* problemPtr_;
};

} // end namespace Dumux

#endif
