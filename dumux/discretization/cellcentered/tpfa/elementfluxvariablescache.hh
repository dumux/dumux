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
 * \brief The global object of flux var caches
 */
#ifndef DUMUX_DISCRETIZATION_CCTPFA_ELEMENT_FLUXVARSCACHE_HH
#define DUMUX_DISCRETIZATION_CCTPFA_ELEMENT_FLUXVARSCACHE_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the stencil local flux variables cache
 */
template<class TypeTag, bool EnableGlobalFluxVariablesCache>
class CCTpfaElementFluxVariablesCache;

/*!
 * \ingroup ImplicitModel
 * \brief Spezialization when caching globally
 */
template<class TypeTag>
class CCTpfaElementFluxVariablesCache<TypeTag, true>
{
    // the local jacobian needs to be able to update the cache during assembly
    friend typename GET_PROP_TYPE(TypeTag, LocalJacobian);

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using GlobalFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GlobalFluxVariablesCache);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

public:
    CCTpfaElementFluxVariablesCache(const GlobalFluxVariablesCache& global)
    : globalFluxVarsCachePtr_(&global) {}

    // Specialization for the global caching being enabled - do nothing here
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars) {}

    // Specialization for the global caching being enabled - do nothing here
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars) {}

    // Specialization for the global caching being enabled - do nothing here
    void bindScvf(const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf) {}

    // access operators in the case of caching
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return (*globalFluxVarsCachePtr_)[scvf.index()]; }

    //! The global object we are a restriction of
    const GlobalFluxVariablesCache& globalFluxVarsCache() const
    {  return *globalFluxVarsCachePtr_; }

private:
    const GlobalFluxVariablesCache* globalFluxVarsCachePtr_;

    // Specialization for the global caching being enabled - do nothing here
    void update(const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars) {}
};

/*!
 * \ingroup ImplicitModel
 * \brief Spezialization when not using global caching
 */
template<class TypeTag>
class CCTpfaElementFluxVariablesCache<TypeTag, false>
{
    // the local jacobian needs to be able to update the cache during assembly
    friend typename GET_PROP_TYPE(TypeTag, LocalJacobian);

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using GlobalFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GlobalFluxVariablesCache);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    static const bool solDependentAdvection = GET_PROP_VALUE(TypeTag, SolutionDependentAdvection);

public:
    CCTpfaElementFluxVariablesCache(const GlobalFluxVariablesCache& global)
    : globalFluxVarsCachePtr_(&global) {}

    // This function has to be called prior to flux calculations on the element.
    // Prepares the transmissibilities of the scv faces in an element. The FvGeometry is assumed to be bound.
    void bindElement(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars)
    {
        // resizing of the cache
        const auto numScvf = fvGeometry.numScvf();
        fluxVarsCache_.resize(numScvf);
        globalScvfIndices_.resize(numScvf);
        IndexType localScvfIdx = 0;

        // fill the containers
        for (auto&& scvf : scvfs(fvGeometry))
        {
            fluxVarsCache_[localScvfIdx].update(globalFluxVarsCache().problem_(), element, fvGeometry, elemVolVars, scvf);
            globalScvfIndices_[localScvfIdx] = scvf.index();
            localScvfIdx++;
        }
    }

    // This function is called by the CCLocalResidual before flux calculations during assembly.
    // Prepares the transmissibilities of the scv faces in the stencil. The FvGeometries are assumed to be bound.
    void bind(const Element& element,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars)
    {
        const auto& problem = globalFluxVarsCache().problem_();
        const auto globalI = problem.elementMapper().index(element);
        const auto& assemblyMapI = problem.model().localJacobian().assemblyMap()[globalI];
        const auto numNeighbors = assemblyMapI.size();

        // find the number of scv faces that need to be prepared
        auto numScvf = fvGeometry.numScvf();
        for (unsigned int localIdxJ = 0; localIdxJ < numNeighbors; ++localIdxJ)
            numScvf += assemblyMapI[localIdxJ].scvfsJ.size();

        // fill the containers with the data on the scv faces inside the actual element
        fluxVarsCache_.resize(numScvf);
        globalScvfIndices_.resize(numScvf);
        unsigned int localScvfIdx = 0;
        for (auto&& scvf : scvfs(fvGeometry))
        {
            fluxVarsCache_[localScvfIdx].update(problem, element, fvGeometry, elemVolVars, scvf);
            globalScvfIndices_[localScvfIdx] = scvf.index();
            localScvfIdx++;
        }

        // add required data on the scv faces in the neighboring elements
        for (unsigned int localIdxJ = 0; localIdxJ < numNeighbors; ++localIdxJ)
        {
            const auto elementJ = fvGeometry.globalFvGeometry().element(assemblyMapI[localIdxJ].globalJ);
            for (auto scvfIdx : assemblyMapI[localIdxJ].scvfsJ)
            {
                auto&& scvfJ = fvGeometry.scvf(scvfIdx);
                fluxVarsCache_[localScvfIdx].update(problem, elementJ, fvGeometry, elemVolVars, scvfJ);
                globalScvfIndices_[localScvfIdx] = scvfJ.index();
                localScvfIdx++;
            }
        }
    }

    void bindScvf(const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const SubControlVolumeFace& scvf)
    {
        fluxVarsCache_.resize(1);
        globalScvfIndices_.resize(1);

        fluxVarsCache_[0].update(globalFluxVarsCache().problem_(), element, fvGeometry, elemVolVars, scvf);
        globalScvfIndices_[0] = scvf.index();
    }

    // access operators in the case of no caching
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    FluxVariablesCache& operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    //! The global object we are a restriction of
    const GlobalFluxVariablesCache& globalFluxVarsCache() const
    {  return *globalFluxVarsCachePtr_; }

private:
    const GlobalFluxVariablesCache* globalFluxVarsCachePtr_;

    // This function updates the transmissibilities after the solution has been deflected during jacobian assembly
    void update(const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars)
    {
        if (solDependentAdvection)
            bind(element, fvGeometry, elemVolVars);
    }

    // get index of scvf in the local container
    int getLocalScvfIdx_(const int scvfIdx) const
    {
        auto it = std::find(globalScvfIndices_.begin(), globalScvfIndices_.end(), scvfIdx);
        assert(it != globalScvfIndices_.end() && "Could not find the flux vars cache for scvfIdx");
        return std::distance(globalScvfIndices_.begin(), it);
    }

    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<IndexType> globalScvfIndices_;
};

} // end namespace

#endif
