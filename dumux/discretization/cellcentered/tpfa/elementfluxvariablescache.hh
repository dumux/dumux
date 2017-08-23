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
#include <dumux/discretization/cellcentered/tpfa/fluxvariablescachefiller.hh>

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

    // Specialization for the global caching being enabled - do nothing here
    void update(const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars) {}

    // access operators in the case of caching
    const FluxVariablesCache& operator [](const SubControlVolumeFace& scvf) const
    { return (*globalFluxVarsCachePtr_)[scvf.index()]; }

    //! The global object we are a restriction of
    const GlobalFluxVariablesCache& globalFluxVarsCache() const
    {  return *globalFluxVarsCachePtr_; }

private:
    const GlobalFluxVariablesCache* globalFluxVarsCachePtr_;
};

/*!
 * \ingroup ImplicitModel
 * \brief Spezialization when not using global caching
 */
template<class TypeTag>
class CCTpfaElementFluxVariablesCache<TypeTag, false>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using GlobalFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GlobalFluxVariablesCache);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FluxVariablesCacheFiller = CCTpfaFluxVariablesCacheFiller<TypeTag>;

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

        // instantiate helper class to fill the caches
        FluxVariablesCacheFiller filler(globalFluxVarsCache().problem_());

        IndexType localScvfIdx = 0;
        // fill the containers
        for (auto&& scvf : scvfs(fvGeometry))
        {
            filler.fill(*this, fluxVarsCache_[localScvfIdx], element, fvGeometry, elemVolVars, scvf);
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

        // instantiate helper class to fill the caches
        FluxVariablesCacheFiller filler(problem);

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
            filler.fill(*this, fluxVarsCache_[localScvfIdx], element, fvGeometry, elemVolVars, scvf);
            globalScvfIndices_[localScvfIdx] = scvf.index();
            localScvfIdx++;
        }

        // add required data on the scv faces in the neighboring elements
        for (unsigned int localIdxJ = 0; localIdxJ < numNeighbors; ++localIdxJ)
        {
            const auto elementJ = fvGeometry.fvGridGeometry().element(assemblyMapI[localIdxJ].globalJ);
            for (auto scvfIdx : assemblyMapI[localIdxJ].scvfsJ)
            {
                auto&& scvfJ = fvGeometry.scvf(scvfIdx);
                filler.fill(*this, fluxVarsCache_[localScvfIdx], elementJ, fvGeometry, elemVolVars, scvfJ);
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

        // instantiate helper class to fill the caches
        FluxVariablesCacheFiller filler(globalFluxVarsCache().problem_());

        filler.fill(*this, fluxVarsCache_[0], element, fvGeometry, elemVolVars, scvf);
        globalScvfIndices_[0] = scvf.index();
    }

    // This function updates the transmissibilities after the solution has been deflected during jacobian assembly
    void update(const Element& element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars)
    {
        static const bool isSolIndependent = FluxVariablesCacheFiller::isSolutionIndependent();

        if (!isSolIndependent)
        {
            const auto& problem = globalFluxVarsCache().problem_();
            const auto globalI = problem.elementMapper().index(element);

            // instantiate filler class
            FluxVariablesCacheFiller filler(problem);

            // let the filler class update the caches
            for (unsigned int localScvfIdx = 0; localScvfIdx < fluxVarsCache_.size(); ++localScvfIdx)
            {
                const auto& scvf = fvGeometry.scvf(globalScvfIndices_[localScvfIdx]);

                const auto scvfInsideScvIdx = scvf.insideScvIdx();
                const auto& insideElement = scvfInsideScvIdx == globalI ?
                                            element :
                                            problem.model().fvGridGeometry().element(scvfInsideScvIdx);

                filler.fill(*this, fluxVarsCache_[localScvfIdx], insideElement, fvGeometry, elemVolVars, scvf);
            }
        }
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
