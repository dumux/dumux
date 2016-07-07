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
 * \brief Base class for the volume variables vector
 */
#ifndef DUMUX_DISCRETIZATION_CC_FLUXVARSCACHEVECTOR_HH
#define DUMUX_DISCRETIZATION_CC_FLUXVARSCACHEVECTOR_HH

#include <dumux/implicit/properties.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the flux variables cache vector, we store one cache per face
 */
template<class TypeTag>
class CCFluxVariablesCacheVector
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Element = typename GridView::template Codim<0>::Entity;
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

public:
    // When global caching is enabled, precompute transmissibilities and stencils for all the scv faces
    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache)>::type update(Problem& problem)
    {
        problemPtr_ = &problem;

        fluxVarsCache_.resize(problem.model().fvGeometries().numScvf());
        for (const auto& element : elements(problem.gridView()))
        {
            // Prepare the geometries within the elements of the stencil
            problem.model().fvGeometries_().bind(element);

            const auto& fvGeometry = problem.model().fvGeometries(element);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                (*this)[scvf].update(problem, element, scvf);
            }
        }
    }

    // When global flux variables caching is disabled, we don't need to update the cache
    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache)>::type update(Problem& problem)
    {
        problemPtr_ = &problem;
    }

    // This function has to be called prior to flux calculations on the element.
    // Prepares the transmissibilities of the scv faces in an element. The FvGeometry is assumed to be bound.
    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache)>::type bindElement(const Element& element)
    {
        const auto& fvGeometry = problem_().model().fvGeometries(element);

        // resizing of the cache
        const auto numScvf = fvGeometry.numScvf();
        fluxVarsCache_.resize(numScvf);
        globalScvfIndices_.resize(numScvf);
        IndexType localScvfIdx = 0;

        // fill the containers
        for (const auto& scvf : scvfs(fvGeometry))
        {
            fluxVarsCache_[localScvfIdx].update(problem_(), element, scvf);
            globalScvfIndices_[localScvfIdx] = scvf.index();
            localScvfIdx++;
        }
    }

    // Specialization for the global caching being enabled - do nothing here
    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache)>::type bindElement(const Element& element) {}

    // This function is called by the CCLocalResidual before flux calculations during assembly.
    // Prepares the transmissibilities of the scv faces in the stencil. The FvGeometries are assumed to be bound.
    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache)>::type bind(const Element& element)
    {
        const auto globalI = problem_().elementMapper().index(element);
        const auto& fvGeometry = problem_().model().fvGeometries(element);
        const auto& neighborStencil = problem_().model().stencils(element).neighborStencil();
        const auto numNeighbors = neighborStencil.size();

        // find the number of scv faces that need to be prepared
        IndexType numScvf = fvGeometry.numScvf();
        for (IndexType localIdxJ = 0; localIdxJ < numNeighbors; ++localIdxJ)
        {
            const auto& fluxVarIndicesJ = problem_().model().localJacobian().assemblyMap()[globalI][localIdxJ];
            numScvf += fluxVarIndicesJ.size();
        }

        // fill the containers with the data on the scv faces inside the actual element
        fluxVarsCache_.resize(numScvf);
        globalScvfIndices_.resize(numScvf);
        IndexType localScvfIdx = 0;
        for (const auto& scvf : scvfs(fvGeometry))
        {
            fluxVarsCache_[localScvfIdx].update(problem_(), element, scvf);
            globalScvfIndices_[localScvfIdx] = scvf.index();
            localScvfIdx++;
        }

        // add required data on the scv faces in the neighboring elements
        for (IndexType localIdxJ = 0; localIdxJ < numNeighbors; ++localIdxJ)
        {
            const auto& fluxVarIndicesJ = problem_().model().localJacobian().assemblyMap()[globalI][localIdxJ];
            const auto elementJ = problem_().model().fvGeometries().element(neighborStencil[localIdxJ]);
            for (auto fluxVarIdx : fluxVarIndicesJ)
            {
                const auto& scvfJ = problem_().model().fvGeometries().subControlVolumeFace(fluxVarIdx);
                fluxVarsCache_[localScvfIdx].update(problem_(), elementJ, scvfJ);
                globalScvfIndices_[localScvfIdx] = scvfJ.index();
                localScvfIdx++;
            }
        }
    }

    // Specialization for the global caching being enabled - do nothing here
    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache)>::type bind(const Element& element) {}

    // access operators in the case of caching
    template <typename T = TypeTag>
    const typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache), FluxVariablesCache>::type&
    operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[scvf.index()]; }

    template <typename T = TypeTag>
    typename std::enable_if<GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache), FluxVariablesCache>::type&
    operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[scvf.index()]; }

    // access operators in the case of no caching
    template <typename T = TypeTag>
    const typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache), FluxVariablesCache>::type&
    operator [](const SubControlVolumeFace& scvf) const
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

    template <typename T = TypeTag>
    typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache), FluxVariablesCache>::type&
    operator [](const SubControlVolumeFace& scvf)
    { return fluxVarsCache_[getLocalScvfIdx_(scvf.index())]; }

private:
    // get index of scvf in the local container
    template <typename T = TypeTag>
    const typename std::enable_if<!GET_PROP_VALUE(T, EnableGlobalFluxVariablesCache), int>::type
    getLocalScvfIdx_(const int scvfIdx) const
    {
        auto it = std::find(globalScvfIndices_.begin(), globalScvfIndices_.end(), scvfIdx);

        if (it != globalScvfIndices_.end())
            return std::distance(globalScvfIndices_.begin(), it);
        else
            DUNE_THROW(Dune::InvalidStateException, "Could not find the flux vars cache for scvfIdx = " << scvfIdx);
    }

    const Problem& problem_() const
    { return *problemPtr_; }

    const Problem* problemPtr_;

    std::vector<FluxVariablesCache> fluxVarsCache_;
    std::vector<IndexType> globalScvfIndices_;
};

} // end namespace

#endif
