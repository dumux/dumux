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
 * \brief Base class for a local finite volume geometry for cell-centered TPFA models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element in the local scope we are restricting to, e.g. stencil or element.
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_CCMPFA_FV_ELEMENT_GEOMETRY_HH

#include <dune/common/iteratorrange.hh>

#include <dumux/discretization/scvandscvfiterators.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>

namespace Dumux
{
/*!
 * \ingroup ImplicitModel
 * \brief Base class for the finite volume geometry vector for cell-centered TPFA models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class TypeTag, bool EnableGlobalFVGeometryCache>
class CCMpfaFVElementGeometry
{};

//! specialization in case the FVElementGeometries are stored globally
//! In this case we just forward internally to the global object
template<class TypeTag>
class CCMpfaFVElementGeometry<TypeTag, true>
{
    using ThisType = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalFVGeometry = typename GET_PROP_TYPE(TypeTag, GlobalFVGeometry);

    using ScvIterator = Dumux::ScvIterator<SubControlVolume, std::vector<IndexType>, ThisType>;
    using ScvfIterator = Dumux::ScvfIterator<SubControlVolumeFace, std::vector<IndexType>, ThisType>;

public:
    //! Constructor
    CCMpfaFVElementGeometry(const GlobalFVGeometry& globalFvGeometry)
    : globalFvGeometryPtr_(&globalFvGeometry) {}

    //! Get an element sub control volume with a global scv index
    const SubControlVolume& scv(IndexType scvIdx) const
    {
        return globalFvGeometry().scv(scvIdx);
    }

    //! Get an element sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(IndexType scvfIdx) const
    {
        return globalFvGeometry().scvf(scvfIdx);
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element (not including neighbor scvs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<ScvIterator>
    scvs(const CCMpfaFVElementGeometry& fvGeometry)
    {
        return Dune::IteratorRange<ScvIterator>(ScvIterator(fvGeometry.scvIndices_.begin(), fvGeometry),
                                                ScvIterator(fvGeometry.scvIndices_.end(), fvGeometry));
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element (not including neighbor scvfs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange<ScvfIterator>
    scvfs(const CCMpfaFVElementGeometry& fvGeometry)
    {
        const auto& g = fvGeometry.globalFvGeometry();
        const auto scvIdx = fvGeometry.scvIndices_[0];
        return Dune::IteratorRange<ScvfIterator>(ScvfIterator(g.scvfIndicesOfScv(scvIdx).begin(), fvGeometry),
                                                 ScvfIterator(g.scvfIndicesOfScv(scvIdx).end(), fvGeometry));
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScv() const
    {
        return scvIndices_.size();
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScvf() const
    {
        return globalFvGeometry().scvfIndicesOfScv(scvIndices_[0]).size();
    }

    //! Binding of an element, called by the local jacobian to prepare element assembly
    void bind(const Element& element)
    {
        this->bindElement(element);
    }

    //! Bind only element-local
    void bindElement(const Element& element)
    {
        elementPtr_ = &element;
        scvIndices_ = std::vector<IndexType>({globalFvGeometry().problem_().elementMapper().index(*elementPtr_)});
    }

    //! The global finite volume geometry we are a restriction of
    const GlobalFVGeometry& globalFvGeometry() const
    { return *globalFvGeometryPtr_; }

private:

    const Element* elementPtr_;
    std::vector<IndexType> scvIndices_;
    const GlobalFVGeometry* globalFvGeometryPtr_;
};

//! specialization in case the FVElementGeometries are not stored
template<class TypeTag>
class CCMpfaFVElementGeometry<TypeTag, false>
{
    using ThisType = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalFVGeometry = typename GET_PROP_TYPE(TypeTag, GlobalFVGeometry);

    using ScvIterator = Dumux::ScvIterator<SubControlVolume, std::vector<IndexType>, ThisType>;
    using ScvfIterator = Dumux::ScvfIterator<SubControlVolumeFace, std::vector<IndexType>, ThisType>;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using CoordScalar = typename GridView::ctype;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    //! Constructor
    CCMpfaFVElementGeometry(const GlobalFVGeometry& globalFvGeometry)
    : globalFvGeometryPtr_(&globalFvGeometry) {}

    //! Get an elment sub control volume with a global scv index
    //! We separate element and neighbor scvs to speed up mapping
    const SubControlVolume& scv(IndexType scvIdx) const
    {
        if (scvIdx == scvIndices_[0])
            return scvs_[0];
        else
            return neighborScvs_[findLocalIndex(scvIdx, neighborScvIndices_)];
    }

    //! Get an element sub control volume face with a global scvf index
    //! We separate element and neighbor scvfs to speed up mapping
    const SubControlVolumeFace& scvf(IndexType scvfIdx) const
    {
        auto it = std::find(scvfIndices_.begin(), scvfIndices_.end(), scvfIdx);
        if (it != scvfIndices_.end())
            return scvfs_[std::distance(scvfIndices_.begin(), it)];
        else
            return neighborScvfs_[findLocalIndex(scvfIdx, neighborScvfIndices_)];
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element (not including neighbor scvs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolume>::const_iterator>
    scvs(const ThisType& g)
    {
        using Iter = typename std::vector<SubControlVolume>::const_iterator;
        return Dune::IteratorRange<Iter>(g.scvs_.begin(), g.scvs_.end());
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element (not including neighbor scvfs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolumeFace>::const_iterator>
    scvfs(const ThisType& g)
    {
        using Iter = typename std::vector<SubControlVolumeFace>::const_iterator;
        return Dune::IteratorRange<Iter>(g.scvfs_.begin(), g.scvfs_.end());
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScv() const
    { return scvs_.size(); }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScvf() const
    { return scvfs_.size(); }

    //! Binding of an element preparing the geometries of the whole stencil
    //! called by the local jacobian to prepare element assembly
    void bind(const Element& element)
    {
        bindElement(element);

        // get some references for convenience
        const auto& problem = globalFvGeometry().problem_();
        const auto& assemblyMap = problem.model().localJacobian().assemblyMap();
        const auto& neighborStencil = problem.model().stencils(element).neighborStencil();
        const auto globalI = problem.elementMapper().index(element);

        // reserve memory
        const auto numNeighbors = neighborStencil.size();
        const auto estimate = numNeighbors*dim*(4*dim-4); // estimate holds for quadrilaterals in 2D/3D, overestimates else
        neighborScvs_.reserve(numNeighbors);
        neighborScvIndices_.reserve(numNeighbors);
        neighborScvfIndices_.reserve(estimate);
        neighborScvfs_.reserve(estimate);

        // build scvfs in neighbors
        // use the assembly map to determine which faces are necessary
        unsigned int j = 0;
        for (auto globalJ : neighborStencil)
        {
            // get data on the neighbor
            auto elementJ = globalFvGeometry().element(globalJ);
            const auto& scvfIndices = assemblyMap[globalI][j];

            // make neighbor geometries
            makeNeighborGeometries(elementJ, globalJ, scvfIndices);

            // increment counter
            j++;
        }

        // free unused memory
        neighborScvs_.shrink_to_fit();
        neighborScvIndices_.shrink_to_fit();
        neighborScvfIndices_.shrink_to_fit();
        neighborScvfs_.shrink_to_fit();
    }

    //! Binding of an element preparing the geometries only inside the element
    void bindElement(const Element& element)
    {
        clear();
        elementPtr_ = &element;
        makeElementGeometries(element);
    }

    //! The global finite volume geometry we are a restriction of
    const GlobalFVGeometry& globalFvGeometry() const
    { return *globalFvGeometryPtr_; }

private:

    //! create scvs and scvfs of the bound element
    void makeElementGeometries(const Element& element)
    {
        // the problem
        const auto& problem = globalFvGeometry().problem_();

        // make the scv
        auto eIdx = problem.elementMapper().index(element);
        scvs_.emplace_back(element.geometry(), eIdx);
        scvIndices_.emplace_back(eIdx);

        // get data on the scv faces
        const auto& scvFaceIndices = globalFvGeometry().scvfIndicesOfScv(eIdx);
        const auto& neighborVolVarIndices = globalFvGeometry().neighborVolVarIndices(eIdx);

        // Instantiate helper class to pass it to scvf constructors
        MpfaHelper helper;

        // the quadrature point to be used on the scvf
        const Scalar q = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Mpfa, Q);

        // reserve memory for the scv faces
        auto numLocalScvf = scvFaceIndices.size();
        scvfIndices_.reserve(numLocalScvf);
        scvfs_.reserve(numLocalScvf);

        // for network grids we only want to do one scvf per half facet
        // this approach assumes conforming grids at branching facets
        std::vector<IndexType> finishedFacets;

        int scvfCounter = 0;
        for (const auto& is : intersections(globalFvGeometry().gridView(), element))
        {
            // if we are dealing with a lower dimensional network
            // only make a new scvf if we haven't handled it yet
            if (dim < dimWorld)
            {
                auto indexInInside = is.indexInInside();
                if(MpfaHelper::contains(finishedFacets, indexInInside))
                    continue;
                else
                    finishedFacets.push_back(indexInInside);
            }

            // get the intersection corners according to generic numbering
            auto numCorners = is.geometry().corners();
            std::vector<GlobalPosition> isCorners(numCorners);

            // if outside level > inside level, use the outside element in the following
            bool useNeighbor = is.neighbor() && is.outside().level() > element.level();
            const auto& e = useNeighbor ? is.outside() : element;
            const auto indexInElement = useNeighbor ? is.indexInOutside() : is.indexInInside();
            const auto eg = e.geometry();
            const auto& refElement = ReferenceElements::general(eg.type());

            for (unsigned int c = 0; c < numCorners; ++c)
                isCorners[c] = eg.global(refElement.position(refElement.subEntity(indexInElement, 1, c, dim), dim));

            // make the scv faces belonging to each corner of the intersection
            for (unsigned int c = 0; c < numCorners; ++c)
            {
                // get the global vertex index the scv face is connected to
                auto vIdxLocal = refElement.subEntity(indexInElement, 1, c, dim);
                auto vIdxGlobal = problem.vertexMapper().subIndex(e, vIdxLocal, dim);

                // do not build scvfs connected to a processor boundary
                if (globalFvGeometry().isGhostVertex(vIdxGlobal))
                    continue;

                scvfs_.emplace_back(helper,
                                    helper.getScvfCorners(isCorners, c),
                                    is.centerUnitOuterNormal(),
                                    vIdxGlobal,
                                    vIdxLocal,
                                    scvFaceIndices[scvfCounter],
                                    eIdx,
                                    neighborVolVarIndices[scvfCounter],
                                    q,
                                    is.boundary()
                                    );

                scvfIndices_.emplace_back(scvFaceIndices[scvfCounter]);
                scvfCounter++;
            }
        }
    }

    //! create the scv and necessary scvfs of a neighboring element
    template<typename IndexVector>
    void makeNeighborGeometries(const Element& element, IndexType eIdxGlobal, const IndexVector& scvfIndices)
    {
        // create the neighbor scv if it doesn't exist yet
        neighborScvs_.emplace_back(element.geometry(), eIdxGlobal);
        neighborScvIndices_.push_back(eIdxGlobal);

        // get data on the scv faces
        const auto& scvFaceIndices = globalFvGeometry().scvfIndicesOfScv(eIdxGlobal);
        const auto& neighborVolVarIndices = globalFvGeometry().neighborVolVarIndices(eIdxGlobal);

        // Instantiate a helper class to pass it to scvf constructors
        MpfaHelper helper;

        // the quadrature point to be used on the scvf
        const Scalar q = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Mpfa, Q);

        // for network grids we only want to do one scvf per half facet
        // this approach assumes conforming grids at branching facets
        std::vector<IndexType> finishedFacets;

        int scvfCounter = 0;
        for (const auto& is : intersections(globalFvGeometry().gridView(), element))
        {
            // if we are dealing with a lower dimensional network
            // only make a new scvf if we haven't handled it yet
            if (dim < dimWorld)
            {
                auto indexInInside = is.indexInInside();
                if(MpfaHelper::contains(finishedFacets, indexInInside))
                    continue;
                else
                    finishedFacets.push_back(indexInInside);
            }

            // get the intersection corners according to generic numbering
            auto numCorners = is.geometry().corners();
            std::vector<GlobalPosition> isCorners(numCorners);

            // if outside level > inside level, use the outside element in the following
            bool useNeighbor = is.neighbor() && is.outside().level() > element.level();
            const auto& e = useNeighbor ? is.outside() : element;
            const auto indexInElement = useNeighbor ? is.indexInOutside() : is.indexInInside();
            const auto eg = e.geometry();
            const auto& refElement = ReferenceElements::general(eg.type());

            for (unsigned int c = 0; c < numCorners; ++c)
                isCorners[c] = eg.global(refElement.position(refElement.subEntity(indexInElement, 1, c, dim), dim));

            // make the scv faces belonging to each corner of the intersection
            for (unsigned int c = 0; c < numCorners; ++c)
            {
                // get the global vertex index the scv face is connected to
                auto vIdxLocal = refElement.subEntity(indexInElement, 1, c, dim);
                auto vIdxGlobal = globalFvGeometry().problem_().vertexMapper().subIndex(e, vIdxLocal, dim);

                // do not build scvfs connected to a processor boundary
                if (globalFvGeometry().isGhostVertex(vIdxGlobal))
                    continue;

                // only build the scvf if it is in the list of necessary indices
                if (MpfaHelper::contains(scvfIndices, scvFaceIndices[scvfCounter]))
                {
                    neighborScvfs_.emplace_back(helper,
                                                helper.getScvfCorners(isCorners, c),
                                                is.centerUnitOuterNormal(),
                                                vIdxGlobal,
                                                vIdxLocal,
                                                scvFaceIndices[scvfCounter],
                                                eIdxGlobal,
                                                neighborVolVarIndices[scvfCounter],
                                                q,
                                                is.boundary()
                                                );

                    neighborScvfIndices_.emplace_back(scvFaceIndices[scvfCounter]);
                }

                // increment counter either way
                scvfCounter++;
            }
        }
    }

    const IndexType findLocalIndex(const IndexType idx,
                                   const std::vector<IndexType>& indices) const
    {
        auto it = std::find(indices.begin(), indices.end(), idx);
        assert(it != indices.end() && "Could not find the scv/scvf! Make sure to properly bind this class!");
        return std::distance(indices.begin(), it);
    }

    void clear()
    {
        scvIndices_.clear();
        scvfIndices_.clear();
        scvs_.clear();
        scvfs_.clear();

        neighborScvIndices_.clear();
        neighborScvfIndices_.clear();
        neighborScvs_.clear();
        neighborScvfs_.clear();
    }

    // the bound element
    const Element* elementPtr_;

    const GlobalFVGeometry* globalFvGeometryPtr_;

    // local storage after binding an element
    std::vector<IndexType> scvIndices_;
    std::vector<IndexType> scvfIndices_;
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;

    std::vector<IndexType> neighborScvIndices_;
    std::vector<IndexType> neighborScvfIndices_;
    std::vector<SubControlVolume> neighborScvs_;
    std::vector<SubControlVolumeFace> neighborScvfs_;
};

} // end namespace

#endif
