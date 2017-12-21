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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredFVElementGeometry
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FV_ELEMENT_GEOMETRY_HH

#include <dune/common/iteratorrange.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/scvandscvfiterators.hh>

namespace Dumux
{

// forward declaration
template<class TypeTag, bool EnableFVGridGeometryCache>
class StaggeredFVGridGeometry;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for the finite volume geometry vector for staggered models
 *        This locally builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class TypeTag, bool EnableFVGridGeometryCache>
class StaggeredFVElementGeometry
{};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for the finite volume geometry vector for staggered models
 *        This locally builds up the sub control volumes and sub control volume faces
 *        for each element. Specialization in case the FVElementGeometries are stored globally.
          In this case we just forward internally to the global object.
 */
template<class TypeTag>
class StaggeredFVElementGeometry<TypeTag, true>
{
    using ThisType = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

    using ScvIterator = Dumux::ScvIterator<SubControlVolume, std::vector<IndexType>, ThisType>;
    using ScvfIterator = Dumux::ScvfIterator<SubControlVolumeFace, std::vector<IndexType>, ThisType>;

public:
    //! Constructor
    StaggeredFVElementGeometry(const FVGridGeometry& fvGridGeometry)
    : fvGridGeometryPtr_(&fvGridGeometry) {}

    //! Get an elment sub control volume with a global scv index
    //! We separate element and neighbor scvs to speed up mapping
    const SubControlVolume& scv(IndexType scvIdx) const
    {
        return fvGridGeometry().scv(scvIdx);
    }

    //! Get an element sub control volume face with a global scvf index
    //! We separate element and neighbor scvfs to speed up mapping
    const SubControlVolumeFace& scvf(IndexType scvfIdx) const
    {
        return fvGridGeometry().scvf(scvfIdx);
    }

    const SubControlVolumeFace& scvf(IndexType eIdx ,IndexType localScvfIdx) const
    {
        return fvGridGeometry().scvf(eIdx, localScvfIdx);
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element (not including neighbor scvs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<ScvIterator>
    scvs(const StaggeredFVElementGeometry& fvGeometry)
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
    scvfs(const StaggeredFVElementGeometry& fvGeometry)
    {
        const auto& g = fvGeometry.fvGridGeometry();
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
        return fvGridGeometry().scvfIndicesOfScv(scvIndices_[0]).size();
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
        scvIndices_ = std::vector<IndexType>({fvGridGeometry().elementMapper().index(*elementPtr_)});
    }

    //! The global finite volume geometry we are a restriction of
    const FVGridGeometry& fvGridGeometry() const
    { return *fvGridGeometryPtr_; }

private:

    const Element* elementPtr_;
    std::vector<IndexType> scvIndices_;
    const FVGridGeometry* fvGridGeometryPtr_;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class for the finite volume geometry vector for staggered models
 *        This locally builds up the sub control volumes and sub control volume faces
 *        for each element. Specialization in case the FVElementGeometries are not stored globally.
 */
template<class TypeTag>
class StaggeredFVElementGeometry<TypeTag, false>
{
    using ThisType = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

    using ScvIterator = Dumux::ScvIterator<SubControlVolume, std::vector<IndexType>, ThisType>;
    using ScvfIterator = Dumux::ScvfIterator<SubControlVolumeFace, std::vector<IndexType>, ThisType>;

public:
    //! Constructor
    StaggeredFVElementGeometry(const FVGridGeometry& fvGridGeometry)
    : fvGridGeometryPtr_(&fvGridGeometry) {}

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
        for (const auto& intersection : intersections(fvGridGeometry().gridView(), element))
        {
            neighborScvs_.reserve(element.subEntities(1));
            neighborScvfIndices_.reserve(element.subEntities(1));
            if (intersection.neighbor())
                makeNeighborGeometries(intersection.outside());
        }
    }

    //! Binding of an element preparing the geometries only inside the element
    void bindElement(const Element& element)
    {
        clear();
        elementPtr_ = &element;
        makeElementGeometries(element);
    }

    //! The global finite volume geometry we are a restriction of
    const FVGridGeometry& fvGridGeometry() const
    { return *fvGridGeometryPtr_; }

private:

    //! create scvs and scvfs of the bound element
    void makeElementGeometries(const Element& element)
    {
        auto eIdx = fvGridGeometry().problem_().elementMapper().index(element);
        scvs_.emplace_back(element.geometry(), eIdx);
        scvIndices_.emplace_back(eIdx);

        const auto& scvFaceIndices = fvGridGeometry().scvfIndicesOfScv(eIdx);
        const auto& neighborVolVarIndices = fvGridGeometry().neighborVolVarIndices(eIdx);

        int scvfCounter = 0;
        for (const auto& intersection : intersections(fvGridGeometry().gridView(), element))
        {
            if (intersection.neighbor() || intersection.boundary())
            {
                scvfs_.emplace_back(intersection,
                                    intersection.geometry(),
                                    scvFaceIndices[scvfCounter],
                                    std::vector<IndexType>({eIdx, neighborVolVarIndices[scvfCounter]})
                                    );

                scvfIndices_.emplace_back(scvFaceIndices[scvfCounter]);
                scvfCounter++;
            }
        }
    }

    //! create the necessary scvs and scvfs of the neighbor elements to the bound elements
    void makeNeighborGeometries(const Element& element)
    {
        // create the neighbor scv
        auto eIdx = fvGridGeometry().problem_().elementMapper().index(element);
        neighborScvs_.emplace_back(element.geometry(), eIdx);
        neighborScvIndices_.push_back(eIdx);

        const auto& scvFaceIndices = fvGridGeometry().scvfIndicesOfScv(eIdx);
        const auto& neighborVolVarIndices = fvGridGeometry().neighborVolVarIndices(eIdx);

        int scvfCounter = 0;
        for (const auto& intersection : intersections(fvGridGeometry().gridView(), element))
        {
            if (intersection.neighbor())
            {
                // only create subcontrol faces where the outside element is the bound element
                if (intersection.outside() == *elementPtr_)
                {
                    neighborScvfs_.emplace_back(intersection,
                                                intersection.geometry(),
                                                scvFaceIndices[scvfCounter],
                                                std::vector<IndexType>({eIdx, neighborVolVarIndices[scvfCounter]})
                                                );

                    neighborScvfIndices_.push_back(scvFaceIndices[scvfCounter]);
                }
                // always increase local counter
                scvfCounter++;
            }
            else if (intersection.boundary())
            {
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

    const FVGridGeometry* fvGridGeometryPtr_;

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
