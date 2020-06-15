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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredFVElementGeometry
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_FV_ELEMENT_GEOMETRY_HH

#include <type_traits>

#include <dune/common/reservedvector.hh>

#include <dumux/common/indextraits.hh>
#include <dune/common/iteratorrange.hh>
#include <dumux/discretization/scvandscvfiterators.hh>
#include <bitset>

namespace Dumux {

template<class GG, bool cachingEnabled>
class FaceCenteredStaggeredFVElementGeometry;

template<class GG>
class FaceCenteredStaggeredFVElementGeometry<GG, /*cachingEnabled*/true>
{
    using ThisType = FaceCenteredStaggeredFVElementGeometry<GG, /*cachingEnabled*/true>;
    using GridView = typename GG::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

public:
    //! export type of subcontrol volume face
    using SubControlVolume = typename GG::SubControlVolume;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridGeometry = GG;

        //! the maximum number of scvs per element
    static constexpr std::size_t maxNumElementScvs = 2*GridView::dimension;
    // //! the maximum number of scvfs per element (use cubes for maximum)
    // static constexpr std::size_t maxNumElementScvfs = 2*GridView::dimension;

    FaceCenteredStaggeredFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometry_(&gridGeometry)
    {}

    //! Get a sub control volume  with a global scv index
    const SubControlVolume& scv(GridIndexType scvIdx) const
    { return gridGeometry().scv(scvIdx); }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const
    { return gridGeometry().scvf(scvfIdx); }

     const SubControlVolumeFace& scvfWithCommonEntity(const SubControlVolumeFace& scvf) const
     {
         assert(scvf.isLateral());
         return gridGeometry().scvf(scvf.scvfIdxWithCommonEntity());
     }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element (not including neighbor scvs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline auto
    scvs(const FaceCenteredStaggeredFVElementGeometry& fvGeometry)
    {
        using IndexContainerType = std::decay_t<decltype(fvGeometry.scvIndices_())>;
        using ScvIterator = Dumux::ScvIterator<SubControlVolume, IndexContainerType, ThisType>;
        return Dune::IteratorRange<ScvIterator>(ScvIterator(fvGeometry.scvIndices_().begin(), fvGeometry),
                                                ScvIterator(fvGeometry.scvIndices_().end(), fvGeometry));
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline auto
    scvfs(const FaceCenteredStaggeredFVElementGeometry& fvGeometry)
    {
        using IndexContainerType = std::decay_t<decltype(fvGeometry.scvfIndices_())>;
        using ScvfIterator = Dumux::ScvfIterator<SubControlVolumeFace, IndexContainerType, ThisType>;
        return Dune::IteratorRange<ScvfIterator>(ScvfIterator(fvGeometry.scvfIndices_().begin(), fvGeometry),
                                                 ScvfIterator(fvGeometry.scvfIndices_().end(), fvGeometry));
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScv() const
    {
        return scvIndices_().size();
    }

    //! number of sub control volumes in this fv element geometry
    std::size_t numScvf() const
    {
        return scvfIndices_().size();
    }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return gridGeometry().hasBoundaryScvf(eIdx_); }

    //! Binding of an element, called by the local jacobian to prepare element assembly
    void bind(const Element& element)
    {
        this->bindElement(element);
    }

    //! Bind only element-local
    void bindElement(const Element& element)
    {
        elementPtr_ = &element;
        eIdx_ = gridGeometry().elementMapper().index(element);
    }

    //! The grid geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    {
        assert(gridGeometry_);
        return *gridGeometry_;
    }

    std::size_t elementIndex() const
    { return eIdx_; }

private:

    const auto& scvIndices_() const
    {
        return gridGeometry().scvIndicesOfElement(eIdx_);
    }

    const auto& scvfIndices_() const
    {
        return gridGeometry().scvfIndicesOfElement(eIdx_);
    }

    const Element* elementPtr_; //TODO maybe remove
    GridIndexType eIdx_;
    const GridGeometry* gridGeometry_;
};




template<class GG>
class FaceCenteredStaggeredFVElementGeometry<GG, /*cachingEnabled*/false>
{
    using GridView = typename GG::GridView;

    static constexpr std::size_t maxNumScvfs = 16; // TODO 3D

    using Scalar = double; // TODO

    //TODO include assert that checks for quad geometry
    static constexpr auto codimIntersection =  1;
    static constexpr auto dim = GridView::Grid::dimension;

    static constexpr auto numElementFaces = dim * 2;
    static constexpr auto numLateralFacesPerElementFace = 2 * (dim - 1);
    static constexpr auto numLateralFaces = numElementFaces*numLateralFacesPerElementFace;
    static constexpr auto numFacesWithoutRearBoundaryFaces = numLateralFaces + numElementFaces;


    static_assert(numLateralFaces == 8); //TODO remove
    static_assert(numFacesWithoutRearBoundaryFaces == 12); //TODO remove

    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;

public:
    //! export type of subcontrol volume face
    using SubControlVolume = typename GG::SubControlVolume;
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridGeometry = GG;

    FaceCenteredStaggeredFVElementGeometry(const GridGeometry& faceGridGeometry)
    : gridGeometry_(&faceGridGeometry)
    {}

    //! Get a sub control volume face with a local scv index
    const SubControlVolumeFace& scvf(LocalIndexType scvfIdx) const // TODO global scvf idx?
    { return scvfs_[scvfIdx]; }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline auto
    scvfs(const FaceCenteredStaggeredFVElementGeometry& fvGeometry)
    {
        // using Iter = typename decltype(fvGeometry.scvfs_)::const_iterator;
        // return Dune::IteratorRange<Iter>(fvGeometry.scvfs_.begin(), fvGeometry.scvfs_.end());
    }

    //! Get a half sub control volume with a global scv index
    const SubControlVolume& scv(const std::size_t scvIdx) const
    { return scvs_[findLocalIndex_(scvIdx, globalToLocalScvIdx_)]; }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element (not including neighbor scvs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline auto
    scvs(const FaceCenteredStaggeredFVElementGeometry& g)
    {
        using IteratorType = typename std::array<SubControlVolume, 1>::const_iterator;
        return Dune::IteratorRange<IteratorType>(g.scvs_.begin(), g.scvs_.end());
    }

    //! Binding of an element preparing the geometries of the whole stencil
    //! called by the local jacobian to prepare element assembly
    void bind(const Element& element)
    {
        // TODO!

    }

    //! The grid geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    {
        assert(gridGeometry_);
        return *gridGeometry_;
    }

private:

    template<class Entry, class Container>
    const LocalIndexType findLocalIndex_(const Entry& entry,
                                         const Container& container) const
    {
        auto it = std::find(container.begin(), container.end(), entry);
        assert(it != container.end() && "Could not find the entry! Make sure to properly bind this class!");
        return std::distance(container.begin(), it);
    }

    Dune::ReservedVector<SubControlVolumeFace, maxNumScvfs> scvfs_;
    std::array<SubControlVolume, numElementFaces> scvs_;
    std::array<std::size_t, numElementFaces> globalToLocalScvIdx_;

    const GridGeometry* gridGeometry_;
};

} // end namespace

#endif
