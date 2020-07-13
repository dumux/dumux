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
 * \ingroup PoreNetworkDiscretization
 * \brief Base class for the local geometry for porenetworks
 */
#ifndef DUMUX_DISCRETIZATION_PNM_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_PNM_FV_ELEMENT_GEOMETRY_HH

#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/scvandscvfiterators.hh>

namespace Dumux {

/*!
 * \ingroup PoreNetworkDiscretization
 * \brief Base class for the local geometry for porenetworks
 * \tparam GG the finite volume grid geometry type
 * \tparam enableFVGridGeometryCache if the grid geometry is cached or not
 */
template<class GG, bool enableFVGridGeometryCache>
class PNMFVElementGeometry;

//! specialization in case the FVElementGeometries are stored
template<class GG>
class PNMFVElementGeometry<GG, true>
{
    using GridView = typename GG::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Element = typename GridView::template Codim<0>::Entity;
    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
public:
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;
    using FVGridGeometry [[deprecated("Use GridGeometry")]] = GG;
    //! the maximum number of scvs per element
    static constexpr std::size_t maxNumElementScvs = 2;

    //! Constructor
    PNMFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometryPtr_(&gridGeometry) {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(LocalIndexType scvIdx) const
    {
        return gridGeometry().scvs(eIdx_)[scvIdx];
    }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(LocalIndexType scvfIdx) const
    {
        return gridGeometry().scvfs(eIdx_)[scvfIdx];
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline auto scvs(const PNMFVElementGeometry& fvGeometry)
    {
        const auto& g = fvGeometry.gridGeometry();
        using Iter = typename std::decay_t<decltype(g.scvs(fvGeometry.eIdx_))>::const_iterator;
        return Dune::IteratorRange<Iter>(g.scvs(fvGeometry.eIdx_).begin(), g.scvs(fvGeometry.eIdx_).end());
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline auto scvfs(const PNMFVElementGeometry& fvGeometry)
    {
        const auto& g = fvGeometry.gridGeometry();
        using Iter = typename std::decay_t<decltype(g.scvfs(fvGeometry.eIdx_))>::const_iterator;
        return Dune::IteratorRange<Iter>(g.scvfs(fvGeometry.eIdx_).begin(), g.scvfs(fvGeometry.eIdx_).end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    {
        return gridGeometry().feCache().get(elementPtr_->geometry().type()).localBasis();
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return gridGeometry().scvs(eIdx_).size();
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return gridGeometry().scvfs(eIdx_).size();
    }

    //! this function is for compatibility reasons with cc methods
    //! The box stencil is always element-local so bind and bindElement
    //! are identical.
    void bind(const Element& element)
    {
        this->bindElement(element);
    }

    //! Binding of an element, has to be called before using the fvgeometries
    //! Prepares all the volume variables within the element
    //! For compatibility reasons with the FVGeometry cache being disabled
    void bindElement(const Element& element)
    {
        elementPtr_ = &element;
        eIdx_ = gridGeometry().elementMapper().index(element);
    }

    //! The global finite volume geometry we are a restriction of
    [[deprecated("Use gridGeometry")]]
    const GridGeometry& fvGridGeometry() const
    { return *gridGeometryPtr_; }

    //! The global finite volume geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return *gridGeometryPtr_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return gridGeometry().hasBoundaryScvf(eIdx_); }

private:
    const Element* elementPtr_;
    const GridGeometry* gridGeometryPtr_;

    GridIndexType eIdx_;
};

//! specialization in case the FVElementGeometries are not stored
template<class GG>
class PNMFVElementGeometry<GG, false>
{
    using GridView = typename GG::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Element = typename GridView::template Codim<0>::Entity;
    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
public:
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;
    using FVGridGeometry [[deprecated("Use GridGeometry")]] = GG;
        //! the maximum number of scvs per element
    static constexpr std::size_t maxNumElementScvs = 2;

    //! Constructor
    PNMFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometryPtr_(&gridGeometry) {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(LocalIndexType scvIdx) const
    {
        return scvs_[scvIdx];
    }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(LocalIndexType scvfIdx) const
    {
        return scvfs_[0];
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline auto scvs(const PNMFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::decay_t<decltype(fvGeometry.scvs_)>::const_iterator;
        return Dune::IteratorRange<Iter>(fvGeometry.scvs_.begin(), fvGeometry.scvs_.end());
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline auto scvfs(const PNMFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::decay_t<decltype(fvGeometry.scvfs_)>::const_iterator;
        return Dune::IteratorRange<Iter>(fvGeometry.scvfs_.begin(), fvGeometry.scvfs_.end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    {
        return gridGeometry().feCache().get(elementPtr_->geometry().type()).localBasis();
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return scvs_.size();
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return scvfs_.size();
    }

    //! this function is for compatibility reasons with cc methods
    //! The box stencil is always element-local so bind and bindElement
    //! are identical.
    void bind(const Element& element)
    {
        this->bindElement(element);
    }

    //! Binding of an element, has to be called before using the fvgeometries
    //! Prepares all the volume variables within the element
    //! For compatibility reasons with the FVGeometry cache being disabled
    void bindElement(const Element& element)
    {
        elementPtr_ = &element;
        eIdx_ = gridGeometry().elementMapper().index(element);
        makeElementGeometries(element);
    }

    //! The global finite volume geometry we are a restriction of
    [[deprecated("Use gridGeometry")]]
    const GridGeometry& fvGridGeometry() const
    { return *gridGeometryPtr_; }

    //! The global finite volume geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return *gridGeometryPtr_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return hasBoundaryScvf_; }

private:

    void makeElementGeometries(const Element& element)
    {
        hasBoundaryScvf_ = false;

        // get the element geometry
        auto elementGeometry = element.geometry();

        // construct the sub control volumes
        for (LocalIndexType scvLocalIdx = 0; scvLocalIdx < elementGeometry.corners(); ++scvLocalIdx)
        {
            // get asssociated dof index
            const auto dofIdxGlobal = gridGeometry().vertexMapper().subIndex(element, scvLocalIdx, dim);

            // get the fractional volume asssociated with the scv
            const auto volume = gridGeometry().poreVolume(dofIdxGlobal) / gridGeometry().coordinationNumber(dofIdxGlobal);

            // add scv to the local container
            scvs_[scvLocalIdx] = SubControlVolume(dofIdxGlobal,
                                                  scvLocalIdx,
                                                  eIdx_,
                                                  elementGeometry.corner(scvLocalIdx),
                                                  volume);

            if (gridGeometry().poreLabel(dofIdxGlobal) > 0)
                hasBoundaryScvf_ = true;
        }

        // construct the inner sub control volume face
        auto unitOuterNormal = elementGeometry.corner(1)-elementGeometry.corner(0);
        unitOuterNormal /= unitOuterNormal.two_norm();
        LocalIndexType scvfLocalIdx = 0;
        scvfs_[0] = SubControlVolumeFace(elementGeometry.center(),
                                         std::move(unitOuterNormal),
                                         scvfLocalIdx++,
                                         std::vector<LocalIndexType>({0, 1}));
    }

    //! The bound element
    const Element* elementPtr_;
    GridIndexType eIdx_;

    //! The global geometry this is a restriction of
    const GridGeometry* gridGeometryPtr_;

    //! vectors to store the geometries locally after binding an element
    std::array<SubControlVolume, 2> scvs_;
    std::array<SubControlVolumeFace, 1> scvfs_;

    bool hasBoundaryScvf_ = false;
};

} // end namespace Dumux

#endif
