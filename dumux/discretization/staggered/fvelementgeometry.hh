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
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FV_ELEMENT_GEOMETRY_HH

#include <dune/common/reservedvector.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/cellcentered/tpfa/fvelementgeometry.hh>
#include <bitset>

namespace Dumux {

template<class GG, bool cachingEnabled>
class FaceStaggeredFVElementGeometry
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
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    using StaggeredSubControlVolumeFace = typename GG::StaggeredSubControlVolumeFace;
    using StaggeredHalfSubControlVolume = typename GG::StaggeredHalfSubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridGeometry = GG;

    FaceStaggeredFVElementGeometry(const GridGeometry& faceGridGeometry)
    : gridGeometry_(&faceGridGeometry)
    {}

    //! Get a sub control volume face with a local scv index
    const StaggeredSubControlVolumeFace& scvf(LocalIndexType scvfIdx) const // TODO global scvf idx?
    { return scvfs_[scvfIdx]; }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline auto
    scvfs(const FaceStaggeredFVElementGeometry& fvGeometry)
    {
        using Iter = typename decltype(fvGeometry.scvfs_)::const_iterator;
        return Dune::IteratorRange<Iter>(fvGeometry.scvfs_.begin(), fvGeometry.scvfs_.end());
    }

    //! Get a half sub control volume with a global scv index
    const StaggeredHalfSubControlVolume& scv(const std::size_t scvIdx) const
    { return scvs_[findLocalIndex_(scvIdx, globalToLocalScvIdx_)]; }

    //! Binding of an element preparing the geometries of the whole stencil
    //! called by the local jacobian to prepare element assembly
    void bind(const Element& element)
    {
        scvfs_.clear();
        scvs_ = {};
        normalDistancesOutsideElement_ = {};
        normalDistancesInsideElement_ = {};
        tangentialDistancesOutsideElement_ = {};
        elementLocalLaterScvfIdxToElementLocalScvfIdx_ = {};

        std::size_t scvfIdx = 0;
        std::size_t localLateralScvfIdx = 0;

        std::bitset<dim> inAxisDistanceHandled;

        typename GridGeometry::GeometryHelper geometryHelper(element, gridGeometry().gridView());

        for (const auto& intersection : intersections(gridGeometry_->gridView(), element))
        {
            geometryHelper.updateLocalFace(gridGeometry().actualGridGeometry().intersectionMapper(), intersection);
            const auto pairData = geometryHelper.pairData();
            const auto axisData = geometryHelper.axisData();

            const auto eIdx = gridGeometry().elementMapper().index(element);
            const auto& scvFaceIndices = gridGeometry().scvfIndicesOfScv(eIdx); // TODO hack, rename

            const auto localScvIdx = intersection.indexInInside();
            const auto insideScvIdx = axisData.selfDof;
            const auto directionIdx = geometryHelper.directionIndex();
            const auto intersectionGeometry = intersection.geometry();
            const auto elementCenter = element.geometry().center();

            globalToLocalScvIdx_[localScvIdx] = insideScvIdx;
            const auto scvCenter = 0.5*(intersectionGeometry.center() + elementCenter);
            scvs_[localScvIdx] = StaggeredHalfSubControlVolume(scvCenter,
                                                               intersectionGeometry.center(),
                                                               insideScvIdx,
                                                               scvFaceIndices[localScvIdx],
                                                               directionIdx,
                                                               sign(intersection.centerUnitOuterNormal()[directionIdx]),
                                                               !intersection.neighbor());


            // TODO write helper
            // frontal face
            scvfs_.push_back(StaggeredSubControlVolumeFace(elementCenter,
                                                           elementCenter,
                                                           std::array{insideScvIdx, axisData.oppositeDof},
                                                           intersectionGeometry.volume(),
                                                           directionIdx,
                                                           sign(intersection.centerUnitOuterNormal()[directionIdx]),
                                                           scvfIdx++,
                                                           StaggeredSubControlVolumeFace::FaceType::frontal,
                                                           false)); // TODO why does reservedvector not have emplace back?

            // rear frontal face on boundary
            if (!intersection.neighbor())
            {
                const auto boundaryCenter = intersectionGeometry.center();
                scvfs_.push_back(StaggeredSubControlVolumeFace(boundaryCenter,
                                                               boundaryCenter,
                                                               std::array{insideScvIdx, insideScvIdx}, // TODO outside boundary
                                                               intersectionGeometry.volume(),
                                                               directionIdx,
                                                               sign(intersection.centerUnitOuterNormal()[directionIdx]),
                                                               scvfIdx++,
                                                               StaggeredSubControlVolumeFace::FaceType::frontal,
                                                               true));
           }

            if (!inAxisDistanceHandled[directionIdx])
            {
                inAxisDistanceHandled.set(directionIdx);

                // Set the Self to Opposite Distance
                const auto oppIdx = localOppositeIdx_(localScvIdx);
                normalDistancesInsideElement_[directionIdx] = (getFacet_(localScvIdx, element).geometry().center() - getFacet_(oppIdx, element).geometry().center()).two_norm();
            }

            // lateral faces
            for (int i = 0; i < numLateralFacesPerElementFace; ++i)
            {
                elementLocalLaterScvfIdxToElementLocalScvfIdx_[localLateralScvfIdx] = scvfIdx;
                localLateralFaceIndexPerElementFace_[localLateralScvfIdx] = std::pair(scvfIdx, i);
                const auto& data = pairData[i];
                const auto& lateralFacetGeometry = getFacet_(data.localLateralFaceIdx, element).geometry();
                const auto& unitOuterNormal = data.unitOuterNormal;
                const auto dirIdx = directionIndex_(unitOuterNormal);
                const auto lateralFaceCenter = 0.5*(lateralFacetGeometry.center() + data.lateralStaggeredFaceCenter);
                scvfs_.push_back(StaggeredSubControlVolumeFace(lateralFaceCenter,
                                                               data.lateralStaggeredFaceCenter,
                                                               std::array{insideScvIdx, data.parallelDofs[0]}, // TODO higher order
                                                               lateralFacetGeometry.volume()*0.5,
                                                               dirIdx,
                                                               sign(unitOuterNormal[dirIdx]),
                                                               scvfIdx++,
                                                               StaggeredSubControlVolumeFace::FaceType::lateral,
                                                               data.lateralFaceOnBoundary));

                normalDistancesOutsideElement_[localLateralScvfIdx] = (faceLength_(i, intersectionGeometry) + data.parallelCellWidths[0]) * 0.5; // higher order

                for (int k = 0; k < dim-1; ++k) // TODO k == 0 for staggered, also in 3D, right?
                    tangentialDistancesOutsideElement_[localLateralScvfIdx][k] = data.lateralDistance;

                ++localLateralScvfIdx;
            }
        }

        assert(localLateralScvfIdx == 8);
    }

    //! The grid geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    {
        assert(gridGeometry_);
        return *gridGeometry_;
    }

    Scalar normalDistanceForGradient(const StaggeredSubControlVolumeFace& scvf) const
    {
        if (scvf.isFrontal())
        {
            assert(!scvf.boundary());
            return normalDistancesInsideElement_[scvf.directionIndex()];
        }
        else
            return normalDistancesOutsideElement_[findLocalIndex_(scvf.localScvfIdx(), elementLocalLaterScvfIdxToElementLocalScvfIdx_)];
    }

    Scalar tangentialDistanceForGradient(const StaggeredSubControlVolumeFace& scvf, LocalIndexType dirIdx = 0) const
    {
        assert(scvf.isLateral());
        return tangentialDistancesOutsideElement_[findLocalIndex_(scvf.localScvfIdx(), elementLocalLaterScvfIdxToElementLocalScvfIdx_)][dirIdx]; // TODO check this
    }

    LocalIndexType localLateralFaceIndex(const StaggeredSubControlVolumeFace& scvf) const
    {
        assert(scvf.isLateral());
        const auto it = std::find_if(localLateralFaceIndexPerElementFace_.begin(), localLateralFaceIndexPerElementFace_.end(),
                                     [&scvf](const auto& x) { return x.first == scvf.localScvfIdx(); });
        assert(it != localLateralFaceIndexPerElementFace_.end() && "Could not find the entry! Make sure to properly bind this class!");
        return it->second;
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


    //! Returns the length of the face in a certain direction (adaptation of area() for 3d)
    template<class Geo>
    Scalar faceLength_([[maybe_unused]] const LocalIndexType localSubFaceIdx, const Geo& geo) const
    {
        if constexpr (dim == 3)
        {
            if (localSubFaceIdx < 2)
                return (geo.corner(1) - geo.corner(0)).two_norm();
            else
                return (geo.corner(2) - geo.corner(0)).two_norm();
        }
        else
            return (geo.corner(1) - geo.corner(0)).two_norm();
    }

    // TODO use from helper
    template<class Vector>
    unsigned int directionIndex_(Vector&& vector) const
    {
        const auto eps = 1e-8;
        return std::find_if(vector.begin(), vector.end(), [eps](const auto& x) { return std::abs(x) > eps; } ) - vector.begin();
    }

    /*!
     * \brief Returns the local opposing intersection index
     *
     * \param idx The local index of the intersection itself
     */
    int localOppositeIdx_(const int idx) const
    {
        return (idx % 2) ? (idx - 1) : (idx + 1);
    }

    auto getFacet_(const int localFacetIdx, const Element& element) const
    {
        return element.template subEntity <1> (localFacetIdx);
    };


    Dune::ReservedVector<StaggeredSubControlVolumeFace, maxNumScvfs> scvfs_;
    std::array<StaggeredHalfSubControlVolume, numElementFaces> scvs_;

    std::array<Scalar, numLateralFaces> normalDistancesOutsideElement_;
    std::array<Scalar, dim> normalDistancesInsideElement_;

    std::array<std::array<Scalar, dim-1>, numLateralFaces> tangentialDistancesOutsideElement_;

    std::array<std::pair<LocalIndexType, LocalIndexType>, numLateralFaces> localLateralFaceIndexPerElementFace_;
    std::array<LocalIndexType, numLateralFaces> elementLocalLaterScvfIdxToElementLocalScvfIdx_;

    std::array<std::size_t, numElementFaces> globalToLocalScvIdx_;


    const GridGeometry* gridGeometry_;
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for staggered models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element in the local scope we are restricting to, e.g. stencil or element.
 * \tparam GG the finite volume grid geometry type
 * \tparam enableGridGeometryCache if the grid geometry is cached or not
 * \note This class is specialized for versions with and without caching the fv geometries on the grid view
 */
template<class GG, bool enableGridGeometryCache>
class StaggeredFVElementGeometry;

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for the finite volume geometry vector for staggered models
 *        This locally builds up the sub control volumes and sub control volume faces
 *        for each element.
 *        Specialization for grid caching enabled
 * \tparam GG the finite volume grid geometry type
 */
template<class GG>
class StaggeredFVElementGeometry<GG, true> : public CCTpfaFVElementGeometry<GG, true>
{
    using ParentType = CCTpfaFVElementGeometry<GG, true>;
    using GridView = typename GG::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
public:
    //! export type of the element
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;

    using ParentType::ParentType;

    //! Constructor getting a auxiliary cell center of face specific FvGridGeometry type.
    //! Needed for the multi-domain framework.


    StaggeredFVElementGeometry(const GG& gridGeometry)
    : ParentType(gridGeometry) {}

    StaggeredFVElementGeometry(const typename GG::CellCenterFVGridGeometryType& gridGeometry)
    : ParentType(gridGeometry.actualGridGeometry()) {}

    StaggeredFVElementGeometry(const typename GG::FaceFVGridGeometryType& gridGeometry)
    : ParentType(gridGeometry.actualGridGeometry()) {}


    // StaggeredFVElementGeometry(const CellCenterOrFaceFVGridGeometry& gridGeometry)
    // : ParentType(gridGeometry.actualGridGeometry()) {}

    //! Get a sub control volume face with an element index and a local scvf index
    using ParentType::scvf;
    const SubControlVolumeFace& scvf(GridIndexType eIdx, LocalIndexType localScvfIdx) const
    {
        return this->gridGeometry().scvf(eIdx, localScvfIdx);
    }
};

/*!
 * \ingroup StaggeredDiscretization
 * \brief Base class for the finite volume geometry vector for staggered models
 *        This locally builds up the sub control volumes and sub control volume faces
 *        for each element.
 *        Specialization for grid caching enabled
 * \tparam GG the finite volume grid geometry type
 */
template<class GG>
class StaggeredFVElementGeometry<GG, false>
{
    using ThisType = StaggeredFVElementGeometry<GG, false>;
    using GridView = typename GG::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
public:
    //! export type of the element
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;

    //! Constructor getting a auxiliary cell center of face specific FvGridGeometry type.
    //! Needed for the multi-domain framework.
    template<class CellCenterOrFaceFVGridGeometry>
    StaggeredFVElementGeometry(const CellCenterOrFaceFVGridGeometry& gridGeometry)
    : gridGeometryPtr_(&gridGeometry.actualGridGeometry()) {}

    //! Constructor
    StaggeredFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometryPtr_(&gridGeometry) {}

    //! Get a sub control volume face with an element index and a local scvf index
    const SubControlVolumeFace& scvf(GridIndexType eIdx, LocalIndexType localScvfIdx) const
    {
        return scvf(this->gridGeometry().localToGlobalScvfIndex(eIdx, localScvfIdx));
    }

    //! Get an elment sub control volume with a global scv index
    //! We separate element and neighbor scvs to speed up mapping
    const SubControlVolume& scv(GridIndexType scvIdx) const
    {
        if (scvIdx == scvIndices_[0])
            return scvs_[0];
        else
            return neighborScvs_[findLocalIndex_(scvIdx, neighborScvIndices_)];
    }

    //! Get an element sub control volume face with a global scvf index
    //! We separate element and neighbor scvfs to speed up mapping
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const
    {
        auto it = std::find(scvfIndices_.begin(), scvfIndices_.end(), scvfIdx);
        if (it != scvfIndices_.end())
            return scvfs_[std::distance(scvfIndices_.begin(), it)];
        else
            return neighborScvfs_[findLocalIndex_(scvfIdx, neighborScvfIndices_)];
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element (not including neighbor scvs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::array<SubControlVolume, 1>::const_iterator>
    scvs(const ThisType& g)
    {
        using IteratorType = typename std::array<SubControlVolume, 1>::const_iterator;
        return Dune::IteratorRange<IteratorType>(g.scvs_.begin(), g.scvs_.end());
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element (not including neighbor scvfs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolumeFace>::const_iterator>
    scvfs(const ThisType& g)
    {
        using IteratorType = typename std::vector<SubControlVolumeFace>::const_iterator;
        return Dune::IteratorRange<IteratorType>(g.scvfs_.begin(), g.scvfs_.end());
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

        neighborScvs_.reserve(element.subEntities(1));
        neighborScvfIndices_.reserve(element.subEntities(1));
        neighborScvfs_.reserve(element.subEntities(1));

        std::vector<GridIndexType> handledNeighbors;
        handledNeighbors.reserve(element.subEntities(1));
        for (const auto& intersection : intersections(gridGeometry().gridView(), element))
        {
            if (intersection.neighbor())
            {
                const auto outside = intersection.outside();
                const auto outsideIdx = gridGeometry().elementMapper().index(outside);

                // make outside geometries only if not done yet (could happen on non-conforming grids)
                if ( std::find(handledNeighbors.begin(), handledNeighbors.end(), outsideIdx) == handledNeighbors.end() )
                {
                    makeNeighborGeometries_(outside, outsideIdx);
                    handledNeighbors.push_back(outsideIdx);
                }
            }
        }
    }

    //! Binding of an element preparing the geometries only inside the element
    void bindElement(const Element& element)
    {
        clear_();
        elementPtr_ = &element;
        scvfs_.reserve(element.subEntities(1));
        scvfIndices_.reserve(element.subEntities(1));
        makeElementGeometries_(element);
    }

    //! The grid finite volume geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return *gridGeometryPtr_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return hasBoundaryScvf_; }

private:

    //! create scvs and scvfs of the bound element
    void makeElementGeometries_(const Element& element)
    {
        const auto eIdx = gridGeometry().elementMapper().index(element);
        scvs_[0] = SubControlVolume(element.geometry(), eIdx);
        scvIndices_[0] = eIdx;

        const auto& scvFaceIndices = gridGeometry().scvfIndicesOfScv(eIdx);
        const auto& neighborVolVarIndices = gridGeometry().neighborVolVarIndices(eIdx);

        typename GridGeometry::GeometryHelper geometryHelper(element, gridGeometry().gridView());

        int scvfCounter = 0;
        for (const auto& intersection : intersections(gridGeometry().gridView(), element))
        {
            const auto& scvfNeighborVolVarIndex = neighborVolVarIndices[scvfCounter];

            if (intersection.neighbor() || intersection.boundary())
            {
                geometryHelper.updateLocalFace(gridGeometry().intersectionMapper(), intersection);
                std::vector<GridIndexType> scvIndices{eIdx, scvfNeighborVolVarIndex};
                scvfs_.emplace_back(intersection,
                                    intersection.geometry(),
                                    scvFaceIndices[scvfCounter],
                                    scvIndices,
                                    geometryHelper);
                scvfIndices_.emplace_back(scvFaceIndices[scvfCounter]);
                scvfCounter++;

                if (intersection.boundary())
                    hasBoundaryScvf_ = true;
            }
        }
    }

    //! create the necessary scvs and scvfs of the neighbor elements to the bound elements
    void makeNeighborGeometries_(const Element& element, const GridIndexType eIdx)
    {
        // using ScvfGridIndexStorage = typename SubControlVolumeFace::Traits::GridIndexStorage;

        // create the neighbor scv
        neighborScvs_.emplace_back(element.geometry(), eIdx);
        neighborScvIndices_.push_back(eIdx);

        const auto& scvFaceIndices = gridGeometry().scvfIndicesOfScv(eIdx);
        const auto& neighborVolVarIndices = gridGeometry().neighborVolVarIndices(eIdx);

        typename GridGeometry::GeometryHelper geometryHelper(element, gridGeometry().gridView());

        int scvfCounter = 0;
        for (const auto& intersection : intersections(gridGeometry().gridView(), element))
        {
            const auto& scvfNeighborVolVarIndex = neighborVolVarIndices[scvfCounter];
            geometryHelper.updateLocalFace(gridGeometry().intersectionMapper(), intersection);

            if (intersection.neighbor())
            {
                // only create subcontrol faces where the outside element is the bound element
                if (intersection.outside() == *elementPtr_)
                {
                    std::vector<GridIndexType> scvIndices{eIdx, scvfNeighborVolVarIndex};
                    neighborScvfs_.emplace_back(intersection,
                                                intersection.geometry(),
                                                scvFaceIndices[scvfCounter],
                                                scvIndices,
                                                geometryHelper);

                    neighborScvfIndices_.push_back(scvFaceIndices[scvfCounter]);
                }
                scvfCounter++;
            }
            else if (intersection.boundary())
                scvfCounter++;
        }
    }

    const LocalIndexType findLocalIndex_(const GridIndexType idx,
                                         const std::vector<GridIndexType>& indices) const
    {
        auto it = std::find(indices.begin(), indices.end(), idx);
        assert(it != indices.end() && "Could not find the scv/scvf! Make sure to properly bind this class!");
        return std::distance(indices.begin(), it);
    }

    //! Clear all local data
    void clear_()
    {
        scvfIndices_.clear();
        scvfs_.clear();

        neighborScvIndices_.clear();
        neighborScvfIndices_.clear();
        neighborScvs_.clear();
        neighborScvfs_.clear();

        hasBoundaryScvf_ = false;
    }

    const Element* elementPtr_; //!< the element to which this fvgeometry is bound
    const GridGeometry* gridGeometryPtr_;  //!< the grid fvgeometry

    // local storage after binding an element
    std::array<GridIndexType, 1> scvIndices_;
    std::array<SubControlVolume, 1> scvs_;

    std::vector<GridIndexType> scvfIndices_;
    std::vector<SubControlVolumeFace> scvfs_;

    std::vector<GridIndexType> neighborScvIndices_;
    std::vector<SubControlVolume> neighborScvs_;

    std::vector<GridIndexType> neighborScvfIndices_;
    std::vector<SubControlVolumeFace> neighborScvfs_;

    bool hasBoundaryScvf_ = false;
};


} // end namespace

#endif
