// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredFVElementGeometry
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FV_ELEMENT_GEOMETRY_HH

#include <optional>
#include <bitset>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/cellcentered/tpfa/fvelementgeometry.hh>
#include <dumux/discretization/facecentered/staggered/normalaxis.hh>

namespace Dumux {

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
    template<class CellCenterOrFaceFVGridGeometry>
    StaggeredFVElementGeometry(const CellCenterOrFaceFVGridGeometry& gridGeometry)
    : ParentType(gridGeometry.actualGridGeometry()) {}

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
 *        Specialization for grid caching disabled
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

    //! Get an element sub control volume with a global scv index
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

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    StaggeredFVElementGeometry bind(const Element& element) &&
    {
        this->bind_(element);
        return std::move(*this);
    }

    //! bind this local view to a specific element (full stencil)
    void bind(const Element& element) &
    { this->bind_(element); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    StaggeredFVElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement_(element);
        return std::move(*this);
    }

    //! bind this local view to a specific element
    void bindElement(const Element& element) &
    { this->bindElement_(element); }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    //! The grid finite volume geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return *gridGeometryPtr_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return hasBoundaryScvf_; }

    //! Create the geometry of a given sub control volume
    typename SubControlVolume::Traits::Geometry geometry(const SubControlVolume& scv) const
    { return gridGeometryPtr_->element(scv.dofIndex()).geometry(); }

    //! Create the geometry of a given sub control volume face
    typename SubControlVolumeFace::Traits::Geometry geometry(const SubControlVolumeFace& scvf) const
    {
        const auto insideElementIndex = scvf.insideScvIdx();
        const auto elementGeometry = (insideElementIndex != scvIndices_[0]) ?
            element_->geometry() :
            gridGeometryPtr_->element(insideElementIndex).geometry();
        const auto center = elementGeometry.center();
        const auto normalAxis = Dumux::normalAxis(scvf.unitOuterNormal());

        auto lowerLeft = elementGeometry.corner(0);
        auto upperRight = elementGeometry.corner(elementGeometry.corners()-1);

        // shift corners to scvf plane and halve lateral faces
        lowerLeft[normalAxis] = center[normalAxis];
        upperRight[normalAxis] = center[normalAxis];

        auto inPlaneAxes = std::move(std::bitset<SubControlVolumeFace::Traits::dimWorld>{}.set());
        inPlaneAxes.set(normalAxis, false);

        return {lowerLeft, upperRight, inPlaneAxes};
    }

private:
    //! Binding of an element preparing the geometries only inside the element
    void bindElement_(const Element& element)
    {
        clear_();
        element_ = element;
        scvfs_.reserve(element.subEntities(1));
        scvfIndices_.reserve(element.subEntities(1));
        makeElementGeometries_();
    }

    //! Binding of an element preparing the geometries of the whole stencil
    //! called by the local jacobian to prepare element assembly
    void bind_(const Element& element)
    {
        bindElement_(element);

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

    //! create scvs and scvfs of the bound element
    void makeElementGeometries_()
    {
        const auto& element = *element_;
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
                if (intersection.outside() == *element_)
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

    std::optional<Element> element_; //!< the element to which this fvgeometry is bound
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
