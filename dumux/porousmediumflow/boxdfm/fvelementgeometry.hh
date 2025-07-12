// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxDFMModel
 * \brief Base class for the local finite volume geometry for the box discrete
 *        fracture model.
 *
 * This builds up the sub control volumes and sub control
 * volume faces for an element.
 */

#ifndef DUMUX_POROUSMEDIUMFLOW_BOXDFM_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_POROUSMEDIUMFLOW_BOXDFM_FV_ELEMENT_GEOMETRY_HH

#include <optional>
#include <utility>
#include <ranges>

#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>


#include <dumux/common/indextraits.hh>
#include <dumux/discretization/scvandscvfiterators.hh>
#include "geometryhelper.hh"

namespace Dumux {

/*!
 * \ingroup BoxDFMModel
 * \brief Base class for the finite volume geometry vector for box discrete fracture model.
 *
 * This builds up the sub control volumes and sub control volume faces for each element.
 *
 * \tparam GG the finite volume grid geometry type
 * \tparam enableGridGeometryCache if the grid geometry is cached or not
 */
template<class GG, bool enableGridGeometryCache>
class BoxDfmFVElementGeometry;

//! Specialization in case the FVElementGeometries are stored
template<class GG>
class BoxDfmFVElementGeometry<GG, true>
{
    using GridView = typename GG::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
    using GGCache = typename GG::Cache;
    using GeometryHelper = typename GGCache::GeometryHelper;
public:
    //! export type of the element
    using Element = typename GridView::template Codim<0>::Entity;
    //! Export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! Export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! Export type of finite volume grid geometry
    using GridGeometry = GG;

    //! The maximum number of scvs per element (2^dim for cubes)
    //! multiplied by 3 for the maximum number of fracture scvs per vertex
    static constexpr std::size_t maxNumElementScvs = (1<<dim)*3;

    //! Constructor
    [[deprecated("This Constructor is deprecated and will be removed after release 3.11. Always use localView(gridGeometry).")]]
    BoxDfmFVElementGeometry(const GridGeometry& gridGeometry)
    : BoxDfmFVElementGeometry(gridGeometry.cache_) {}

    //! Constructor
    BoxDfmFVElementGeometry(const GGCache& ggCache)
    : ggCache_(&ggCache) {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(LocalIndexType scvIdx) const
    { return ggCache_->scvs(eIdx_)[scvIdx]; }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(LocalIndexType scvfIdx) const
    { return ggCache_->scvfs(eIdx_)[scvfIdx]; }

   /*!
    * \brief Iterator range for sub control volumes.
    *
    * Iterates over all scvs of the bound element.
    * This is a free function found by means of ADL.
    * To iterate over all sub control volumes of this FVElementGeometry use
    * for (auto&& scv : scvs(fvGeometry)).
    */
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolume>::const_iterator>
    scvs(const BoxDfmFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolume>::const_iterator;
        const auto& s = fvGeometry.ggCache_->scvs(fvGeometry.eIdx_);
        return Dune::IteratorRange<Iter>(s.begin(), s.end());
    }

    //! range over sub control volumes related to a local dof.
    template<class LocalDof>
    friend inline std::ranges::range auto
    scvs(const BoxDfmFVElementGeometry& fvGeometry, const LocalDof& localDof)
    {
       return std::views::iota(0u, fvGeometry.numScv())
            | std::views::filter([&](std::size_t i) { return fvGeometry.scv(i).localDofIndex() == localDof.index(); })
            | std::views::transform([&](std::size_t i) -> const SubControlVolume& { return fvGeometry.scv(i); });
    }

    //! range over local dofs
    friend inline auto localDofs(const BoxDfmFVElementGeometry& fvGeometry)
    {
        return Dune::transformedRangeView(
            Dune::range(fvGeometry.numLocalDofs()),
            [&](const auto i) { return CVFE::LocalDof
            {
                static_cast<LocalIndexType>(fvGeometry.scv(i).localDofIndex()),
                static_cast<GridIndexType>(fvGeometry.scv(i).dofIndex()),
                static_cast<GridIndexType>(fvGeometry.scv(i).elementIndex())
            }; }
        );
    }

    //! range over control-volume local dofs
    friend inline auto cvLocalDofs(const BoxDfmFVElementGeometry& fvGeometry)
    {
        // Default it that all dofs are cv dofs
        return localDofs(fvGeometry);
    }

   /*!
    * \brief Iterator range for sub control volumes faces.
    *
    * Iterates over all scvfs of the bound element.
    * This is a free function found by means of ADL.
    * To iterate over all sub control volume faces of this FVElementGeometry use
    * for (auto&& scvf : scvfs(fvGeometry)).
    */
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolumeFace>::const_iterator>
    scvfs(const BoxDfmFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolumeFace>::const_iterator;
        const auto& s = fvGeometry.ggCache_->scvfs(fvGeometry.eIdx_);
        return Dune::IteratorRange<Iter>(s.begin(), s.end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    { return gridGeometry().feCache().get(element_->type()).localBasis(); }

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return ggCache_->scvs(eIdx_).size(); }

    //! The total number of element-local dofs
    std::size_t numLocalDofs() const
    { return element().geometry().corners(); }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    { return ggCache_->scvfs(eIdx_).size(); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    BoxDfmFVElementGeometry bind(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! This function is for compatibility reasons with cc methods
    //! The box stencil is always element-local so bind and bindElement are identical.
    void bind(const Element& element) &
    { this->bindElement(element); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    BoxDfmFVElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

   /*!
    * \brief Binding of an element, has to be called before using the fvgeometries
    *
    * Prepares all the volume variables within the element.
    * For compatibility reasons with the FVGeometry cache being disabled.
    */
    void bindElement(const Element& element) &
    {
        element_ = element;
        eIdx_ = gridGeometry().elementMapper().index(element);
    }

    //! The global finite volume geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return ggCache_->gridGeometry(); }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    //! Create the geometry of a given sub control volume
    typename SubControlVolume::Traits::Geometry geometry(const SubControlVolume& scv) const
    {
        if (scv.isOnFracture())
            DUNE_THROW(Dune::InvalidStateException, "The geometry object cannot be defined for fracture scvs "
                                                    "because the number of known corners is insufficient. "
                                                    "You can do this manually by extracting the corners from this scv "
                                                    "and extruding them by the corresponding aperture. ");

        const GeometryHelper geometryHelper(element().geometry());
        const auto corners = geometryHelper.getScvCorners(scv.index());
        using ScvGeometry = typename SubControlVolume::Traits::Geometry;
        return { Dune::GeometryTypes::cube(ScvGeometry::mydimension), corners };
    }

    //! Create the geometry of a given sub control volume face
    typename SubControlVolumeFace::Traits::Geometry geometry(const SubControlVolumeFace& scvf) const
    {
        if (scvf.isOnFracture())
            DUNE_THROW(Dune::InvalidStateException, "The geometry object cannot be defined for fracture scvs "
                                                    "because the number of known corners is insufficient. "
                                                    "You can do this manually by extracting the corners from this scv "
                                                    "and extruding them by the corresponding aperture. ");
        const GeometryHelper geometryHelper(element().geometry());
        const auto corners = geometryHelper.getScvfCorners(scvf.indexInElement());
        using ScvfGeometry = typename SubControlVolumeFace::Traits::Geometry;
        return { Dune::GeometryTypes::cube(ScvfGeometry::mydimension), corners };
    }

private:
    const GGCache* ggCache_;

    std::optional<Element> element_;
    GridIndexType eIdx_;
};

//! Specialization in case the FVElementGeometries are not stored
template<class GG>
class BoxDfmFVElementGeometry<GG, false>
{
    using GridView = typename GG::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;

    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
    using GGCache = typename GG::Cache;
    using GeometryHelper = typename GGCache::GeometryHelper;
public:
    //! export type of the element
    using Element = typename GridView::template Codim<0>::Entity;
    //! Export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! Export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! Export type of finite volume grid geometry
    using GridGeometry = GG;
    //! The maximum number of scvs per element (2^dim for cubes)
    //! multiplied by 3 for the maximum number of fracture scvs per vertex
    static constexpr std::size_t maxNumElementScvs = (1<<dim)*3;

    //! Constructor
    [[deprecated("This Constructor is deprecated and will be removed after release 3.11. Always use localView(gridGeometry).")]]
    BoxDfmFVElementGeometry(const GridGeometry& gridGeometry)
    : BoxDfmFVElementGeometry(gridGeometry.cache_) {}

    //! Constructor
    BoxDfmFVElementGeometry(const GGCache& ggCache)
    : ggCache_(&ggCache) {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(LocalIndexType scvIdx) const
    { return scvs_[scvIdx]; }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(LocalIndexType scvfIdx) const
    { return scvfs_[scvfIdx]; }

   /*!
    * \brief Iterator range for sub control volumes.
    *
    * Iterates over all scvs of the bound element.
    * This is a free function found by means of ADL.
    * To iterate over all sub control volumes of this FVElementGeometry use
    * for (auto&& scv : scvs(fvGeometry)).
    */
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolume>::const_iterator>
    scvs(const BoxDfmFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolume>::const_iterator;
        return Dune::IteratorRange<Iter>(fvGeometry.scvs_.begin(), fvGeometry.scvs_.end());
    }

    //! range over sub control volumes related to a local dof.
    template<class LocalDof>
    friend inline std::ranges::range auto
    scvs(const BoxDfmFVElementGeometry& fvGeometry, const LocalDof& localDof)
    {
       return std::views::iota(0u, fvGeometry.numScv())
            | std::views::filter([&](std::size_t i) { return fvGeometry.scv(i).localDofIndex() == localDof.index(); })
            | std::views::transform([&](std::size_t i) -> const SubControlVolume& { return fvGeometry.scv(i); });
    }

    //! range over local dofs
    friend inline auto localDofs(const BoxDfmFVElementGeometry& fvGeometry)
    {
        return Dune::transformedRangeView(
            Dune::range(fvGeometry.numLocalDofs()),
            [&](const auto i) { return CVFE::LocalDof
            {
                static_cast<LocalIndexType>(fvGeometry.scv(i).localDofIndex()),
                static_cast<GridIndexType>(fvGeometry.scv(i).dofIndex()),
                static_cast<GridIndexType>(fvGeometry.scv(i).elementIndex())
            }; }
        );
    }

    //! range over control-volume local dofs
    friend inline auto cvLocalDofs(const BoxDfmFVElementGeometry& fvGeometry)
    {
        // Default it that all dofs are cv dofs
        return localDofs(fvGeometry);
    }

   /*!
    * \brief Iterator range for sub control volumes faces.
    *
    * Iterates over all scvfs of the bound element.
    * This is a free function found by means of ADL.
    * To iterate over all sub control volume faces of this FVElementGeometry use
    * for (auto&& scvf : scvfs(fvGeometry)).
    */
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolumeFace>::const_iterator>
    scvfs(const BoxDfmFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolumeFace>::const_iterator;
        return Dune::IteratorRange<Iter>(fvGeometry.scvfs_.begin(), fvGeometry.scvfs_.end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    { return gridGeometry().feCache().get(element_->type()).localBasis(); }

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return scvs_.size(); }

    //! The total number of element-local dofs
    std::size_t numLocalDofs() const
    { return element().geometry().corners(); }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    { return scvfs_.size(); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    BoxDfmFVElementGeometry bind(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    /*!
     * \brief Binding of an element, has to be called before using the fvgeometries
     *        Prepares all the volume variables within the element.
     * \note For the box scheme, bind() and bindElement() are identical, but the
     *       distinction is here for the sake of compatibility with cc schemes.
     */
    void bind(const Element& element) &
    { this->bindElement(element); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    BoxDfmFVElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

   /*!
    * \brief Binding of an element, has to be called before using the fvgeometries
    *        Prepares all the volume variables within the element.
    */
    void bindElement(const Element& element) &
    {
        element_ = element;
        eIdx_ = gridGeometry().elementMapper().index(element);
        makeElementGeometries_();
    }

    //! The global finite volume geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return ggCache_->gridGeometry(); }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    //! Create the geometry of a given sub control volume
    typename SubControlVolume::Traits::Geometry geometry(const SubControlVolume& scv) const
    {
        if (scv.isOnFracture())
            DUNE_THROW(Dune::InvalidStateException, "The geometry object cannot be defined for fracture scvs "
                                                    "because the number of known corners is insufficient. "
                                                    "You can do this manually by extracting the corners from this scv "
                                                    "and extruding them by the corresponding aperture. ");

        const GeometryHelper geometryHelper(element().geometry());
        const auto corners = geometryHelper.getScvCorners(scv.index());
        using ScvGeometry = typename SubControlVolume::Traits::Geometry;
        return { Dune::GeometryTypes::cube(ScvGeometry::mydimension), corners };
    }

    //! Create the geometry of a given sub control volume face
    typename SubControlVolumeFace::Traits::Geometry geometry(const SubControlVolumeFace& scvf) const
    {
        if (scvf.isOnFracture())
            DUNE_THROW(Dune::InvalidStateException, "The geometry object cannot be defined for fracture scvs "
                                                    "because the number of known corners is insufficient. "
                                                    "You can do this manually by extracting the corners from this scv "
                                                    "and extruding them by the corresponding aperture. ");
        const GeometryHelper geometryHelper(element().geometry());
        const auto corners = geometryHelper.getScvfCorners(scvf.indexInElement());
        using ScvfGeometry = typename SubControlVolumeFace::Traits::Geometry;
        return { Dune::GeometryTypes::cube(ScvfGeometry::mydimension), corners };
    }

private:

    void makeElementGeometries_()
    {
        // get the element geometry
        const auto& element = *element_;
        const auto elementGeometry = element.geometry();
        const auto refElement = referenceElement(element);

        // get the sub control volume geometries of this element
        GeometryHelper geometryHelper(elementGeometry);

        // construct the sub control volumes
        scvs_.resize(elementGeometry.corners());
        for (LocalIndexType scvLocalIdx = 0; scvLocalIdx < elementGeometry.corners(); ++scvLocalIdx)
        {
            // get associated dof index
            const auto dofIdxGlobal = gridGeometry().vertexMapper().subIndex(element, scvLocalIdx, dim);

            // add scv to the local container
            scvs_[scvLocalIdx] = SubControlVolume(geometryHelper,
                                                  scvLocalIdx,
                                                  eIdx_,
                                                  dofIdxGlobal);
        }

        // construct the sub control volume faces
        const auto numInnerScvf = (dim==1) ? 1 : element.subEntities(dim-1);
        scvfs_.resize(numInnerScvf);

        LocalIndexType scvfLocalIdx = 0;
        for (; scvfLocalIdx < numInnerScvf; ++scvfLocalIdx)
        {
            // find the local scv indices this scvf is connected to
            std::vector<LocalIndexType> localScvIndices({static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 0, dim)),
                                                         static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 1, dim))});

            scvfs_[scvfLocalIdx] = SubControlVolumeFace(geometryHelper,
                                                        element,
                                                        scvfLocalIdx,
                                                        std::move(localScvIndices));
        }

        // construct the ...
        // ... sub-control volume faces on the domain boundary
        // ... sub-control volumes on fracture facets
        // ... sub-control volume faces on fracture facets
        // NOTE We do not construct fracture scvfs on boundaries here!
        //      That means specifying Neumann fluxes on only the fractures is not possible
        //      However, it is difficult to interpret this here as a fracture ending on the boundary
        //      could also be connected to a facet which is both boundary and fracture at the same time!
        //      In that case, the fracture boundary scvf wouldn't make sense. In order to do it properly
        //      we would have to find only those fractures that are at the boundary and aren't connected
        //      to a fracture which is a boundary.
        LocalIndexType scvLocalIdx = element.subEntities(dim);
        for (const auto& intersection : intersections(gridGeometry().gridView(), element))
        {
            // first, obtain all vertex indices on this intersection
            const auto& isGeometry = intersection.geometry();
            const auto numCorners = isGeometry.corners();
            const auto idxInInside = intersection.indexInInside();

            std::vector<GridIndexType> isVertexIndices(numCorners);
            for (unsigned int vIdxLocal = 0; vIdxLocal < numCorners; ++vIdxLocal)
                isVertexIndices[vIdxLocal] = gridGeometry().vertexMapper().subIndex(element,
                                                                                      refElement.subEntity(idxInInside, 1, vIdxLocal, dim),
                                                                                      dim);

            if (intersection.boundary())
            {
                for (unsigned int isScvfLocalIdx = 0; isScvfLocalIdx < numCorners; ++isScvfLocalIdx)
                {
                    // find the scv this scvf is connected to
                    const LocalIndexType insideScvIdx = static_cast<LocalIndexType>(refElement.subEntity(idxInInside, 1, isScvfLocalIdx, dim));
                    std::vector<LocalIndexType> localScvIndices = {insideScvIdx, insideScvIdx};

                    scvfs_.emplace_back(geometryHelper,
                                        intersection,
                                        isScvfLocalIdx,
                                        scvfLocalIdx,
                                        std::move(localScvIndices));

                    // increment local counter
                    scvfLocalIdx++;
                }
            }

            // maybe add fracture scvs & scvfs
            if (this->gridGeometry().isOnFracture(element, intersection))
            {
                // add fracture scv for each vertex of intersection
                const auto curNumScvs = scvs_.size();
                scvs_.reserve(curNumScvs+numCorners);
                for (unsigned int vIdxLocal = 0; vIdxLocal < numCorners; ++vIdxLocal)
                    scvs_.emplace_back(geometryHelper,
                                       intersection,
                                       isGeometry,
                                       vIdxLocal,
                                       static_cast<LocalIndexType>(refElement.subEntity(idxInInside, 1, vIdxLocal, dim)),
                                       scvLocalIdx++,
                                       idxInInside,
                                       eIdx_,
                                       isVertexIndices[vIdxLocal]);

                // add fracture scvf for each edge of the intersection in 3d
                if (dim == 3)
                {
                    const auto& faceRefElement = referenceElement(isGeometry);
                    for (unsigned int edgeIdx = 0; edgeIdx < faceRefElement.size(1); ++edgeIdx)
                    {
                        // inside/outside scv indices in face local node numbering
                        std::vector<LocalIndexType> localScvIndices({static_cast<LocalIndexType>(faceRefElement.subEntity(edgeIdx, 1, 0, dim-1)),
                                                                     static_cast<LocalIndexType>(faceRefElement.subEntity(edgeIdx, 1, 1, dim-1))});

                        // add offset to get the right scv indices
                        std::for_each( localScvIndices.begin(),
                                       localScvIndices.end(),
                                       [curNumScvs] (auto& elemLocalIdx) { elemLocalIdx += curNumScvs; } );

                        // add scvf
                        scvfs_.emplace_back(geometryHelper,
                                            intersection,
                                            edgeIdx,
                                            scvfLocalIdx++,
                                            std::move(localScvIndices),
                                            intersection.boundary());
                    }
                }

                // dim == 2, intersection is an edge, make 1 scvf
                else
                {
                    // inside/outside scv indices in face local node numbering
                    std::vector<LocalIndexType> localScvIndices({0, 1});

                    // add offset such that the fracture scvs above are addressed
                    std::for_each( localScvIndices.begin(),
                                   localScvIndices.end(),
                                   [curNumScvs] (auto& elemLocalIdx) { elemLocalIdx += curNumScvs; } );

                    // add scvf
                    scvfs_.emplace_back(geometryHelper,
                                        intersection,
                                        /*idxOnIntersection*/0,
                                        scvfLocalIdx++,
                                        std::move(localScvIndices),
                                        intersection.boundary());
                }
            }
        }
    }

    //! The bound element
    std::optional<Element> element_;
    GridIndexType eIdx_;

    //! The global geometry cache
    const GGCache* ggCache_;

    //! The element-local scvs & scvfs
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;
};

} // end namespace Dumux

#endif
