// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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

#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

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
    using GridIndexType = typename GridView::IndexSet::IndexType;
    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
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
    BoxDfmFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometryPtr_(&gridGeometry) {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(std::size_t scvIdx) const
    { return gridGeometry().scvs(eIdx_)[scvIdx]; }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(std::size_t scvfIdx) const
    { return gridGeometry().scvfs(eIdx_)[scvfIdx]; }

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
        const auto& g = fvGeometry.gridGeometry();
        using Iter = typename std::vector<SubControlVolume>::const_iterator;
        return Dune::IteratorRange<Iter>(g.scvs(fvGeometry.eIdx_).begin(), g.scvs(fvGeometry.eIdx_).end());
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
        const auto& g = fvGeometry.gridGeometry();
        using Iter = typename std::vector<SubControlVolumeFace>::const_iterator;
        return Dune::IteratorRange<Iter>(g.scvfs(fvGeometry.eIdx_).begin(), g.scvfs(fvGeometry.eIdx_).end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    { return gridGeometry().feCache().get(element_->type()).localBasis(); }

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return gridGeometry().scvs(eIdx_).size(); }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    { return gridGeometry().scvfs(eIdx_).size(); }

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
    { return *gridGeometryPtr_; }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    // suppress warnings due to current implementation
    // these interfaces should be used!
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"

    //! Create the geometry of a given sub control volume
    typename SubControlVolume::Traits::Geometry geometry(const SubControlVolume& scv) const
    { return scv.geometry(); }

    //! Create the geometry of a given sub control volume face
    typename SubControlVolumeFace::Traits::Geometry geometry(const SubControlVolumeFace& scvf) const
    { return scvf.geometry(); }

    #pragma GCC diagnostic pop

private:
    const GridGeometry* gridGeometryPtr_;

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

    using GridIndexType = typename GridView::IndexSet::IndexType;

    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
    using GeometryHelper = typename GG::GeometryHelper;
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
    BoxDfmFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometryPtr_(&gridGeometry) {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(std::size_t scvIdx) const
    { return scvs_[scvIdx]; }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(std::size_t scvfIdx) const
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
    { return *gridGeometryPtr_; }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    // suppress warnings due to current implementation
    // these interfaces should be used!
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"

    //! Create the geometry of a given sub control volume
    typename SubControlVolume::Traits::Geometry geometry(const SubControlVolume& scv) const
    { return scv.geometry(); }

    //! Create the geometry of a given sub control volume face
    typename SubControlVolumeFace::Traits::Geometry geometry(const SubControlVolumeFace& scvf) const
    { return scvf.geometry(); }

    #pragma GCC diagnostic pop

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
        using LocalIndexType = typename SubControlVolumeFace::Traits::LocalIndexType;
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

        unsigned int scvfLocalIdx = 0;
        for (; scvfLocalIdx < numInnerScvf; ++scvfLocalIdx)
        {
            // find the local scv indices this scvf is connected to
            std::vector<LocalIndexType> localScvIndices({static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 0, dim)),
                                                         static_cast<LocalIndexType>(refElement.subEntity(scvfLocalIdx, dim-1, 1, dim))});

            scvfs_[scvfLocalIdx] = SubControlVolumeFace(geometryHelper,
                                                        element,
                                                        elementGeometry,
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
                                        isGeometry,
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
                                            isGeometry,
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
                                        isGeometry,
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

    //! The global geometry this is a restriction of
    const GridGeometry* gridGeometryPtr_;

    //! The element-local scvs & scvfs
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;
};

} // end namespace Dumux

#endif
