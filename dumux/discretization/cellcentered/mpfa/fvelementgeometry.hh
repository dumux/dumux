// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCMpfaDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for cell-centered mpfa models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element in the local scope we are restricting to, e.g. stencil or element.
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_CCMPFA_FV_ELEMENT_GEOMETRY_HH

#include <optional>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/geometry/type.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/scvandscvfiterators.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail::Mpfa {

template<typename GridGeometry, typename SubControlVolumeFace>
typename SubControlVolumeFace::Traits::Geometry makeScvfGeometry(const GridGeometry& gridGeometry,
                                                                 const SubControlVolumeFace& scvf)
{
    static constexpr int dim = GridGeometry::GridView::dimension;

    const auto& facetInfo = scvf.facetInfo();
    const auto element = gridGeometry.element(facetInfo.elementIndex);
    const auto elemGeo = element.geometry();
    const auto refElement = referenceElement(elemGeo);
    for (const auto& is : intersections(gridGeometry.gridView(), element))
    {
        if (is.indexInInside() == facetInfo.facetIndex)
        {
            const auto numCorners = is.geometry().corners();
            const auto isPositions = GridGeometry::MpfaHelper::computeScvfCornersOnIntersection(
                elemGeo, refElement, facetInfo.facetIndex, numCorners
            );
            return {
                Dune::GeometryTypes::cube(dim-1),
                GridGeometry::MpfaHelper::getScvfCorners(
                    isPositions, numCorners, facetInfo.facetCornerIndex
                )
            };
        }
    }
    DUNE_THROW(Dune::InvalidStateException, "Could not construct scvf geometry");
}

template<typename GridGeometry, typename SubControlVolumeFace>
auto getVertexCorner(const GridGeometry& gridGeometry, const SubControlVolumeFace& scvf)
{
    static constexpr int dim = GridGeometry::GridView::dimension;

    const auto& facetInfo = scvf.facetInfo();
    const auto element = gridGeometry.element(facetInfo.elementIndex);
    const auto elemGeo = element.geometry();
    const auto refElement = referenceElement(elemGeo);
    return elemGeo.global(refElement.position(
        refElement.subEntity(facetInfo.facetIndex, 1, facetInfo.facetCornerIndex, dim),
        dim
    ));
}

template<typename GridGeometry, typename SubControlVolumeFace>
auto getFacetCorner(const GridGeometry& gridGeometry, const SubControlVolumeFace& scvf)
{
    const auto& facetInfo = scvf.facetInfo();
    const auto element = gridGeometry.element(facetInfo.elementIndex);
    const auto elemGeo = element.geometry();
    const auto refElement = referenceElement(elemGeo);
    return elemGeo.global(refElement.position(facetInfo.facetIndex, 1));
}

} // namespace Detail::Mpfa
#endif // DOXYGEN

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for cell-centered mpfa models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element in the local scope we are restricting to, e.g. stencil or element.
 * \tparam GG the finite volume grid geometry type
 * \tparam enableGridGeometryCache if the grid geometry is cached or not
 * \note This class is specialized for versions with and without caching the fv geometries on the grid view
 */
template<class GG, bool enableGridGeometryCache>
class CCMpfaFVElementGeometry;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for cell-centered mpfa models
 *        Specialization for grid caching enabled
 * \note The finite volume geometries are stored in the corresponding FVGridGeometry
 */
template<class GG>
class CCMpfaFVElementGeometry<GG, true>
{
    using ThisType = CCMpfaFVElementGeometry<GG, true>;
    using GridView = typename GG::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    //! export type of the element
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;
    //! the maximum number of scvs per element
    static constexpr std::size_t maxNumElementScvs = 1;
    //! the maximum number of scvfs per element (use cubes for maximum)
    static constexpr std::size_t maxNumElementScvfs = dim == 3 ? 24 : 8;

    //! Constructor
    CCMpfaFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometryPtr_(&gridGeometry) {}

    //! Get an element sub control volume with a global scv index
    const SubControlVolume& scv(GridIndexType scvIdx) const
    {
        return gridGeometry().scv(scvIdx);
    }

    //! Get an element sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const
    {
        return gridGeometry().scvf(scvfIdx);
    }

    //! Get the scvf on the same face but from the other side
    //! Note that e.g. the normals might be different in the case of surface grids
    const SubControlVolumeFace& flipScvf(GridIndexType scvfIdx, unsigned int outsideScvIdx = 0) const
    {
        return gridGeometry().flipScvf(scvfIdx, outsideScvIdx);
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element (not including neighbor scvs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange< ScvIterator<SubControlVolume, std::array<GridIndexType, 1>, ThisType> >
    scvs(const CCMpfaFVElementGeometry& fvGeometry)
    {
        using ScvIterator = Dumux::ScvIterator<SubControlVolume, std::array<GridIndexType, 1>, ThisType>;
        return Dune::IteratorRange<ScvIterator>(ScvIterator(fvGeometry.scvIndices_.begin(), fvGeometry),
                                                ScvIterator(fvGeometry.scvIndices_.end(), fvGeometry));
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element (not including neighbor scvfs)
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange< ScvfIterator<SubControlVolumeFace, std::vector<GridIndexType>, ThisType> >
    scvfs(const CCMpfaFVElementGeometry& fvGeometry)
    {
        const auto& g = fvGeometry.gridGeometry();
        const auto scvIdx = fvGeometry.scvIndices_[0];
        using ScvfIterator = Dumux::ScvfIterator<SubControlVolumeFace, std::vector<GridIndexType>, ThisType>;
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
        return gridGeometry().scvfIndicesOfScv(scvIndices_[0]).size();
    }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    CCMpfaFVElementGeometry bind(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    void bind(const Element& element) &
    {
        this->bindElement(element);
    }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    CCMpfaFVElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! Bind only element-local
    void bindElement(const Element& element) &
    {
        element_ = element;
        scvIndices_[0] = gridGeometry().elementMapper().index(element);
    }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    //! The global finite volume geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return *gridGeometryPtr_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return gridGeometry().hasBoundaryScvf(scvIndices_[0]); }

    //! Create the geometry of a given sub control volume
    typename Element::Geometry geometry(const SubControlVolume& scv) const
    { return gridGeometryPtr_->element(scv.dofIndex()).geometry(); }

    //! Create the geometry of a given sub control volume face
    typename SubControlVolumeFace::Traits::Geometry geometry(const SubControlVolumeFace& scvf) const
    { return Detail::Mpfa::makeScvfGeometry(gridGeometry(), scvf); }

    //! Return the position of the scvf corner that coincides with an element vertex
    typename SubControlVolumeFace::Traits::GlobalPosition vertexCorner(const SubControlVolumeFace& scvf) const
    { return Detail::Mpfa::getVertexCorner(gridGeometry(), scvf); }

    //! Return the corner of the scvf that is inside the facet the scvf is embedded in
    typename SubControlVolumeFace::Traits::GlobalPosition facetCorner(const SubControlVolumeFace& scvf) const
    { return Detail::Mpfa::getFacetCorner(gridGeometry(), scvf); }

private:

    std::optional<Element> element_;
    std::array<GridIndexType, 1> scvIndices_;
    const GridGeometry* gridGeometryPtr_;
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for cell-centered TPFA models
 *        Specialization for grid caching disabled
 */
template<class GG>
class CCMpfaFVElementGeometry<GG, false>
{
    using ThisType = CCMpfaFVElementGeometry<GG, false>;
    using GridView = typename GG::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using MpfaHelper = typename GG::MpfaHelper;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using CoordScalar = typename GridView::ctype;

public:
    //! export type of the element
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometries
    using GridGeometry = GG;
    //! the maximum number of scvs per element
    static constexpr std::size_t maxNumElementScvs = 1;
    //! the maximum number of scvfs per element (use cubes for maximum)
    static constexpr std::size_t maxNumElementScvfs = dim == 3 ? 24 : 8;

    //! Constructor
    CCMpfaFVElementGeometry(const GridGeometry& gridGeometry)
    : gridGeometryPtr_(&gridGeometry) {}

    //! Get an element sub control volume with a global scv index
    //! We separate element and neighbor scvs to speed up mapping
    const SubControlVolume& scv(GridIndexType scvIdx) const
    {
        if (scvIdx == scvIndices_[0])
            return scvs_[0];
        else
            return neighborScvs_[findLocalIndex(scvIdx, neighborScvIndices_)];
    }

    //! Get an element sub control volume face with a global scvf index
    //! We separate element and neighbor scvfs to speed up mapping
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const
    {
        auto it = std::find(scvfIndices_.begin(), scvfIndices_.end(), scvfIdx);
        if (it != scvfIndices_.end())
            return scvfs_[std::distance(scvfIndices_.begin(), it)];
        else
            return neighborScvfs_[findLocalIndex(scvfIdx, neighborScvfIndices_)];
    }

    //! Get the scvf on the same face but from the other side
    //! Note that e.g. the normals might be different in the case of surface grids
    const SubControlVolumeFace& flipScvf(GridIndexType scvfIdx, unsigned int outsideScvIdx = 0) const
    {
        return scvf( gridGeometry().flipScvfIdx(scvfIdx, outsideScvIdx) );
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
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    CCMpfaFVElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement_(element);
        return std::move(*this);
    }

    void bindElement(const Element& element) &
    {
        this->bindElement_(element);
    }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    CCMpfaFVElementGeometry bind(const Element& element) &&
    {
        this->bind_(element);
        return std::move(*this);
    }

    void bind(const Element& element) &
    {
        this->bind_(element);
    }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    //! The global finite volume geometry we are a restriction of
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
    { return Detail::Mpfa::makeScvfGeometry(gridGeometry(), scvf); }

    //! Return the position of the scvf corner that coincides with an element vertex
    typename SubControlVolumeFace::Traits::GlobalPosition vertexCorner(const SubControlVolumeFace& scvf) const
    { return Detail::Mpfa::getVertexCorner(gridGeometry(), scvf); }

    //! Return the corner of the scvf that is inside the facet the scvf is embedded in
    typename SubControlVolumeFace::Traits::GlobalPosition facetCorner(const SubControlVolumeFace& scvf) const
    { return Detail::Mpfa::getFacetCorner(gridGeometry(), scvf); }

private:

    //! Binding of an element preparing the geometries of the whole stencil
    //! called by the local assembler to prepare element assembly
    void bind_(const Element& element)
    {
        // make inside geometries
        bindElement_(element);

        // get some references for convenience
        const auto globalI = gridGeometry().elementMapper().index(element);
        const auto& assemblyMapI = gridGeometry().connectivityMap()[globalI];

        // reserve memory
        const auto numNeighbors = assemblyMapI.size();
        const auto numNeighborScvfs = numNeighborScvfs_(assemblyMapI);
        neighborScvs_.reserve(numNeighbors);
        neighborScvIndices_.reserve(numNeighbors);
        neighborScvfIndices_.reserve(numNeighborScvfs);
        neighborScvfs_.reserve(numNeighborScvfs);

        // make neighbor geometries
        // use the assembly map to determine which faces are necessary
        for (const auto& dataJ : assemblyMapI)
            makeNeighborGeometries(gridGeometry().element(dataJ.globalJ),
                                   dataJ.globalJ,
                                   dataJ.scvfsJ,
                                   dataJ.additionalScvfs);

        // //! TODO Check if user added additional DOF dependencies, i.e. the residual of DOF globalI depends
        // //! on additional DOFs not included in the discretization schemes' occupation pattern
        // const auto& additionalDofDependencies = problem.getAdditionalDofDependencies(globalI);
        // if (!additionalDofDependencies.empty())
        // {
        //     const auto newNumNeighbors = neighborScvs_.size() + additionalDofDependencies.size();
        //     neighborScvs_.reserve(newNumNeighbors);
        //     neighborScvIndices_.reserve(newNumNeighbors);
        //     for (auto globalJ : additionalDofDependencies)
        //     {
        //         neighborScvs_.emplace_back(gridGeometry().element(globalJ).geometry(), globalJ);
        //         neighborScvIndices_.emplace_back(globalJ);
        //     }
        // }
    }

    //! Binding of an element preparing the geometries only inside the element
    void bindElement_(const Element& element)
    {
        clear();
        element_ = element;
        makeElementGeometries(element);
    }

    //! Computes the number of neighboring scvfs that have to be prepared
    template<class DataJContainer>
    std::size_t numNeighborScvfs_(const DataJContainer& dataJContainer)
    {
        std::size_t numNeighborScvfs = 0;
        for (const auto& dataJ : dataJContainer)
            numNeighborScvfs += dataJ.scvfsJ.size() + dataJ.additionalScvfs.size();
        return numNeighborScvfs;
    }

    //! create scvs and scvfs of the bound element
    void makeElementGeometries(const Element& element)
    {
        // make the scv
        const auto eIdx = gridGeometry().elementMapper().index(element);
        scvs_[0] = SubControlVolume(element.geometry(), eIdx);
        scvIndices_[0] = eIdx;

        // get data on the scv faces
        const auto& scvFaceIndices = gridGeometry().scvfIndicesOfScv(eIdx);
        const auto& neighborVolVarIndices = gridGeometry().neighborVolVarIndices(eIdx);

        // the quadrature point parameterizaion to be used on scvfs
        static const auto q = getParam<CoordScalar>("MPFA.Q");

        // reserve memory for the scv faces
        const auto numLocalScvf = scvFaceIndices.size();
        scvfIndices_.reserve(numLocalScvf);
        scvfs_.reserve(numLocalScvf);

        // for network grids we only want to do one scvf per half facet
        // this approach assumes conforming grids at branching facets
        std::vector<bool> finishedFacets;
        if (dim < dimWorld)
            finishedFacets.resize(element.subEntities(1), false);

        int scvfCounter = 0;
        for (const auto& is : intersections(gridGeometry().gridView(), element))
        {
            // if we are dealing with a lower dimensional network
            // only make a new scvf if we haven't handled it yet
            if (dim < dimWorld)
            {
                const auto indexInInside = is.indexInInside();
                if (finishedFacets[indexInInside])
                    continue;
                else
                    finishedFacets[indexInInside] = true;
            }

            // if outside level > inside level, use the outside element in the following
            bool useNeighbor = is.neighbor() && is.outside().level() > element.level();
            const auto& e = useNeighbor ? is.outside() : element;
            const auto indexInElement = useNeighbor ? is.indexInOutside() : is.indexInInside();
            const auto eg = e.geometry();
            const auto refElement = referenceElement(eg);

            // Set up a container with all relevant positions for scvf corner computation
            const auto numCorners = is.geometry().corners();
            const auto isPositions = MpfaHelper::computeScvfCornersOnIntersection(eg,
                                                                                  refElement,
                                                                                  indexInElement,
                                                                                  numCorners);

            // make the scv faces belonging to each corner of the intersection
            for (int c = 0; c < numCorners; ++c)
            {
                // get the global vertex index the scv face is connected to
                auto vIdxLocal = refElement.subEntity(indexInElement, 1, c, dim);
                auto vIdxGlobal = gridGeometry().vertexMapper().subIndex(e, vIdxLocal, dim);

                // do not build scvfs connected to a processor boundary
                if (gridGeometry().isGhostVertex(vIdxGlobal))
                    continue;

                hasBoundaryScvf_ = (hasBoundaryScvf_ || is.boundary());
                typename SubControlVolumeFace::FacetInfo facetInfo{
                    gridGeometry().elementMapper().index(e),
                    useNeighbor ? is.indexInOutside() : is.indexInInside(),
                    c
                };
                scvfs_.emplace_back(MpfaHelper(),
                                    MpfaHelper::getScvfCorners(isPositions, numCorners, c),
                                    is,
                                    std::move(facetInfo),
                                    vIdxGlobal,
                                    vIdxLocal,
                                    scvFaceIndices[scvfCounter],
                                    eIdx,
                                    neighborVolVarIndices[scvfCounter],
                                    q,
                                    is.boundary());

                scvfIndices_.emplace_back(scvFaceIndices[scvfCounter]);
                scvfCounter++;
            }
        }
    }

    //! create the scv and necessary scvfs of a neighboring element
    template<typename IndexVector>
    void makeNeighborGeometries(const Element& element,
                                GridIndexType eIdxGlobal,
                                const IndexVector& scvfIndices,
                                const IndexVector& additionalScvfs)
    {
        // create the neighbor scv if it doesn't exist yet
        neighborScvs_.emplace_back(element.geometry(), eIdxGlobal);
        neighborScvIndices_.push_back(eIdxGlobal);

        // get data on the scv faces
        const auto& scvFaceIndices = gridGeometry().scvfIndicesOfScv(eIdxGlobal);
        const auto& neighborVolVarIndices = gridGeometry().neighborVolVarIndices(eIdxGlobal);

        // the quadrature point parameterizaion to be used on scvfs
        static const auto q = getParam<CoordScalar>("MPFA.Q");

        // for network grids we only want to do one scvf per half facet
        // this approach assumes conforming grids at branching facets
        std::vector<bool> finishedFacets;
        if (dim < dimWorld)
            finishedFacets.resize(element.subEntities(1), false);

        int scvfCounter = 0;
        for (const auto& is : intersections(gridGeometry().gridView(), element))
        {
            // if we are dealing with a lower dimensional network
            // only make a new scvf if we haven't handled it yet
            if (dim < dimWorld)
            {
                auto indexInInside = is.indexInInside();
                if(finishedFacets[indexInInside])
                    continue;
                else
                    finishedFacets[indexInInside] = true;
            }

            // if outside level > inside level, use the outside element in the following
            bool useNeighbor = is.neighbor() && is.outside().level() > element.level();
            const auto& e = useNeighbor ? is.outside() : element;
            const auto indexInElement = useNeighbor ? is.indexInOutside() : is.indexInInside();
            const auto eg = e.geometry();
            const auto refElement = referenceElement(eg);

            // Set up a container with all relevant positions for scvf corner computation
            const auto numCorners = is.geometry().corners();
            const auto isPositions = MpfaHelper::computeScvfCornersOnIntersection(eg,
                                                                                  refElement,
                                                                                  indexInElement,
                                                                                  numCorners);

            // make the scv faces belonging to each corner of the intersection
            for (int c = 0; c < numCorners; ++c)
            {
                // get the global vertex index the scv face is connected to
                auto vIdxLocal = refElement.subEntity(indexInElement, 1, c, dim);
                auto vIdxGlobal = gridGeometry().vertexMapper().subIndex(e, vIdxLocal, dim);

                // do not build scvfs connected to a processor boundary
                if (gridGeometry().isGhostVertex(vIdxGlobal))
                    continue;

                // only build the scvf if it is in the list of necessary indices
                if (!MpfaHelper::vectorContainsValue(scvfIndices, scvFaceIndices[scvfCounter])
                    && !MpfaHelper::vectorContainsValue(additionalScvfs, scvFaceIndices[scvfCounter]))
                {
                    // increment counter either way
                    scvfCounter++;
                    continue;
                }

                // build scvf
                typename SubControlVolumeFace::FacetInfo facetInfo{
                    gridGeometry().elementMapper().index(e),
                    useNeighbor ? is.indexInOutside() : is.indexInInside(),
                    c
                };
                neighborScvfs_.emplace_back(MpfaHelper(),
                                            MpfaHelper::getScvfCorners(isPositions, numCorners, c),
                                            is,
                                            std::move(facetInfo),
                                            vIdxGlobal,
                                            vIdxLocal,
                                            scvFaceIndices[scvfCounter],
                                            eIdxGlobal,
                                            neighborVolVarIndices[scvfCounter],
                                            q,
                                            is.boundary());

                neighborScvfIndices_.emplace_back(scvFaceIndices[scvfCounter]);

                // increment counter
                scvfCounter++;
            }
        }
    }

    //! map a global index to the local storage index
    unsigned int findLocalIndex(const GridIndexType idx,
                                const std::vector<GridIndexType>& indices) const
    {
        auto it = std::find(indices.begin(), indices.end(), idx);
        assert(it != indices.end() && "Could not find the scv/scvf! Make sure to properly bind this class!");
        return std::distance(indices.begin(), it);
    }

    //! clear all containers
    void clear()
    {
        scvfIndices_.clear();
        scvfs_.clear();

        neighborScvIndices_.clear();
        neighborScvfIndices_.clear();
        neighborScvs_.clear();
        neighborScvfs_.clear();

        hasBoundaryScvf_ = false;
    }

    const GridGeometry* gridGeometryPtr_;
    std::optional<Element> element_;

    // local storage after binding an element
    std::array<GridIndexType, 1> scvIndices_;
    std::vector<GridIndexType> scvfIndices_;
    std::array<SubControlVolume, 1> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;

    std::vector<GridIndexType> neighborScvIndices_;
    std::vector<GridIndexType> neighborScvfIndices_;
    std::vector<SubControlVolume> neighborScvs_;
    std::vector<SubControlVolumeFace> neighborScvfs_;

    bool hasBoundaryScvf_ = false;
};

} // end namespace

#endif
