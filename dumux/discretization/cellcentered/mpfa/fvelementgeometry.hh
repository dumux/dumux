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
 * \ingroup CCMpfaDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for cell-centered mpfa models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element in the local scope we are restricting to, e.g. stencil or element.
 */
#ifndef DUMUX_DISCRETIZATION_CCMPFA_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_CCMPFA_FV_ELEMENT_GEOMETRY_HH

#include <dune/common/exceptions.hh>
#include <dune/common/iteratorrange.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/scvandscvfiterators.hh>

namespace Dumux
{
/*!
 * \ingroup CCMpfaDiscretization
 * \brief Stencil-local finite volume geometry (scvs and scvfs) for cell-centered mpfa models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element in the local scope we are restricting to, e.g. stencil or element.
 * \tparam GG the finite volume grid geometry type
 * \tparam enableFVGridGeometryCache if the grid geometry is cached or not
 * \note This class is specialized for versions with and without caching the fv geometries on the grid view
 */
template<class GG, bool EnableFVGridGeometryCache>
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
    using Element = typename GridView::template Codim<0>::Entity;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

    static constexpr int dim = GridView::dimension;

public:
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using FVGridGeometry = GG;
    //! the maximum number of scvs per element
    static constexpr std::size_t maxNumElementScvs = 1;
    //! the maximum number of scvfs per element (use cubes for maximum)
    static constexpr std::size_t maxNumElementScvfs = dim == 3 ? 24 : 8;

    //! Constructor
    CCMpfaFVElementGeometry(const FVGridGeometry& fvGridGeometry)
    : fvGridGeometryPtr_(&fvGridGeometry) {}

    //! Get an element sub control volume with a global scv index
    const SubControlVolume& scv(GridIndexType scvIdx) const
    {
        return fvGridGeometry().scv(scvIdx);
    }

    //! Get an element sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const
    {
        return fvGridGeometry().scvf(scvfIdx);
    }

    //! Get the scvf on the same face but from the other side
    //! Note that e.g. the normals might be different in the case of surface grids
    const SubControlVolumeFace& flipScvf(GridIndexType scvfIdx, unsigned int outsideScvIdx = 0) const
    {
        return fvGridGeometry().flipScvf(scvfIdx, outsideScvIdx);
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
        const auto& g = fvGeometry.fvGridGeometry();
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
        return fvGridGeometry().scvfIndicesOfScv(scvIndices_[0]).size();
    }

    //! Binding of an element, called by the local assembler to prepare element assembly
    void bind(const Element& element)
    {
        this->bindElement(element);
    }

    //! Bind only element-local
    void bindElement(const Element& element)
    {
        scvIndices_[0] = fvGridGeometry().elementMapper().index(element);
    }

    //! The global finite volume geometry we are a restriction of
    const FVGridGeometry& fvGridGeometry() const
    { return *fvGridGeometryPtr_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return fvGridGeometry().hasBoundaryScvf(scvIndices_[0]); }

private:

    std::array<GridIndexType, 1> scvIndices_;
    const FVGridGeometry* fvGridGeometryPtr_;
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
    using Element = typename GridView::template Codim<0>::Entity;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using MpfaHelper = typename GG::MpfaHelper;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using CoordScalar = typename GridView::ctype;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

public:
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using FVGridGeometry = GG;
    //! the maximum number of scvs per element
    static constexpr std::size_t maxNumElementScvs = 1;
    //! the maximum number of scvfs per element (use cubes for maximum)
    static constexpr std::size_t maxNumElementScvfs = dim == 3 ? 24 : 8;

    //! Constructor
    CCMpfaFVElementGeometry(const FVGridGeometry& fvGridGeometry)
    : fvGridGeometryPtr_(&fvGridGeometry) {}

    //! Get an elment sub control volume with a global scv index
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
        return scvf( fvGridGeometry().flipScvfIdx(scvfIdx, outsideScvIdx) );
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
    //! called by the local assembler to prepare element assembly
    void bind(const Element& element)
    {
        // make inside geometries
        bindElement(element);

        // get some references for convenience
        const auto globalI = fvGridGeometry().elementMapper().index(element);
        const auto& assemblyMapI = fvGridGeometry().connectivityMap()[globalI];

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
            makeNeighborGeometries(fvGridGeometry().element(dataJ.globalJ),
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
        //         neighborScvs_.emplace_back(fvGridGeometry().element(globalJ).geometry(), globalJ);
        //         neighborScvIndices_.emplace_back(globalJ);
        //     }
        // }
    }

    //! Binding of an element preparing the geometries only inside the element
    void bindElement(const Element& element)
    {
        clear();
        makeElementGeometries(element);
    }

    //! The global finite volume geometry we are a restriction of
    const FVGridGeometry& fvGridGeometry() const
    { return *fvGridGeometryPtr_; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return hasBoundaryScvf_; }

private:

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
        const auto eIdx = fvGridGeometry().elementMapper().index(element);
        scvs_[0] = SubControlVolume(element.geometry(), eIdx);
        scvIndices_[0] = eIdx;

        // get data on the scv faces
        const auto& scvFaceIndices = fvGridGeometry().scvfIndicesOfScv(eIdx);
        const auto& neighborVolVarIndices = fvGridGeometry().neighborVolVarIndices(eIdx);

        // the quadrature point parameterizaion to be used on scvfs
        static const auto q = getParam<CoordScalar>("Mpfa.Q");

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
        for (const auto& is : intersections(fvGridGeometry().gridView(), element))
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
            const auto refElement = ReferenceElements::general(eg.type());

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
                auto vIdxGlobal = fvGridGeometry().vertexMapper().subIndex(e, vIdxLocal, dim);

                // do not build scvfs connected to a processor boundary
                if (fvGridGeometry().isGhostVertex(vIdxGlobal))
                    continue;

                hasBoundaryScvf_ = (hasBoundaryScvf_ || is.boundary());

                scvfs_.emplace_back(MpfaHelper(),
                                    MpfaHelper::getScvfCorners(isPositions, numCorners, c),
                                    is.centerUnitOuterNormal(),
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
        const auto& scvFaceIndices = fvGridGeometry().scvfIndicesOfScv(eIdxGlobal);
        const auto& neighborVolVarIndices = fvGridGeometry().neighborVolVarIndices(eIdxGlobal);

        // the quadrature point parameterizaion to be used on scvfs
        static const auto q = getParam<CoordScalar>("Mpfa.Q");

        // for network grids we only want to do one scvf per half facet
        // this approach assumes conforming grids at branching facets
        std::vector<bool> finishedFacets;
        if (dim < dimWorld)
            finishedFacets.resize(element.subEntities(1), false);

        int scvfCounter = 0;
        for (const auto& is : intersections(fvGridGeometry().gridView(), element))
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
            const auto refElement = ReferenceElements::general(eg.type());

            // Set up a container with all relevant positions for scvf corner computation
            const auto numCorners = is.geometry().corners();
            const auto isPositions = MpfaHelper::computeScvfCornersOnIntersection(eg,
                                                                                  refElement,
                                                                                  indexInElement,
                                                                                  numCorners);

            // make the scv faces belonging to each corner of the intersection
            for (int c = 0; c < numCorners; ++c)
            {
                // only build the scvf if it is in the list of necessary indices
                if (!MpfaHelper::vectorContainsValue(scvfIndices, scvFaceIndices[scvfCounter])
                    && !MpfaHelper::vectorContainsValue(additionalScvfs, scvFaceIndices[scvfCounter]))
                {
                    // increment counter either way
                    scvfCounter++;
                    continue;
                }

                // get the global vertex index the scv face is connected to
                auto vIdxLocal = refElement.subEntity(indexInElement, 1, c, dim);
                auto vIdxGlobal = fvGridGeometry().vertexMapper().subIndex(e, vIdxLocal, dim);

                // do not build scvfs connected to a processor boundary
                if (fvGridGeometry().isGhostVertex(vIdxGlobal))
                    continue;

                // build scvf
                neighborScvfs_.emplace_back(MpfaHelper(),
                                            MpfaHelper::getScvfCorners(isPositions, numCorners, c),
                                            is.centerUnitOuterNormal(),
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

    const FVGridGeometry* fvGridGeometryPtr_;

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
