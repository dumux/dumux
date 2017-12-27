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
 * \brief The finite volume geometry (scvs and scvfs) for cell-centered mpfa models on a grid view
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element of the grid partition.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FV_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FV_GRID_GEOMETRY_HH

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/discretization/basefvgridgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/fvelementgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/connectivitymap.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>
#include <dumux/discretization/cellcentered/mpfa/gridinteractionvolumeindexsets.hh>

namespace Dumux
{
/*!
 * \ingroup CCMpfaDiscretization
 * \brief The finite volume geometry (scvs and scvfs) for cell-centered mpfa models on a grid view
 *        This builds up the sub control volumes and sub control volume faces
 * \note This class is specialized for versions with and without caching the fv geometries on the grid view
 */
template<class TypeTag, bool EnableFVElementGeometryCache>
class CCMpfaFVGridGeometry;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The finite volume geometry (scvs and scvfs) for cell-centered mpfa models on a grid view
 *        This builds up the sub control volumes and sub control volume faces
 * \note For caching enabled we store the fv geometries for the whole grid view which is memory intensive but faster
 */
template<class TypeTag>
class CCMpfaFVGridGeometry<TypeTag, true> : public BaseFVGridGeometry<TypeTag>
{
    using ParentType = BaseFVGridGeometry<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementMapper = typename GET_PROP_TYPE(TypeTag, ElementMapper);

    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);

    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Intersection = typename GridView::Intersection;
    using CoordScalar = typename GridView::ctype;
    using GridIndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = typename PrimaryInteractionVolume::Traits::LocalIndexType;

    using GridIVIndexSets = CCMpfaGridInteractionVolumeIndexSets<TypeTag>;
    using ConnectivityMap = CCMpfaConnectivityMap<TypeTag>;

    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

public:
    using SecondaryIvIndicatorType = std::function<bool(const Element&, const Intersection&, bool)>;

    //! Constructor without indicator function for secondary interaction volumes
    //! Per default, we use the secondary IVs at branching points & boundaries
    CCMpfaFVGridGeometry(const GridView& gridView)
    : ParentType(gridView)
    , secondaryIvIndicator_([] (const Element& e, const Intersection& is, bool isBranching)
                               { return is.boundary() || isBranching; } )
    {}

    //! Constructor with user-defined indicator function for secondary interaction volumes
    CCMpfaFVGridGeometry(const GridView& gridView, const SecondaryIvIndicatorType& indicator)
    : ParentType(gridView)
    , secondaryIvIndicator_(indicator)
    {}

    //! the element mapper is the dofMapper
    //! this is convenience to have better chance to have the same main files for box/tpfa/mpfa...
    const ElementMapper& dofMapper() const { return this->elementMapper(); }

    //! The total number of sub control volumes
    std::size_t numScv() const { return scvs_.size(); }

    //! The total number of sub control volume faces
    std::size_t numScvf() const { return scvfs_.size(); }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const { return numBoundaryScvf_; }

    //! The total number of degrees of freedom
    std::size_t numDofs() const { return this->gridView().size(0); }

    //! Get an element from a global element index
    Element element(GridIndexType eIdx) const { return this->elementMap()[eIdx]; }

    //! Get an element from a sub control volume contained in it
    Element element(const SubControlVolume& scv) const { return this->elementMap()[scv.elementIndex()]; }

    //! Returns true if secondary interaction volumes are used around a given vertex (index).
    //! This specialization is enabled if the use of secondary interaction volumes is active.
    template<bool useSecondary = MpfaHelper::considerSecondaryIVs(), std::enable_if_t<useSecondary, int> = 0>
    bool vertexUsesSecondaryInteractionVolume(GridIndexType vIdxGlobal) const
    { return secondaryInteractionVolumeVertices_[vIdxGlobal]; }

    //! Returns true if secondary interaction volumes are used around a given vertex (index).
    //! If the use of secondary interaction volumes is disabled, this can be evaluated at compile time.
    template<bool useSecondary = MpfaHelper::considerSecondaryIVs(), std::enable_if_t<!useSecondary, int> = 0>
    constexpr bool vertexUsesSecondaryInteractionVolume(GridIndexType vIdxGlobal) const { return false; }

    //! update all fvElementGeometries (do this again after grid adaption)
    void update()
    {
        ParentType::update();

        // stop the time required for the update
        Dune::Timer timer;

        // clear scvfs container
        scvfs_.clear();

        // determine the number of geometric entities
        const auto numVert = this->gridView().size(dim);
        const auto numScvs = numDofs();
        std::size_t numScvf = MpfaHelper::getGlobalNumScvf(this->gridView());

        // resize containers
        scvs_.resize(numScvs);
        scvfs_.reserve(numScvf);
        scvfIndicesOfScv_.resize(numScvs);

        // Some methods require to use a second type of interaction volume, e.g.
        // around vertices on the boundary or branching points (surface grids)
        secondaryInteractionVolumeVertices_.resize(numVert, false);

        // find vertices on processor boundaries
        const auto isGhostVertex = MpfaHelper::findGhostVertices(this->gridView(), this->vertexMapper());

        // instantiate the dual grid index set (to be used for construction of interaction volumes)
        CCMpfaDualGridIndexSet<GridIndexType, LocalIndexType, dim> dualIdSet(this->gridView());

        // Build the SCVs and SCV faces
        GridIndexType scvfIdx = 0;
        numBoundaryScvf_ = 0;
        for (const auto& element : elements(this->gridView()))
        {
            const auto eIdx = this->elementMapper().index(element);

            // the element geometry
            auto elementGeometry = element.geometry();

            // The local scvf index set
            std::vector<GridIndexType> scvfIndexSet;
            scvfIndexSet.reserve(MpfaHelper::getNumLocalScvfs(elementGeometry.type()));

            // for network grids there might be multiple intersections with the same geometryInInside
            // we indentify those by the indexInInside for now (assumes conforming grids at branching facets)
            std::vector<std::vector<GridIndexType>> outsideIndices;
            if (dim < dimWorld)
            {
                outsideIndices.resize(element.subEntities(1));
                for (const auto& intersection : intersections(this->gridView(), element))
                {
                    if (intersection.neighbor())
                    {
                        const auto& outside = intersection.outside();
                        auto nIdx = this->elementMapper().index(outside);
                        outsideIndices[intersection.indexInInside()].push_back(nIdx);
                    }
                }
            }

            // construct the sub control volume faces
            for (const auto& is : intersections(this->gridView(), element))
            {
                const auto indexInInside = is.indexInInside();
                const bool boundary = is.boundary();
                const bool neighbor = is.neighbor();

                // for surface grids, skip the rest if handled already
                if (dim < dimWorld && neighbor && outsideIndices[indexInInside].empty())
                    continue;

                // if outside level > inside level, use the outside element in the following
                const bool useNeighbor = neighbor && is.outside().level() > element.level();
                const auto& e = useNeighbor ? is.outside() : element;
                const auto indexInElement = useNeighbor ? is.indexInOutside() : indexInInside;
                const auto eg = e.geometry();
                const auto& refElement = ReferenceElements::general(eg.type());

                // Set up a container with all relevant positions for scvf corner computation
                const auto numCorners = is.geometry().corners();
                const auto isPositions = MpfaHelper::computeScvfCornersOnIntersection(eg,
                                                                                      refElement,
                                                                                      indexInElement,
                                                                                      numCorners);

                // evaluate if vertices on this intersection use primary/secondary IVs
                const bool isBranchingPoint = dim < dimWorld ? outsideIndices[indexInInside].size() > 1 : false;
                const bool usesSecondaryIV = secondaryIvIndicator_(element, is, isBranchingPoint);

                // make the scv faces belonging to each corner of the intersection
                for (std::size_t c = 0; c < numCorners; ++c)
                {
                    // get the global vertex index the scv face is connected to
                    const auto vIdxLocal = refElement.subEntity(indexInElement, 1, c, dim);
                    const auto vIdxGlobal = this->vertexMapper().subIndex(e, vIdxLocal, dim);

                    // do not build scvfs connected to a processor boundary
                    if (isGhostVertex[vIdxGlobal])
                        continue;

                    // if this vertex is tagged to use the secondary IVs, store info
                    if (usesSecondaryIV)
                        secondaryInteractionVolumeVertices_[vIdxGlobal] = true;

                    // the quadrature point parameterizarion to be used on scvfs
                    static const Scalar q = getParamFromGroup<Scalar>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Mpfa.Q");

                    // make the scv face (for non-boundary scvfs on network grids, use precalculated outside indices)
                    const auto& outsideScvIndices = [&] ()
                                                    {
                                                        if (!boundary)
                                                        {
                                                            return dim == dimWorld ?
                                                                   std::vector<GridIndexType>({this->elementMapper().index(is.outside())}) :
                                                                   outsideIndices[indexInInside];
                                                        }
                                                        else
                                                        {
                                                            // compute boundary scv idx and increment counter
                                                            const GridIndexType bIdx = numScvs + numBoundaryScvf_;
                                                            numBoundaryScvf_++;
                                                            return std::vector<GridIndexType>(1, bIdx);
                                                        }
                                                    } ();

                    scvfIndexSet.push_back(scvfIdx);
                    scvfs_.emplace_back(MpfaHelper(),
                                        MpfaHelper::getScvfCorners(isPositions, numCorners, c),
                                        is.centerUnitOuterNormal(),
                                        vIdxGlobal,
                                        vIdxLocal,
                                        scvfIdx,
                                        eIdx,
                                        outsideScvIndices,
                                        q,
                                        boundary);

                    // insert the scvf data into the dual grid index set
                    dualIdSet[vIdxGlobal].insert(scvfs_.back());

                    // increment scvf counter
                    scvfIdx++;
                }

                // for network grids, clear outside indices to not make a second scvf on that facet
                if (dim < dimWorld)
                    outsideIndices[indexInInside].clear();
            }

            // create sub control volume for this element
            scvs_[eIdx] = SubControlVolume(std::move(elementGeometry), eIdx);

            // Save the scvf indices belonging to this scv to build up fv element geometries fast
            scvfIndicesOfScv_[eIdx] = scvfIndexSet;
        }

        // Make the flip index set for network and surface grids
        if (dim < dimWorld)
        {
            flipScvfIndices_.resize(scvfs_.size());
            for (const auto& scvf : scvfs_)
            {
                if (scvf.boundary())
                    continue;

                const auto numOutsideScvs = scvf.numOutsideScvs();
                const auto vIdxGlobal = scvf.vertexIndex();
                const auto insideScvIdx = scvf.insideScvIdx();

                flipScvfIndices_[scvf.index()].resize(numOutsideScvs);
                for (std::size_t i = 0; i < numOutsideScvs; ++i)
                {
                    const auto outsideScvIdx = scvf.outsideScvIdx(i);
                    for (auto outsideScvfIndex : scvfIndicesOfScv_[outsideScvIdx])
                    {
                        const auto& outsideScvf = this->scvf(outsideScvfIndex);
                        if (outsideScvf.vertexIndex() == vIdxGlobal &&
                            MpfaHelper::vectorContainsValue(outsideScvf.outsideScvIndices(), insideScvIdx))
                        {
                            flipScvfIndices_[scvf.index()][i] = outsideScvfIndex;
                            // there is always only one flip face in an outside element
                            break;
                        }
                    }

                }
            }
        }

        // building the geometries has finished
        std::cout << "Initializing of the grid finite volume geometry took " << timer.elapsed() << " seconds." << std::endl;

        // Initialize the grid interaction volume seeds
        timer.reset();
        ivIndexSets_.update(*this, std::move(dualIdSet));
        std::cout << "Initializing of the grid interaction volume index sets took " << timer.elapsed() << " seconds." << std::endl;

        // build the connectivity map for an efficient assembly
        timer.reset();
        connectivityMap_.update(*this);
        std::cout << "Initializing of the connectivity map took " << timer.elapsed() << " seconds." << std::endl;
    }

    //! Get a sub control volume with a global scv index
    const SubControlVolume& scv(GridIndexType scvIdx) const { return scvs_[scvIdx]; }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const { return scvfs_[scvfIdx]; }

    //! Returns the connectivity map of which dofs
    //! have derivatives with respect to a given dof.
    const ConnectivityMap& connectivityMap() const { return connectivityMap_; }

    //! Returns the grid interaction volume seeds class.
    const GridIVIndexSets& gridInteractionVolumeIndexSets() const { return ivIndexSets_; }

    //! Get the scvf on the same face but from the other side
    //! Note that e.g. the normals might be different in the case of surface grids
    const SubControlVolumeFace& flipScvf(GridIndexType scvfIdx, unsigned int outsideScvfIdx = 0) const
    { return scvfs_[flipScvfIndices_[scvfIdx][outsideScvfIdx]]; }

    //! Get the sub control volume face indices of an scv by global index
    const std::vector<GridIndexType>& scvfIndicesOfScv(GridIndexType scvIdx) const
    { return scvfIndicesOfScv_[scvIdx]; }

private:
    // connectivity map for efficient assembly
    ConnectivityMap connectivityMap_;

    // the finite volume grid geometries
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;

    // containers storing the global data
    std::vector<std::vector<GridIndexType>> scvfIndicesOfScv_;
    std::vector<bool> secondaryInteractionVolumeVertices_;
    std::size_t numBoundaryScvf_;

    // needed for embedded surface and network grids (dim < dimWorld)
    std::vector<std::vector<GridIndexType>> flipScvfIndices_;

    // The grid interaction volume index set
    GridIVIndexSets ivIndexSets_;

    // Indicator function on where to use the secondary IVs
    SecondaryIvIndicatorType secondaryIvIndicator_;
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The finite volume geometry (scvs and scvfs) for cell-centered mpfa models on a grid view
 *        This builds up the sub control volumes and sub control volume faces
 * \note For caching disabled we store only some essential index maps to build up local systems on-demand in
 *       the corresponding FVElementGeometry
 */
template<class TypeTag>
class CCMpfaFVGridGeometry<TypeTag, false> : public BaseFVGridGeometry<TypeTag>
{
    using ParentType = BaseFVGridGeometry<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementMapper = typename GET_PROP_TYPE(TypeTag, ElementMapper);

    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using PrimaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);

    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Intersection = typename GridView::Intersection;
    using CoordScalar = typename GridView::ctype;
    using GridIndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = typename PrimaryInteractionVolume::Traits::LocalIndexType;

    using GridIVIndexSets = CCMpfaGridInteractionVolumeIndexSets<TypeTag>;
    using ConnectivityMap = CCMpfaConnectivityMap<TypeTag>;

    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

public:
    using SecondaryIvIndicator = std::function<bool(const Element&, const Intersection&, bool)>;

    //! Constructor without indicator function for secondary interaction volumes
    //! Per default, we use the secondary IVs at branching points & boundaries
    CCMpfaFVGridGeometry(const GridView& gridView)
    : ParentType(gridView)
    , secondaryIvIndicator_([] (const Element& e, const Intersection& is, bool isBranching)
                               { return is.boundary() || isBranching; } )
    {}

    //! Constructor with user-defined indicator function for secondary interaction volumes
    CCMpfaFVGridGeometry(const GridView& gridView, const SecondaryIvIndicator& indicator)
    : ParentType(gridView)
    , secondaryIvIndicator_(indicator)
    {}

    //! the element mapper is the dofMapper
    //! this is convenience to have better chance to have the same main files for box/tpfa/mpfa...
    const ElementMapper& dofMapper() const
    { return this->elementMapper(); }

    //! Returns the total number of sub control volumes.
    std::size_t numScv() const { return numScvs_; }

    //! Returns the total number of sub control volume faces.
    std::size_t numScvf() const { return numScvf_; }

    //! Returns the number of scvfs on the domain boundary.
    std::size_t numBoundaryScvf() const { return numBoundaryScvf_; }

    //! Returns the total number of degrees of freedom.
    std::size_t numDofs() const { return this->gridView().size(0); }

    //! Gets an element from a global element index.
    Element element(GridIndexType eIdx) const { return this->elementMap()[eIdx]; }

    //! Gets an element from a sub control volume contained in it.
    Element element(const SubControlVolume& scv) const { return this->elementMap()[scv.elementIndex()]; }

    //! Returns true if secondary interaction volumes are used around a given vertex (index).
    //! This specialization is enabled if the use of secondary interaction volumes is active.
    template<bool useSecondary = MpfaHelper::considerSecondaryIVs(), std::enable_if_t<useSecondary, int> = 0>
    bool vertexUsesSecondaryInteractionVolume(GridIndexType vIdxGlobal) const
    { return secondaryInteractionVolumeVertices_[vIdxGlobal]; }

    //! Returns true if secondary interaction volumes are used around a given vertex (index).
    //! If the use of secondary interaction volumes is disabled, this can be evaluated at compile time.
    template<bool useSecondary = MpfaHelper::considerSecondaryIVs(), std::enable_if_t<!useSecondary, int> = 0>
    constexpr bool vertexUsesSecondaryInteractionVolume(GridIndexType vIdxGlobal) const { return false; }

    //! Returns true if a given vertex lies on a processor boundary inside a ghost element.
    bool isGhostVertex(const Vertex& v) const { return isGhostVertex_[this->vertexMapper().index(v)]; }

    //! Returns true if the vertex (index) lies on a processor boundary inside a ghost element.
    bool isGhostVertex(GridIndexType vIdxGlobal) const { return isGhostVertex_[vIdxGlobal]; }

    //! Updates all finite volume geometries of the grid. Has to be called again after grid adaption.
    void update()
    {
        ParentType::update();

        // stop the time required for the update
        Dune::Timer timer;

        // resize containers
        numScvs_ = numDofs();
        scvfIndicesOfScv_.resize(numScvs_);
        neighborVolVarIndices_.resize(numScvs_);

        // Some methods require to use a second type of interaction volume, e.g.
        // around vertices on the boundary or branching points (surface grids)
        const auto numVert = this->gridView().size(dim);
        secondaryInteractionVolumeVertices_.resize(numVert, false);

        // find vertices on processor boundaries HERE!!
        isGhostVertex_ = MpfaHelper::findGhostVertices(this->gridView(), this->vertexMapper());

        // instantiate the dual grid index set (to be used for construction of interaction volumes)
        CCMpfaDualGridIndexSet<GridIndexType, LocalIndexType, dim> dualIdSet(this->gridView());

        // Build the SCVs and SCV faces
        numScvf_ = 0;
        numBoundaryScvf_ = 0;
        for (const auto& element : elements(this->gridView()))
        {
            const auto eIdx = this->elementMapper().index(element);

            // the element geometry
            auto elementGeometry = element.geometry();

            // the element-wise index sets for finite volume geometry
            const auto numLocalFaces = MpfaHelper::getNumLocalScvfs(elementGeometry.type());
            std::vector<GridIndexType> scvfsIndexSet;
            std::vector< std::vector<GridIndexType> > neighborVolVarIndexSet;
            scvfsIndexSet.reserve(numLocalFaces);
            neighborVolVarIndexSet.reserve(numLocalFaces);

            // for network grids there might be multiple intersections with the same geometryInInside
            // we indentify those by the indexInInside for now (assumes conforming grids at branching facets)
            std::vector<std::vector<GridIndexType>> outsideIndices;
            if (dim < dimWorld)
            {
                outsideIndices.resize(element.subEntities(1));
                for (const auto& intersection : intersections(this->gridView(), element))
                {
                    if (intersection.neighbor())
                    {
                        const auto& outside = intersection.outside();
                        auto nIdx = this->elementMapper().index(outside);
                        outsideIndices[intersection.indexInInside()].push_back(nIdx);
                    }
                }
            }

            // construct the sub control volume faces
            for (const auto& is : intersections(this->gridView(), element))
            {
                const auto indexInInside = is.indexInInside();
                const bool boundary = is.boundary();
                const bool neighbor = is.neighbor();

                // for surface grids, skip the rest if handled already
                if (dim < dimWorld && neighbor && outsideIndices[indexInInside].empty())
                    continue;

                // if outside level > inside level, use the outside element in the following
                const bool useNeighbor = neighbor && is.outside().level() > element.level();
                const auto& e = useNeighbor ? is.outside() : element;
                const auto indexInElement = useNeighbor ? is.indexInOutside() : indexInInside;
                const auto eg = e.geometry();
                const auto& refElement = ReferenceElements::general(eg.type());

                // evaluate if vertices on this intersection use primary/secondary IVs
                const bool isBranchingPoint = dim < dimWorld ? outsideIndices[indexInInside].size() > 1 : false;
                const bool usesSecondaryIV = secondaryIvIndicator_(element, is, isBranchingPoint);

                // make the scv faces belonging to each corner of the intersection
                for (std::size_t c = 0; c < is.geometry().corners(); ++c)
                {
                    // get the global vertex index the scv face is connected to
                    const auto vIdxLocal = refElement.subEntity(indexInElement, 1, c, dim);
                    const auto vIdxGlobal = this->vertexMapper().subIndex(e, vIdxLocal, dim);

                    // do not build scvfs connected to a processor boundary
                    if (isGhostVertex_[vIdxGlobal])
                        continue;

                    // if this vertex is tagged to use the secondary IVs, store info
                    if (usesSecondaryIV)
                        secondaryInteractionVolumeVertices_[vIdxGlobal] = true;

                    // make the scv face (for non-boundary scvfs on network grids, use precalculated outside indices)
                    const auto& outsideScvIndices = [&] ()
                                                    {
                                                        if (!boundary)
                                                        {
                                                            return dim == dimWorld ?
                                                                   std::vector<GridIndexType>({this->elementMapper().index(is.outside())}) :
                                                                   outsideIndices[indexInInside];
                                                        }
                                                        else
                                                        {
                                                            // compute boundary scv idx and increment counter
                                                            const GridIndexType bIdx = numScvs_ + numBoundaryScvf_;
                                                            numBoundaryScvf_++;
                                                            return std::vector<GridIndexType>(1, bIdx);
                                                        }
                                                    } ();

                    // insert the scvf data into the dual grid index set
                    dualIdSet[vIdxGlobal].insert(boundary, numScvf_, eIdx, outsideScvIndices);

                    // store information on the scv face
                    scvfsIndexSet.push_back(numScvf_++);
                    neighborVolVarIndexSet.emplace_back(std::move(outsideScvIndices));
                }

                // for network grids, clear outside indices to not make a second scvf on that facet
                if (dim < dimWorld)
                    outsideIndices[indexInInside].clear();
            }

            // store the sets of indices in the data container
            scvfIndicesOfScv_[eIdx] = scvfsIndexSet;
            neighborVolVarIndices_[eIdx] = neighborVolVarIndexSet;
        }

        // building the geometries has finished
        std::cout << "Initializing of the grid finite volume geometry took " << timer.elapsed() << " seconds." << std::endl;

        // Initialize the grid interaction volume seeds
        timer.reset();
        ivIndexSets_.update(*this, std::move(dualIdSet));
        std::cout << "Initializing of the grid interaction volume index sets took " << timer.elapsed() << " seconds." << std::endl;

        // build the connectivity map for an effecient assembly
        timer.reset();
        connectivityMap_.update(*this);
        std::cout << "Initializing of the connectivity map took " << timer.elapsed() << " seconds." << std::endl;
    }

    //! Returns the sub control volume face indices of an scv by global index.
    const std::vector<GridIndexType>& scvfIndicesOfScv(GridIndexType scvIdx) const
    { return scvfIndicesOfScv_[scvIdx]; }

    //! Returns the neighboring vol var indices for each scvf contained in an scv.
    const std::vector< std::vector<GridIndexType> >& neighborVolVarIndices(GridIndexType scvIdx) const
    { return neighborVolVarIndices_[scvIdx]; }

    //! Returns the connectivity map of which dofs
    //! have derivatives with respect to a given dof.
    const ConnectivityMap& connectivityMap() const { return connectivityMap_; }

    //! Returns the grid interaction volume seeds class.
    const GridIVIndexSets& gridInteractionVolumeIndexSets() const { return ivIndexSets_; }

private:
    // connectivity map for efficient assembly
    ConnectivityMap connectivityMap_;

    // containers storing the global data
    std::vector<std::vector<GridIndexType>> scvfIndicesOfScv_;
    std::vector< std::vector< std::vector<GridIndexType> > > neighborVolVarIndices_;
    std::vector<bool> secondaryInteractionVolumeVertices_;
    std::vector<bool> isGhostVertex_;
    std::size_t numScvs_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;

    // The grid interaction volume index set
    GridIVIndexSets ivIndexSets_;

    // Indicator function on where to use the secondary IVs
    SecondaryIvIndicator secondaryIvIndicator_;
};

} // end namespace Dumux

#endif
