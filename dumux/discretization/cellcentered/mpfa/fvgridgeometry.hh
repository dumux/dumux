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
 * \ingroup CCMpfaDiscretization
 * \brief The finite volume geometry (scvs and scvfs) for cell-centered mpfa models on a grid view
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element of the grid partition.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FV_GRID_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FV_GRID_GEOMETRY_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/basegridgeometry.hh>
#include <dumux/discretization/checkoverlapsize.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The finite volume geometry (scvs and scvfs) for cell-centered mpfa models on a grid view
 *        This builds up the sub control volumes and sub control volume faces
 * \note This class is specialized for versions with and without caching the fv geometries on the grid view
 *
 * \tparam GridView the grid view
 * \tparam GridIvIs the interaction volume index sets
 * \tparam Traits traits class
 * \tparam enableCache
 */
template<class GridView, class Traits, bool enableCache>
class CCMpfaFVGridGeometry;

//! check the overlap size for parallel computations
template<class GridView>
void checkOverlapSizeCCMpfa(const GridView& gridView)
{
    // Check if the overlap size is what we expect
    if (!CheckOverlapSize<DiscretizationMethod::ccmpfa>::isValid(gridView))
        DUNE_THROW(Dune::InvalidStateException, "The ccmpfa discretization method needs at least an overlap of 1 for parallel computations. "
                                                 << " Set the parameter \"Grid.Overlap\" in the input file.");
}

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The finite volume geometry (scvs and scvfs) for cell-centered mpfa models on a grid view
 *        This builds up the sub control volumes and sub control volume faces
 * \note For caching enabled we store the fv geometries for the whole grid view which is memory intensive but faster
 */
template<class GV, class Traits>
class CCMpfaFVGridGeometry<GV, Traits, true>
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = CCMpfaFVGridGeometry<GV, Traits, true>;
    using ParentType = BaseGridGeometry<GV, Traits>;

    static constexpr int dim = GV::dimension;
    static constexpr int dimWorld = GV::dimensionworld;

    using Element = typename GV::template Codim<0>::Entity;
    using Vertex = typename GV::template Codim<dim>::Entity;
    using Intersection = typename GV::Intersection;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using CoordScalar = typename GV::ctype;

    using ScvfOutsideGridIndexStorage = typename Traits::SubControlVolumeFace::Traits::OutsideGridIndexStorage;

public:
    //! export the flip scvf index set type
    using FlipScvfIndexSet = std::vector<ScvfOutsideGridIndexStorage>;
    //! export the grid interaction volume index set type
    using GridIVIndexSets = typename Traits::template GridIvIndexSets<ThisType>;
    //! export the type to be used for indicators where to use the secondary ivs
    using SecondaryIvIndicatorType = std::function<bool(const Element&, const Intersection&, bool)>;

    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, true>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;
    //! export the connectivity map type
    using ConnectivityMap = typename Traits::template ConnectivityMap<ThisType>;
    //! export dof mapper type
    using DofMapper = typename Traits::ElementMapper;
    //! export the grid view type
    using GridView = GV;
    //! export the mpfa helper type
    using MpfaHelper = typename Traits::template MpfaHelper<ThisType>;

    //! export the discretization method this geometry belongs to
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::ccmpfa;

    //! The maximum admissible stencil size (used for static memory allocation during assembly)
    static constexpr int maxElementStencilSize = Traits::maxElementStencilSize;

    //! State if only a single type is used for interaction volumes
    static constexpr bool hasSingleInteractionVolumeType = !MpfaHelper::considerSecondaryIVs();

    //! Constructor without indicator function for secondary interaction volumes
    //! Per default, we use the secondary IVs at branching points & boundaries
    CCMpfaFVGridGeometry(const GridView& gridView)
    : ParentType(gridView)
    , secondaryIvIndicator_([] (const Element& e, const Intersection& is, bool isBranching)
                               { return is.boundary() || isBranching; } )
    {
        checkOverlapSizeCCMpfa(gridView);
    }

    //! Constructor with user-defined indicator function for secondary interaction volumes
    CCMpfaFVGridGeometry(const GridView& gridView, const SecondaryIvIndicatorType& indicator)
    : ParentType(gridView)
    , secondaryIvIndicator_(indicator)
    {
        checkOverlapSizeCCMpfa(gridView);
    }

    //! the element mapper is the dofMapper
    //! this is convenience to have better chance to have the same main files for box/tpfa/mpfa...
    const DofMapper& dofMapper() const
    { return this->elementMapper(); }

    //! The total number of sub control volumes
    std::size_t numScv() const
    { return scvs_.size(); }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    { return scvfs_.size(); }

    //! The total number of boundary sub control volume faces
    std::size_t numBoundaryScvf() const
    { return numBoundaryScvf_; }

    //! The total number of degrees of freedom
    std::size_t numDofs() const
    { return this->gridView().size(0); }

    //! Returns true if secondary interaction volumes are used around a given vertex (index).
    //! This specialization is enabled if the use of secondary interaction volumes is active.
    template<bool useSecondary = !hasSingleInteractionVolumeType, std::enable_if_t<useSecondary, bool> = 0>
    bool vertexUsesSecondaryInteractionVolume(GridIndexType vIdxGlobal) const
    { return secondaryInteractionVolumeVertices_[vIdxGlobal]; }

    //! Returns true if secondary interaction volumes are used around a given vertex (index).
    //! If the use of secondary interaction volumes is disabled, this can be evaluated at compile time.
    template<bool useSecondary = !hasSingleInteractionVolumeType, std::enable_if_t<!useSecondary, bool> = 0>
    constexpr bool vertexUsesSecondaryInteractionVolume(GridIndexType vIdxGlobal) const
    { return false; }

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
        hasBoundaryScvf_.resize(numScvs, false);

        // Some methods require to use a second type of interaction volume, e.g.
        // around vertices on the boundary or branching points (surface grids)
        secondaryInteractionVolumeVertices_.resize(numVert, false);

        // find vertices on processor boundaries
        const auto isGhostVertex = MpfaHelper::findGhostVertices(this->gridView(), this->vertexMapper());

        // instantiate the dual grid index set (to be used for construction of interaction volumes)
        typename GridIVIndexSets::DualGridIndexSet dualIdSet(this->gridView());

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

            // for network grids there might be multiple intersection with the same geometryInInside
            // we indentify those by the indexInInside for now (assumes conforming grids at branching facets)
            std::vector<ScvfOutsideGridIndexStorage> outsideIndices;
            if (dim < dimWorld)
            {
                outsideIndices.resize(element.subEntities(1));
                for (const auto& intersection : intersections(this->gridView(), element))
                {
                    if (intersection.neighbor())
                    {
                        const auto nIdx = this->elementMapper().index( intersection.outside() );
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

                if (boundary)
                    hasBoundaryScvf_[eIdx] = true;

                // for surface grids, skip the rest if handled already
                if (dim < dimWorld && neighbor && outsideIndices[indexInInside].empty())
                    continue;

                // if outside level > inside level, use the outside element in the following
                const bool useNeighbor = neighbor && is.outside().level() > element.level();
                const auto& e = useNeighbor ? is.outside() : element;
                const auto indexInElement = useNeighbor ? is.indexInOutside() : indexInInside;
                const auto eg = e.geometry();
                const auto refElement = referenceElement(eg);

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
                    static const auto q = getParam<CoordScalar>("MPFA.Q");

                    // make the scv face (for non-boundary scvfs on network grids, use precalculated outside indices)
                    const auto& outsideScvIndices = [&] ()
                                                    {
                                                        if (!boundary)
                                                            return dim == dimWorld ?
                                                                   ScvfOutsideGridIndexStorage({this->elementMapper().index(is.outside())}) :
                                                                   outsideIndices[indexInInside];
                                                        else
                                                            return ScvfOutsideGridIndexStorage({GridIndexType(numScvs) + numBoundaryScvf_++});
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

        // building the geometries has finished
        std::cout << "Initializing of the grid finite volume geometry took " << timer.elapsed() << " seconds." << std::endl;

        // Initialize the grid interaction volume index sets
        timer.reset();
        ivIndexSets_.update(*this, std::move(dualIdSet));
        std::cout << "Initializing of the grid interaction volume index sets took " << timer.elapsed() << " seconds." << std::endl;

        // build the connectivity map for an efficient assembly
        timer.reset();
        connectivityMap_.update(*this);
        std::cout << "Initializing of the connectivity map took " << timer.elapsed() << " seconds." << std::endl;
    }

    //! Returns instance of the mpfa helper type
    MpfaHelper mpfaHelper() const
    { return MpfaHelper(); }

    //! Get a sub control volume with a global scv index
    const SubControlVolume& scv(GridIndexType scvIdx) const
    { return scvs_[scvIdx]; }

    //! Get a sub control volume face with a global scvf index
    const SubControlVolumeFace& scvf(GridIndexType scvfIdx) const
    { return scvfs_[scvfIdx]; }

    //! Returns the connectivity map of which dofs
    //! have derivatives with respect to a given dof.
    const ConnectivityMap& connectivityMap() const
    { return connectivityMap_; }

    //! Returns the grid interaction volume index set class.
    const GridIVIndexSets& gridInteractionVolumeIndexSets() const
    { return ivIndexSets_; }

    //! Get the sub control volume face indices of an scv by global index
    const std::vector<GridIndexType>& scvfIndicesOfScv(GridIndexType scvIdx) const
    { return scvfIndicesOfScv_[scvIdx]; }

    //! Returns the flip scvf index set
    const FlipScvfIndexSet& flipScvfIndexSet() const
    { return flipScvfIndices_; }

    //! Get the scvf on the same face but from the other side
    //! Note that e.g. the normals might be different in the case of surface grids
    const SubControlVolumeFace& flipScvf(GridIndexType scvfIdx, unsigned int outsideScvfIdx = 0) const
    { return scvfs_[flipScvfIndices_[scvfIdx][outsideScvfIdx]]; }

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf(GridIndexType eIdx) const
    { return hasBoundaryScvf_[eIdx]; }

private:
    // connectivity map for efficient assembly
    ConnectivityMap connectivityMap_;

    // the finite volume grid geometries
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;

    // containers storing the global data
    std::vector<std::vector<GridIndexType>> scvfIndicesOfScv_;
    std::vector<bool> secondaryInteractionVolumeVertices_;
    GridIndexType numBoundaryScvf_;
    std::vector<bool> hasBoundaryScvf_;

    // needed for embedded surface and network grids (dim < dimWorld)
    FlipScvfIndexSet flipScvfIndices_;

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
template<class GV, class Traits>
class CCMpfaFVGridGeometry<GV, Traits, false>
: public BaseGridGeometry<GV, Traits>
{
    using ThisType = CCMpfaFVGridGeometry<GV, Traits, false>;
    using ParentType = BaseGridGeometry<GV, Traits>;

    static constexpr int dim = GV::dimension;
    static constexpr int dimWorld = GV::dimensionworld;

    using Element = typename GV::template Codim<0>::Entity;
    using Vertex = typename GV::template Codim<dim>::Entity;
    using Intersection = typename GV::Intersection;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using CoordScalar = typename GV::ctype;

    using ScvfOutsideGridIndexStorage = typename Traits::SubControlVolumeFace::Traits::OutsideGridIndexStorage;

public:
    //! export the flip scvf index set type
    using FlipScvfIndexSet = std::vector<ScvfOutsideGridIndexStorage>;
    //! export the grid interaction volume index set type
    using GridIVIndexSets = typename Traits::template GridIvIndexSets<ThisType>;
    //! export the type to be used for indicators where to use the secondary ivs
    using SecondaryIvIndicatorType = std::function<bool(const Element&, const Intersection&, bool)>;

    //! export the type of the fv element geometry (the local view type)
    using LocalView = typename Traits::template LocalView<ThisType, false>;
    //! export the type of sub control volume
    using SubControlVolume = typename Traits::SubControlVolume;
    //! export the type of sub control volume
    using SubControlVolumeFace = typename Traits::SubControlVolumeFace;
    //! export the type of extrusion
    using Extrusion = Extrusion_t<Traits>;
    //! export the connectivity map type
    using ConnectivityMap = typename Traits::template ConnectivityMap<ThisType>;
    //! export dof mapper type
    using DofMapper = typename Traits::ElementMapper;
    //! export the grid view type
    using GridView = GV;
    //! export the mpfa helper type
    using MpfaHelper = typename Traits::template MpfaHelper<ThisType>;

    //! export the discretization method this geometry belongs to
    static constexpr DiscretizationMethod discMethod = DiscretizationMethod::ccmpfa;

    //! The maximum admissible stencil size (used for static memory allocation during assembly)
    static constexpr int maxElementStencilSize = Traits::maxElementStencilSize;

    //! State if only a single type is used for interaction volumes
    static constexpr bool hasSingleInteractionVolumeType = !MpfaHelper::considerSecondaryIVs();

    //! Constructor without indicator function for secondary interaction volumes
    //! Per default, we use the secondary IVs at branching points & boundaries
    CCMpfaFVGridGeometry(const GridView& gridView)
    : ParentType(gridView)
    , secondaryIvIndicator_([] (const Element& e, const Intersection& is, bool isBranching)
                               { return is.boundary() || isBranching; } )
    {
        checkOverlapSizeCCMpfa(gridView);
    }

    //! Constructor with user-defined indicator function for secondary interaction volumes
    CCMpfaFVGridGeometry(const GridView& gridView, const SecondaryIvIndicatorType& indicator)
    : ParentType(gridView)
    , secondaryIvIndicator_(indicator)
    {
        checkOverlapSizeCCMpfa(gridView);
    }

    //! the element mapper is the dofMapper
    //! this is convenience to have better chance to have the same main files for box/tpfa/mpfa...
    const DofMapper& dofMapper() const
    { return this->elementMapper(); }

    //! Returns the total number of sub control volumes.
    std::size_t numScv() const
    { return numScvs_; }

    //! Returns the total number of sub control volume faces.
    std::size_t numScvf() const
    { return numScvf_; }

    //! Returns the number of scvfs on the domain boundary.
    std::size_t numBoundaryScvf() const
    { return numBoundaryScvf_; }

    //! Returns the total number of degrees of freedom.
    std::size_t numDofs() const
    { return this->gridView().size(0); }

    //! Returns true if secondary interaction volumes are used around a given vertex (index).
    //! This specialization is enabled if the use of secondary interaction volumes is active.
    template<bool useSecondary = !hasSingleInteractionVolumeType, std::enable_if_t<useSecondary, bool> = 0>
    bool vertexUsesSecondaryInteractionVolume(GridIndexType vIdxGlobal) const
    { return secondaryInteractionVolumeVertices_[vIdxGlobal]; }

    //! Returns true if secondary interaction volumes are used around a given vertex (index).
    //! If the use of secondary interaction volumes is disabled, this can be evaluated at compile time.
    template<bool useSecondary = !hasSingleInteractionVolumeType, std::enable_if_t<!useSecondary, bool> = 0>
    constexpr bool vertexUsesSecondaryInteractionVolume(GridIndexType vIdxGlobal) const
    { return false; }

    //! Returns true if a given vertex lies on a processor boundary inside a ghost element.
    bool isGhostVertex(const Vertex& v) const
    { return isGhostVertex_[this->vertexMapper().index(v)]; }

    //! Returns true if the vertex (index) lies on a processor boundary inside a ghost element.
    bool isGhostVertex(GridIndexType vIdxGlobal) const
    { return isGhostVertex_[vIdxGlobal]; }

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
        typename GridIVIndexSets::DualGridIndexSet dualIdSet(this->gridView());

        // keep track of boundary scvfs and scvf vertex indices in order to set up flip scvf index set
        const auto maxNumScvfs = numScvs_*LocalView::maxNumElementScvfs;
        std::vector<bool> scvfIsOnBoundary;
        std::vector<GridIndexType> scvfVertexIndex;
        scvfIsOnBoundary.reserve(maxNumScvfs);
        scvfVertexIndex.reserve(maxNumScvfs);

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
            std::vector<ScvfOutsideGridIndexStorage> neighborVolVarIndexSet;
            scvfsIndexSet.reserve(numLocalFaces);
            neighborVolVarIndexSet.reserve(numLocalFaces);

            // for network grids there might be multiple intersections with the same geometryInInside
            // we indentify those by the indexInInside for now (assumes conforming grids at branching facets)
            std::vector<ScvfOutsideGridIndexStorage> outsideIndices;
            if (dim < dimWorld)
            {
                outsideIndices.resize(element.subEntities(1));
                for (const auto& intersection : intersections(this->gridView(), element))
                {
                    if (intersection.neighbor())
                    {
                        auto nIdx = this->elementMapper().index( intersection.outside() );
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
                const auto refElement = referenceElement(eg);

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
                                                            return dim == dimWorld ?
                                                                   ScvfOutsideGridIndexStorage({this->elementMapper().index(is.outside())}) :
                                                                   outsideIndices[indexInInside];
                                                        else
                                                            return ScvfOutsideGridIndexStorage({GridIndexType(numScvs_) + numBoundaryScvf_++});
                                                    } ();

                    // insert the scvf data into the dual grid index set
                    dualIdSet[vIdxGlobal].insert(numScvf_, eIdx, boundary);

                    // store information on the scv face
                    scvfsIndexSet.push_back(numScvf_++);
                    scvfIsOnBoundary.push_back(boundary);
                    scvfVertexIndex.push_back(vIdxGlobal);
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

        // Make the flip scvf index set
        flipScvfIndices_.resize(numScvf_);
        for (std::size_t scvIdx = 0; scvIdx < numScvs_; ++scvIdx)
        {
            const auto& scvfIndices = scvfIndicesOfScv_[scvIdx];
            for (unsigned int i = 0; i < scvfIndices.size(); ++i)
            {
                // boundary scvf have no flip scvfs
                if (scvfIsOnBoundary[ scvfIndices[i] ])
                    continue;

                const auto scvfIdx = scvfIndices[i];
                const auto vIdxGlobal = scvfVertexIndex[scvfIdx];
                const auto numOutsideScvs = neighborVolVarIndices_[scvIdx][i].size();

                flipScvfIndices_[scvfIdx].resize(numOutsideScvs);
                for (unsigned int j = 0; j < numOutsideScvs; ++j)
                {
                    const auto outsideScvIdx = neighborVolVarIndices_[scvIdx][i][j];
                    const auto& outsideScvfIndices = scvfIndicesOfScv_[outsideScvIdx];
                    for (unsigned int k = 0; k < outsideScvfIndices.size(); ++k)
                    {
                        const auto outsideScvfIndex = outsideScvfIndices[k];
                        const auto outsideScvfVertexIndex = scvfVertexIndex[outsideScvfIndex];
                        const auto& outsideScvfNeighborIndices = neighborVolVarIndices_[outsideScvIdx][k];
                        if (outsideScvfVertexIndex == vIdxGlobal &&
                            MpfaHelper::vectorContainsValue(outsideScvfNeighborIndices, scvIdx))
                        {
                            flipScvfIndices_[scvfIdx][j] = outsideScvfIndex;
                            // there is always only one flip face in an outside element
                            break;
                        }
                    }
                }
            }
        }

        // building the geometries has finished
        std::cout << "Initializing of the grid finite volume geometry took " << timer.elapsed() << " seconds." << std::endl;

        // Initialize the grid interaction volume index sets
        timer.reset();
        ivIndexSets_.update(*this, std::move(dualIdSet));
        std::cout << "Initializing of the grid interaction volume index sets took " << timer.elapsed() << " seconds." << std::endl;

        // build the connectivity map for an efficient assembly
        timer.reset();
        connectivityMap_.update(*this);
        std::cout << "Initializing of the connectivity map took " << timer.elapsed() << " seconds." << std::endl;
    }

    //! Returns instance of the mpfa helper type
    MpfaHelper mpfaHelper() const
    { return MpfaHelper(); }

    //! Returns the sub control volume face indices of an scv by global index.
    const std::vector<GridIndexType>& scvfIndicesOfScv(GridIndexType scvIdx) const
    { return scvfIndicesOfScv_[scvIdx]; }

    //! Returns the neighboring vol var indices for each scvf contained in an scv.
    const std::vector<ScvfOutsideGridIndexStorage>& neighborVolVarIndices(GridIndexType scvIdx) const
    { return neighborVolVarIndices_[scvIdx]; }

    //! Get the index scvf on the same face but from the other side
    //! Note that e.g. the normals might be different in the case of surface grids
    const GridIndexType flipScvfIdx(GridIndexType scvfIdx, unsigned int outsideScvfIdx = 0) const
    { return flipScvfIndices_[scvfIdx][outsideScvfIdx]; }

    //! Returns the flip scvf index set
    const FlipScvfIndexSet& flipScvfIndexSet() const
    { return flipScvfIndices_; }

    //! Returns the connectivity map of which dofs
    //! have derivatives with respect to a given dof.
    const ConnectivityMap& connectivityMap() const
    { return connectivityMap_; }

    //! Returns the grid interaction volume seeds class.
    const GridIVIndexSets& gridInteractionVolumeIndexSets() const
    { return ivIndexSets_; }

private:
    // connectivity map for efficient assembly
    ConnectivityMap connectivityMap_;

    // containers storing the global data
    std::vector<std::vector<GridIndexType>> scvfIndicesOfScv_;
    std::vector<std::vector<ScvfOutsideGridIndexStorage>> neighborVolVarIndices_;
    std::vector<bool> secondaryInteractionVolumeVertices_;
    std::vector<bool> isGhostVertex_;
    GridIndexType numScvs_;
    GridIndexType numScvf_;
    GridIndexType numBoundaryScvf_;

    // needed for embedded surface and network grids (dim < dimWorld)
    FlipScvfIndexSet flipScvfIndices_;

    // The grid interaction volume index set
    GridIVIndexSets ivIndexSets_;

    // Indicator function on where to use the secondary IVs
    SecondaryIvIndicatorType secondaryIvIndicator_;
};

} // end namespace Dumux

#endif
