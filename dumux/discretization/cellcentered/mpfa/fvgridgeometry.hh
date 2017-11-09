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
 * \brief Base class for the finite volume geometry vector for mpfa models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_FV_GRID_GEOMETRY_BASE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_FV_GRID_GEOMETRY_BASE_HH

#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/elementmap.hh>
#include <dumux/common/boundingboxtree.hh>

#include <dumux/discretization/basefvgridgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/fvelementgeometry.hh>
#include <dumux/discretization/cellcentered/mpfa/connectivitymap.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>
#include <dumux/discretization/cellcentered/mpfa/gridinteractionvolumeindexsets.hh>

namespace Dumux
{
/*!
 * \ingroup ImplicitModel
 * \brief Base class for the finite volume geometry vector for mpfa models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 */
template<class TypeTag, bool EnableFVElementGeometryCache>
class CCMpfaFVGridGeometry
{};

// specialization in case the finite volume grid geometries are stored
template<class TypeTag>
class CCMpfaFVGridGeometry<TypeTag, true> : public BaseFVGridGeometry<TypeTag>
{
    using ParentType = BaseFVGridGeometry<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridIVIndexSets = CCMpfaGridInteractionVolumeIndexSets<TypeTag>;
    using ConnectivityMap = CCMpfaConnectivityMap<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Intersection = typename GridView::Intersection;
    using CoordScalar = typename GridView::ctype;
    using IndexType = typename GridView::IndexSet::IndexType;

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

    //! The local class needs access to the scv, scvfs and the fv element geometry
    //! as they are globally cached
    friend typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

public:
    //! Constructor
    CCMpfaFVGridGeometry(const GridView gridView)
    : ParentType(gridView), elementMap_(gridView)
    {}

    /*!
     * \brief Returns the total number of sub control volumes.
     */
    std::size_t numScv() const
    { return scvs_.size(); }

    /*!
     * \brief Returns the total number of sub control volume faces.
     */
    std::size_t numScvf() const
    { return scvfs_.size(); }

    /*!
     * \brief Returns the number of scvfs on the domain boundary.
     */
    std::size_t numBoundaryScvf() const
    { return numBoundaryScvf_; }


    /*!
     * \brief Returns the total number of degrees of freedom.
     */
    std::size_t numDofs() const
    { return this->gridView().size(0); }

    /*!
     * \brief Gets an element from a sub control volume contained in it.
     */
    Element element(const SubControlVolume& scv) const
    { return elementMap_.element(scv.elementIndex()); }

    /*!
     * \brief Gets an element from a global element index.
     */
    Element element(IndexType eIdx) const
    { return elementMap_.element(eIdx); }

    /*!
     * \brief Returns true if primary interaction volumes are used around a given vertex,
     *        false otherwise.
     */
    bool vertexUsesPrimaryInteractionVolume(const Vertex& v) const
    { return primaryInteractionVolumeVertices_[this->vertexMapper().index(v)]; }

    /*!
     * \brief Returns true if primary interaction volumes are used around a given vertex index,
     *        false otherwise.
     */
    bool vertexUsesPrimaryInteractionVolume(IndexType vIdxGlobal) const
    { return primaryInteractionVolumeVertices_[vIdxGlobal]; }

    /*!
     * \brief Returns if primary interaction volumes are used around a given vertex.
     */
    bool vertexUsesSecondaryInteractionVolume(const Vertex& v) const
    { return secondaryInteractionVolumeVertices_[this->vertexMapper().index(v)]; }

    /*!
     * \brief Returns true if primary interaction volumes are used around a given vertex index,
     *        false otherwise.
     */
    bool vertexUsesSecondaryInteractionVolume(IndexType vIdxGlobal) const
    { return secondaryInteractionVolumeVertices_[vIdxGlobal]; }

    /*!
     * \brief Updates all finite volume geometries of the grid. This has to be called again
     *        after grid adaption. A function can be passed to this method specifying where the secondary
     *        interaction volume type should be used. Per default we use it on boundaries
     *        and branching points.
     *
     * \param useSecondaryIV Indicator function at which vertices to apply the secondary interaction volume
     */
    void update( std::function<bool(const Element&, const Intersection&, bool)> useSecondaryIV
                              = [] (const Element& e, const Intersection& is, bool isBranching)
                                   { return is.boundary() || isBranching; } )
    {
        // stop the time required for the update
        Dune::Timer timer;

        // clear containers (necessary after grid refinement)
        scvs_.clear();
        scvfs_.clear();
        scvfIndicesOfScv_.clear();
        elementMap_.clear();

        // determine the number of geometric entities
        const auto numVert = this->gridView().size(dim);
        const auto numScvs = numDofs();
        std::size_t numScvf = MpfaHelper::estimateNumScvf(this->gridView());

        // resize containers
        scvs_.resize(numScvs);
        scvfs_.reserve(numScvf);
        scvfIndicesOfScv_.resize(numScvs);
        elementMap_.resize(numScvs);

        // Some methods require to use a second type of interaction volume, e.g.
        // around vertices on the boundary or branching points (surface grids)
        primaryInteractionVolumeVertices_.resize(numVert, true);
        secondaryInteractionVolumeVertices_.resize(numVert, false);

        // find vertices on processor boundaries
        const auto isGhostVertex = MpfaHelper::findGhostVertices(this->gridView(), this->vertexMapper());

        // instantiate the dual grid index set (to be used for construction of interaction volumes)
        CCMpfaDualGridIndexSet<TypeTag> dualIdSet(this->gridView());

        // Build the SCVs and SCV faces
        IndexType scvfIdx = 0;
        numBoundaryScvf_ = 0;
        for (const auto& element : elements(this->gridView()))
        {
            const auto eIdx = this->elementMapper().index(element);

            // the element geometry
            auto elementGeometry = element.geometry();

            // fill the element map with seeds
            elementMap_[eIdx] = element.seed();

            // The local scvf index set
            std::vector<IndexType> scvfIndexSet;
            scvfIndexSet.reserve(MpfaHelper::getNumLocalScvfs(elementGeometry.type()));

            // for network grids there might be multiple intersections with the same geometryInInside
            // we indentify those by the indexInInside for now (assumes conforming grids at branching facets)
            std::vector<std::vector<IndexType>> outsideIndices;
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
                const bool usesSecondaryIV = useSecondaryIV(element, is, isBranchingPoint);

                // make the scv faces belonging to each corner of the intersection
                for (int c = 0; c < numCorners; ++c)
                {
                    // get the global vertex index the scv face is connected to
                    const auto vIdxLocal = refElement.subEntity(indexInElement, 1, c, dim);
                    const auto vIdxGlobal = this->vertexMapper().subIndex(e, vIdxLocal, dim);

                    // do not build scvfs connected to a processor boundary
                    if (isGhostVertex[vIdxGlobal])
                        continue;

                    // if this vertex is tagged to use the secondary IVs, store info
                    if (usesSecondaryIV)
                    {
                        primaryInteractionVolumeVertices_[vIdxGlobal] = false;
                        secondaryInteractionVolumeVertices_[vIdxGlobal] = true;
                    }

                    // the quadrature point to be used on the scvfs
                    static const Scalar q = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Mpfa, Q);

                    // make the scv face (for non-boundary scvfs on network grids, use precalculated outside indices)
                    const auto& outsideScvIndices = [&] ()
                                                    {
                                                        if (!boundary)
                                                        {
                                                            return dim == dimWorld ?
                                                                   std::vector<IndexType>({this->elementMapper().index(is.outside())}) :
                                                                   outsideIndices[indexInInside];
                                                        }
                                                        else
                                                        {
                                                            // compute boundary scv idx and increment counter
                                                            const IndexType bIdx = numScvs + numBoundaryScvf_;
                                                            numBoundaryScvf_++;
                                                            return std::vector<IndexType>(1, bIdx);
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
                if (dim < dimWorld) outsideIndices[indexInInside].clear();
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
            for (auto&& scvf : scvfs_)
            {
                if (scvf.boundary())
                    continue;

                const auto numOutsideScvs = scvf.numOutsideScvs();
                const auto vIdxGlobal = scvf.vertexIndex();
                const auto insideScvIdx = scvf.insideScvIdx();

                flipScvfIndices_[scvf.index()].resize(numOutsideScvs);
                for (unsigned int i = 0; i < numOutsideScvs; ++i)
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

        // build the connectivity map for an effecient assembly
        timer.reset();
        connectivityMap_.update(*this);
        std::cout << "Initializing of the connectivity map took " << timer.elapsed() << " seconds." << std::endl;
    }

    /*!
     * \brief Returns a sub control volume for a given global scv index.
     */
    const SubControlVolume& scv(IndexType scvIdx) const
    { return scvs_[scvIdx]; }

    /*!
     * \brief Returns a sub control volume face for a given global scvf index.
     */
    const SubControlVolumeFace& scvf(IndexType scvfIdx) const
    { return scvfs_[scvfIdx]; }

    /*!
     * \brief Returns the "flipped" scvf, i.e. the correspongin scvf in an outside element.
     *        The second argument is optional and only comes into play on network/surface grids.
     */
    const SubControlVolumeFace& flipScvf(IndexType scvfIdx, unsigned int outsideScvfIdx = 0) const
    { return scvfs_[flipScvfIndices_[scvfIdx][outsideScvfIdx]]; }

    /*!
     * \brief Returns the sub control volume face indices of an scv by global index.
     */
    const std::vector<IndexType>& scvfIndicesOfScv(IndexType scvIdx) const
    { return scvfIndicesOfScv_[scvIdx]; }

    /*!
     * \brief Returns the connectivity map of which dofs have derivatives with respect
     *        to a given dof.
     */
    const ConnectivityMap& connectivityMap() const
    { return connectivityMap_; }

    /*!
     * \brief Returns the grit interaction volume seeds class.
     */
    const GridIVIndexSets& gridInteractionVolumeIndexSets() const
    { return ivIndexSets_; }

private:
    // mappers
    Dumux::ElementMap<GridView> elementMap_;
    ConnectivityMap connectivityMap_;

    // the finite volume grid geometries
    std::vector<SubControlVolume> scvs_;
    std::vector<SubControlVolumeFace> scvfs_;

    // containers storing the global data
    std::vector<std::vector<IndexType>> scvfIndicesOfScv_;
    std::vector<bool> primaryInteractionVolumeVertices_;
    std::vector<bool> secondaryInteractionVolumeVertices_;
    std::size_t numBoundaryScvf_;

    // needed for embedded surface and network grids (dim < dimWorld)
    std::vector<std::vector<IndexType>> flipScvfIndices_;

    // The grid interaction volume index set
    GridIVIndexSets ivIndexSets_;
};

// specialization in case the FVElementGeometries are not stored
template<class TypeTag>
class CCMpfaFVGridGeometry<TypeTag, false> : public BaseFVGridGeometry<TypeTag>
{
    using ParentType = BaseFVGridGeometry<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using GridIVIndexSets = CCMpfaGridInteractionVolumeIndexSets<TypeTag>;
    using ConnectivityMap = CCMpfaConnectivityMap<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Intersection = typename GridView::Intersection;
    using CoordScalar = typename GridView::ctype;
    using IndexType = typename GridView::IndexSet::IndexType;

    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dim>;

    //! Todo is this necessary?
    friend typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

public:
    //! Constructor
    CCMpfaFVGridGeometry(const GridView gridView)
    : ParentType(gridView), elementMap_(gridView)
    {}

    /*!
     * \brief Returns the total number of sub control volumes.
     */
    std::size_t numScv() const
    { return numScvs_; }

    /*!
     * \brief Returns the total number of sub control volume faces.
     */
    std::size_t numScvf() const
    { return numScvf_; }

    /*!
     * \brief Returns the number of scvfs on the domain boundary.
     */
    std::size_t numBoundaryScvf() const
    { return numBoundaryScvf_; }

    /*!
     * \brief Returns the total number of degrees of freedom.
     */
    std::size_t numDofs() const
    { return this->gridView().size(0); }

    /*!
     * \brief Gets an element from a sub control volume contained in it.
     */
    Element element(const SubControlVolume& scv) const
    { return elementMap_.element(scv.elementIndex()); }

    /*!
     * \brief Gets an element from a global element index.
     */
    Element element(IndexType eIdx) const
    { return elementMap_.element(eIdx); }

    /*!
     * \brief Returns true if primary interaction volumes are used around a given vertex,
     *        false otherwise.
     */
    bool vertexUsesPrimaryInteractionVolume(const Vertex& v) const
    { return primaryInteractionVolumeVertices_[this->vertexMapper().index(v)]; }

    /*!
     * \brief Returns true if primary interaction volumes are used around a given vertex index,
     *        false otherwise.
     */
    bool vertexUsesPrimaryInteractionVolume(IndexType vIdxGlobal) const
    { return primaryInteractionVolumeVertices_[vIdxGlobal]; }

    /*!
     * \brief Returns if primary interaction volumes are used around a given vertex.
     */
    bool vertexUsesSecondaryInteractionVolume(const Vertex& v) const
    { return secondaryInteractionVolumeVertices_[this->vertexMapper().index(v)]; }

    /*!
     * \brief Returns true if primary interaction volumes are used around a given vertex index,
     *        false otherwise.
     */
    bool vertexUsesSecondaryInteractionVolume(IndexType vIdxGlobal) const
    { return secondaryInteractionVolumeVertices_[vIdxGlobal]; }

    /*!
     * \brief Returns true if a given vertex lies on a processor boundary inside a ghost element.
     */
    bool isGhostVertex(const Vertex& v) const
    { return isGhostVertex_[this->vertexMapper().index(v)]; }

    /*!
     * \brief Returns true if the vertex corresponding to a given vertex index lies on a
     *        processor boundary inside a ghost element.
     */
    bool isGhostVertex(IndexType vIdxGlobal) const
    { return isGhostVertex_[vIdxGlobal]; }

    /*!
     * \brief Updates all finite volume geometries of the grid. This has to be called again
     *        after grid adaption. A function can be passed to this method specifying where the secondary
     *        interaction volume type should be used. Per default we use it on boundaries
     *        and branching points.
     *
     * \param useSecondaryIV Indicator function at which vertices to apply the secondary interaction volume
     */
    void update( std::function<bool(const Element&, const Intersection&, bool)> useSecondaryIV
                              = [] (const Element& e, const Intersection& is, bool isBranching)
                                   { return is.boundary() || isBranching; } )
    {
        // stop the time required for the update
        Dune::Timer timer;

        // clear containers (necessary after grid refinement)
        scvfIndicesOfScv_.clear();
        neighborVolVarIndices_.clear();
        elementMap_.clear();

        // resize containers
        numScvs_ = numDofs();
        scvfIndicesOfScv_.resize(numScvs_);
        neighborVolVarIndices_.resize(numScvs_);
        elementMap_.resize(numScvs_);

        // Some methods require to use a second type of interaction volume, e.g.
        // around vertices on the boundary or branching points (surface grids)
        const auto numVert = this->gridView().size(dim);
        primaryInteractionVolumeVertices_.resize(numVert, true);
        secondaryInteractionVolumeVertices_.resize(numVert, false);

        // find vertices on processor boundaries HERE!!
        isGhostVertex_ = MpfaHelper::findGhostVertices(this->gridView(), this->vertexMapper());

        // instantiate the dual grid index set (to be used for construction of interaction volumes)
        CCMpfaDualGridIndexSet<TypeTag> dualIdSet(this->gridView());

        // Build the SCVs and SCV faces
        numScvf_ = 0;
        numBoundaryScvf_ = 0;
        for (const auto& element : elements(this->gridView()))
        {
            const auto eIdx = this->elementMapper().index(element);

            // the element geometry
            auto elementGeometry = element.geometry();

            // fill the element map with seeds
            elementMap_[eIdx] = element.seed();

            // the element-wise index sets for finite volume geometry
            const auto numLocalFaces = MpfaHelper::getNumLocalScvfs(elementGeometry.type());
            std::vector<IndexType> scvfsIndexSet;
            std::vector< std::vector<IndexType> > neighborVolVarIndexSet;
            scvfsIndexSet.reserve(numLocalFaces);
            neighborVolVarIndexSet.reserve(numLocalFaces);

            // for network grids there might be multiple intersections with the same geometryInInside
            // we indentify those by the indexInInside for now (assumes conforming grids at branching facets)
            std::vector<std::vector<IndexType>> outsideIndices;
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

            unsigned int localFaceIdx = 0;
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
                const bool usesSecondaryIV = useSecondaryIV(element, is, isBranchingPoint);

                // make the scv faces belonging to each corner of the intersection
                for (int c = 0; c < is.geometry().corners(); ++c)
                {
                    // get the global vertex index the scv face is connected to
                    const auto vIdxLocal = refElement.subEntity(indexInElement, 1, c, dim);
                    const auto vIdxGlobal = this->vertexMapper().subIndex(e, vIdxLocal, dim);

                    // do not build scvfs connected to a processor boundary
                    if (isGhostVertex_[vIdxGlobal])
                        continue;

                    // if this vertex is tagged to use the secondary IVs, store info
                    if (usesSecondaryIV)
                    {
                        primaryInteractionVolumeVertices_[vIdxGlobal] = false;
                        secondaryInteractionVolumeVertices_[vIdxGlobal] = true;
                    }

                    // make the scv face (for non-boundary scvfs on network grids, use precalculated outside indices)
                    const auto& outsideScvIndices = [&] ()
                                                    {
                                                        if (!boundary)
                                                        {
                                                            return dim == dimWorld ?
                                                                   std::vector<IndexType>({this->elementMapper().index(is.outside())}) :
                                                                   outsideIndices[indexInInside];
                                                        }
                                                        else
                                                        {
                                                            // compute boundary scv idx and increment counter
                                                            const IndexType bIdx = numScvs_ + numBoundaryScvf_;
                                                            numBoundaryScvf_++;
                                                            return std::vector<IndexType>(1, bIdx);
                                                        }
                                                    } ();

                    // insert the scvf data into the dual grid index set
                    dualIdSet[vIdxGlobal].insert(boundary, numScvf_, eIdx, outsideScvIndices);

                    // store information on the scv face
                    scvfsIndexSet.push_back(numScvf_++);
                    neighborVolVarIndexSet.emplace_back(std::move(outsideScvIndices));

                    // increment counter
                    localFaceIdx++;
                }

                // for network grids, clear outside indices to not make a second scvf on that facet
                if (dim < dimWorld) outsideIndices[indexInInside].clear();
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

    /*!
     * \brief Returns the sub control volume face indices of an scv by global index.
     */
    const std::vector<IndexType>& scvfIndicesOfScv(IndexType scvIdx) const
    { return scvfIndicesOfScv_[scvIdx]; }

    /*!
     * \brief Returns the neighboring vol var indices for each scvf contained in an scv.
     */
    const std::vector< std::vector<IndexType> >& neighborVolVarIndices(IndexType scvIdx) const
    { return neighborVolVarIndices_[scvIdx]; }

    /*!
     * \brief Returns the connectivity map of which dofs have derivatives with respect
     *        to a given dof.
     */
    const ConnectivityMap& connectivityMap() const
    { return connectivityMap_; }

    /*!
     * \brief Returns the grit interaction volume seeds class.
     */
    const GridIVIndexSets& gridInteractionVolumeIndexSets() const
    { return ivIndexSets_; }

private:
    // mappers
    Dumux::ElementMap<GridView> elementMap_;
    ConnectivityMap connectivityMap_;

    // containers storing the global data
    std::vector<std::vector<IndexType>> scvfIndicesOfScv_;
    std::vector< std::vector< std::vector<IndexType> > > neighborVolVarIndices_;
    std::vector<bool> primaryInteractionVolumeVertices_;
    std::vector<bool> secondaryInteractionVolumeVertices_;
    std::vector<bool> isGhostVertex_;
    std::size_t numScvs_;
    std::size_t numScvf_;
    std::size_t numBoundaryScvf_;

    // The grid interaction volume index set
    GridIVIndexSets ivIndexSets_;
};

} // end namespace

#endif
