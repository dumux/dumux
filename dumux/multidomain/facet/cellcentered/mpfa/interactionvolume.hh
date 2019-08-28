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
 * \ingroup FacetCoupling
 * \copydoc Dumux::CCMpfaOFacetCouplingInteractionVolume
 */
#ifndef DUMUX_MULTIDOMAIN_FACET_CC_MPFA_O_INTERACTIONVOLUME_HH
#define DUMUX_MULTIDOMAIN_FACET_CC_MPFA_O_INTERACTIONVOLUME_HH

#include <dumux/discretization/cellcentered/mpfa/omethod/interactionvolume.hh>
#include <dumux/discretization/cellcentered/mpfa/omethod/scvgeometryhelper.hh>
#include <dumux/discretization/cellcentered/mpfa/interactionvolumebase.hh>
#include <dumux/discretization/cellcentered/mpfa/localfacedata.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

#include "localassembler.hh"
#include "localsubcontrolentities.hh"

namespace Dumux {

//! Forward declaration of the facet coupling MPFA-O interaction volume
template< class Traits > class CCMpfaOFacetCouplingInteractionVolume;

/*!
 * \ingroup FacetCoupling
 * \brief The default interaction volume traits class for the mpfa-o method
 *        in the context of facet coupling. This uses dynamic types types for
 *        matrices/vectors in order to work on general grids.
 *
 * \tparam NodalIndexSet The type used for the dual grid's nodal index sets
 * \tparam Scalar The Type used for scalar values
 */
template< class NodalIndexSet, class Scalar >
struct CCMpfaOFacetCouplingDefaultInteractionVolumeTraits
: public CCMpfaODefaultInteractionVolumeTraits< NodalIndexSet, Scalar >
{
private:
    //! export the type for the interaction volume index set
    using IVIndexSet = CCMpfaOInteractionVolumeIndexSet< NodalIndexSet >;

    static constexpr int dim = NodalIndexSet::Traits::GridView::dimension;
    static constexpr int dimWorld = NodalIndexSet::Traits::GridView::dimensionworld;

public:
    //! export the type of interaction-volume local scvs
    using LocalScvType = CCMpfaOFacetCouplingInteractionVolumeLocalScv< IVIndexSet, Scalar, dim, dimWorld >;
    //! export the type of interaction-volume local scvfs
    using LocalScvfType = CCMpfaOFacetCouplingInteractionVolumeLocalScvf< IVIndexSet >;

    //! Use the assembler that considers the coupled domain
    template<class Problem, class FVElementGeometry, class ElemVolVars>
    using LocalAssembler = MpfaOFacetCouplingInteractionVolumeAssembler<Problem, FVElementGeometry, ElemVolVars>;
};

/*!
 * \ingroup FacetCoupling
 * \brief Class for the interaction volume of the mpfa-o scheme in the
 *        context of models involving coupling to a lower-dimensional
 *        domain across the element facets.
 */
template< class Traits >
class CCMpfaOFacetCouplingInteractionVolume
: public CCMpfaInteractionVolumeBase< Traits >
{
    using GridView = typename Traits::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using IndexSet = typename Traits::IndexSet;
    using GridIndexType = typename IndexSet::GridIndexType;
    using LocalIndexType = typename IndexSet::LocalIndexType;
    using Stencil = typename IndexSet::NodalGridStencilType;

    using LocalScvType = typename Traits::LocalScvType;
    using LocalScvfType = typename Traits::LocalScvfType;
    using LocalFaceData = typename Traits::LocalFaceData;

public:
    //! Reuse standard o-scheme's Dirichlet Data class
    using DirichletData = typename CCMpfaOInteractionVolume<Traits>::DirichletData;

    //! Define data structure to store which scvfs lie on interior boundaries
    class InteriorBoundaryData
    {
        GridIndexType scvfIdx_;

    public:
        //! Constructor
        InteriorBoundaryData(GridIndexType scvfIdx) : scvfIdx_(scvfIdx) {}

        //! Return corresponding scvf index
        GridIndexType scvfIndex() const { return scvfIdx_; }
    };

    //! publicly state the mpfa-scheme this interaction volume is associated with
    static constexpr MpfaMethods MpfaMethod = MpfaMethods::oMethod;

    //! Sets up the local scope for a given iv index set
    template< class Problem, class FVElementGeometry >
    void bind(const IndexSet& indexSet,
              const Problem& problem,
              const FVElementGeometry& fvGeometry)
    {
        // for the o-scheme, the stencil is equal to the scv
        // index set of the dual grid's nodal index set
        stencil_ = &indexSet.nodalIndexSet().gridScvIndices();

        // find out how many facet elements appear in this iv
        std::size_t numFacetElems = 0;
        std::size_t numOutsideFaces = 0;
        std::vector<bool> isOnInteriorBoundary(indexSet.numFaces(), false);
        for (LocalIndexType fIdx = 0; fIdx < indexSet.numFaces(); ++fIdx)
        {
            const auto& scvf = fvGeometry.scvf(indexSet.gridScvfIndex(fIdx));
            const auto element = fvGeometry.gridGeometry().element(scvf.insideScvIdx());
            if (problem.couplingManager().isOnInteriorBoundary(element, scvf))
            {
                numFacetElems++;
                numOutsideFaces += scvf.numOutsideScvs();
                isOnInteriorBoundary[fIdx] = true;
                interiorBoundaryData_.emplace_back( scvf.index() );
            }
        }

        // number of interaction-volume-local scvs(=node-local for o-scheme) and scvfs
        numFaces_ = indexSet.numFaces() + numOutsideFaces;
        const auto numLocalScvs = indexSet.numScvs();
        const auto numGlobalScvfs = indexSet.nodalIndexSet().numScvfs();

        // reserve memory for local entities
        elements_.clear();      elements_.reserve(numLocalScvs);
        scvs_.clear();          scvs_.reserve(numLocalScvs);
        scvfs_.clear();         scvfs_.reserve(numFaces_);
        localFaceData_.clear(); localFaceData_.reserve(numGlobalScvfs);
        dirichletData_.clear(); dirichletData_.reserve(numFaces_);

        // keep track of the number of unknowns etc
        numUnknowns_ = 0;
        numKnowns_ = numLocalScvs + numFacetElems;

        // index map from grid scvf index to local scvf index
        std::unordered_map<GridIndexType, LocalIndexType> scvfIndexMap;

        // set up objects related to sub-control volume faces
        LocalIndexType facetElementCounter = 0;
        for (LocalIndexType faceIdxLocal = 0; faceIdxLocal < indexSet.numFaces(); ++faceIdxLocal)
        {
            const auto gridScvfIdx = indexSet.gridScvfIndex(faceIdxLocal);
            const auto& flipScvfIdxSet = fvGeometry.gridGeometry().flipScvfIndexSet()[gridScvfIdx];
            const auto& scvf = fvGeometry.scvf(gridScvfIdx);
            const auto element = fvGeometry.gridGeometry().element(scvf.insideScvIdx());

            // the neighboring scvs in local indices (order: 0 - inside scv, 1..n - outside scvs)
            const auto& neighborScvIndicesLocal = indexSet.neighboringLocalScvIndices(faceIdxLocal);
            const auto numNeighborScvs = neighborScvIndicesLocal.size();

            // the ív-local scvf index of the face about to be created
            const auto curLocalScvfIdx = scvfs_.size();
            scvfIndexMap[gridScvfIdx] = curLocalScvfIdx;
            localFaceData_.emplace_back(curLocalScvfIdx, neighborScvIndicesLocal[0], scvf.index());

            // on interior boundaries, create local scvfs for inside AND all outside scvfs
            if (isOnInteriorBoundary[faceIdxLocal])
            {
                const LocalIndexType facetLocalDofIdx = numLocalScvs + facetElementCounter++;
                const bool isDirichlet = problem.interiorBoundaryTypes(element, scvf).hasOnlyDirichlet();

                // create local scvf
                if (isDirichlet)
                    scvfs_.emplace_back(scvf, neighborScvIndicesLocal, facetLocalDofIdx, /*isDirichlet*/true, facetLocalDofIdx);
                else
                    scvfs_.emplace_back(scvf, neighborScvIndicesLocal, numUnknowns_++, /*isDirichlet*/false, facetLocalDofIdx);

                // create "outside" local scvfs
                for (LocalIndexType i = 1; i < numNeighborScvs; ++i)
                {
                    const auto outsideGridScvfIdx = flipScvfIdxSet[i-1];
                    const auto& flipScvf = fvGeometry.scvf(outsideGridScvfIdx);
                    const auto& outsideFlipScvfIdxSet = fvGeometry.gridGeometry().flipScvfIndexSet()[outsideGridScvfIdx];

                    // rearrange the neighbor scv index vector corresponding to this scvfs flip scvf index set
                    using std::swap;
                    auto outsideNeighborScvIdxSet = neighborScvIndicesLocal;
                    outsideNeighborScvIdxSet[0] = outsideNeighborScvIdxSet[i];
                    for (LocalIndexType j = 0; j < outsideFlipScvfIdxSet.size(); ++j)
                    {
                        const auto flipScvfIdx = outsideFlipScvfIdxSet[j];
                        auto it = std::find(flipScvfIdxSet.begin(), flipScvfIdxSet.end(), flipScvfIdx);

                        // if we found the index, use corresponding local scv index
                        if (it != flipScvfIdxSet.end())
                            outsideNeighborScvIdxSet[j+1] = neighborScvIndicesLocal[std::distance(flipScvfIdxSet.begin(), it)+1];

                        // otherwise this must be the "inside" scvf again
                        else
                        {
                            assert(flipScvfIdx == gridScvfIdx);
                            outsideNeighborScvIdxSet[j+1] = neighborScvIndicesLocal[0];
                        }
                    }

                    scvfIndexMap[outsideGridScvfIdx] = curLocalScvfIdx+i;
                    localFaceData_.emplace_back(curLocalScvfIdx+i, outsideNeighborScvIdxSet[0], flipScvf.index());
                    if (isDirichlet)
                        scvfs_.emplace_back(flipScvf, outsideNeighborScvIdxSet, facetLocalDofIdx, /*isDirichlet*/true, facetLocalDofIdx);
                    else
                        scvfs_.emplace_back(flipScvf, outsideNeighborScvIdxSet, numUnknowns_++, /*isDirichlet*/false, facetLocalDofIdx);
                }
            }

            // otherwise crate boundary scvf ...
            else if (scvf.boundary())
            {
                if (problem.boundaryTypes(element, scvf).hasOnlyDirichlet())
                {
                    scvfs_.emplace_back(scvf, neighborScvIndicesLocal, numKnowns_++, /*isDirichlet*/true);
                    dirichletData_.emplace_back(scvf.outsideScvIdx());
                }
                else
                    scvfs_.emplace_back(scvf, neighborScvIndicesLocal, numUnknowns_++, /*isDirichlet*/false);
            }

            // ... or interior scvf
            else
            {
                scvfs_.emplace_back(scvf, neighborScvIndicesLocal, numUnknowns_++, /*isDirichlet*/false);

                // add local face data objects for the outside faces
                for (LocalIndexType i = 1; i < numNeighborScvs; ++i)
                {
                    // loop over scvfs in outside scv until we find the one coinciding with current scvf
                    const auto outsideLocalScvIdx = neighborScvIndicesLocal[i];
                    const auto& flipScvfIndex = fvGeometry.gridGeometry().flipScvfIndexSet()[scvf.index()][i-1];
                    const auto& flipScvf = fvGeometry.scvf(flipScvfIndex);
                    scvfIndexMap[flipScvfIndex] = curLocalScvfIdx;
                    localFaceData_.emplace_back(curLocalScvfIdx,    // iv-local scvf idx
                                                outsideLocalScvIdx, // iv-local scv index
                                                i-1,                // scvf-local index in outside faces
                                                flipScvf.index());  // global scvf index
                }
            }
        }

        // set up stuff related to sub-control volumes
        for (LocalIndexType scvIdxLocal = 0; scvIdxLocal < numLocalScvs; scvIdxLocal++)
        {
            elements_.emplace_back(fvGeometry.gridGeometry().element( stencil()[scvIdxLocal] ));
            scvs_.emplace_back(fvGeometry.gridGeometry().mpfaHelper(),
                               fvGeometry,
                               fvGeometry.scv( stencil()[scvIdxLocal] ),
                               scvIdxLocal,
                               indexSet,
                               scvfIndexMap);
        }
    }

    //! returns the number of primary scvfs of this interaction volume
    std::size_t numFaces() const
    { return numFaces_; }

    //! returns the number of intermediate unknowns within this interaction volume
    std::size_t numUnknowns() const
    { return numUnknowns_; }

    //! returns the number of (in this context) known solution values within this interaction volume
    std::size_t numKnowns() const
    { return numKnowns_; }

    //! returns the number of scvs embedded in this interaction volume
    std::size_t numScvs() const
    { return scvs_.size(); }

    //! returns the cell-stencil of this interaction volume
    const Stencil& stencil() const
    { return *stencil_; }

    //! returns the grid element corresponding to a given iv-local scv idx
    const Element& element(LocalIndexType ivLocalScvIdx) const
    { return elements_[ivLocalScvIdx]; }

    //! returns the local scvf entity corresponding to a given iv-local scvf idx
    const LocalScvfType& localScvf(LocalIndexType ivLocalScvfIdx) const
    { return scvfs_[ivLocalScvfIdx]; }

    //! returns the local scv entity corresponding to a given iv-local scv idx
    const LocalScvType& localScv(LocalIndexType ivLocalScvIdx) const
    { return scvs_[ivLocalScvIdx]; }

    //! returns a reference to the container with the local face data
    const std::vector<LocalFaceData>& localFaceData() const
    { return localFaceData_; }

    //! returns a reference to the information container on Dirichlet BCs within this iv
    const std::vector<DirichletData>& dirichletData() const
    { return dirichletData_; }

    //! returns a reference to the data container on interior boundaries
    const std::vector<InteriorBoundaryData>& interiorBoundaryData() const
    { return interiorBoundaryData_; }

    //! returns the geometry of the i-th local scv
    template< class FVElementGeometry >
    auto getScvGeometry(LocalIndexType ivLocalScvIdx, const FVElementGeometry& fvGeometry) const
    { return CCMpfaOScvGeometryHelper<LocalScvType>::computeScvGeometry(ivLocalScvIdx, *this, fvGeometry); }

    //! returns the number of interaction volumes living around a vertex
    template< class NI >
    static constexpr std::size_t numIVAtVertex(const NI& nodalIndexSet)
    { return 1; }

    //! adds the iv index sets living around a vertex to a given container
    //! and stores the the corresponding index in a map for each scvf
    template< class IvIndexSetContainer,
              class ScvfIndexMap,
              class NodalIndexSet,
              class FlipScvfIndexSet >
    static void addIVIndexSets(IvIndexSetContainer& ivIndexSetContainer,
                               ScvfIndexMap& scvfIndexMap,
                               const NodalIndexSet& nodalIndexSet,
                               const FlipScvfIndexSet& flipScvfIndexSet)
    {
        // reuse the function of the standard mpfa-o interaction volume
        CCMpfaOInteractionVolume<Traits>::addIVIndexSets(ivIndexSetContainer,
                                                         scvfIndexMap,
                                                         nodalIndexSet,
                                                         flipScvfIndexSet);
    }

private:
    // pointer to cell stencil (in iv index set)
    const Stencil* stencil_;

    // Variables defining the local scope
    std::vector<Element> elements_;
    std::vector<LocalScvType> scvs_;
    std::vector<LocalScvfType> scvfs_;
    std::vector<LocalFaceData> localFaceData_;
    std::vector<DirichletData> dirichletData_;
    std::vector<InteriorBoundaryData> interiorBoundaryData_;

    // sizes involved in the local system equations
    std::size_t numFaces_;
    std::size_t numUnknowns_;
    std::size_t numKnowns_;
};

} // end namespace

#endif
