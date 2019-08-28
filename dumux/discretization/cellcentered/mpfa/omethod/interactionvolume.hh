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
 * \brief Class for the interaction volume of the mpfa-o scheme.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_INTERACTIONVOLUME_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_INTERACTIONVOLUME_HH

#include <type_traits>

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>

#include <dumux/discretization/cellcentered/mpfa/interactionvolumebase.hh>
#include <dumux/discretization/cellcentered/mpfa/localfacedata.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

#include "localassembler.hh"
#include "localsubcontrolentities.hh"
#include "interactionvolumeindexset.hh"
#include "scvgeometryhelper.hh"

namespace Dumux {

//! Forward declaration of the o-method's interaction volume
template< class Traits > class CCMpfaOInteractionVolume;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The default interaction volume traits class for the mpfa-o method.
 *        This uses dynamic types types for matrices/vectors in order to work
 *        on general grids. For interaction volumes known at compile time use
 *        the static interaction volume implementation.
 *
 * \tparam NodalIndexSet The type used for the dual grid's nodal index sets
 * \tparam Scalar The Type used for scalar values
 */
template< class NodalIndexSet, class Scalar >
struct CCMpfaODefaultInteractionVolumeTraits
{
private:
    using GridIndexType = typename NodalIndexSet::GridIndexType;
    using LocalIndexType = typename NodalIndexSet::LocalIndexType;

    static constexpr int dim = NodalIndexSet::Traits::GridView::dimension;
    static constexpr int dimWorld = NodalIndexSet::Traits::GridView::dimensionworld;

    using DimVector = Dune::FieldVector<Scalar, dim>;
    using FaceOmegas = typename std::conditional< (dim<dimWorld),
                                                  std::vector<DimVector>,
                                                  Dune::ReservedVector<DimVector, 2> >::type;

    //! Matrix/Vector traits to be used by the data handle
    struct MVTraits
    {
        using OmegaStorage = std::vector< FaceOmegas >;

        using AMatrix = Dune::DynamicMatrix< Scalar >;
        using BMatrix = Dune::DynamicMatrix< Scalar >;
        using CMatrix = Dune::DynamicMatrix< Scalar >;
        using DMatrix = Dune::DynamicMatrix< Scalar >;
        using TMatrix = Dune::DynamicMatrix< Scalar >;
        using CellVector = Dune::DynamicVector< Scalar >;
        using FaceVector = Dune::DynamicVector< Scalar >;
    };

public:
    //! export the type of grid view
    using GridView = typename NodalIndexSet::Traits::GridView;
    //! export the type for the interaction volume index set
    using IndexSet = CCMpfaOInteractionVolumeIndexSet< NodalIndexSet >;
    //! export the type of interaction-volume local scvs
    using LocalScvType = CCMpfaOInteractionVolumeLocalScv< IndexSet, Scalar, dim, dimWorld >;
    //! export the type of interaction-volume local scvfs
    using LocalScvfType = CCMpfaOInteractionVolumeLocalScvf< IndexSet >;
    //! export the type of used for the iv-local face data
    using LocalFaceData = InteractionVolumeLocalFaceData<GridIndexType, LocalIndexType>;
    //! export the matrix/vector traits to be used by the iv
    using MatVecTraits = MVTraits;

    //! the type of assembler used for the o-method's iv-local eq systems
    template<class Problem, class FVElementGeometry, class ElemVolVars>
    using LocalAssembler = MpfaOInteractionVolumeAssembler<Problem, FVElementGeometry, ElemVolVars>;
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Class for the interaction volume of the mpfa-o method.
 *        This implementation creates dynamic objects of the local geometries
 *        and can be used at boundaries and on unstructured grids.
 */
template< class Traits >
class CCMpfaOInteractionVolume
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
    //! Data attached to scvf touching Dirichlet boundaries.
    //! For the default o-scheme, we only store the corresponding vol vars index.
    class DirichletData
    {
        GridIndexType volVarIndex_;
    public:
        //! Constructor
        DirichletData(const GridIndexType index) : volVarIndex_(index) {}

        //! Return corresponding vol var index
        GridIndexType volVarIndex() const { return volVarIndex_; }
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

        // number of interaction-volume-local scvs(=node-local for o-scheme) and scvfs
        numFaces_ = indexSet.numFaces();
        const auto numLocalScvs = indexSet.numScvs();
        const auto numGlobalScvfs = indexSet.nodalIndexSet().numScvfs();

        // reserve memory for local entities
        elements_.clear();      elements_.reserve(numLocalScvs);
        scvs_.clear();          scvs_.reserve(numLocalScvs);
        scvfs_.clear();         scvfs_.reserve(numFaces_);
        localFaceData_.clear(); localFaceData_.reserve(numGlobalScvfs);
        dirichletData_.clear(); dirichletData_.reserve(numFaces_);

        // set up stuff related to sub-control volumes
        for (LocalIndexType scvIdxLocal = 0; scvIdxLocal < numLocalScvs; scvIdxLocal++)
        {
            elements_.emplace_back(fvGeometry.gridGeometry().element( stencil()[scvIdxLocal] ));
            scvs_.emplace_back(fvGeometry.gridGeometry().mpfaHelper(),
                               fvGeometry,
                               fvGeometry.scv( stencil()[scvIdxLocal] ),
                               scvIdxLocal,
                               indexSet);
        }

        // keep track of the number of unknowns etc
        numUnknowns_ = 0;
        numKnowns_ = numLocalScvs;

        // set up quantitites related to sub-control volume faces
        for (LocalIndexType faceIdxLocal = 0; faceIdxLocal < numFaces_; ++faceIdxLocal)
        {
            const auto& scvf = fvGeometry.scvf(indexSet.gridScvfIndex(faceIdxLocal));

            // the neighboring scvs in local indices (order: 0 - inside scv, 1..n - outside scvs)
            const auto& neighborScvIndicesLocal = indexSet.neighboringLocalScvIndices(faceIdxLocal);
            const auto numNeighborScvs = neighborScvIndicesLocal.size();
            localFaceData_.emplace_back(faceIdxLocal, neighborScvIndicesLocal[0], scvf.index());

            // create iv-local scvf object
            if (scvf.boundary())
            {
                if (problem.boundaryTypes(elements_[neighborScvIndicesLocal[0]], scvf).hasOnlyDirichlet())
                {
                    scvfs_.emplace_back(scvf, neighborScvIndicesLocal, numKnowns_++, /*isDirichlet*/true);
                    dirichletData_.emplace_back(scvf.outsideScvIdx());
                }
                else
                    scvfs_.emplace_back(scvf, neighborScvIndicesLocal, numUnknowns_++, /*isDirichlet*/false);
            }
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
                    localFaceData_.emplace_back(faceIdxLocal,       // iv-local scvf idx
                                                outsideLocalScvIdx, // iv-local scv index
                                                i-1,                // scvf-local index in outside faces
                                                flipScvf.index());  // global scvf index
                }
            }
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
        // the global index of the iv index set that is about to be created
        const auto curGlobalIndex = ivIndexSetContainer.size();

        // make the one index set for this node
        ivIndexSetContainer.emplace_back(nodalIndexSet, flipScvfIndexSet);

        // store the index mapping
        for (const auto scvfIdx : nodalIndexSet.gridScvfIndices())
            scvfIndexMap[scvfIdx] = curGlobalIndex;
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

    // sizes involved in the local system equations
    std::size_t numFaces_;
    std::size_t numUnknowns_;
    std::size_t numKnowns_;
};

} // end namespace Dumux

#endif
