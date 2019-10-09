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
 * \brief Class for the interaction volume of the mpfa-o scheme to be used when the
 *        sizes are known at compile time, e.g. for non-boundary interaction volumes
 *        on structured grids.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_STATIC_INTERACTIONVOLUME_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_STATIC_INTERACTIONVOLUME_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/math.hh>

#include <dumux/discretization/cellcentered/mpfa/interactionvolumebase.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>
#include <dumux/discretization/cellcentered/mpfa/localfacedata.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

#include "localassembler.hh"
#include "localsubcontrolentities.hh"
#include "interactionvolumeindexset.hh"
#include "interactionvolume.hh"
#include "scvgeometryhelper.hh"

namespace Dumux {

//! Forward declaration of the o-method's static interaction volume
template< class Traits > class CCMpfaOStaticInteractionVolume;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief The default interaction volume traits class for the mpfa-o method
 *        with known size of the interaction volumes at compile time. It uses
 *        statically sized containers for the iv-local data structures and static
 *        matrices and vectors.
 *
 * \tparam NI The type used for the dual grid's nodal index sets
 * \tparam S The Type used for scalar values
 * \tparam C The number of sub-control volumes (cells) in the ivs
 * \tparam F The number of sub-control volume faces in the ivs
 */
template< class NI, class S, int C, int F >
struct CCMpfaODefaultStaticInteractionVolumeTraits
{
private:
    using GridIndexType = typename NI::GridIndexType;
    using LocalIndexType = typename NI::LocalIndexType;

    static constexpr int dim = NI::Traits::GridView::dimension;
    static constexpr int dimWorld = NI::Traits::GridView::dimensionworld;

    using DimVector = Dune::FieldVector<S, dim>;
    using FaceOmegas = Dune::ReservedVector<DimVector, 2>;

    //! Matrix/Vector traits to be used by the data handle
    struct MVTraits
    {
        using OmegaStorage = std::array< FaceOmegas, F >;

        using AMatrix = Dune::FieldMatrix< S, F, F >;
        using BMatrix = Dune::FieldMatrix< S, F, C >;
        using CMatrix = Dune::FieldMatrix< S, F, F >;
        using DMatrix = Dune::FieldMatrix< S, F, C >;
        using TMatrix = Dune::FieldMatrix< S, F, C >;
        using CellVector = Dune::FieldVector< S, C >;
        using FaceVector = Dune::FieldVector< S, F >;
    };

public:
    //! export the type of grid view
    using GridView = typename NI::Traits::GridView;
    //! export the type for the interaction volume index set
    using IndexSet = CCMpfaOInteractionVolumeIndexSet< NI >;
    //! export the type of interaction-volume local scvs
    using LocalScvType = CCMpfaOInteractionVolumeLocalScv< IndexSet, S, dim, dimWorld >;
    //! export the type of interaction-volume local scvfs
    using LocalScvfType = CCMpfaOInteractionVolumeLocalScvf< IndexSet >;
    //! export the type of used for the iv-local face data
    using LocalFaceData = InteractionVolumeLocalFaceData<GridIndexType, LocalIndexType>;
    //! export the matrix/vector traits to be used by the iv
    using MatVecTraits = MVTraits;
    //! export the number of scvs in the interaction volumes
    static constexpr int numScvs = C;
    //! export the number of scvfs in the interaction volumes
    static constexpr int numScvfs = F;

    //! the type of assembler used for the o-method's iv-local eq systems
    template<class Problem, class FVElementGeometry, class ElemVolVars>
    using LocalAssembler = MpfaOInteractionVolumeAssembler<Problem, FVElementGeometry, ElemVolVars>;
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Class for the interaction volume of the mpfa-o method.
 *        This implementation creates static objects of the local geometries
 *        and can only be used for interior interaction volumes with a static
 *        size known at compile time. This size has to match the sizes of the
 *        vectors/matrices defined in the traits class in case static types are used.
 *
 * \tparam Traits The type traits class to be used
 */
template< class Traits >
class CCMpfaOStaticInteractionVolume
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

    static constexpr int numScvf = Traits::numScvfs;
    static constexpr int numScv = Traits::numScvs;

public:
    //! This does not work on surface grids
    static_assert(int(GridView::dimension)==int(GridView::dimensionworld), "static iv does not work on surface grids");

    //! export the standard o-methods dirichlet data
    //! this is fake Dirichlet data as this implementation cannot be used on boundaries
    struct DirichletData
    {
        GridIndexType volVarIndex() const
        { DUNE_THROW(Dune::InvalidStateException, "Static interaction volume cannot be used on bounaries!"); }
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
        assert(indexSet.numScvs() == numScv);
        stencil_ = &indexSet.nodalIndexSet().gridScvIndices();

        // set up stuff related to sub-control volumes
        for (LocalIndexType scvIdxLocal = 0; scvIdxLocal < numScv; scvIdxLocal++)
        {
            elements_[scvIdxLocal] = fvGeometry.gridGeometry().element( stencil()[scvIdxLocal] );
            scvs_[scvIdxLocal] = LocalScvType(fvGeometry.gridGeometry().mpfaHelper(),
                                              fvGeometry,
                                              fvGeometry.scv( stencil()[scvIdxLocal] ),
                                              scvIdxLocal,
                                              indexSet);
        }

        // set up quantitites related to sub-control volume faces
        for (LocalIndexType faceIdxLocal = 0; faceIdxLocal < numScvf; ++faceIdxLocal)
        {
            const auto& scvf = fvGeometry.scvf(indexSet.gridScvfIndex(faceIdxLocal));
            assert(!scvf.boundary());

            // the neighboring scvs in local indices (order: 0 - inside scv, 1..n - outside scvs)
            const auto& neighborScvIndicesLocal = indexSet.neighboringLocalScvIndices(faceIdxLocal);

            // create iv-local scvf objects
            scvfs_[faceIdxLocal] = LocalScvfType(scvf, neighborScvIndicesLocal, faceIdxLocal, /*isDirichlet*/false);
            localFaceData_[faceIdxLocal*2] = LocalFaceData(faceIdxLocal, neighborScvIndicesLocal[0], scvf.index());

            // add local face data objects for the outside face
            const auto outsideLocalScvIdx = neighborScvIndicesLocal[1];
            for (int coord = 0; coord < GridView::dimension; ++coord)
            {
                if (indexSet.localScvfIndex(outsideLocalScvIdx, coord) == faceIdxLocal)
                {
                    const auto globalScvfIdx = indexSet.nodalIndexSet().gridScvfIndex(outsideLocalScvIdx, coord);
                    const auto& flipScvf = fvGeometry.scvf(globalScvfIdx);
                    localFaceData_[faceIdxLocal*2+1] = LocalFaceData(faceIdxLocal,       // iv-local scvf idx
                                                                     outsideLocalScvIdx, // iv-local scv index
                                                                     0,                  // scvf-local index in outside faces
                                                                     flipScvf.index());  // global scvf index
                    break; // go to next outside face
                }
            }

            // make sure we found it
            assert(localFaceData_[faceIdxLocal*2+1].ivLocalInsideScvIndex() == outsideLocalScvIdx);
        }
    }

    //! returns the number of primary scvfs of this interaction volume
    static constexpr std::size_t numFaces()
    { return numScvf; }

    //! returns the number of intermediate unknowns within this interaction volume
    static constexpr std::size_t numUnknowns()
    { return numScvf; }

    //! returns the number of (in this context) known solution values within this interaction volume
    static constexpr std::size_t numKnowns()
    { return numScv; }

    //! returns the number of scvs embedded in this interaction volume
    static constexpr std::size_t numScvs()
    { return numScv; }

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
    const std::array<LocalFaceData, numScvf*2>& localFaceData() const
    { return localFaceData_; }

    //! Returns a reference to the information container on Dirichlet BCs within this iv.
    //! Here, we return an empty container as this implementation cannot be used on boundaries.
    const std::array<DirichletData, 0>& dirichletData() const
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
        // reuse the standard o-method's implementation of this
        CCMpfaOInteractionVolume<Traits>::addIVIndexSets(ivIndexSetContainer,
                                                         scvfIndexMap,
                                                         nodalIndexSet,
                                                         flipScvfIndexSet);
    }

private:
    // pointer to cell stencil (in iv index set)
    const Stencil * stencil_ = nullptr;

    // Variables defining the local scope
    std::array<Element, numScv> elements_;
    std::array<LocalScvType, numScv> scvs_;
    std::array<LocalScvfType, numScvf> scvfs_;
    std::array<LocalFaceData, numScvf*2> localFaceData_;

    // Dummy dirichlet data container (compatibility with dynamic o-iv)
    std::array<DirichletData, 0> dirichletData_;
};

} // end namespace Dumux

#endif
