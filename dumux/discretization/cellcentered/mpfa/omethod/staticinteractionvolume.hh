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
 * \brief Class for the interaction volume of the mpfa-o scheme to be used when the
 *        sizes are known at compile time, e.g. for non-boundary interaction volumes
 *        on structured grids.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_STATIC_INTERACTIONVOLUME_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_STATIC_INTERACTIONVOLUME_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/matrixvectorhelper.hh>

#include <dumux/discretization/cellcentered/mpfa/interactionvolumedatahandle.hh>
#include <dumux/discretization/cellcentered/mpfa/interactionvolumebase.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>
#include <dumux/discretization/cellcentered/mpfa/localfacedata.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

#include "localsubcontrolentities.hh"
#include "interactionvolumeindexset.hh"

namespace Dumux
{
//! Forward declaration of the o-method's static interaction volume
template< class Traits, int localSize > class CCMpfaOStaticInteractionVolume;

//! Specialization of the default interaction volume traits class for the mpfa-o method
template< class TypeTag, int localSize >
struct CCMpfaODefaultStaticInteractionVolumeTraits
{
private:
    using GridIndexType = typename GET_PROP_TYPE(TypeTag, GridView)::IndexSet::IndexType;
    static constexpr int dim = GET_PROP_TYPE(TypeTag, GridView)::dimension;
    static constexpr int dimWorld = GET_PROP_TYPE(TypeTag, GridView)::dimensionworld;

public:
    //! export the problem type (needed for iv-local assembly)
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    //! export the type of the local view on the finite volume grid geometry
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    //! export the type of the local view on the grid volume variables
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    //! export the type of grid view
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    //! export the type used for scalar values
    using ScalarType = typename GET_PROP_TYPE(TypeTag, Scalar);
    //! export the type of the mpfa helper class
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    //! export the type used for local indices
    using LocalIndexType = std::uint8_t;
    //! export the type for the interaction volume index set
    using IndexSet = CCMpfaOInteractionVolumeIndexSet< DualGridNodalIndexSet<GridIndexType, LocalIndexType, dim> >;
    //! export the type of interaction-volume local scvs
    using LocalScvType = CCMpfaOInteractionVolumeLocalScv< IndexSet, ScalarType, dim, dimWorld >;
    //! export the type of interaction-volume local scvfs
    using LocalScvfType = CCMpfaOInteractionVolumeLocalScvf< IndexSet >;
    //! export the type of used for the iv-local face data
    using LocalFaceData = InteractionVolumeLocalFaceData<GridIndexType, LocalIndexType>;
    //! export the type used for iv-local matrices
    using Matrix = Dune::FieldMatrix< ScalarType, localSize, localSize >;
    //! export the type used for iv-local vectors
    using Vector = Dune::FieldVector< ScalarType, localSize >;
    //! export the type used for the iv-stencils
    using Stencil = std::vector< GridIndexType >;
    //! export the data handle type for this iv
    using DataHandle = InteractionVolumeDataHandle< TypeTag, Matrix, Vector, LocalIndexType >;
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
 * \tparam localSize The size of the local eq system
 *                   This also is the number of local scvs/scvfs
 */
template< class Traits, int localSize >
class CCMpfaOStaticInteractionVolume : public CCMpfaInteractionVolumeBase< CCMpfaOInteractionVolume<Traits>, Traits >
{
    using ThisType = CCMpfaOInteractionVolume< Traits >;
    using ParentType = CCMpfaInteractionVolumeBase< CCMpfaOInteractionVolume<Traits>, Traits >;

    using Scalar = typename Traits::ScalarType;
    using Helper = typename Traits::MpfaHelper;
    using Problem = typename Traits::Problem;
    using FVElementGeometry = typename Traits::FVElementGeometry;

    using GridView = typename Traits::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    using DimVector = Dune::FieldVector<Scalar, dim>;

    using Matrix = typename Traits::Matrix;
    using LocalScvType = typename Traits::LocalScvType;
    using LocalScvfType = typename Traits::LocalScvfType;
    using LocalFaceData = typename Traits::LocalFaceData;

    using IndexSet = typename Traits::IndexSet;
    using GridIndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = typename Traits::LocalIndexType;
    using Stencil = typename Traits::Stencil;

    //! For the o method, the interaction volume stencil can be taken directly
    //! from the nodal index set, which always uses dynamic types to be compatible
    //! on the boundaries and unstructured grids. Thus, we have to make sure that
    //! the type set for the stencils in the traits is castable.
    static_assert( std::is_convertible<Stencil*, typename IndexSet::GridIndexContainer*>::value,
                   "The o-method uses the (dynamic) nodal index set's stencil as the interaction volume stencil. "
                   "Using a different type is not permissive here." );

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

public:
    //! publicly state the mpfa-scheme this interaction volume is associated with
    static constexpr MpfaMethods MpfaMethod = MpfaMethods::oMethod;

    //! Sets up the local scope for a given iv index set
    void setUpLocalScope(const IndexSet& indexSet,
                         const Problem& problem,
                         const FVElementGeometry& fvGeometry)
    {
        // for the o-scheme, the stencil is equal to the scv
        // index set of the dual grid's nodal index set
        stencil_ = &indexSet.nodalIndexSet().globalScvIndices();

        // set up local geometry containers
        const auto& scvIndices = indexSet.globalScvIndices();
        for (LocalIndexType localIdx = 0; localIdx < localSize; localIdx++)
        {
            // scv-related quantities
            scvs_[localIdx] = LocalScvType(Helper(),
                                              fvGeometry,
                                              fvGeometry.scv( scvIndices[localIdx] ),
                                              localIdx,
                                              indexSet);
            elements_[localIdx] = fvGeometry.fvGridGeometry().element( scvIndices[localIdx] );

            // scvf-related quantities
            const auto& scvf = fvGeometry.scvf(indexSet.scvfIdxGlobal(localIdx));

            // this interaction volume implementation does not work on boundaries
            assert(!scvf.boundary() && "static mpfa-o interaction volume cannot be used on boundaries");

            const auto& neighborScvIndicesLocal = indexSet.neighboringLocalScvIndices(localIdx);
            scvfs_[localIdx] = LocalScvfType(scvf, neighborScvIndicesLocal, localIdx, /*isDirichlet*/false);
            localFaceData_[localIdx*2] = LocalFaceData(localIdx, neighborScvIndicesLocal[0], scvf.index());

            const auto outsideLocalScvIdx = neighborScvIndicesLocal[1];
            // loop over scvfs in outside scv until we find the one coinciding with current scvf
            for (int coord = 0; coord < dim; ++coord)
            {
                if (indexSet.scvfIdxLocal(outsideLocalScvIdx, coord) == localIdx)
                {
                    const auto globalScvfIdx = indexSet.nodalIndexSet().scvfIdxGlobal(outsideLocalScvIdx, coord);
                    const auto& flipScvf = fvGeometry.scvf(globalScvfIdx);
                    localFaceData_[localIdx*2 + 1] = LocalFaceData(localIdx,           // iv-local scvf idx
                                                                   outsideLocalScvIdx, // iv-local scv index
                                                                   0,                  // scvf-local index in outside faces
                                                                   flipScvf.index());  // global scvf index
                }
            }
        }

        // Maybe resize local matrices if dynamic types are used
        resizeMatrix(A_, localSize, localSize);
        resizeMatrix(B_, localSize, localSize);
        resizeMatrix(C_, localSize, localSize);
    }

    //! returns the number of primary scvfs of this interaction volume
    static constexpr std::size_t numFaces() { return localSize; }

    //! returns the number of intermediate unknowns within this interaction volume
    static constexpr std::size_t numUnknowns() { return localSize; }

    //! returns the number of (in this context) known solution values within this interaction volume
    static constexpr std::size_t numKnowns() { return localSize; }

    //! returns the number of scvs embedded in this interaction volume
    static constexpr std::size_t numScvs() { return localSize; }

    //! returns the cell-stencil of this interaction volume
    const Stencil& stencil() const { return *stencil_; }

    //! returns the grid element corresponding to a given iv-local scv idx
    const Element& element(const LocalIndexType ivLocalScvIdx) const
    {
        assert(ivLocalScvIdx < localSize);
        return elements_[ivLocalScvIdx];
    }

    //! returns the local scvf entity corresponding to a given iv-local scvf idx
    const LocalScvfType& localScvf(const LocalIndexType ivLocalScvfIdx) const
    {
        assert(ivLocalScvfIdx < localSize);
        return scvfs_[ivLocalScvfIdx];
    }

    //! returns the local scv entity corresponding to a given iv-local scv idx
    const LocalScvType& localScv(const LocalIndexType ivLocalScvIdx) const
    {
        assert(ivLocalScvIdx < localSize);
        return scvs_[ivLocalScvIdx];
    }

    //! returns a reference to the container with the local face data
    const std::array<LocalFaceData, 2*localSize>& localFaceData() const { return localFaceData_; }

    //! Returns a reference to the information container on Dirichlet BCs within this iv.
    //! Here, we return an empty container as this implementation cannot be used on boundaries.
    const std::array<DirichletData, 0>& dirichletData() const { return dirichletData_; }

    //! returns the matrix associated with face unknowns in local equation system
    const Matrix& A() const { return A_; }
    Matrix& A() { return A_; }

    //! returns the matrix associated with cell unknowns in local equation system
    const Matrix& B() const { return B_; }
    Matrix& B() { return B_; }

    //! returns the matrix associated with face unknowns in flux expressions
    const Matrix& C() const { return C_; }
    Matrix& C() { return C_; }

    //! returns container storing the transmissibilities for each face & coordinate
    const std::array< std::array<DimVector, 2>, localSize >& omegas() const { return wijk_; }
    std::array< std::array<DimVector, 2>, localSize >& omegas() { return wijk_; }

    //! returns the number of interaction volumes living around a vertex
    //! the mpfa-o scheme always constructs one iv per vertex
    template< class NodalIndexSet >
    static std::size_t numInteractionVolumesAtVertex(const NodalIndexSet& nodalIndexSet) { return 1; }

    //! adds the iv index sets living around a vertex to a given container
    //! and stores the the corresponding index in a map for each scvf
    template<class IvIndexSetContainer, class ScvfIndexMap, class NodalIndexSet>
    static void addInteractionVolumeIndexSets(IvIndexSetContainer& ivIndexSetContainer,
                                              ScvfIndexMap& scvfIndexMap,
                                              const NodalIndexSet& nodalIndexSet)
    {
        // the global index of the iv index set that is about to be created
        const auto curGlobalIndex = ivIndexSetContainer.size();

        // make the one index set for this node
        ivIndexSetContainer.emplace_back(nodalIndexSet);

        // store the index mapping
        for (const auto scvfIdx : nodalIndexSet.globalScvfIndices())
            scvfIndexMap[scvfIdx] = curGlobalIndex;
    }

private:
    // pointer to cell stencil (in iv index set)
    const Stencil* stencil_;

    // Variables defining the local scope
    std::array<Element, localSize> elements_;
    std::array<LocalScvType, localSize> scvs_;
    std::array<LocalScvfType, localSize> scvfs_;
    std::array<LocalFaceData, 2*localSize> localFaceData_;

    // Empty Dirichlet container to be compatible with dynamic assembly
    std::array<DirichletData, 0> dirichletData_;

    // Matrices needed for computation of transmissibilities
    Matrix A_;
    Matrix B_;
    Matrix C_;

    // The omega factors are stored during assembly of local system
    // we assume one "outside" face per scvf (does not work on surface grids)
    std::array< std::array<DimVector, 2>, localSize > wijk_;
};

} // end namespace

#endif
