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
 * \brief Classes for sub control entities of the mpfa-o method.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_LOCAL_SUBCONTROLENTITIES_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_LOCAL_SUBCONTROLENTITIES_HH

#include <dune/common/fvector.hh>

namespace Dumux
{

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Class for the interaction volume-local sub-control volume used
 *        in the mpfa-o scheme.
 *
 * \tparam IvIndexSet The type used for index sets within interaction volumes
 * \tparam dim The dimensionality of the grid
 * \tparam dimWorld The dimension of the world the grid is embedded in
 */
template< class IvIndexSet, class Scalar, int dim, int dimWorld>
class CCMpfaOInteractionVolumeLocalScv
{

public:
    // export some types
    using GridIndexType = typename IvIndexSet::GridIndexType;
    using LocalIndexType = typename IvIndexSet::LocalIndexType;
    using GlobalCoordinate = Dune::FieldVector<Scalar, dimWorld>;
    using ctype = typename GlobalCoordinate::value_type;
    using LocalBasis = std::array< GlobalCoordinate, dim >;

    static constexpr int myDimension = dim;
    static constexpr int worldDimension = dimWorld;

    //! The default constructor
    CCMpfaOInteractionVolumeLocalScv() = default;

    /*!
     * \brief The constructor
     *
     * \param helper Helper class for mpfa schemes
     * \param fvGeometry The element finite volume geometry
     * \param scv The grid sub-control volume
     * \param localIndex The iv-local index of this scvIdx
     * \param indexSet The interaction volume index set
     */
    template<class MpfaHelper, class FVElementGeometry, class SubControlVolume>
    CCMpfaOInteractionVolumeLocalScv(const MpfaHelper& helper,
                                     const FVElementGeometry& fvGeometry,
                                     const SubControlVolume& scv,
                                     const LocalIndexType localIndex,
                                     const IvIndexSet& indexSet)
    : indexSet_(&indexSet)
    , globalScvIndex_(scv.dofIndex())
    , localDofIndex_(localIndex)
    {
        // center of the global scv
        const auto& center = scv.center();

        // set up local basis
        LocalBasis localBasis;
        for (unsigned int coordIdx = 0; coordIdx < myDimension; ++coordIdx)
        {
            const auto scvfIdx = indexSet.nodalIndexSet().gridScvfIndex(localDofIndex_, coordIdx);
            const auto& scvf = fvGeometry.scvf(scvfIdx);
            localBasis[coordIdx] = scvf.ipGlobal();
            localBasis[coordIdx] -= center;
        }

        nus_ = helper.calculateInnerNormals(localBasis);
        detX_ = helper.calculateDetX(localBasis);
    }

    //! detX is needed for setting up the omegas in the interaction volumes
    ctype detX() const
    { return detX_; }

    //! grid index related to this scv
    GridIndexType gridScvIndex() const
    { return globalScvIndex_; }

    //! returns the index in the set of cell unknowns of the iv
    LocalIndexType localDofIndex() const
    { return localDofIndex_; }

    //! iv-local index of the coordir's scvf in this scv
    LocalIndexType localScvfIndex(unsigned int coordDir) const
    {
        assert(coordDir < myDimension);
        return indexSet_->localScvfIndex(localDofIndex_, coordDir);
    }

    //! the nu vectors are needed for setting up the omegas of the iv
    const GlobalCoordinate& nu(unsigned int coordDir) const
    {
        assert(coordDir < myDimension);
        return nus_[coordDir];
    }

private:
    const IvIndexSet* indexSet_;
    GridIndexType globalScvIndex_;
    LocalIndexType localDofIndex_;
    LocalBasis nus_;
    ctype detX_;
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Class for the interaction volume-local sub-control volume face
 *        used in the mpfa-o scheme.
 *
 * \tparam IvIndexSet The type used for index sets within interaction volumes
 */
template< class IvIndexSet >
struct CCMpfaOInteractionVolumeLocalScvf
{
  using ScvfNeighborLocalIndexSet = typename IvIndexSet::ScvfNeighborLocalIndexSet;

public:
    // export index types
    using GridIndexType = typename IvIndexSet::GridIndexType;
    using LocalIndexType = typename IvIndexSet::LocalIndexType;

    //! The default constructor
    CCMpfaOInteractionVolumeLocalScvf() = default;

    /*!
     * \brief The constructor
     *
     * \param scvf The grid sub-control volume face
     * \param localScvIndices The iv-local neighboring scv indices
     * \param localDofIdx This scvf's interaction volume-local dof index
     * \param isDirichlet Specifies if this scv is on a Dirichlet boundary
     */
    template< class SubControlVolumeFace >
    CCMpfaOInteractionVolumeLocalScvf(const SubControlVolumeFace& scvf,
                                      const ScvfNeighborLocalIndexSet& localScvIndices,
                                      const LocalIndexType localDofIdx,
                                      const bool isDirichlet)
    : isDirichlet_(isDirichlet)
    , scvfIdxGlobal_(scvf.index())
    , localDofIndex_(localDofIdx)
    , neighborScvIndicesLocal_(&localScvIndices)
    {}

    //! This is either the iv-local index of the intermediate unknown (interior/Neumann face)
    //! or the index of the Dirichlet boundary within the vol vars (Dirichlet faces)
    LocalIndexType localDofIndex() const { return localDofIndex_; }

    //! returns the grid view-global index of this scvf
    GridIndexType gridScvfIndex() const { return scvfIdxGlobal_; }

    //! Returns the local indices of the scvs neighboring this scvf
    const ScvfNeighborLocalIndexSet& neighboringLocalScvIndices() const { return *neighborScvIndicesLocal_; }

    //! states if this is scvf is on a Dirichlet boundary
    bool isDirichlet() const { return isDirichlet_; }

private:
    bool isDirichlet_;
    GridIndexType scvfIdxGlobal_;
    LocalIndexType localDofIndex_;
    const ScvfNeighborLocalIndexSet* neighborScvIndicesLocal_;
};

} // end namespace Dumux

#endif
