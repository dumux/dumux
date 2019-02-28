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
 * \brief Classes for sub control entities of the
 *        mpfa-o method in the context of facet coupling
 */
#ifndef DUMUX_MULTDOMAIN_FACET_CC_MPFA_O_LOCAL_SUBCONTROLENTITIES_HH
#define DUMUX_MULTDOMAIN_FACET_CC_MPFA_O_LOCAL_SUBCONTROLENTITIES_HH

#include <array>

#include <dune/common/fvector.hh>
#include <dumux/discretization/cellcentered/mpfa/omethod/localsubcontrolentities.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Class for the interaction volume-local sub-control volume used
 *        in the mpfa-o scheme in the context of facet coupling.
 *
 * \tparam IvIndexSet The type used for index sets within interaction volumes
 * \tparam dim The dimensionality of the grid
 * \tparam dimWorld The dimension of the world the grid is embedded in
 */
template< class IvIndexSet, class Scalar, int dim, int dimWorld>
class CCMpfaOFacetCouplingInteractionVolumeLocalScv
: public CCMpfaOInteractionVolumeLocalScv<IvIndexSet, Scalar, dim, dimWorld>
{
    using ParentType = CCMpfaOInteractionVolumeLocalScv<IvIndexSet, Scalar, dim, dimWorld>;

public:
    // export some types
    using typename ParentType::LocalIndexType;

    // export dimension
    static constexpr int myDimension = dim;

    //! The default constructor
    CCMpfaOFacetCouplingInteractionVolumeLocalScv() = default;

    /*!
     * \brief The constructor
     *
     * \param helper Helper class for mpfa schemes
     * \param fvGeometry The element finite volume geometry
     * \param scv The grid sub-control volume
     * \param localIndex The iv-local index of this scvIdx
     * \param indexSet The interaction volume index set
     * \param scvfGridToLocalIndexMap maps to grid scvf indices the iv-local scvf index
     */
    template<class MpfaHelper, class FVElementGeometry, class SubControlVolume, class IndexMap>
    CCMpfaOFacetCouplingInteractionVolumeLocalScv(const MpfaHelper& helper,
                                                  const FVElementGeometry& fvGeometry,
                                                  const SubControlVolume& scv,
                                                  const LocalIndexType localIndex,
                                                  const IvIndexSet& indexSet,
                                                  const IndexMap& scvfGridToLocalIndexMap)
    : ParentType(helper, fvGeometry, scv, localIndex, indexSet)
    {
        // set up local scvf indices
        const auto& nis = indexSet.nodalIndexSet();
        for (unsigned int dir = 0; dir < myDimension; ++dir)
            localScvfIndices_[dir] = scvfGridToLocalIndexMap.at(nis.gridScvfIndex(this->localDofIndex(), dir));
    }

    //! iv-local index of the coordir's scvf in this scv
    LocalIndexType localScvfIndex(unsigned int coordDir) const
    {
        assert(coordDir < myDimension);
        return localScvfIndices_[coordDir];
    }

private:
    std::array<LocalIndexType, dim> localScvfIndices_;
};

/*!
 * \ingroup FacetCoupling
 * \brief Class for the interaction volume-local sub-control volume face
 *        used in the mpfa-o scheme in the context of facet coupling.
 *
 * \tparam IvIndexSet The type used for index sets within interaction volumes
 */
template< class IvIndexSet >
struct CCMpfaOFacetCouplingInteractionVolumeLocalScvf
: public CCMpfaOInteractionVolumeLocalScvf< IvIndexSet >
{
  using ParentType = CCMpfaOInteractionVolumeLocalScvf< IvIndexSet >;
  using ScvfNeighborLocalIndexSet = typename IvIndexSet::ScvfNeighborLocalIndexSet;

public:
    // pull up index types
    using typename ParentType::LocalIndexType;
    using typename ParentType::GridIndexType;

    //! pull up parent's constructors
    using ParentType::ParentType;

    /*!
     * \brief The constructor for interior boundary faces.
     *
     * \param scvf The grid sub-control volume face
     * \param localScvIndices The iv-local neighboring scv indices
     * \param localDofIdx This scvf's interaction volume-local dof index
     * \param isDirichlet Specifies if this scv is on a Dirichlet boundary
     * \param coupledFacetLocalDof The local index of the coupled facet element
     *                             in the set of cell&Dirichlet values.
     */
    template< class SubControlVolumeFace >
    CCMpfaOFacetCouplingInteractionVolumeLocalScvf(const SubControlVolumeFace& scvf,
                                                   const ScvfNeighborLocalIndexSet& localScvIndices,
                                                   const LocalIndexType localDofIdx,
                                                   const bool isDirichlet,
                                                   const LocalIndexType coupledFacetLocalDof)
    : ParentType(scvf, neighborLocalScvIndices_, localDofIdx, isDirichlet)
    , isInteriorBoundary_(true)
    , coupledFacetLocalDofIndex_(coupledFacetLocalDof)
    , neighborLocalScvIndices_(localScvIndices)
    {}

    //! Returns the iv-local dof index of the coupled facet element
    LocalIndexType coupledFacetLocalDofIndex() const
    { assert(isInteriorBoundary_); return coupledFacetLocalDofIndex_; }

    //! Returns the local indices of the scvs neighboring this scvf
    const ScvfNeighborLocalIndexSet& neighboringLocalScvIndices() const
    { return isOnInteriorBoundary() ? neighborLocalScvIndices_ : ParentType::neighboringLocalScvIndices(); }

    //! Returns true if this face is on an interior boundary
    bool isOnInteriorBoundary() const { return isInteriorBoundary_; }

private:
    bool isInteriorBoundary_{false};
    LocalIndexType coupledFacetLocalDofIndex_{0};
    ScvfNeighborLocalIndexSet neighborLocalScvIndices_;
};

} // end namespace Dumux

#endif
