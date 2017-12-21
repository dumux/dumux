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
 * \brief Base class for sub control entities of the mpfa-o method.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_LOCALSUBCONTROLENTITIES_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_LOCALSUBCONTROLENTITIES_HH

#include <dune/common/fvector.hh>
#include <dumux/common/properties.hh>

namespace Dumux
{
template<class TypeTag, class IvIndexSet>
class CCMpfaOInteractionVolumeLocalScv
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using GlobalIndexType = typename GridView::IndexSet::IndexType;
    using Helper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);

    //! The mpfa-o method always uses the dynamic types
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
    using LocalIndexType = typename InteractionVolume::Traits::LocalIndexType;
    using LocalBasis = typename InteractionVolume::Traits::ScvBasis;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    explicit CCMpfaOInteractionVolumeLocalScv(const FVElementGeometry& fvGeometry,
                                              const SubControlVolume& scv,
                                              const LocalIndexType localIndex,
                                              const IvIndexSet& indexSet)
             : indexSet_(indexSet)
             , globalScvIndex_(scv.dofIndex())
             , localDofIndex_(localIndex)
    {
        // center of the global scv
        const auto center = scv.center();

        // set up local basis
        LocalBasis localBasis;
        for (unsigned int coordIdx = 0; coordIdx < dim; ++coordIdx)
        {
            const auto scvfIdx = indexSet.nodalIndexSet().scvfIdxGlobal(localDofIndex_, coordIdx);
            const auto& scvf = fvGeometry.scvf(scvfIdx);
            localBasis[coordIdx] = scvf.ipGlobal();
            localBasis[coordIdx] -= center;
        }

        innerNormals_ = Helper::calculateInnerNormals(localBasis);
        detX_ = Helper::calculateDetX(localBasis);
    }

    //! detX is needed for setting up the omegas in the interaction volumes
    Scalar detX() const { return detX_; }

    //! the inner normals are needed for setting up the omegas in the interaction volumes
    const GlobalPosition& innerNormal(const LocalIndexType coordDir) const
    {
        assert(coordDir < dim);
        return innerNormals_[coordDir];
    }

    //! grid view-global index related to this scv
    GlobalIndexType globalScvIndex() const { return globalScvIndex_; }

    //! returns the index in the set of cell unknowns of the iv
    LocalIndexType localDofIndex() const { return localDofIndex_; }

    //! iv-local index of the coordir's scvf in this scv
    LocalIndexType scvfIdxLocal(const unsigned int coordDir) const
    {
        assert(coordDir < dim);
        return indexSet_.scvfIdxLocal(localDofIndex_, coordDir);
    }

private:
    const IvIndexSet& indexSet_;
    GlobalIndexType globalScvIndex_;
    LocalIndexType localDofIndex_;
    LocalBasis innerNormals_;
    Scalar detX_;
};


template<class TypeTag>
struct CCMpfaOInteractionVolumeLocalScvf
{
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    //! The mpfa-o method always uses the dynamic types
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using GlobalIndexType = typename GridView::IndexSet::IndexType;
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
    using LocalIndexContainer = typename InteractionVolume::Traits::DynamicLocalIndexContainer;
    using LocalIndexType = typename LocalIndexContainer::value_type;

public:
    CCMpfaOInteractionVolumeLocalScvf(const SubControlVolumeFace& scvf,
                                      const LocalIndexContainer& localScvIndices,
                                      const bool isDirichlet,
                                      const LocalIndexType localDofIdx)
    : scvfIdxGlobal_(scvf.index())
    , neighborScvIndicesLocal_(localScvIndices)
    , isDirichlet_(isDirichlet)
    , localDofIndex_(localDofIdx)
    {}

    //! returns the grid view-global index of this scvf
    GlobalIndexType globalScvfIndex() const { return scvfIdxGlobal_; }

    //! Returns the local indices of the scvs neighboring this scvf
    const LocalIndexContainer& neighboringLocalScvIndices() const { return neighborScvIndicesLocal_; }

    //! states if this is scvf is on a Dirichlet boundary
    bool isDirichlet() const { return isDirichlet_; }

    //! This is either the iv-local index of the intermediate unknown (interior/Neumann face)
    //! or the index of the Dirichlet boundary within the vol vars (Dirichlet faces)
    LocalIndexType localDofIndex() const { return localDofIndex_; }

private:
    GlobalIndexType scvfIdxGlobal_;
    const LocalIndexContainer& neighborScvIndicesLocal_;

    bool isDirichlet_;
    LocalIndexType localDofIndex_;
};
} // end namespace

#endif
