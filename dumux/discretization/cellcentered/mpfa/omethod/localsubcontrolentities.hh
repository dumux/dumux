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

#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>

namespace Dumux
{
template<class TypeTag>
class CCMpfaOLocalScv
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Helper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    using Element = typename GridView::template Codim<0>::Entity;

    // we use the seed types of the boundary interaction volume to be compatible with other mpfa
    // methods that use o-type interaction volumes on the boundary but differing ones inside the domain.
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using LocalScvSeed = typename InteractionVolume::Seed::LocalScvSeed;
    using GlobalIndexType = typename InteractionVolume::GlobalIndexType;
    using LocalIndexType = typename InteractionVolume::LocalIndexType;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using LocalBasis = std::array<GlobalPosition, dim>;

public:
    // constructor has the same signature as the LocalScv entity
    CCMpfaOLocalScv(const Problem& problem,
                    const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const LocalScvSeed& scvSeed)
    : seed_(scvSeed)
    {
        // set up local basis
        center_ = element.geometry().center();
        LocalBasis localBasis;

        LocalIndexType coordIdx = 0;
        for (auto globalScvfIdx : scvSeed.globalScvfIndices())
        {
            const auto& scvf = fvGeometry.scvf(globalScvfIdx);
            localBasis[coordIdx] = scvf.ipGlobal();
            localBasis[coordIdx] -= center_;
            coordIdx++;
        }

        innerNormals_ = Helper::calculateInnerNormals(localBasis);
        detX_ = Helper::calculateDetX(localBasis);
    }

    GlobalIndexType globalIndex() const
    { return scvSeed_().globalIndex(); }

    GlobalIndexType localScvfIndex(const LocalIndexType coordDir) const
    {
        assert(coordDir < dim);
        return scvSeed_().localScvfIndices()[coordDir];
    }

    GlobalPosition center() const
    { return center_; }

    GlobalPosition innerNormal(const LocalIndexType coordDir) const
    {
        assert(coordDir < dim);
        return innerNormals_[coordDir];
    }

    Scalar detX() const
    { return detX_; }

private:
    const LocalScvSeed& scvSeed_() const
    { return seed_; }

    const LocalScvSeed& seed_;
    GlobalPosition center_;
    LocalBasis innerNormals_;
    Scalar detX_;
};


template<class TypeTag>
struct CCMpfaOLocalScvf
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Helper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    using Element = typename GridView::template Codim<0>::Entity;

    // we use the seed types of the boundary interaction volume to be compatible with other mpfa
    // methods that use o-type interaction volumes on the boundary but differing ones inside the domain.
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    using LocalScvfSeed = typename InteractionVolume::Seed::LocalScvfSeed;
    using GlobalIndexType = typename InteractionVolume::GlobalIndexType;
    using GlobalIndexSet = typename InteractionVolume::GlobalIndexSet;
    using LocalIndexType = typename InteractionVolume::LocalIndexType;
    using LocalIndexSet = typename InteractionVolume::LocalIndexSet;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    CCMpfaOLocalScvf(const LocalScvfSeed& scvfSeed,
                     const SubControlVolumeFace& scvf)
    : seedPtr_(&scvfSeed),
      scvfPtr_(&scvf)
    {}

    GlobalIndexType insideGlobalScvfIndex() const
    { return scvfSeed_().insideGlobalScvfIndex(); }

    GlobalIndexType insideGlobalScvIndex() const
    { return scvfSeed_().insideGlobalScvIndex(); }

    LocalIndexType insideLocalScvIndex() const
    { return scvfSeed_().insideLocalScvIndex(); }

    const GlobalIndexSet& outsideGlobalScvfIndices() const
    { return scvfSeed_().outsideGlobalScvfIndices(); }

    const GlobalIndexSet& outsideGlobalScvIndices() const
    { return scvfSeed_().outsideGlobalScvIndices(); }

    const LocalIndexSet& outsideLocalScvIndices() const
    { return scvfSeed_().outsideLocalScvIndices(); }

    GlobalIndexType outsideGlobalScvfIndex(unsigned int outsideIdx = 0) const
    { return scvfSeed_().outsideGlobalScvfIndex(outsideIdx); }

    GlobalIndexType outsideGlobalScvIndex(unsigned int outsideIdx = 0) const
    { return scvfSeed_().outsideGlobalScvIndex(outsideIdx); }

    LocalIndexType outsideLocalScvIndex(unsigned int outsideIdx = 0) const
    { return scvfSeed_().outsideLocalScvIndex(outsideIdx); }

    MpfaFaceTypes faceType() const
    { return scvfSeed_().faceType(); }

    GlobalPosition ip() const
    { return globalScvf().ipGlobal(); }

    GlobalPosition unitOuterNormal() const
    { return globalScvf().unitOuterNormal(); }

    Scalar area() const
    { return globalScvf().area(); }

    bool boundary() const
    { return scvfSeed_().boundary(); }

    const SubControlVolumeFace& globalScvf() const
    { return *scvfPtr_; }

private:
    const LocalScvfSeed& scvfSeed_() const
    { return *seedPtr_; }

    const LocalScvfSeed* seedPtr_;
    const SubControlVolumeFace* scvfPtr_;
};
} // end namespace

#endif
