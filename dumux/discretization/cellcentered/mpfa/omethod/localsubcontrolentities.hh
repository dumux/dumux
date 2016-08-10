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
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Helper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);

    using Element = typename GridView::template Codim<0>::Entity;
    using LocalScvSeed = typename InteractionVolume::Seed::LocalScvSeed;

    using GlobalIndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = typename InteractionVolume::LocalIndexType;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using LocalBasis = std::array<DimVector, dim>;

public:
    // constructor has the same signature as the LocalScv entity
    CCMpfaOLocalScv(const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const LocalScvSeed& scvSeed)
    : seed_(scvSeed)
    {

        // et up local basis
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

    DimVector center() const
    { return center_; }

    DimVector innerNormal(const LocalIndexType coordDir) const
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
    DimVector center_;
    LocalBasis innerNormals_;
    Scalar detX_;
};


template<class TypeTag>
struct CCMpfaOLocalScvf
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Helper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using InteractionVolume = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    using Element = typename GridView::template Codim<0>::Entity;
    using LocalScvfSeed = typename InteractionVolume::Seed::LocalScvfSeed;

    using GlobalIndexType = typename GridView::IndexSet::IndexType;
    using LocalIndexType = typename InteractionVolume::LocalIndexType;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;

public:
    CCMpfaOLocalScvf(const LocalScvfSeed& scvfSeed,
                     const SubControlVolumeFace& scvf)
    : seed_(scvfSeed),
      center_(scvf.center()),
      ip_(scvf.ipGlobal()),
      normal_(scvf.unitOuterNormal()),
      area_(scvf.area())
    {}

    GlobalIndexType insideGlobalScvfIndex() const
    { return scvfSeed_().insideGlobalScvfIndex(); }

    GlobalIndexType outsideGlobalScvfIndex() const
    {
        assert(!boundary());
        return scvfSeed_().globalScvfIndices()[1];
    }

    LocalIndexType insideLocalScvIndex() const
    { return scvfSeed_().insideLocalScvIndex(); }

    LocalIndexType outsideLocalScvIndex() const
    {
        assert(faceType() != MpfaFaceTypes::neumann && faceType() != MpfaFaceTypes::dirichlet);
        return scvfSeed_().localScvIndices()[1];
    }

    GlobalIndexType insideGlobalScvIndex() const
    { return scvfSeed_().globalScvIndices()[0]; }

    GlobalIndexType outsideGlobalScvIndex() const
    { return scvfSeed_().globalScvIndices()[1];}

    MpfaFaceTypes faceType() const
    { return scvfSeed_().faceType(); }

    DimVector center() const
    { return center_; }

    DimVector ip() const
    { return ip_; }

    DimVector unitOuterNormal() const
    { return normal_; }

    Scalar area() const
    { return area_; }

    bool boundary() const
    { return faceType() != MpfaFaceTypes::interior; }

private:
    const LocalScvfSeed& scvfSeed_() const
    { return seed_; }

    const LocalScvfSeed& seed_;
    DimVector center_;
    DimVector ip_;
    DimVector normal_;
    Scalar area_;
};
} // end namespace

#endif
