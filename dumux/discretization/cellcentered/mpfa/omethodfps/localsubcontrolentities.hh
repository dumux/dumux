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
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_FPS_LOCALSUBCONTROLENTITIES_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_FPS_LOCALSUBCONTROLENTITIES_HH

#include "mpfafpsgeometryhelper.hh"

namespace Dumux
{
template<class TypeTag>
class CCMpfaOFpsLocalScv
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
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
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using GeometryHelper = Dumux::MpfaFpsGeometryHelper<GridView, dim>;

public:
    using Geometry = Dune::MultiLinearGeometry<Scalar, dim, dimWorld>;

    //! The constructor
    //! Initialization in the initializer list is desired to not have to use the optional (see scv).
    //! Here, we use the first scvf in the seed to extract an scvf and the local vertex index it is connected to.
    //! We know that in the fps scheme all the scvfs are connected to the same vertex. The MpfaFpsGeometry Helper
    //! is then used to extract the corners of the given scv (in the geometry constructor).
    CCMpfaOFpsLocalScv(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const LocalScvSeed& scvSeed)
    : seedPtr_(&scvSeed),
      geometry_(Dune::GeometryType(Dune::GeometryType::cube, dim),
                GeometryHelper(element.geometry()).getScvCorners(fvGeometry.scvf(scvSeed.globalScvfIndices()[0]).vertexIndexInElement()))
    {}

    GlobalIndexType globalIndex() const
    { return scvSeed_().globalIndex(); }

    GlobalIndexType localScvfIndex(LocalIndexType coordDir) const
    {
        assert(coordDir < dim);
        return scvSeed_().localScvfIndices()[coordDir];
    }

    LocalIndexType getScvfIdxInScv(LocalIndexType localScvfIndex) const
    {
        auto it = std::find(scvSeed_().localScvfIndices().begin(), scvSeed_().localScvfIndices().end(), localScvfIndex);
        assert(it != scvSeed_().localScvfIndices().end() && "Could not find the local coordinate of the local scvf");
        return std::distance(scvSeed_().localScvfIndices().begin(), it);
    }

    //! Returns the center of the element the scv is embedded in
    GlobalPosition center() const
    { return geometry().corner(0); }

    const Geometry& geometry() const
    { return geometry_; }

private:
    const LocalScvSeed& scvSeed_() const
    { return *seedPtr_; }

    const LocalScvSeed* seedPtr_;
    Geometry geometry_;
};
} // end namespace

#endif
