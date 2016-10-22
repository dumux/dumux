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

#include <dumux/implicit/cellcentered/mpfa/properties.hh>

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

    // the default constructor
    CCMpfaOFpsLocalScv() = default;

    // the constructor
    CCMpfaOFpsLocalScv(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const LocalScvSeed& scvSeed)
    : seedPtr_(&scvSeed)
    {
        // the geometry helper that will give us the scv corners
        GeometryHelper geomHelper(element.geometry());

        // extract the vertex index from the first scvf
        auto vIdxGlobal = fvGeometry.scvf(scvSeed.globalScvfIndices()[0]).vertexIndex();

        // find local index of the vertex in the element
        int vIdxLocal = -1;
        for (unsigned int localIdx = 0; localIdx < element.subEntities(dim); ++localIdx)
        {
            if (problem.vertexMapper().subIndex(element, localIdx, dim) == vIdxGlobal)
            {
                vIdxLocal = localIdx;
                break;
            }
        }

        assert(vIdxLocal != -1 && "could not find the local index of the scv in the element!");

        // construct the geometry of the scv
        geometry_ = Geometry(Dune::GeometryType(Dune::GeometryType::cube, dim), geomHelper.getScvCorners(vIdxLocal));
    }

    //! The copy constrcutor
    CCMpfaOFpsLocalScv(const CCMpfaOFpsLocalScv& other) = default;

    //! The move constrcutor
    CCMpfaOFpsLocalScv(CCMpfaOFpsLocalScv&& other) = default;

    //! The copy assignment operator
    CCMpfaOFpsLocalScv& operator=(const CCMpfaOFpsLocalScv& other)
    {
        // We want to use the default copy/move assignment.
        // But since geometry is not copy assignable :( we
        // have to construct it again
        geometry_.release();
        geometry_.emplace(other.geometry_.value());
        return *this;
    }

    //! The move assignment operator
    CCMpfaOFpsLocalScv& operator=(CCMpfaOFpsLocalScv&& other)
    {
        // We want to use the default copy/move assignment.
        // But since geometry is not copy assignable :( we
        // have to construct it again
        geometry_.release();
        geometry_.emplace(other.geometry_.value());
        return *this;
    }

    GlobalIndexType globalIndex() const
    { return scvSeed_().globalIndex(); }

    GlobalIndexType localScvfIndex(const LocalIndexType coordDir) const
    {
        assert(coordDir < dim);
        return scvSeed_().localScvfIndices()[coordDir];
    }

    LocalIndexType getScvfIdxInScv(const LocalIndexType localScvfIndex) const
    {
        auto it = std::find(scvSeed_().localScvfIndices().begin(), scvSeed_().localScvfIndices().end(), localScvfIndex);
        assert(it != scvSeed_().localScvfIndices().end() && "Could not find the local coordinate of the local scvf");
        return std::distance(scvSeed_().localScvfIndices().begin(), it);
    }

    //! Returns the center of the element the scv is embedded in
    GlobalPosition center() const
    { return geometry().corner(0); }

    const Geometry& geometry() const
    { return geometry_.value(); }

private:
    const LocalScvSeed& scvSeed_() const
    { return *seedPtr_; }

    const LocalScvSeed* seedPtr_;
    Optional<Geometry> geometry_;
};
} // end namespace

#endif
