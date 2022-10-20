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
#ifndef DUMUX_DISCRETIZATION_CCTPFA_DETAIL_HH
#define DUMUX_DISCRETIZATION_CCTPFA_DETAIL_HH
#ifndef DOXYGEN

#include <dumux/experimental/new_assembly/dumux/common/size.hh>
#include <dumux/experimental/new_assembly/dumux/discretization/fvgridgeometries.hh>

namespace Dumux::CCTpfa::Detail {

template<typename Traits>
inline constexpr auto maxNumElementFaces = Traits::maxNumScvfsPerElement;

template<typename Traits>
inline constexpr auto maxNumFaceNeighbors = Traits::maxNumBranchesPerScvf + 1;

template<typename Traits>
inline constexpr auto maxNumStencilScvfs = maxNumElementFaces<Traits>*maxNumFaceNeighbors<Traits>;

template<typename Traits>
inline constexpr auto maxElementStencilSize = maxNumElementFaces<Traits>*Traits::maxNumBranchesPerScvf + 1;

// wrapper around entities to attach an id used for identification by the storing class
template<typename BaseEntity, typename Friend>
class WrappedEntity : public BaseEntity
{
public:
    using BaseEntity::BaseEntity;

    operator const BaseEntity&() const
    { return static_cast<const BaseEntity&>(*this); }

private:
    friend Friend;
    std::size_t id;
};

template<typename GV> using Coordinate = typename GV::template Codim<0>::Entity::Geometry::GlobalCoordinate;
template<typename GV> using Face = Dumux::Face<Coordinate<GV>>;
template<typename GV> using Scvf = SubControlVolumeFace<Coordinate<GV>>;
template<typename GV> using Scv = CCSubControlVolume<typename GV::IndexSet::IndexType, Coordinate<GV>>;
template<typename GV, typename Friend> using ScvWithId = WrappedEntity<Scv<GV>, Friend>;
template<typename GV, typename Friend> using ScvfWithId = WrappedEntity<Scvf<GV>, Friend>;

} // namespace Dumux::CCTpfa::Detail

#endif // DOXYGEN
#endif
