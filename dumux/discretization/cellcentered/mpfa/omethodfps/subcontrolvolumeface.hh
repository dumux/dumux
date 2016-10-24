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
 * \brief Class for an mpfao-fps sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_FPS_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_FPS_SUBCONTROLVOLUMEFACE_HH

#include <dumux/discretization/cellcentered/mpfa/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/omethod/subcontrolvolumeface.hh>

namespace Dumux
{

/*!
 * \ingroup Discretization
 * \brief Class for a sub control volume face in the mpfao-fps method.
 */
template<class G, typename I>
class CCMpfaSubControlVolumeFace<MpfaMethods::oMethodFps, G, I> : public CCMpfaSubControlVolumeFace<MpfaMethods::oMethod, G, I>
{
    using ParentType = CCMpfaSubControlVolumeFace<MpfaMethods::oMethod, G, I>;
    using Geometry = G;
    using IndexType = I;

    using Scalar = typename Geometry::ctype;
    static const int dim = Geometry::mydimension;
    static const int dimworld = Geometry::coorddimension;

    using GlobalPosition = Dune::FieldVector<Scalar, dimworld>;

public:
    template<class MpfaGeometryHelper>
    CCMpfaSubControlVolumeFace(const MpfaGeometryHelper& geomHelper,
                               std::vector<GlobalPosition>&& corners,
                               GlobalPosition&& unitOuterNormal,
                               IndexType vIdxGlobal,
                               unsigned int vIdxLocal,
                               IndexType scvfIndex,
                               std::array<IndexType, 2>&& scvIndices,
                               Scalar q,
                               bool boundary)
    : ParentType(geomHelper,
                 std::move(corners),
                 std::move(unitOuterNormal),
                 vIdxGlobal,
                 vIdxLocal,
                 scvfIndex,
                 std::move(scvIndices),
                 q,
                 boundary),
      vIdxInElement_(vIdxLocal)
      {}

    //! Returns the local index inside the element of the vertex the scvf is connected to
    IndexType vertexIndexInElement() const
    { return vIdxInElement_; }

private:
    unsigned int vIdxInElement_;
};

} // end namespace

#endif
