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
 * \brief Class for an mpfa-o sub control volume face
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_L_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_L_SUBCONTROLVOLUMEFACE_HH

#include <dumux/discretization/cellcentered/mpfa/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/subcontrolvolumefacebase.hh>

namespace Dumux
{

/*!
 * \ingroup Mpfa
 * \brief Class for a sub control volume face in the mpfa-o method. We simply inherit from the base class here.
 */
template<class G, class GT, typename I>
class CCMpfaSubControlVolumeFaceImplementation<MpfaMethods::lMethod, G, GT, I> : public CCMpfaSubControlVolumeFaceBase<G, GT, I>
{
    using ParentType = CCMpfaSubControlVolumeFaceBase<G, GT, I>;
    using Geometry = G;
    using IndexType = I;

    using Scalar = typename Geometry::ctype;
    static const int dim = Geometry::mydimension;
    static const int dimWorld = Geometry::coorddimension;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Corners = typename GT::template CornerStorage<dim, dimWorld>;

public:
    //! We do not use the localIndex variable here.
    //! It is here to satisfy the general mpfa scvf interface.
    template<class MpfaGeometryHelper>
    CCMpfaSubControlVolumeFaceImplementation(const MpfaGeometryHelper& geomHelper,
                                             Corners&& corners,
                                             GlobalPosition&& unitOuterNormal,
                                             IndexType vertexIndex,
                                             unsigned int localIndex,
                                             IndexType scvfIndex,
                                             IndexType insideScvIdx,
                                             const std::vector<IndexType>& outsideScvIndices,
                                             Scalar q,
                                             bool boundary)
    : ParentType(geomHelper,
                 std::forward<Corners>(corners),
                 std::forward<GlobalPosition>(unitOuterNormal),
                 vertexIndex,
                 localIndex,
                 scvfIndex,
                 insideScvIdx,
                 outsideScvIndices,
                 0.0, // q should always be zero for the mpfa-l method
                 boundary) {}
};

} // end namespace

#endif
