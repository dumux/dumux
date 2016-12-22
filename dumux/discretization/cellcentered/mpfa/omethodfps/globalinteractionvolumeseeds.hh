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
 * \brief Base class for the global interaction volumes of the mpfa-o-fps method.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_O_FPS_GLOBALINTERACTIONVOLUMESEEDS_HH
#define DUMUX_DISCRETIZATION_MPFA_O_FPS_GLOBALINTERACTIONVOLUMESEEDS_HH

#include <dumux/discretization/cellcentered/mpfa/globalinteractionvolumeseedsbase.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

namespace Dumux
{
/*!
 * \ingroup Mpfa
 * \brief Specialization of the class for the mpfa-o-fps method.
 */
template<class TypeTag>
class CCMpfaGlobalInteractionVolumeSeedsImplementation<TypeTag, MpfaMethods::oMethodFps>
      : public CCMpfaGlobalInteractionVolumeSeedsImplementation<TypeTag, MpfaMethods::oMethod>
{
    using ParentType = CCMpfaGlobalInteractionVolumeSeedsImplementation<TypeTag, MpfaMethods::oMethod>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
public:
    CCMpfaGlobalInteractionVolumeSeedsImplementation(const GridView& gridView) : ParentType(gridView) {}
};
} // end namespace


#endif
