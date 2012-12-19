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
 *
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 */
#ifndef DUMUX_BOX_VOLUME_VARIABLES_HH
#define DUMUX_BOX_VOLUME_VARIABLES_HH

#include <dumux/implicit/common/implicitvolumevariables.hh>

namespace Dumux
{

/*!
 * \ingroup BoxModel
 * \ingroup BoxVolumeVariables
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 */
template <class TypeTag>
class BoxVolumeVariables : public ImplicitVolumeVariables<TypeTag>
{
    typedef ImplicitVolumeVariables<TypeTag> ParentType;

public:
    DUNE_DEPRECATED_MSG("Use ImplicitVolumeVariables from "
                        "dumux/implicit/common/implicitvolumevariables.hh.")
    BoxVolumeVariables() : ParentType() {}

    DUNE_DEPRECATED_MSG("Use ImplicitVolumeVariables from "
                        "dumux/implicit/common/implicitvolumevariables.hh.")
    BoxVolumeVariables(const BoxVolumeVariables &v)
    : ParentType(v) {}
};

} // end namepace

#endif
