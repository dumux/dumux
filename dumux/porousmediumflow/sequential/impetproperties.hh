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
#ifndef DUMUX_IMPET_PROPERTIES_HH
#define DUMUX_IMPET_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/sequential/properties.hh>
#include <dumux/porousmediumflow/sequential/pressureproperties.hh>
#include <dumux/porousmediumflow/sequential/transportproperties.hh>

/*!
 * \ingroup IMPET
 * \ingroup IMPETProperties
 */
/*!
 * \file
 * \brief Base file for properties related to sequential IMPET algorithms
 */
namespace Dumux
{

template<class TypeTag>
class IMPET;

namespace Properties
{
/*!
 *
 * \brief General properties for sequential IMPET algorithms
 *
 * This class holds properties necessary for the sequential IMPET solution.
 */

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the diffusion-scheme
// Create new type tags
namespace TTag {
struct IMPET { using InheritsFrom = std::tuple<SequentialModel>; };
} // end namespace TTag
}
}

#include <dumux/porousmediumflow/sequential/impet.hh>

namespace Dumux
{
namespace Properties
{
//set impet model
template<class TypeTag>
struct Model<TypeTag, TTag::IMPET> { using type = IMPET<TypeTag>; };
}
}

#endif
