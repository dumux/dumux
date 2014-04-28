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
#ifndef DUMUX_SUBDOMAIN_PROPERTIES_HH
#define DUMUX_SUBDOMAIN_PROPERTIES_HH

#include <dumux/common/propertysystem.hh>

namespace Dumux
{
namespace Properties
{
/*!
 * \addtogroup ModelCoupling
 */
// \{

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag from which sub-problems of coupling models inherit
NEW_TYPE_TAG(SubDomain);

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
//! Specifies the host grid
NEW_PROP_TAG(Grid);

//! Specifies the scalar grid function space used for sub-problems
NEW_PROP_TAG(ScalarGridFunctionSpace);

//! Specifies the grid function space used for sub-problems
NEW_PROP_TAG(GridFunctionSpace);

//! Specifies the grid operator used for sub-problems
NEW_PROP_TAG(GridOperator);

//! Specifies the grid operator space used for sub-problems
NEW_PROP_TAG(GridOperatorSpace);

//! Specifies the type of the constraints
NEW_PROP_TAG(Constraints);

//! Specifies the type of the constraints transformation
NEW_PROP_TAG(ConstraintsTrafo);

//! Specifies the local finite element space
NEW_PROP_TAG(LocalFEMSpace);

//! Specifies the local operator
NEW_PROP_TAG(LocalOperator);

} // namespace Properties
} // namespace Dumux
#endif // DUMUX_SUBDOMAIN_PROPERTIES_HH
