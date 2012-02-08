// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup Properties
 * \ingroup BoxProperties
 * \ingroup BoxStokes2cModel
 *
 * \file
 *
 * \brief Defines the supplementary properties required for the compositional
 * stokes BOX model.
 */
#ifndef DUMUX_STOKES2CPROPERTIES_HH
#define DUMUX_STOKES2CPROPERTIES_HH

#include <dumux/freeflow/stokes/stokesproperties.hh>

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the compositional stokes problems
NEW_TYPE_TAG(BoxStokes2c, INHERITS_FROM(BoxStokes));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(Stokes2cIndices); //!< Enumerations for the compositional stokes models
NEW_PROP_TAG(NumComponents); //!< Number of components

}

}

#endif
