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
 * \ingroup ZeroEqModel
 *
 * \brief A single-phase, isothermal Reynolds-Averaged Navier-Stokes 0-Eq. model
 *
 * \copydoc RANSModel
 *
 * These models calculate the eddy viscosity without solving additional PDEs,
 * only based on the wall distance and the velocity gradient.
 * The following models are available:
 *  -# Prandtl's mixing length, e.g. \cite Oertel2012a
 *  -# Van-Driest modification, \cite vanDriest1956a and \cite Hanna1981a
 */

#ifndef DUMUX_ZEROEQ_MODEL_HH
#define DUMUX_ZEROEQ_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/rans/model.hh>

#include "indices.hh"
#include "volumevariables.hh"

namespace Dumux
{

// \{
///////////////////////////////////////////////////////////////////////////
// properties for the single-phase Reynolds-Averaged Navier-Stokes 0-Eq. model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal Reynolds-Averaged Navier-Stokes 0-Eq. model
NEW_TYPE_TAG(ZeroEq, INHERITS_FROM(RANS));

//! The type tag for the corresponding non-isothermal model
NEW_TYPE_TAG(ZeroEqNI, INHERITS_FROM(ZeroEq, RANSNI));

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////

//! The indices
SET_PROP(ZeroEq, Indices)
{
private:
    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int dim = GET_PROP_TYPE(TypeTag, GridView)::dimension;
public:
    using type = ZeroEqIndices<dim, numEq>;
};

//! The volume variables
SET_TYPE_PROP(ZeroEq, VolumeVariables, ZeroEqVolumeVariables<TypeTag>);

// \}
}

} // end namespace

#endif // DUMUX_ZEROEQ_MODEL_HH
