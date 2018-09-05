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
 *
 * The following models are available:
 *  -# Prandtl's mixing length, e.g. \cite Oertel2012a
 *  -# Van-Driest modification, \cite vanDriest1956a and \cite Hanna1981a
 *  -# Baldwin-Lomax, \cite Baldwin1978a
 */

#ifndef DUMUX_ZEROEQ_MODEL_HH
#define DUMUX_ZEROEQ_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/rans/model.hh>
#include <dumux/freeflow/navierstokes/volumevariables.hh>

#include "volumevariables.hh"
#include "gridvariables.hh"

namespace Dumux {
namespace Properties {

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal RANS 0-Eq. model
///////////////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal Reynolds-Averaged Navier-Stokes 0-Eq. model
NEW_TYPE_TAG(ZeroEq, INHERITS_FROM(RANS));

//! Set the volume variables property
SET_PROP(ZeroEq, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    static_assert(FSY::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = ZeroEqVolumeVariables<Traits, NSVolVars>;
};

//! Set the grid variables (volume, flux and face variables)
SET_PROP(ZeroEq, GridVariables)
{
private:
    using GG = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GVV = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using GFVC = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache);
    using GFV = typename GET_PROP_TYPE(TypeTag, GridFaceVariables);
    using NavierStokesGridVariables = StaggeredGridVariables<GG, GVV, GFVC, GFV>;
public:
    using type = ZeroEqGridVariables<NavierStokesGridVariables>;
};

//////////////////////////////////////////////////////////////////
// default property values for the non-isothermal RANS 0-Eq. model
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal Reynolds-Averaged Navier-Stokes model
NEW_TYPE_TAG(ZeroEqNI, INHERITS_FROM(RANSNI));

//! Set the volume variables property
SET_PROP(ZeroEqNI, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    static_assert(FSY::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = ZeroEqVolumeVariables<Traits, NSVolVars>;
};

//! Set the grid variables (volume, flux and face variables)
SET_PROP(ZeroEqNI, GridVariables)
{
private:
    using GG = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GVV = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables);
    using GFVC = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache);
    using GFV = typename GET_PROP_TYPE(TypeTag, GridFaceVariables);
    using NavierStokesGridVariables = StaggeredGridVariables<GG, GVV, GFVC, GFV>;
public:
    using type = ZeroEqGridVariables<NavierStokesGridVariables>;
};

// \}
}

} // end namespace

#endif // DUMUX_ZEROEQ_MODEL_HH
