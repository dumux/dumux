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
 * \ingroup RANSModel
 *
 * \brief A single-phase, multi-component isothermal Navier-Stokes model
 *
 * \copydoc Dumux::RANSModel
 *
 * The system is closed by a <B> component mass/mole balance equation </B> for each component \f$\kappa\f$:
 * \f[
 *    \frac{\partial \left(\varrho X^\kappa\right)}{\partial t}
 *    + \nabla \cdot \left( \varrho {\boldsymbol{v}} X^\kappa
 *    - (D^\kappa + D_\text{t}) \varrho \frac{M^\kappa}{M} \textbf{grad}\, x^\kappa \right)
 *    - q^\kappa = 0
 * \f]
 *
 * Alternatively, one component balance equation can be replace by a <B> total mass/mole balance equation </B>:
 * \f[
 *    \frac{\partial \varrho_g}{\partial t}
 *    + \nabla \cdot \left(
 *        \varrho {\boldsymbol{v}}
 *        - \sum_\kappa (D^\kappa + D_\text{t}) \varrho \frac{M^\kappa}{M} \textbf{grad}\, x^\kappa
 *      \right)
 *    - q = 0
 * \f]
 *
 * The eddy diffusivity \f$ D_\text{t} \f$ is related to the eddy viscosity \f$ \nu_\text{t} \f$
 * by the turbulent Schmidt number:
 * \f[ D_\text{t} = \frac{\nu_\text{t}}{\mathrm{Sc}_\text{t}} \f]
 *
 * So far, only the staggered grid spatial discretization (for structured grids) is available.
 */

#ifndef DUMUX_RANS_NC_MODEL_HH
#define DUMUX_RANS_NC_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/freeflow/nonisothermal/ransvtkoutputfields.hh>
#include <dumux/freeflow/rans/zeroeq/volumevariables.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux {

///////////////////////////////////////////////////////////////////////////
// properties for the single-phase, multi-component RANS model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the single-phase, multi-component isothermal RANS model
NEW_TYPE_TAG(RANSNC, INHERITS_FROM(NavierStokesNC));
NEW_TYPE_TAG(ZeroEqNC, INHERITS_FROM(RANSNC));

///////////////////////////////////////////////////////////////////////////
// default property values
///////////////////////////////////////////////////////////////////////////

//! Set the volume variables property
SET_PROP(ZeroEqNC, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using CompositionalVolVars = FreeflowNCVolumeVariables<Traits>;
    using RANSVolVars = ZeroEqVolumeVariables<Traits, CompositionalVolVars>;
public:
    using type = RANSNCVolumeVariables<Traits, RANSVolVars>;
};

//! The specific vtk output fields
SET_PROP(ZeroEqNC, VtkOutputFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    using SinglePhaseVtkOutputFields = RANSVtkOutputFields<FVGridGeometry>;
public:
    using type = RANSNCVtkOutputFields<FVGridGeometry, FluidSystem, phaseIdx, SinglePhaseVtkOutputFields>;
};

//////////////////////////////////////////////////////////////////////////
// Property values for non-isothermal multi-component RANS model
//////////////////////////////////////////////////////////////////////////

//! The type tags for the single-phase, multi-component non-isothermal RANS models
NEW_TYPE_TAG(RANSNCNI, INHERITS_FROM(NavierStokesNCNI));
NEW_TYPE_TAG(ZeroEqNCNI, INHERITS_FROM(RANSNCNI));

//! Set the volume variables property
SET_PROP(ZeroEqNCNI, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = FreeflowNCVolumeVariables<Traits>;
    using RANSVolVars = ZeroEqVolumeVariables<Traits, NSVolVars>;
public:
    using type = RANSNCVolumeVariables<Traits, RANSVolVars>;
};

//! The specific vtk output fields
SET_PROP(ZeroEqNCNI, VtkOutputFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    using NavierStokesFields = NavierStokesVtkOutputFields<FVGridGeometry>;
    using RANSFields = RANSVtkOutputFields<FVGridGeometry>;
    using SinglePhaseVtkOutputFields = RANSNonIsothermalVtkOutputFields<NavierStokesFields, RANSFields>;
public:
    using type = RANSNCVtkOutputFields<FVGridGeometry, FluidSystem, phaseIdx, SinglePhaseVtkOutputFields>;
};

// \}
} // end namespace Properties
} // end namespace Dumux

#endif
