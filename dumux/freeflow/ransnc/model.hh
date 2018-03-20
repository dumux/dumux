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
 * \copydoc RANSModel
 * *
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
 * The eddy diffusivity \$[ D_\text{t} \$] is related to the eddy viscosity by the turbulent
 * Schmidt number:
 * \f[ D_\text{t} = \frac{\nu_\text{t}}{\mathrm{Sc}_\text{t}} \f]
 *
 * So far, only the staggered grid spatial discretization (for structured grids) is available.
 */

#ifndef DUMUX_RANS_NC_MODEL_HH
#define DUMUX_RANS_NC_MODEL_HH


// TODO: remove unused headers
#include <dumux/common/properties.hh>

#include <dumux/freeflow/rans/model.hh>
#include <dumux/freeflow/rans/zeroeq/model.hh>
#include <dumux/freeflow/navierstokesnc/model.hh>
#include <dumux/freeflow/nonisothermal/model.hh>
#include <dumux/freeflow/nonisothermal/indices.hh>
#include <dumux/freeflow/nonisothermal/vtkoutputfields.hh>
#include <dumux/discretization/fickslaw.hh>
#include <dumux/discretization/fourierslaw.hh>

#include "volumevariables.hh"
#include "vtkoutputfields.hh"

#include <dumux/assembly/staggeredlocalresidual.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include <dumux/material/fluidstates/compositional.hh>

namespace Dumux {

/*!
 * \ingroup RANSModel
 * \brief Traits for the Reynolds-averaged Navier-Stokes multi-component model
 */
template<int dim, int nComp>
struct RANSNCModelTraits : NavierStokesNCModelTraits<dim, nComp>
{ };

///////////////////////////////////////////////////////////////////////////
// properties for the single-phase, multi-component Reynolds-averaged Navier-Stokes model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, multi-component isothermal Reynolds-averaged Navier-Stokes model
NEW_TYPE_TAG(RANSNC, INHERITS_FROM(RANS, NavierStokesNC));
NEW_TYPE_TAG(ZeroEqNC, INHERITS_FROM(RANSNC, NavierStokesNC));

// //! The type tag for the single-phase, multi-component non-isothermal Reynolds-averaged Navier-Stokes model
// NEW_TYPE_TAG(RANSNCNI, INHERITS_FROM(RANS, NavierStokesNC));

///////////////////////////////////////////////////////////////////////////
// default property values
///////////////////////////////////////////////////////////////////////////

//! The volume variables
SET_TYPE_PROP(RANSNC, VolumeVariables, RANSNCVolumeVariables<TypeTag>);

//! The specific vtk output fields
SET_PROP(RANSNC, VtkOutputFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
public:
     using type = RANSNCVtkOutputFields<FVGridGeometry, FluidSystem, phaseIdx>;
};

//////////////////////////////////////////////////////////////////////////
// Property values for non-isothermal multi-component Reynolds-averaged Navier-Stokes model
//////////////////////////////////////////////////////////////////////////

// //! The non-isothermal vtk output fields
// SET_PROP(RANSNCNI, VtkOutputFields)
// {
// private:
//     using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
//     using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
//     static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
//     using IsothermalFields = RANSNCVtkOutputFields<FVGridGeometry, FluidSystem, phaseIdx>;
// public:
//      using type = NavierStokesNonIsothermalVtkOutputFields<IsothermalFields>;
// };
//
// //! Use Fourier's Law as default heat conduction type
// SET_TYPE_PROP(RANSNCNI, HeatConductionType, FouriersLaw<TypeTag>);

// \}
} // end namespace Properties
} // end namespace Dumux


#endif
