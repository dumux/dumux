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
 * \ingroup NavierStokesNCModel
 *
 * \brief A single-phase, multi-component isothermal Navier-Stokes model
 *
 * This model implements a single-phase, multi-component isothermal Navier-Stokes model, solving the <B> momentum balance equation </B>
 * \f[
 \frac{\partial (\varrho \textbf{v})}{\partial t} + \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\textup{T}}) = \nabla \cdot (\mu (\nabla \textbf{v} + \nabla \textbf{v}^{\textup{T}}))
     - \nabla p + \varrho \textbf{g} - \textbf{f}
 * \f]
 * By setting the property <code>EnableInertiaTerms</code> to <code>false</code> the Stokes
 * equation can be solved. In this case the term
 * \f[
 *    \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\textup{T}})
 * \f]
 * is neglected.
 *
 *
 * The system is closed by a <B> component mass/mole balance equation </B> for each component \f$\kappa\f$:
 * \f[
 *    \frac{\partial \left(\varrho X^\kappa\right)}{\partial t}
 *    + \nabla \cdot \left( \varrho {\boldsymbol{v}} X^\kappa
 *    - D^\kappa \varrho \frac{M^\kappa}{M} \textbf{grad}\, x^\kappa \right)
 *    - q^\kappa = 0
 * \f]
 *
 * Alternatively, one component balance equation can be replace by a <B> total mass/mole balance equation </B>:
 *
 * \f[
 *    \frac{\partial \varrho_g}{\partial t}
 *    + \nabla \cdot \left(
 *        \varrho {\boldsymbol{v}}
 *        - \sum_\kappa D^\kappa \varrho \frac{M^\kappa}{M} \textbf{grad}\, x^\kappa
 *      \right)
 *    - q = 0
 * \f]
 *
 * So far, only the staggered grid spatial discretization (for structured grids) is available.
 */

#ifndef DUMUX_NAVIERSTOKES_NC_MODEL_HH
#define DUMUX_NAVIERSTOKES_NC_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/nonisothermal/model.hh>
#include <dumux/discretization/fickslaw.hh>

#include "volumevariables.hh"
#include "indices.hh"
#include "localresidual.hh"
#include "fluxvariables.hh"
#include "vtkoutputfields.hh"

#include <dumux/assembly/staggeredlocalresidual.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include <dumux/material/fluidstates/compositional.hh>

namespace Dumux
{

///////////////////////////////////////////////////////////////////////////
// properties for the single-phase, multi-component Navier-Stokes model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, multi-component isothermal Navier-Stokes model
NEW_TYPE_TAG(NavierStokesNC, INHERITS_FROM(NavierStokes));

//! The type tag for the single-phase, multi-component non-isothermal Navier-Stokes model
NEW_TYPE_TAG(NavierStokesNCNI, INHERITS_FROM(NavierStokesNC, NavierStokesNonIsothermal));

///////////////////////////////////////////////////////////////////////////
// default property values
///////////////////////////////////////////////////////////////////////////

/*!
* \brief Set the property for the number of components.
*
* We just forward the number from the fluid system
*
*/
SET_PROP(NavierStokesNC, NumComponents)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    static constexpr int value = FluidSystem::numComponents;
};

/*!
* \brief The number of equations.
*         There are as many momentum balance equations as dimensions
*         and as many balance equations as components.
*/
SET_PROP(NavierStokesNC, NumEq)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr auto dim = GridView::dimension;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    static constexpr int value = dim + FluidSystem::numComponents;
};

SET_INT_PROP(NavierStokesNC, PhaseIdx, 0); //!< Defines the phaseIdx
SET_BOOL_PROP(NavierStokesNC, EnableMolecularDiffusion, true); //!< Enable molecular diffusion
SET_BOOL_PROP(NavierStokesNC, UseMoles, false); //!< Defines whether molar (true) or mass (false) density is used
SET_INT_PROP(NavierStokesNC, ReplaceCompEqIdx, 0); //<! Set the ReplaceCompEqIdx to 0 by default

//! The local residual
SET_TYPE_PROP(NavierStokesNC, LocalResidual, NavierStokesNCResidual<TypeTag>);

//! The volume variables
SET_TYPE_PROP(NavierStokesNC, VolumeVariables, NavierStokesNCVolumeVariables<TypeTag>);

//! The flux variables
SET_TYPE_PROP(NavierStokesNC, FluxVariables, NavierStokesNCFluxVariables<TypeTag>);

//! The indices
SET_PROP(NavierStokesNC, Indices)
{
private:
    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int dim = GET_PROP_TYPE(TypeTag, GridView)::dimension;
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
public:
    using type = NavierStokesNCIndices<dim, numEq, phaseIdx, replaceCompEqIdx>;
};

//! The vtk output fields
SET_TYPE_PROP(NavierStokesNC, VtkOutputFields, NavierStokesNCVtkOutputFields<TypeTag>);

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(NavierStokesNC, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};

//! Use Fick's law for molecular diffusion per default
SET_TYPE_PROP(NavierStokesNC, MolecularDiffusionType, FicksLaw<TypeTag>);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////
//! The isothermal indices
SET_PROP(NavierStokesNCNI, IsothermalIndices)
{
private:
    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int dim = GET_PROP_TYPE(TypeTag, GridView)::dimension;
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
public:
    using type = NavierStokesNCIndices<dim, numEq, phaseIdx, replaceCompEqIdx>;
};

//! The isothermal vtk output fields
SET_TYPE_PROP(NavierStokesNCNI, IsothermalVtkOutputFields, NavierStokesNCVtkOutputFields<TypeTag>);

//! The number of equations for the isothermal model
SET_PROP(NavierStokesNCNI, IsothermalNumEq)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr auto dim = GridView::dimension;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    static constexpr int value = dim + FluidSystem::numComponents;
};

// \}
} // end namespace Properties
} // end namespace Dumux


#endif
