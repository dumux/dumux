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
 * \ingroup NavierStokesModel
 *
 * \brief A single-phase, isothermal Navier-Stokes model
 *
 * This model implements a single-phase, isothermal Navier-Stokes model, solving the <B> momentum balance equation </B>
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
 * The <B> mass balance equation </B>
 * \f[
       \frac{\partial \varrho}{\partial t} + \nabla \cdot (\varrho \textbf{v}) - q = 0
 * \f]
 *
 * closes the system.
 *
 *
 * So far, only the staggered grid spatial discretization (for structured grids) is available.
 */

#ifndef DUMUX_NAVIERSTOKES_MODEL_HH
#define DUMUX_NAVIERSTOKES_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/nonisothermal/model.hh>

#include "localresidual.hh"
#include "volumevariables.hh"
#include "fluxvariables.hh"
#include "fluxvariablescache.hh"
#include "indices.hh"
#include "vtkoutputfields.hh"

#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

// \{
///////////////////////////////////////////////////////////////////////////
// properties for the single-phase Navier-Stokes model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal Navier-Stokes model
NEW_TYPE_TAG(NavierStokes, INHERITS_FROM(FreeFlow));

//! The type tag for the corresponding non-isothermal model
NEW_TYPE_TAG(NavierStokesNI, INHERITS_FROM(NavierStokes, NavierStokesNonIsothermal));

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
SET_INT_PROP(NavierStokes, NumPhases, 1); //!< The number of phases in the 1p model is 1
SET_INT_PROP(NavierStokes, NumComponents, 1); //!< The number of components in the 1p model is 1
SET_INT_PROP(NavierStokes, PhaseIdx, 0); //!< The default phase index

SET_BOOL_PROP(NavierStokes, EnableAdvection, true); //!< Enable advection
SET_BOOL_PROP(NavierStokes, EnableMolecularDiffusion, false); //!< The one-phase model has no molecular diffusion
SET_BOOL_PROP(NavierStokes, EnableEnergyBalance, false); //!< The model is isothermal
SET_BOOL_PROP(NavierStokes, EnableInertiaTerms, true); //!< Consider inertia terms by default
SET_BOOL_PROP(NavierStokes, NormalizePressure, true); //!< Normalize the pressure term in the momentum balance by default

/*!
* \brief The number of equations.
*         There are as many momentum balance equations as dimensions
*         and one mass balance equation.
*/
SET_PROP(NavierStokes, NumEq)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr auto dim = GridView::dimension;
public:
    static constexpr int value = dim + 1;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(NavierStokes, FluidState){
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = Dumux::ImmiscibleFluidState<Scalar, FluidSystem>;
};

//! The local residual
SET_TYPE_PROP(NavierStokes, LocalResidual, NavierStokesResidual<TypeTag>);

//! The volume variables
SET_TYPE_PROP(NavierStokes, VolumeVariables, NavierStokesVolumeVariables<TypeTag>);

//! The flux variables
SET_TYPE_PROP(NavierStokes, FluxVariables, NavierStokesFluxVariables<TypeTag>);

//! The flux variables cache class, by default the one for free flow
SET_TYPE_PROP(NavierStokes, FluxVariablesCache, FreeFlowFluxVariablesCache<TypeTag>);

//! The indices required by the isothermal single-phase model
SET_PROP(NavierStokes, Indices)
{
private:
    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int dim = GET_PROP_TYPE(TypeTag, GridView)::dimension;
public:
    using type = NavierStokesIndices<dim, numEq>;
};

//! The specific vtk output fields
SET_TYPE_PROP(NavierStokes, VtkOutputFields, NavierStokesVtkOutputFields<TypeTag>);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////
//! The indices required by the isothermal single-phase model
SET_PROP(NavierStokesNI, IsothermalIndices)
{
private:
    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int dim = GET_PROP_TYPE(TypeTag, GridView)::dimension;
public:
    using type = NavierStokesIndices<dim, numEq>;
};

//! The specific isothermal vtk output fields
SET_TYPE_PROP(NavierStokesNI, IsothermalVtkOutputFields, NavierStokesVtkOutputFields<TypeTag>);

//! The number of equations for the isothermal model
SET_PROP(NavierStokesNI, IsothermalNumEq)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr auto dim = GridView::dimension;
public:
    static constexpr int value = dim + 1;
};
 // \}
}

} // end namespace

#endif // DUMUX_NAVIERSTOKES_MODEL_HH
