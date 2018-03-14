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
 * \brief A single-phase, isothermal Reynolds-Averaged Navier-Stokes model
 *
 * This model implements a single-phase, isothermal Reynolds-Averaged
 * Navier-Stokes model, solving the <B> momentum balance equation </B>
 * \f[
 \frac{\partial (\varrho \textbf{v})}{\partial t} + \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\textup{T}}) = \nabla \cdot (\mu_\textrm{eff} (\nabla \textbf{v} + \nabla \textbf{v}^{\textup{T}}))
   - \nabla p + \varrho \textbf{g} - \textbf{f}
 * \f]
 * The effective viscosity is composed of the fluid and the eddy viscosity:
 * \f[
 *    \mu_\textrm{eff} = \mu + \mu_\textrm{t}
 * \f].
 *
 */

#ifndef DUMUX_RANS_MODEL_HH
#define DUMUX_RANS_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/nonisothermal/model.hh>

// #include <dumux/freeflow/navierstokes/fluxvariables.hh>
// #include <dumux/freeflow/navierstokes/fluxvariablescache.hh>
// #include <dumux/freeflow/navierstokes/indices.hh>
// #include <dumux/freeflow/navierstokes/localresidual.hh>
#include "volumevariables.hh"
#include "vtkoutputfields.hh"

#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

// \{
///////////////////////////////////////////////////////////////////////////
// properties for the single-phase Reynolds-Averaged Navier-Stokes model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal Reynolds-Averaged Navier-Stokes model
NEW_TYPE_TAG(RANS, INHERITS_FROM(NavierStokes));

//! The type tag for the corresponding non-isothermal model
NEW_TYPE_TAG(RANSNI, INHERITS_FROM(RANS, NavierStokesNonIsothermal));

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
// SET_INT_PROP(RANS, NumPhases, 1); //!< The number of phases in the 1p model is 1
// SET_INT_PROP(RANS, NumComponents, 1); //!< The number of components in the 1p model is 1
// SET_INT_PROP(RANS, PhaseIdx, 0); //!< The default phase index
//
// SET_BOOL_PROP(RANS, EnableAdvection, true); //!< Enable advection
// SET_BOOL_PROP(RANS, EnableMolecularDiffusion, false); //!< The one-phase model has no molecular diffusion
// SET_BOOL_PROP(RANS, EnableEnergyBalance, false); //!< The model is isothermal
SET_BOOL_PROP(RANS, EnableInertiaTerms, true); //!< Explicitly force the consideration of inertia terms by default
// SET_BOOL_PROP(RANS, NormalizePressure, true); //!< Normalize the pressure term in the momentum balance by default
//
// /*!
// * \brief The number of equations.
// *         There are as many momentum balance equations as dimensions
// *         and one mass balance equation.
// */
// SET_PROP(RANS, NumEq)
// {
// private:
//     using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
//     static constexpr auto dim = GridView::dimension;
// public:
//     static constexpr int value = dim + 1;
// };
//
// /*!
//  * \brief The fluid state which is used by the volume variables to
//  *        store the thermodynamic state. This should be chosen
//  *        appropriately for the model ((non-)isothermal, equilibrium, ...).
//  *        This can be done in the problem.
//  */
// SET_PROP(RANS, FluidState){
// private:
//     using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
//     using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
// public:
//     using type = Dumux::ImmiscibleFluidState<Scalar, FluidSystem>;
// };
//
// //! The local residual
// SET_TYPE_PROP(RANS, LocalResidual, RANSResidual<TypeTag>);

//! The volume variables
SET_TYPE_PROP(RANS, VolumeVariables, RANSVolumeVariables<TypeTag>);

// //! The flux variables
// SET_TYPE_PROP(RANS, FluxVariables, RANSFluxVariables<TypeTag>);
//
// //! The flux variables cache class, by default the one for free flow
// SET_TYPE_PROP(RANS, FluxVariablesCache, FreeFlowFluxVariablesCache<TypeTag>);
//
// //! The indices required by the isothermal single-phase model
// SET_TYPE_PROP(RANS, Indices, RANSIndices<TypeTag>);

//! The specific vtk output fields
SET_PROP(RANS, VtkOutputFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
public:
     using type = RANSVtkOutputFields<FVGridGeometry>;
};
//
// //////////////////////////////////////////////////////////////////
// // Property values for isothermal model required for the general non-isothermal model
// //////////////////////////////////////////////////////////////////
// //! The indices required by the isothermal single-phase model
// SET_TYPE_PROP(RANSNI, IsothermalIndices, RANSIndices<TypeTag>);
//
// //! The specific isothermal vtk output fields
// SET_TYPE_PROP(RANSNI, IsothermalVtkOutputFields, RANSVtkOutputFields<TypeTag>);
//
// //! The number of equations for the isothermal model
// SET_PROP(RANSNI, IsothermalNumEq)
// {
// private:
//     using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
//     static constexpr auto dim = GridView::dimension;
// public:
//     static constexpr int value = dim + 1;
// };
 // \}
}

} // end namespace

#endif // DUMUX_RANS_MODEL_HH
