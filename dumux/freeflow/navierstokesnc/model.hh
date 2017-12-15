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
 *
 * \brief Base class for all models which use the one-phase,
 *        fully implicit model.
 *        Adaption of the fully implicit scheme to the one-phase flow model.
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


/*!
 * \ingroup NavierStokesModel TODO: doc me properly!
 * \brief A single-phase, isothermal flow model using the fully implicit scheme.
 *
 * Single-phase, isothermal flow model, which uses a standard Darcy approach as the
 * equation for the conservation of momentum:
 * \f[
 v = - \frac{\textbf K}{\mu}
 \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 * \f]
 *
 * and solves the mass continuity equation:
 * \f[
 \phi \frac{\partial \varrho}{\partial t} + \text{div} \left\lbrace
 - \varrho \frac{\textbf K}{\mu} \left( \textbf{grad}\, p -\varrho {\textbf g} \right) \right\rbrace = q,
 * \f]
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 * The model supports compressible as well as incompressible fluids.
 */

 namespace Dumux
 {
 ///////////////////////////////////////////////////////////////////////////
 // properties for the isothermal Navier-Stokes model
 ///////////////////////////////////////////////////////////////////////////
 namespace Properties {

 //////////////////////////////////////////////////////////////////
 // Type tags
 //////////////////////////////////////////////////////////////////

 //! The type tags for the implicit single-phase problems
 NEW_TYPE_TAG(NavierStokesNC, INHERITS_FROM(NavierStokes));
 NEW_TYPE_TAG(NavierStokesNCNI, INHERITS_FROM(NavierStokesNC, NavierStokesNonIsothermal));

 ///////////////////////////////////////////////////////////////////////////
 // default property values for the isothermal single phase model
 ///////////////////////////////////////////////////////////////////////////

 //! The number of equations
 SET_PROP(NavierStokesNC, NumEq)
 {
 private:
     using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
     static constexpr auto dim = GridView::dimension;
     using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
 public:
     static constexpr int value = dim + FluidSystem::numComponents;
 };

 SET_INT_PROP(NavierStokesNC, ReplaceCompEqIdx, 0);

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

 //! the VolumeVariables property
 SET_TYPE_PROP(NavierStokesNC, VolumeVariables, NavierStokesNCVolumeVariables<TypeTag>);
 SET_TYPE_PROP(NavierStokesNC, Indices, NavierStokesNCIndices<TypeTag>);

 SET_TYPE_PROP(NavierStokesNC, FluxVariables, NavierStokesNCFluxVariables<TypeTag>);

 SET_TYPE_PROP(NavierStokesNC, LocalResidual, NavierStokesNCResidual<TypeTag>);

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

 SET_BOOL_PROP(NavierStokesNC, EnableComponentTransport, true);

 //! The one-phase model has no molecular diffusion
 SET_BOOL_PROP(NavierStokesNC, EnableMolecularDiffusion, true);

 SET_TYPE_PROP(NavierStokesNC, MolecularDiffusionType, FicksLaw<TypeTag>);

 SET_BOOL_PROP(NavierStokesNC, UseMoles, false); //!< Defines whether molar (true) or mass (false) density is used

 SET_INT_PROP(NavierStokesNC, PhaseIdx, 0); //!< Defines the phaseIdx

 SET_TYPE_PROP(NavierStokesNC, VtkOutputFields, NavierStokesNCVtkOutputFields<TypeTag>); //! the vtk output fields

 // non-isothermal properties
 SET_TYPE_PROP(NavierStokesNCNI, IsothermalIndices, NavierStokesNCIndices<TypeTag>); //! the isothermal indices
 SET_TYPE_PROP(NavierStokesNCNI, IsothermalVtkOutputFields, NavierStokesNCVtkOutputFields<TypeTag>); //! the isothermal vtk output fields

 //! The number of equations
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
