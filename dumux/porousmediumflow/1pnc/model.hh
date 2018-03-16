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
 * \ingroup OnePNCModel
 * \brief  Adaption of the fully implicit model to the one-phase n-component flow model.
 *
 * This model implements a one-phase flow of a compressible fluid, that consists
 * of n components, using a standard Darcy approach as the equation for the
 * conservation of momentum:
 \f[
 v = - \frac{\textbf K}{\mu}
 \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 \f]
 *
 * Gravity can be enabled or disabled via the property system.
 * By inserting this into the continuity equation, one gets
 \f[
 \phi\frac{\partial \varrho}{\partial t} - \text{div} \left\{
   \varrho \frac{\textbf K}{\mu}  \left(\textbf{grad}\, p - \varrho {\textbf g} \right)
 \right\} = q \;,
 \f]
 *
 * The transport of the components \f$\kappa \in \{ w, a, ... \}\f$ is described by the following equation:
 \f[
 \phi \frac{ \partial \varrho X^\kappa}{\partial t}
 - \text{div} \left\lbrace \varrho X^\kappa \frac{{\textbf K}}{\mu} \left( \textbf{grad}\, p -
 \varrho {\textbf g} \right)
 + \varrho D^\kappa_\text{pm} \frac{M^\kappa}{M_\alpha} \textbf{grad} x^\kappa \right\rbrace = q.
 \f]
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. useMoles is set to true by default.
 *
 * The primary variables are the pressure \f$p\f$ and the mole fraction of dissolved components \f$x^\kappa\f$.
 */

#ifndef DUMUX_1PNC_MODEL_HH
#define DUMUX_1PNC_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/fluidstates/compositional.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/compositional/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux
{
/*!
 * \ingroup OnePNCModel
 * \brief Specifies a number properties of models that
 *        consider a single-phase with multiple components.
 *
 * \ţparam nComp the number of components to be considered.
 */
template<int nComp>
struct OnePNCModelTraits
{
    static constexpr int numEq() { return nComp; }
    static constexpr int numPhases() { return 1; }
    static constexpr int numComponents() { return nComp; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return true; }
    static constexpr bool enableEnergyBalance() { return false; }
};

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the implicit the isothermal & non-isothermal one phase n component problems
NEW_TYPE_TAG(OnePNC, INHERITS_FROM(PorousMediumFlow));
NEW_TYPE_TAG(OnePNCNI, INHERITS_FROM(OnePNC, NonIsothermal));

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////

//! The model traits. Per default, we use the number of components of the fluid system.
SET_PROP(OnePNC, ModelTraits)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem));
public:
    using type = OnePNCModelTraits<FluidSystem::numComponents>;
};

//! Set as default that no component mass balance is replaced by the total mass balance
SET_PROP(OnePNC, ReplaceCompEqIdx)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem));
public:
    static constexpr auto value = FluidSystem::numComponents;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(OnePNC, FluidState){
    private:
        using Scalar =  typename GET_PROP_TYPE(TypeTag, Scalar);
        using FluidSystem =  typename GET_PROP_TYPE(TypeTag, FluidSystem);
    public:
        using type = Dumux::CompositionalFluidState<Scalar, FluidSystem>;
};

//! Use the model after Millington (1961) for the effective diffusivity
SET_TYPE_PROP(OnePNC, EffectiveDiffusivityModel,
              DiffusivityMillingtonQuirk<typename GET_PROP_TYPE(TypeTag, Scalar)>);

SET_INT_PROP(OnePNC, PhaseIdx, 0); //!< The default phase index
SET_TYPE_PROP(OnePNC, LocalResidual, CompositionalLocalResidual<TypeTag>); //!< The local residual function
SET_TYPE_PROP(OnePNC, VolumeVariables, OnePNCVolumeVariables<TypeTag>);   //!< the VolumeVariables property

//! Set the vtk output fields specific to this model
SET_PROP(OnePNC, VtkOutputFields)
{
private:
   using FluidSystem =  typename GET_PROP_TYPE(TypeTag, FluidSystem);
   static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
public:
    using type = OnePNCVtkOutputFields<FluidSystem, phaseIdx>;
};

//! The indices required by the isothermal single-phase model
SET_PROP(OnePNC, Indices)
{
private:
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
public:
    using type = OnePNCIndices<phaseIdx>;
};


///////////////////////////////////////////////////////////////////////////
// properties for the non-isothermal single phase model
///////////////////////////////////////////////////////////////////////////

//! the isothermal vtk output fields
SET_PROP(OnePNCNI, IsothermalVtkOutputFields)
{
private:
   using FluidSystem =  typename GET_PROP_TYPE(TypeTag, FluidSystem);
   static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
public:
    using type = OnePNCVtkOutputFields<FluidSystem, phaseIdx>;
};

SET_TYPE_PROP(OnePNCNI, IsothermalVolumeVariables, OnePNCVolumeVariables<TypeTag>);     //!< Vol vars of the isothermal model
SET_TYPE_PROP(OnePNCNI, IsothermalLocalResidual, CompositionalLocalResidual<TypeTag>);   //!< Local residual of the isothermal model
SET_TYPE_PROP(OnePNCNI,
              ThermalConductivityModel,
              ThermalConductivityAverage<typename GET_PROP_TYPE(TypeTag, Scalar)>); //!< Use the average for effective conductivities

//! Indices of the isothermal model
SET_PROP(OnePNCNI, IsothermalIndices)
{
private:
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
public:
    using type = OnePNCIndices<phaseIdx>;
};

//! model traits of the isothermal model.
SET_PROP(OnePNCNI, IsothermalModelTraits)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem));
public:
    using type = OnePNCModelTraits<FluidSystem::numComponents>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
