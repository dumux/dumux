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
 * \ingroup TwoPModel
 * \brief Adaption of the fully implicit scheme to the two-phase flow model.
 *
 * This model implements two-phase flow of two immiscible fluids
 * \f$\alpha \in \{ w, n \}\f$ using a standard multiphase Darcy
 * approach as the equation for the conservation of momentum, i.e.
 \f[
 v_\alpha = - \frac{k_{r\alpha}}{\mu_\alpha} \textbf{K}
 \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} {\textbf g} \right)
 \f]
 *
 * By inserting this into the equation for the conservation of the
 * phase mass, one gets
 \f[
 \phi \frac{\partial \varrho_\alpha S_\alpha}{\partial t}
 -
 \text{div} \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\textbf{grad}\, p_\alpha - \varrho_{\alpha} \mbox{\bf g} \right)
 \right\} - q_\alpha = 0 \;,
 \f]
 *
 * All equations are discretized using a vertex-centered finite volume (box)
 * or cell-centered finite volume scheme as spatial
 * and the implicit Euler method as time discretization.
 *
 * By using constitutive relations for the capillary pressure \f$p_c =
 * p_n - p_w\f$ and relative permeability \f$k_{r\alpha}\f$ and taking
 * advantage of the fact that \f$S_w + S_n = 1\f$, the number of
 * unknowns can be reduced to two. Currently the model supports
 * choosing either \f$p_w\f$ and \f$S_n\f$ or \f$p_n\f$ and \f$S_w\f$
 * as primary variables. The formulation which ought to be used can be
 * specified by setting the <tt>Formulation</tt> property to either
 * <tt>TwoPFormulation::pwsn</tt> or <tt>TwoPFormulation::pnsw</tt>. By
 * default, the model uses \f$p_w\f$ and \f$S_n\f$.
 */

#ifndef DUMUX_TWOP_MODEL_HH
#define DUMUX_TWOP_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/spatialparams/fv.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux
{

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the isothermal two-phase model
NEW_TYPE_TAG(TwoP, INHERITS_FROM(PorousMediumFlow));
//! The type tag for the non-isothermal two-phase model
NEW_TYPE_TAG(TwoPNI, INHERITS_FROM(TwoP, NonIsothermal));

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal two-phase model
///////////////////////////////////////////////////////////////////////////
SET_INT_PROP(TwoP, NumEq, 2);                                                 //!< Set the number of equations to 2
SET_INT_PROP(TwoP, NumPhases, 2);                                             //!< The number of phases in the 2p model is 2
SET_INT_PROP(TwoP, NumComponents, 2);                                         //!< The number of components in the 2p model is 2
SET_INT_PROP(TwoP, Formulation, TwoPFormulation::pwsn);                       //!< Set the default formulation to pWsN
SET_BOOL_PROP(TwoP, EnableAdvection, true);                                   //!< Enable advection
SET_BOOL_PROP(TwoP, EnableMolecularDiffusion, false);                         //!< The two-phase model has no molecular diffusion
SET_BOOL_PROP(TwoP, EnableEnergyBalance, false);                              //!< Isothermal model (non-isothermal type tag is below)
SET_TYPE_PROP(TwoP, LocalResidual, ImmiscibleLocalResidual<TypeTag>);         //!< Use the immiscible local residual operator for the 2p model
SET_TYPE_PROP(TwoP, VolumeVariables, TwoPVolumeVariables<TypeTag>);           //!< the VolumeVariables property
SET_TYPE_PROP(TwoP, SpatialParams, FVSpatialParams<TypeTag>);                 //!< The spatial parameters. Use FVSpatialParams by default.

//! Set the vtk output fields specific to the twop model
SET_TYPE_PROP(TwoP, VtkOutputFields, TwoPVtkOutputFields<typename GET_PROP_TYPE(TypeTag, Indices)>);

//! The indices required by the isothermal 2p model
SET_PROP(TwoP, Indices)
{
private:
    using Fluidsystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int formulation = GET_PROP_VALUE(TypeTag, Formulation);
public:
    using type = TwoPIndices<Fluidsystem, formulation, 0>;
};

//! The two-phase model uses the immiscible fluid state
SET_PROP(TwoP, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = ImmiscibleFluidState<Scalar, FluidSystem>;
};

////////////////////////////////////////////////////////
// properties for the non-isothermal two-phase model
////////////////////////////////////////////////////////
SET_INT_PROP(TwoPNI, IsothermalNumEq, 2);                                         //!< set isothermal NumEq
SET_TYPE_PROP(TwoPNI, IsothermalVolumeVariables, TwoPVolumeVariables<TypeTag>);   //!< set isothermal VolumeVariables
SET_TYPE_PROP(TwoPNI, IsothermalLocalResidual, ImmiscibleLocalResidual<TypeTag>); //!< set isothermal LocalResidual

//! Set the vtk output fields specific to the isothermal twop model
SET_TYPE_PROP(TwoPNI, IsothermalVtkOutputFields, TwoPVtkOutputFields<typename GET_PROP_TYPE(TypeTag, Indices)>);

//! Set isothermal Indices
SET_PROP(TwoPNI, IsothermalIndices)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int formulation = GET_PROP_VALUE(TypeTag, Formulation);
public:
    using type = TwoPIndices<FluidSystem, formulation, 0>;
};

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(TwoPNI, ThermalConductivityModel)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
public:
    using type = ThermalConductivitySomerton<Scalar, Indices>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
