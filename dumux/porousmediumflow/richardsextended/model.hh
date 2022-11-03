// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup ExtendedRichardsModel
 * \brief This model implements a variant of the extended Richards'
 *        equation for quasi-twophase flow (see e.g. Vanderborght et al. 2017).
 *
 * In the unsaturated zone, Richards' equation
 \f[
 \frac{\partial\;\phi S_w \varrho_w}{\partial t}
 +
 \frac{\partial\;\phi (1-S_w)\varrho_n X_n^w}{\partial t}
 -
 \text{div} \left\lbrace
 \varrho_w \frac{k_{rw}}{\mu_w} \; \mathbf{K} \;
 \left( \text{\textbf{grad}}
 p_w - \varrho_w \textbf{g}
 \right)
 +
 {\bf D_{n, pm}^w} \varrho_n \text{grad}\, X^w_n
 \right\rbrace
 =
 q_w,
 \f]
 * is frequently used to
 * approximate the water distribution above the groundwater level.
 *
 * It can be derived from the two-phase equations, i.e.
 \f[
 \phi\frac{\partial S_\alpha \varrho_\alpha}{\partial t}
 -
 \text{div} \left\lbrace
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha}\; \mathbf{K} \;
 \left( \text{\textbf{grad}}
 p_\alpha - \varrho_\alpha \textbf{g}
 \right)
 \right\rbrace
 =
 q_\alpha,
 \f]
 * where \f$\alpha \in \{w, n\}\f$ is the fluid phase,
 * \f$\kappa \in \{ w, a \}\f$ are the components,
 * \f$\rho_\alpha\f$ is the fluid density, \f$S_\alpha\f$ is the fluid
 * saturation, \f$\phi\f$ is the porosity of the soil,
 * \f$k_{r\alpha}\f$ is the relative permeability for the fluid,
 * \f$\mu_\alpha\f$ is the fluid's dynamic viscosity, \f$\mathbf{K}\f$ is the
 * intrinsic permeability, \f$p_\alpha\f$ is the fluid pressure and
 * \f$g\f$ is the potential of the gravity field.
 *
 * In contrast to the full two-phase model, the Richards model assumes
 * gas as the nonwetting fluid and that it exhibits a much lower
 * viscosity than the (liquid) wetting phase. (For example at
 * atmospheric pressure and at room temperature, the viscosity of air
 * is only about \f$1\%\f$ of the viscosity of liquid water.) As a
 * consequence, the \f$\frac{k_{r\alpha}}{\mu_\alpha}\f$ term
 * typically is much larger for the gas phase than for the wetting
 * phase. For this reason, the Richards model assumes that
 * \f$\frac{k_{rn}}{\mu_n}\f$ is infinitely large. This implies that
 * the pressure of the gas phase is equivalent to the static pressure
 * distribution and that therefore, mass conservation only needs to be
 * considered for the wetting phase.
 *
 * The model thus chooses the absolute pressure of the wetting phase
 * \f$p_w\f$ as its only primary variable. The wetting phase
 * saturation is calculated using the inverse of the capillary
 * pressure, i.e.
 \f[
 S_w = p_c^{-1}(p_n - p_w)
 \f]
 * holds, where \f$p_n\f$ is a given reference pressure. Nota bene,
 * that the last step is assumes that the capillary
 * pressure-saturation curve can be uniquely inverted, so it is not
 * possible to set the capillary pressure to zero when using the
 * Richards model!
 */

#ifndef DUMUX_RICHARDSEXTENDED_MODEL_HH
#define DUMUX_RICHARDSEXTENDED_MODEL_HH

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>

#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivity/somerton.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidstates/immiscible.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>
#include <dumux/porousmediumflow/richards/balanceequationopts.hh>
#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/porousmediumflow/richards/velocityoutput.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "iofields.hh"
#include "localresidual.hh"

namespace Dumux {

/*!
 * \ingroup ExtendedRichardsModel
 * \brief Specifies a number properties of the extended Richards model.
 */
struct ExtendedRichardsModelTraits : public RichardsModelTraits
{
    using Indices = ExtendedRichardsIndices;

    static constexpr bool enableMolecularDiffusion() { return true; }
};

/*!
 * \ingroup RichardsModel
 * \brief Traits class for the Richards model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT, class DT, class EDM>
struct ExtendedRichardsVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using SolidSystem = SSY;
    using SolidState = SST;
    using PermeabilityType = PT;
    using ModelTraits = MT;
    using DiffusionType = DT;
    using EffectiveDiffusivityModel = EDM;
};

// \{
///////////////////////////////////////////////////////////////////////////
// properties for the isothermal Richards model.
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit isothermal one-phase two-component problems
// Create new type tags
namespace TTag {
struct ExtendedRichards { using InheritsFrom = std::tuple<Richards>; };
struct ExtendedRichardsNI { using InheritsFrom = std::tuple<ExtendedRichards>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Properties values
//////////////////////////////////////////////////////////////////

//! The local residual operator
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::ExtendedRichards> { using type = ExtendedRichardsLocalResidual<TypeTag>; };

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::ExtendedRichards>
{
    using type = ExtendedRichardsIOFields;
};

//! The model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ExtendedRichards> { using type = ExtendedRichardsModelTraits; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ExtendedRichards>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    using Traits = ExtendedRichardsVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, DT, EDM>;
public:
    using type = ExtendedRichardsVolumeVariables<Traits>;
};

//! Use the model after Millington (1961) for the effective diffusivity
template<class TypeTag>
struct EffectiveDiffusivityModel<TypeTag, TTag::ExtendedRichards>
{ using type = DiffusivityMillingtonQuirk<GetPropType<TypeTag, Properties::Scalar>>; };

//! The primary variables vector for the richards model
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::ExtendedRichards>
{
private:
    using PrimaryVariablesVector = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>,
                                                     GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
public:
    using type = SwitchablePrimaryVariables<PrimaryVariablesVector, int>;
};

/////////////////////////////////////////////////////
// Property values for non-isothermal Richars model
/////////////////////////////////////////////////////

//! set non-isothermal model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::ExtendedRichardsNI>
{
private:
    using IsothermalTraits = ExtendedRichardsModelTraits;
public:
    using type = PorousMediumFlowNIModelTraits<IsothermalTraits>;
};

//! Set the vtk output fields specific to th non-isothermal model
template<class TypeTag>
struct IOFields<TypeTag, TTag::ExtendedRichardsNI>
{
    using type = EnergyIOFields<ExtendedRichardsIOFields>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::ExtendedRichardsNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using EDM = GetPropType<TypeTag, Properties::EffectiveDiffusivityModel>;
    using BaseTraits = ExtendedRichardsVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, DT, EDM>;

    using ETCM = GetPropType< TypeTag, Properties::ThermalConductivityModel>;
    template<class BaseTraits, class ETCM>
    struct NITraits : public BaseTraits { using EffectiveThermalConductivityModel = ETCM; };

public:
    using type = ExtendedRichardsVolumeVariables<NITraits<BaseTraits, ETCM>>;
};

//! Somerton is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::ExtendedRichardsNI>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = ThermalConductivitySomerton<Scalar>;
};

// \}
} // end namespace Properties
} // end namespace Dumux

#endif
