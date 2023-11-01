// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPModel
 * \brief Adaption of the fully implicit scheme to the two-phase flow model.
 *
 * This model implements two-phase flow of two immiscible fluids
 * \f$\alpha \in \{ w, n \}\f$ using a standard multi-phase Darcy
 * approach as the equation for the conservation of momentum.
 * For details on Darcy's law see dumux/flux/darcyslaw.hh.
 *
 * By inserting Darcy's law into the equations for the conservation of the
 * phase mass, one gets
 \f[
 \frac{\partial (\phi \varrho_\alpha S_\alpha) }{\partial t}
 -
 \nabla \cdot \left\{
 \varrho_\alpha \frac{k_{r\alpha}}{\mu_\alpha} \mathbf{K} \left(\nabla  p_\alpha - \varrho_{\alpha} \mathbf{g} \right)
 \right\} - q_\alpha = 0,
 \f]
 *where:
 * * \f$ \phi \f$ is the porosity of the porous medium,
 * * \f$ S_\alpha \f$ represents the saturation of phase \f$ \alpha \f$,
 * * \f$ \varrho_\alpha \f$ is the mass density of phase \f$ \alpha \f$,
 * * \f$ k_{r\alpha} \f$ is the relative permeability of phase \f$ \alpha \f$,
 * * \f$ \mu_\alpha \f$ is the dynamic viscosity of phase \f$ \alpha \f$,
 * * \f$ \mathbf{K} \f$ is the intrinsic permeability tensor,
 * * \f$ p_\alpha \f$ is the pressure of phase \f$ \alpha \f$,
 * * \f$ \mathbf{g} \f$ is the gravitational acceleration vector,
 * * \f$ q_\alpha \f$ is a source or sink term.
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

#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivity/somerton.hh>
#include <dumux/material/fluidstates/immiscible.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/model.hh>
#include <dumux/porousmediumflow/nonisothermal/indices.hh>
#include <dumux/porousmediumflow/nonisothermal/iofields.hh>

#include "formulation.hh"
#include "indices.hh"
#include "volumevariables.hh"
#include "iofields.hh"
#include "saturationreconstruction.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPModel
 * \brief Specifies a number properties of two-phase models.
 */
template<TwoPFormulation formulation>
struct TwoPModelTraits
{
    using Indices = TwoPIndices;

    static constexpr TwoPFormulation priVarFormulation()
    { return formulation; }

    static constexpr int numEq() { return 2; }
    static constexpr int numFluidPhases() { return 2; }
    static constexpr int numFluidComponents() { return 2; }

    static constexpr bool enableAdvection() { return true; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return false; }
};

/*!
 * \ingroup TwoPModel
 * \brief Traits class for the two-phase model.
 *
 * \tparam PV The type used for primary variables
 * \tparam FSY The fluid system type
 * \tparam FST The fluid state type
 * \tparam SSY The solid system type
 * \tparam SST The solid state type
 * \tparam PT The type used for permeabilities
 * \tparam MT The model traits
 * \tparam SR The class used for reconstruction of
 *            nonwetting phase saturations in scvs
 */
template<class PV, class FSY, class FST, class SSY, class SST, class PT, class MT, class SR>
struct TwoPVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using FluidSystem = FSY;
    using FluidState = FST;
    using SolidSystem = SSY;
    using SolidState = SST;
    using PermeabilityType = PT;
    using ModelTraits = MT;
    using SaturationReconstruction = SR;
};

// necessary for models derived from 2p
class TwoPIOFields;

////////////////////////////////
// properties
////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tag for the isothermal two-phase model
struct TwoP { using InheritsFrom = std::tuple<PorousMediumFlow>; };

//! The type tag for the non-isothermal two-phase model
struct TwoPNI { using InheritsFrom = std::tuple<TwoP>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// properties for the isothermal two-phase model
///////////////////////////////////////////////////////////////////////////
 //!< Set the default formulation to pwsn
template<class TypeTag>
struct Formulation<TypeTag, TTag::TwoP>
{ static constexpr auto value = TwoPFormulation::p0s1; };

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::TwoP> { using type = ImmiscibleLocalResidual<TypeTag>; };         //!< Use the immiscible local residual operator for the 2p model

//! The base model traits class
template<class TypeTag>
struct BaseModelTraits<TypeTag, TTag::TwoP> { using type = TwoPModelTraits<getPropValue<TypeTag, Properties::Formulation>()>; };
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TwoP> { using type = GetPropType<TypeTag, Properties::BaseModelTraits>; }; //!< default the actually used traits to the base traits

//! Set the vtk output fields specific to the twop model
template<class TypeTag>
struct IOFields<TypeTag, TTag::TwoP> { using type = TwoPIOFields; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoP>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using DM = typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of nonwetting phase saturations
    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;

    using Traits = TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;
public:
    using type = TwoPVolumeVariables<Traits>;
};

//! The two-phase model uses the immiscible fluid state
template<class TypeTag>
struct FluidState<TypeTag, TTag::TwoP>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
    using type = ImmiscibleFluidState<Scalar, FluidSystem>;
};

////////////////////////////////////////////////////////
// properties for the non-isothermal two-phase model
////////////////////////////////////////////////////////

//! The non-isothermal model traits class
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::TwoPNI> { using type = PorousMediumFlowNIModelTraits<GetPropType<TypeTag, Properties::BaseModelTraits>>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TwoPNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using PT = typename GetPropType<TypeTag, Properties::SpatialParams>::PermeabilityType;
    using DM = typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod;
    static constexpr bool enableIS = getPropValue<TypeTag, Properties::EnableBoxInterfaceSolver>();
    // class used for scv-wise reconstruction of nonwetting phase saturations
    using SR = TwoPScvSaturationReconstruction<DM, enableIS>;
    using BaseTraits = TwoPVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT, SR>;

    using ETCM = GetPropType< TypeTag, Properties:: ThermalConductivityModel>;

    template<class BaseTraits, class ETCM>
    struct NITraits : public BaseTraits { using EffectiveThermalConductivityModel = ETCM; };

public:
    using type = TwoPVolumeVariables<NITraits<BaseTraits, ETCM>>;
};

//! Set the vtk output fields specific to the non-isothermal twop model
template<class TypeTag>
struct IOFields<TypeTag, TTag::TwoPNI> { using type = EnergyIOFields<TwoPIOFields>; };

//! Somerton is used as default model to compute the effective thermal heat conductivity
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::TwoPNI>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = ThermalConductivitySomerton<Scalar>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
