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
 * \ingroup NonEquilibriumModel
 * \brief This specifies models which are able to capture non-equilibrium mass and / or energy transfer.
 */

#ifndef DUMUX_NONEQUILIBRIUM_MODEL_HH
#define DUMUX_NONEQUILIBRIUM_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/dimensionlessnumbers.hh>
#include <dumux/material/fluidstates/nonequilibrium.hh>

#include <dumux/flux/fourierslawnonequilibrium.hh>
#include <dumux/porousmediumflow/nonequilibrium/thermal/localresidual.hh>
#include <dumux/porousmediumflow/nonequilibrium/localresidual.hh>

#include "localresidual.hh"
#include "indices.hh"
#include "gridvariables.hh"
#include "iofields.hh"

namespace Dumux {

/*!
 * \ingroup NonEquilibriumModel
 * \brief Specifies a number properties of porous-medium flow non-equilibrium models.
 *
 * \tparam ET The model traits of the underlying model assuming equilibrium
 * \tparam chem Boolean to indicate if chemical non-equilibrium is to be considered
 * \tparam therm Boolean to indicate if thermal non-equilibrium is to be considered
 * \tparam numEF Number of energy balance equations to be solved for the fluids
 * \tparam numES Number of energy balance equations to be solved for the solids
 * \tparam nf Formulation used for the computation of the nusselt number
 * \tparam sf Formulation used for the computation of the sherwood number
 */
template<class ET, bool chem, bool therm, int numEF, int numES, NusseltFormulation nf, SherwoodFormulation sf>
struct NonEquilibriumModelTraits : public ET
{
    static constexpr int numEq() { return numEnergyEqFluid()+numEnergyEqSolid()+numTransportEq()+ET::numConstraintEq(); }
    static constexpr int numTransportEq() { return chem ? ET::numFluidPhases()*ET::numFluidComponents() : ET::numFluidComponents(); }

    static constexpr int numEnergyEqFluid() { return therm ? numEF : 0; }
    static constexpr int numEnergyEqSolid() { return therm ? numES : 0; }
    static constexpr int numEnergyEq() { return numEnergyEqFluid()+numEnergyEqSolid(); }

    static constexpr bool enableEnergyBalance() { return ET::enableEnergyBalance() || therm; }
    static constexpr bool enableThermalNonEquilibrium() { return therm; }
    static constexpr bool enableChemicalNonEquilibrium() { return chem; }

    static constexpr NusseltFormulation nusseltFormulation() { return nf; }
    static constexpr SherwoodFormulation sherwoodFormulation() { return sf; }

    static_assert(!(ET::enableEnergyBalance() && therm), "It is not possible to use a nonisothermal model assuming local thermal equilibrium in combination with a model using thermal non-equilibrium");

    using Indices = NonEquilbriumIndices<typename ET::Indices, numEnergyEqFluid(), numEnergyEqSolid(), numEq()>;
};

namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
namespace TTag {
struct NonEquilibrium {};
}

/////////////////////////////////////////////////
// Properties for the non-equilibrium mpnc model
/////////////////////////////////////////////////

//! Set the model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::NonEquilibrium>
{
private:
    using EquiTraits = GetPropType<TypeTag, Properties::EquilibriumModelTraits>;
    static constexpr bool enableTNE = getPropValue<TypeTag, Properties::EnableThermalNonEquilibrium>();
    static constexpr bool enableCNE = getPropValue<TypeTag, Properties::EnableChemicalNonEquilibrium>();
    static constexpr int numEF = getPropValue<TypeTag, Properties::NumEnergyEqFluid>();
    static constexpr int numES = getPropValue<TypeTag, Properties::NumEnergyEqSolid>();
    static constexpr auto nf = getPropValue<TypeTag, Properties::NusseltFormulation>();
    static constexpr auto ns = getPropValue<TypeTag, Properties::SherwoodFormulation>();
public:
    using type = NonEquilibriumModelTraits<EquiTraits, enableCNE, enableTNE, numEF, numES, nf, ns>;
};

//! Per default, we consider both thermal and chemical non-equilibrium
template<class TypeTag>
struct EnableThermalNonEquilibrium<TypeTag, TTag::NonEquilibrium> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableChemicalNonEquilibrium<TypeTag, TTag::NonEquilibrium> { static constexpr bool value = true; };

//! Default values for the number of energy balance equations
template<class TypeTag>
struct NumEnergyEqSolid<TypeTag, TTag::NonEquilibrium> { static constexpr int value = 1; };
template<class TypeTag>
struct NumEnergyEqFluid<TypeTag, TTag::NonEquilibrium> { static constexpr int value = GetPropType<TypeTag, Properties::EquilibriumModelTraits>::numFluidPhases(); };

template<class TypeTag>
struct EnergyLocalResidual<TypeTag, TTag::NonEquilibrium> { using type = EnergyLocalResidualNonEquilibrium<TypeTag, getPropValue<TypeTag, Properties::NumEnergyEqFluid>()>; };
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::NonEquilibrium> { using type = NonEquilibriumLocalResidual<TypeTag>; };
template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::NonEquilibrium> { using type = FouriersLawNonEquilibrium<TypeTag>; };

template<class TypeTag>
struct FluidState<TypeTag, TTag::NonEquilibrium>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
public:
     using type = NonEquilibriumFluidState<Scalar, FluidSystem>;
};

//! The grid variables
template<class TypeTag>
struct GridVariables<TypeTag, TTag::NonEquilibrium> { using type = NonEquilibriumGridVariables<TypeTag>; };

//! indices for non-isothermal models
template<class TypeTag>
struct IOFields<TypeTag, TTag::NonEquilibrium>
{
private:
    using EquilibriumIOFields = GetPropType<TypeTag, Properties::EquilibriumIOFields>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
public:
     using type = NonEquilibriumIOFields<ModelTraits, EquilibriumIOFields>;
};

template<class TypeTag>
struct NusseltFormulation<TypeTag, TTag::NonEquilibrium>
{
public:
    static constexpr Dumux::NusseltFormulation value = Dumux::NusseltFormulation::WakaoKaguei;
};

/*!
 * \brief Set the default formulation for the sherwood correlation
 *        Other possible parametrizations can be found in the dimensionlessnumbers
 */
template<class TypeTag>
struct SherwoodFormulation<TypeTag, TTag::NonEquilibrium>
{
public:
    static constexpr Dumux::SherwoodFormulation value = Dumux::SherwoodFormulation::WakaoKaguei;
};

} //end namespace Properties
} //end namespace Dumux

#endif
