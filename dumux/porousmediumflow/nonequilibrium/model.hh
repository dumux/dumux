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
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief This is the specialization that is able to capture non-equilibrium mass and / or energy transfer.
 * \todo DocMe
 */
#ifndef DUMUX_NONEQUILIBRIUM_MODEL_HH
#define DUMUX_NONEQUILIBRIUM_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/dimensionlessnumbers.hh>
#include <dumux/material/fluidstates/nonequilibrium.hh>

#include <dumux/discretization/fourierslawnonequilibrium.hh>
#include <dumux/discretization/fourierslawnonequilibrium.hh>
#include <dumux/porousmediumflow/nonequilibrium/thermal/localresidual.hh>
#include <dumux/porousmediumflow/nonequilibrium/localresidual.hh>

#include "localresidual.hh"
#include "indices.hh"
#include "gridvariables.hh"
#include "vtkoutputfields.hh"

/*!
 * \ingroup \ingroup PorousmediumNonEquilibriumModel
 * \brief Defines the properties required for non-equilibrium models
 */
namespace Dumux
{

/*!
 * \ingroup PorousmediumNonEquilibriumModel
 * \brief Specifies a number properties of porous-medium flow non-equilibrium models.
 *
 * \tparam ET The model traits of the underlying model assuming equilibrium
 * \tparam chem Boolean to indicate if chemical non-equilibrium is to be considered
 * \tparam therm Boolean to indicate if thermal non-equilibrium is to be considered
 * \tparam numEF Number of energy balance equations to be solved for the fluids
 * \tparam numES Number of energy balance equations to be solved for the solids
 */
template<class ET, bool chem, bool therm, int numEF, int numES>
struct NonEquilibriumModelTraits : public ET
{
    static constexpr int numEq() { return numEnergyEqFluid()+numEnergyEqSolid()+numTransportEq()+ET::numConstraintEq(); }
    static constexpr int numTransportEq() { return chem ? ET::numPhases()*ET::numComponents() : ET::numComponents(); }

    static constexpr int numEnergyEqFluid() { return therm ? numEF : 0; }
    static constexpr int numEnergyEqSolid() { return therm ? numES : 0; }

    static constexpr bool enableEnergyBalance() { return ET::enableEnergyBalance() || therm; }
    static constexpr bool enableThermalNonEquilibrium() { return therm; }
    static constexpr bool enableChemicalNonEquilibrium() { return chem; }
};

namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
NEW_TYPE_TAG(NonEquilibrium);

/////////////////////////////////////////////////
// Properties for the non-equilibrium mpnc model
/////////////////////////////////////////////////

//! Set the model traits
SET_PROP(NonEquilibrium, ModelTraits)
{
private:
    using EquiTraits = typename GET_PROP_TYPE(TypeTag, EquilibriumModelTraits);
    static constexpr bool enableTNE = GET_PROP_VALUE(TypeTag, EnableThermalNonEquilibrium);
    static constexpr bool enableCNE = GET_PROP_VALUE(TypeTag, EnableChemicalNonEquilibrium);
    static constexpr int numEF = GET_PROP_VALUE(TypeTag, NumEnergyEqFluid);
    static constexpr int numES = GET_PROP_VALUE(TypeTag, NumEnergyEqSolid);
public:
    using type = NonEquilibriumModelTraits<EquiTraits, enableCNE, enableTNE, numEF, numES>;
};

//! Per default, we consider both thermal and chemical non-equilibrium
SET_BOOL_PROP(NonEquilibrium, EnableThermalNonEquilibrium, true);
SET_BOOL_PROP(NonEquilibrium, EnableChemicalNonEquilibrium, true);

//! Default values for the number of energy balance equations
SET_INT_PROP(NonEquilibrium, NumEnergyEqSolid, 1);
SET_PROP(NonEquilibrium, NumEnergyEqFluid)
{
private:
    using EquiTraits = typename GET_PROP_TYPE(TypeTag, EquilibriumModelTraits);
public:
    static const int value = EquiTraits::numPhases();
};

SET_TYPE_PROP(NonEquilibrium, EnergyLocalResidual, EnergyLocalResidualNonEquilibrium<TypeTag, GET_PROP_VALUE(TypeTag, NumEnergyEqFluid)>);
SET_TYPE_PROP(NonEquilibrium, LocalResidual, NonEquilibriumLocalResidual<TypeTag>);
SET_TYPE_PROP(NonEquilibrium, HeatConductionType, FouriersLawNonEquilibrium<TypeTag>);

//! indices for non-isothermal models
SET_PROP(NonEquilibrium, Indices)
{
private:
    using EquilibriumIndices = typename GET_PROP_TYPE(TypeTag, EquilibriumIndices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    static constexpr int numEF = ModelTraits::numEnergyEqFluid();
    static constexpr int numES = ModelTraits::numEnergyEqSolid();
    static constexpr int numEq = ModelTraits::numEq();
public:
    using type = NonEquilbriumIndices<EquilibriumIndices, FluidSystem, numEF, numES, numEq, 0>;
};

SET_PROP(NonEquilibrium, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
     using type = NonEquilibriumFluidState<Scalar, FluidSystem>;
};

//! The grid variables
SET_TYPE_PROP(NonEquilibrium, GridVariables, NonEquilibriumGridVariables<TypeTag>);

//! indices for non-isothermal models
SET_PROP(NonEquilibrium, VtkOutputFields)
{
private:
    using EquilibriumVtkOutputFields = typename GET_PROP_TYPE(TypeTag, EquilibriumVtkOutputFields);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    static constexpr int numEF = GET_PROP_TYPE(TypeTag, ModelTraits)::numEnergyEqFluid();
    static constexpr int numES = GET_PROP_TYPE(TypeTag, ModelTraits)::numEnergyEqSolid();
public:
     using type = NonEquilibriumVtkOutputFields<EquilibriumVtkOutputFields, FluidSystem, numEF, numES>;
};

SET_PROP(NonEquilibrium, NusseltFormulation)
{
    private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using DimLessNum =  DimensionlessNumbers<Scalar>;
    public:
    static constexpr int value = DimLessNum::NusseltFormulation::WakaoKaguei;
};

/*!
 * \brief Set the default formulation for the sherwood correlation
 *        Other possible parametrizations can be found in the dimensionlessnumbers
 */
SET_PROP(NonEquilibrium, SherwoodFormulation )
{
    private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using DimLessNum =  DimensionlessNumbers<Scalar>;
    public:
    static constexpr int value = DimLessNum::SherwoodFormulation::WakaoKaguei;
};

} //end namespace Properties
} //end namespace Dumux

#endif
