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
 * \ingroup SolidEnergyModel
 * \brief The energy balance equation for a porous solid
 *
 * The energy balance is described by the following equation:
 \f[
   \frac{ \partial n c_p \varrho T}{\partial t}
   - \text{div} \left\lbrace \lambda_\text{pm} \textbf{grad} T \right\rbrace = q,
 \f]
 * where \f$n\f$ is the volume fraction of the conducting material, \f$c_p\f$ its specific heat capacity,
 * \f$\varrho\f$ its density, \f$T\f$ the temperature, and \f$\lambda\f$ the heat conductivity of the porous solid.
*/

#ifndef DUMUX_SOLID_ENERGY_MODEL_HH
#define DUMUX_SOLID_ENERGY_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/properties.hh>

#include <dumux/porousmediumflow/nonisothermal/iofields.hh>
#include <dumux/porousmediumflow/solidenergy/volumevariables.hh>
#include <dumux/porousmediumflow/solidenergy/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup SolidEnergyModel
 * \brief The indices
 */
struct SolidEnergyIndices
{
    static constexpr int temperatureIdx = 0;
    static constexpr int energyEqIdx = 0;
};

/*!
 * \ingroup SolidEnergyModel
 * \brief The energy balance equation for a porous solid
 */
struct SolidEnergyModelTraits
{
    using Indices = SolidEnergyIndices;
    static constexpr int numEq() { return 1; }
    static constexpr int numFluidPhases() { return 0; }
    static constexpr int numFluidComponents() { return 0; }
    static constexpr int numSolidComponents() { return 1; }

    static constexpr bool enableAdvection() { return false; }
    static constexpr bool enableMolecularDiffusion() { return false; }
    static constexpr bool enableEnergyBalance() { return true; }
};

/*!
 * \ingroup SolidEnergyModel
 * \brief The volume variable traits
 */
template<class PV, class SSY, class SST, class MT>
struct SolidEnergyVolumeVariablesTraits
{
    using PrimaryVariables = PV;
    using SolidSystem = SSY;
    using SolidState = SST;
    using ModelTraits = MT;
};

struct ThermalConductivitySolidEnergy
{
    /*!
     * \brief Relation for a simple effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$
     *
     * \param volVars volume variables
     * \param spatialParams spatial parameters
     * \param element element (to be passed to spatialParams)
     * \param fvGeometry fvGeometry (to be passed to spatialParams)
     * \param scv scv (to be passed to spatialParams)
     *
     * \return effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$
     */
    template<class VolumeVariables, class SpatialParams, class Element, class FVGeometry>
    static typename VolumeVariables::PrimaryVariables::value_type
    effectiveThermalConductivity(const VolumeVariables& volVars,
                                 const SpatialParams& spatialParams,
                                 const Element& element,
                                 const FVGeometry& fvGeometry,
                                 const typename FVGeometry::SubControlVolume& scv)
    {
        return volVars.solidThermalConductivity()*(1.0-volVars.porosity());
    }
};

// \{
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the fully implicit tracer model.
// Create new type tags
namespace TTag {
struct SolidEnergy { using InheritsFrom = std::tuple<PorousMediumFlow>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// Properties
///////////////////////////////////////////////////////////////////////////

//! set the model traits
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::SolidEnergy> { using type = SolidEnergyModelTraits; };

//! set the local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::SolidEnergy> { using type = SolidEnergyLocalResidual<TypeTag>; };

//! Set the vtk output fields specific to this model
template<class TypeTag>
struct IOFields<TypeTag, TTag::SolidEnergy> { using type = EnergyIOFields<>; };

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::SolidEnergy>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
    using SST = GetPropType<TypeTag, Properties::SolidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    using Traits = SolidEnergyVolumeVariablesTraits<PV, SSY, SST, MT>;
public:
    using type = SolidEnergyVolumeVariables<Traits>;
};

//! Use the average for effective conductivities
template<class TypeTag>
struct ThermalConductivityModel<TypeTag, TTag::SolidEnergy>
{ using type = ThermalConductivitySolidEnergy; };

} // end namespace Properties
// \}
} // end namespace Dumux

#endif
