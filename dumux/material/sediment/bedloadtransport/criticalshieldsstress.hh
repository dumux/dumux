// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTransport
 * \brief Formulas to calculate the critical Shields stress and the dimensionless grain diameter
 */
#ifndef DUMUX_MATERIAL_CRITICAL_SHIELDS_STRESS_HH
#define DUMUX_MATERIAL_CRITICAL_SHIELDS_STRESS_HH

namespace Dumux {

/*!
 * \ingroup BedloadTransport
 * \brief Calculate the dimensionless grain diameter.
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \param grainDiameter grain diameter \f$[m]\f$
 * \param gravity acceleration due to gravity \f$[m/s^2]\f$
 * \param kinematicViscosity kinematic viscosity \f$[m^2/s]\f$
 * \param specificDensity the grain density divided by the water density \f$[-]\f$
 *
 * \return dimensionless grain diameter \f$[-]\f$
 */
template<class Scalar>
const Scalar calculateDimensionlessGrainDiameter(const Scalar grainDiameter, const Scalar gravity, const Scalar kinematicViscosity, const Scalar specificDensity)
{
    using std::pow;

    return grainDiameter * pow(((specificDensity - 1) * gravity / (kinematicViscosity * kinematicViscosity)), 1/3.0);
}

/*!
 * \ingroup BedloadTransport
 * \brief Calculate the critical Shields stress.
 *
 * The critical Shields stress is also known as critical dimensionless bed shear stress.
 * At the moment there are two approaches implemented.
 * First, the one proposed by van Rijn in "Sediment transport by currents and waves" (1989).
 * Second, the approach proposed by Yalin and da Silva in "Fluvial Processes" (2001)
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \param dimensionlessGrainDiameter dimensionless grain diameter \f$[-]\f$
 *
 * \return critical Shields stress \f$[-]\f$
 */
template<class Scalar>
const Scalar calculateCriticalShieldsStress(const Scalar dimensionlessGrainDiameter)
{
    using std::pow;
    using std::exp;

    static const std::string formula = getParam<std::string>("Sediment.CriticalShieldsStressCalculation");
    if (formula=="vanRijn") {
        if (dimensionlessGrainDiameter <= 4.0) return 0.24 / dimensionlessGrainDiameter;
        else if (dimensionlessGrainDiameter <= 10.0) return 0.14 * pow(dimensionlessGrainDiameter, -0.64);
        else if (dimensionlessGrainDiameter <= 20.0) return 0.04 * pow(dimensionlessGrainDiameter, -0.1);
        else if (dimensionlessGrainDiameter <= 150.0) return 0.013 * pow(dimensionlessGrainDiameter, 0.29);
        else return 0.055;
    }
    else if (formula == "Yalin") {
        return 0.13 * pow(dimensionlessGrainDiameter, -0.392) * exp(-0.015*dimensionlessGrainDiameter*dimensionlessGrainDiameter)
               + 0.045 * (1 - exp(-0.068*dimensionlessGrainDiameter));
    }
    else {
        DUNE_THROW(Dune::InvalidStateException, "The parameter 'Sediment.CriticalShieldsStressCalculation' is set to '"<<formula
                   <<"', which is not a valid value. Valid values 'vanRijn' and 'Yalin'!");
    }
}
} // end namespace Dumux

#endif
