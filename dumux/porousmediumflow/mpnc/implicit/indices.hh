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
 * \brief The primary variable and equation indices for the MpNc model.
 */
#ifndef DUMUX_MPNC_INDICES_HH
#define DUMUX_MPNC_INDICES_HH

#include "properties.hh"

#include "mass/indices.hh"
#include "energy/indices.hh"

namespace Dumux
{

/*!
 * \ingroup MPNCModel
 * \ingroup ImplicitIndices
 * \brief Enumerates the formulations which the MpNc model accepts.
 */
struct MpNcPressureFormulation
{
    enum {
        mostWettingFirst,
        leastWettingFirst
    };
};

/*!
 * \ingroup MPNCModel
 * \ingroup ImplicitIndices
 * \brief The primary variable and equation indices for the MpNc model.
 */
template <class TypeTag, int BasePVOffset = 0>
struct MPNCIndices :
        public MPNCMassIndices<BasePVOffset,
                               TypeTag,
                               GET_PROP_VALUE(TypeTag, EnableKinetic) >,
        public MPNCEnergyIndices<BasePVOffset +
                                 MPNCMassIndices<0, TypeTag, GET_PROP_VALUE(TypeTag, EnableKinetic) >::numPrimaryVars,
                                 GET_PROP_VALUE(TypeTag, EnableEnergy),
                                 GET_PROP_VALUE(TypeTag, NumEnergyEquations)>
{
private:
            using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
            enum { enableEnergy         = GET_PROP_VALUE(TypeTag, EnableEnergy) };
            enum { enableKinetic        = GET_PROP_VALUE(TypeTag, EnableKinetic) }; //mass transfer
            enum { numEnergyEquations  = GET_PROP_VALUE(TypeTag, NumEnergyEquations) }; // energy transfer
            enum { numPhases = FluidSystem::numPhases };

            using MassIndices = MPNCMassIndices<BasePVOffset, TypeTag, enableKinetic>;
            using EnergyIndices = MPNCEnergyIndices<BasePVOffset + MassIndices::numPrimaryVars, enableEnergy, numEnergyEquations>;

public:
    /*!
     * \brief The number of primary variables / equations.
     */
    // temperature + Mass Balance  + constraints for switch stuff
    static const unsigned int numPrimaryVars =
        MassIndices::numPrimaryVars +
        EnergyIndices::numPrimaryVars +
        numPhases;

    /*!
     * \brief The number of primary variables / equations of the energy module.
     */
    static const unsigned int numPrimaryEnergyVars =
        EnergyIndices::numPrimaryVars ;

    /*!
     * \brief Index of the saturation of the first phase in a vector
     *        of primary variables.
     *
     * The following (numPhases - 1) primary variables represent the
     * saturations for the phases [1, ..., numPhases - 1]
     */
    static const unsigned int s0Idx =
        MassIndices::numPrimaryVars +
        EnergyIndices::numPrimaryVars;

    /*!
     * \brief Index of the first phase' pressure in a vector of
     *        primary variables.
     */
    static const unsigned int p0Idx =
        MassIndices::numPrimaryVars +
        EnergyIndices::numPrimaryVars +
        numPhases - 1;

    /*!
     * \brief Index of the first phase NCP equation.
     *
     * The index for the remaining phases are consecutive.
     */
    static const unsigned int phase0NcpIdx =
        MassIndices::numPrimaryVars +
        EnergyIndices::numPrimaryVars;
};

}

#endif
