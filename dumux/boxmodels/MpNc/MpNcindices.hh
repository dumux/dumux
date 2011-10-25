/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_MPNC_INDICES_HH
#define DUMUX_MPNC_INDICES_HH

#include "MpNcproperties.hh"

#include "mass/MpNcindicesmass.hh"
#include "energy/MpNcindicesenergy.hh"

namespace Dumux
{
/*!
 * \brief The primary variable and equation indices for the MpNc
 *        model.
 */
template <class TypeTag, int BasePVOffset = 0>
struct MPNCIndices :
        public MPNCMassIndices<BasePVOffset,
                               TypeTag,
                               GET_PROP_VALUE(TypeTag, PTAG(EnableKinetic)) >,
        public MPNCEnergyIndices<BasePVOffset +
                                 MPNCMassIndices<0, TypeTag, GET_PROP_VALUE(TypeTag, PTAG(EnableKinetic)) >::NumPrimaryVars,
                                 GET_PROP_VALUE(TypeTag, PTAG(EnableEnergy)),
                                 GET_PROP_VALUE(TypeTag, PTAG(EnableKineticEnergy))>
{
private:
    enum { enableEnergy         = GET_PROP_VALUE(TypeTag, PTAG(EnableEnergy)) };
    enum { enableDiffusion      = GET_PROP_VALUE(TypeTag, PTAG(EnableDiffusion)) };
    enum { enableKinetic        = GET_PROP_VALUE(TypeTag, PTAG(EnableKinetic)) }; //mass transfer
    enum { enableKineticEnergy  = GET_PROP_VALUE(TypeTag, PTAG(EnableKineticEnergy)) }; // energy transfer

    typedef MPNCMassIndices<BasePVOffset, TypeTag, enableKinetic> MassIndices;
    typedef MPNCEnergyIndices<BasePVOffset + MassIndices::NumPrimaryVars, enableEnergy, enableKineticEnergy> EnergyIndices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    enum { numComponents = FluidSystem::numComponents };
    enum { numPhases = FluidSystem::numPhases };

public:
    /*!
     * \brief The number of primary variables / equations.
     */
    // temperature + Mass Balance  + constraints for switch stuff
    static const int NumPrimaryVars =
        MassIndices::NumPrimaryVars +
        EnergyIndices::NumPrimaryVars +
        numPhases;

    /*!
     * \brief The number of primary variables / equations of the energy module.
     */
    static const int NumPrimaryEnergyVars =
        EnergyIndices::NumPrimaryVars ;

    /*!
     * \brief Index of the saturation of the first phase in a vector
     *        of primary variables.
     *
     * The following (numPhases - 1) primary variables represent the
     * saturations for the phases [1, ..., numPhases - 1]
     */
    static const int S0Idx =
        MassIndices::NumPrimaryVars +
        EnergyIndices::NumPrimaryVars;

    /*!
     * \brief Index of the first phase' pressure in a vector of
     *        primary variables.
     */
    static const int p0Idx =
        MassIndices::NumPrimaryVars +
        EnergyIndices::NumPrimaryVars +
        numPhases - 1;

    /*!
     * \brief Index of the first phase NCP equation.
     *
     * The index for the remaining phases are consecutive.
     */
    static const int phase0NcpIdx =
        MassIndices::NumPrimaryVars +
        EnergyIndices::NumPrimaryVars;
};

}

#endif
