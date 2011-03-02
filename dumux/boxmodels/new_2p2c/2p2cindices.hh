// $Id: 2p2cproperties.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
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

/*!
 * \file
 *
 * \brief Defines the indices required for the 2p2c BOX model.
 */
#ifndef DUMUX_2P2C_INDICES_HH
#define DUMUX_2P2C_INDICES_HH

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCModel
 */
// \{

/*!
 * \brief Enumerates the formulations which the 2p2c model accepts.
 */
struct TwoPTwoCFormulation
{
    enum {
        plSg,
        pgSl
    };
};

/*!
 * \brief The indices for the isothermal TwoPTwoC model.
 *
 * \tparam formulation The formulation, either pwSn or pnSw.
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag,
          int formulation = TwoPTwoCFormulation::plSg,
          int PVOffset = 0>
class TwoPTwoCIndices
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    // Phase indices
    static const int lPhaseIdx = FluidSystem::lPhaseIdx; //!< Index of the liquid phase
    static const int gPhaseIdx = FluidSystem::gPhaseIdx; //!< Index of the gas phase

    // Component indices
    static const int lCompIdx = 0; //!< Index of the liquid's primary component
    static const int gCompIdx = 1; //!< Index of the gas' primary component

    // present phases (-> 'pseudo' primary variable)
    static const int lPhaseOnly = 1; //!< Only the non-wetting phase is present
    static const int gPhaseOnly = 0; //!< Only the wetting phase is present
    static const int bothPhases = 2; //!< Both phases are present

    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int switchIdx = PVOffset + 1; //!< Index of the either the saturation or the mass fraction of the non-wetting/wetting phase

    static const int plIdx = pressureIdx; //!< Index for liquid phase pressure in a solution vector
    static const int SgOrXIdx = switchIdx; //!< Index of the either the saturation of the gas phase or the mass fraction secondary component in the only phase

    // equation indices
    static const int contiLEqIdx = PVOffset + lCompIdx; //!< Index of the mass conservation equation for the liquid's primary component
    static const int contiGEqIdx = PVOffset + gCompIdx; //!< Index of the mass conservation equation for the gas' primary component
};

/*!
 * \brief The indices for the isothermal TwoPTwoC model in the pg-Sl
 *        formulation.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <class TypeTag, int PVOffset>
class TwoPTwoCIndices<TypeTag, TwoPTwoCFormulation::pgSl, PVOffset>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    // Phase indices
    static const int lPhaseIdx = FluidSystem::lPhaseIdx; //!< Index of the liquid phase
    static const int gPhaseIdx = FluidSystem::gPhaseIdx; //!< Index of the gas phase

    // Component indices
    static const int lCompIdx = 0; //!< Index of the liquid's primary component
    static const int gCompIdx = 1; //!< Index of the gas' primary component

    // present phases (-> 'pseudo' primary variable)
    static const int lPhaseOnly = 1; //!< Only the non-wetting phase is present
    static const int gPhaseOnly = 2; //!< Only the wetting phase is present
    static const int bothPhases = 3; //!< Both phases are present

    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int switchIdx = PVOffset + 1; //!< Index of the either the saturation or the mass fraction of the non-wetting/wetting phase

    static const int pgIdx = pressureIdx; //!< Index for gas phase pressure in a solution vector
    static const int SlOrXIdx = switchIdx; //!< Index of the either the saturation of the liquid phase or the mass fraction secondary component in the only phase

    // Equation indices
    static const int contiLEqIdx = PVOffset + lCompIdx; //!< Index of the mass conservation equation for the liquid's primary component
    static const int contiGEqIdx = PVOffset + gCompIdx; //!< Index of the mass conservation equation for the gas' primary component
};

// \}

}

#endif
