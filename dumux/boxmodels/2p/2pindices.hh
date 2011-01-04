// $Id: 2pproperties.hh 3357 2010-03-25 13:02:05Z lauser $
/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \brief Defines the indices required for the two-phase box model.
 */
#ifndef DUMUX_BOX_2P_INDICES_HH
#define DUMUX_BOX_2P_INDICES_HH

namespace Dumux
{
/*!
 * \ingroup TwoPBoxModel
 */
// \{

/*!
 * \brief The common indices for the isothermal two-phase model.
 */
struct TwoPCommonIndices
{
    // Formulations
    static const int pwSn = 0; //!< Pw and Sn as primary variables
    static const int pnSw = 1; //!< Pn and Sw as primary variables

    // Phase indices
    static const int wPhaseIdx = 0; //!< Index of the wetting phase in a phase vector
    static const int nPhaseIdx = 1; //!< Index of the non-wetting phase in a phase vector
};

/*!
 * \brief The indices for the \f$p_w-S_n\f$ formulation of the
 *        isothermal two-phase model.
 *
 * \tparam formulation The formulation, either pwSn or pnSw
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int formulation = TwoPCommonIndices::pwSn, int PVOffset = 0>
struct TwoPIndices : public TwoPCommonIndices
{
    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int saturationIdx = PVOffset + 1; //!< Index of the saturation of the non-wetting/wetting phase

    // indices of the primary variables
    static const int pwIdx = PVOffset + 0; //!< Pressure index of the wetting phase
    static const int SnIdx = PVOffset + 1; //!< Saturation index of the wetting phase

    // indices of the equations
    static const int contiWEqIdx = PVOffset + 0; //!< Index of the continuity equation of the wetting phase
    static const int contiNEqIdx = PVOffset + 1; //!< Index of the continuity equation of the non-wetting phase
};

/*!
 * \brief The indices for the \f$p_w-S_n\f$ formulation of the
 *        isothermal two-phase model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int PVOffset>
struct TwoPIndices<TwoPCommonIndices::pnSw, PVOffset>
    : public TwoPCommonIndices
{
    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!< Index for wetting/non-wetting phase pressure (depending on formulation) in a solution vector
    static const int saturationIdx = PVOffset + 1; //!< Index of the saturation of the non-wetting/wetting phase

    // indices of the primary variables
    static const int pnIdx = PVOffset + 0; //!< Pressure index of the wetting phase
    static const int SwIdx = PVOffset + 1; //!< Saturation index of the wetting phase

    // indices of the equations
    static const int contiNEqIdx = PVOffset + 0; //!< Index of the continuity equation of the non-wetting phase
    static const int contiWEqIdx = PVOffset + 1; //!< Index of the continuity equation of the wetting phase
};

// \}
} // namespace Dumux


#endif
