// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
#ifndef DUMUX_DECOUPLED_2P_INDICES_HH
#define DUMUX_DECOUPLED_2P_INDICES_HH

namespace Dumux
{
/*!
 * \ingroup IMPES
 */
// \{

/*!
 * \brief The common indices for the isothermal two-phase model.
 */
struct DecoupledTwoPCommonIndices
{
    // Formulations
    static const int pwSn = 0; //!< Pw and Sn as primary variables
    static const int pnSw = 1; //!< Pn and Sw as primary variables
    static const int pwSw = 2; //!< Pw and Sw as primary variables
    static const int pnSn = 3; //!< Pn and Sn as primary variables
    static const int pGlobalSw = 4; //!< PGlobal and Sw as primary variables
    static const int pGlobalSn = 5; //!< PGlobal and Sn as primary variables

    // Phase indices
    static const int wPhaseIdx = 0; //!< Index of the wetting phase in a phase vector
    static const int nPhaseIdx = 1; //!< Index of the non-wetting phase in a phase vector
    static const int totalPhaseIdx = 2; //!< Index of the total phase (wetting + nonwetting)

    //saturation flags
    static const int saturationW = 0; //!< Indicates wetting phase saturation
    static const int saturationNW = 1; //!<  Indicates non-wetting phase saturation
    //pressure flags
    static const int pressureW = 0; //!< Indicates wetting phase pressure
    static const int pressureNW = 1; //!<  Indicates non-wetting phase pressure
    static const int pressureGlobal = 2; //!<  Indicates global-pressure

    //velocity flags
    static const int velocityW = 0; //!< Indicates wetting phase velocity
    static const int velocityNW = 1; //!<  Indicates non-wetting phase velocity
    static const int velocityTotal = 2; //!<  Indicates total velocity
};

/*!
 * \brief The indices for the \f$p_w-S_n\f$ formulation of the
 *        isothermal two-phase model.
 *
 * \tparam formulation Index of the formulation
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int formulation = DecoupledTwoPCommonIndices::pwSn, int PVOffset = 0>
struct DecoupledTwoPIndices : public DecoupledTwoPCommonIndices
{
    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!<  Index for the primary pressure vairable in a solution vector
    static const int saturationIdx = PVOffset + 1; //!< Index for the primary saturation variable in a solution vector

    // indices of the primary variables
    static const int pwIdx = PVOffset + 0; //!< Pressure index of the wetting phase
    static const int SnIdx = PVOffset + 1; //!< Saturation index of the non-wetting phase

    //! \cond \private
    //Set the types of the single models depending on the formulation
    static const int pressureType = pressureW;
    static const int saturationType = saturationNW;

    static const int velocityDefault = velocityNW;
    //! \endcond

    // indices of the equations
    static const int contiWEqIdx = PVOffset + 0; //!< Index of the continuity equation of the wetting phase
    static const int pressEqIdx = contiWEqIdx; //!< Index of the pressure equation (total mass balance)
    static const int contiNEqIdx = PVOffset + 1; //!< Index of the continuity equation of the non-wetting phase
    static const int satEqIdx = contiNEqIdx; //!< Index of the continuity equation of the non-wetting phase (saturation equation)
    static const int transportEqIdx = satEqIdx; //!< Index of the saturation transport equation
};

/*!
 * \brief The indices for the \f$p_n-S_w\f$ formulation of the
 *        isothermal two-phase model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int PVOffset>
struct DecoupledTwoPIndices<DecoupledTwoPCommonIndices::pnSw, PVOffset>
    : public DecoupledTwoPCommonIndices
{
    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!<  Index for the primary pressure vairable in a solution vector
    static const int saturationIdx = PVOffset + 1; //!< Index for the primary saturation variable in a solution vector

    // indices of the primary variables
    static const int pnIdx = PVOffset + 0; //!< Pressure index of the non-wetting phase
    static const int SwIdx = PVOffset + 1; //!< Saturation index of the wetting phase

    //! \cond \private
    //Set the types of the single models depending on the formulation
    static const int pressureType = pressureNW;
    static const int saturationType = saturationW;

    static const int velocityDefault = velocityW;
    //! \endcond

    // indices of the equations
    static const int contiNEqIdx = PVOffset + 0; //!< Index of the continuity equation of the non-wetting phase
    static const int pressEqIdx = contiNEqIdx; //!< Index of the pressure equation (total mass balance)
    static const int contiWEqIdx = PVOffset + 1; //!< Index of the continuity equation of the wetting phase
    static const int satEqIdx = contiWEqIdx; //!< Index of the continuity equation of the non-wetting phase (saturation equation)
    static const int transportEqIdx = satEqIdx; //!< Index of the saturation transport equation
};


/*!
 * \brief The indices for the \f$p_w-S_w\f$ formulation of the
 *        isothermal two-phase model.
 *
 * \tparam formulation Index of the formulation
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int PVOffset>
struct DecoupledTwoPIndices<DecoupledTwoPCommonIndices::pwSw, PVOffset>
    : public DecoupledTwoPCommonIndices
{
    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!<  Index for the primary pressure vairable in a solution vector
    static const int saturationIdx = PVOffset + 1; //!< Index for the primary saturation variable in a solution vector

    // indices of the primary variables
    static const int pwIdx = PVOffset + 0; //!< Pressure index of the wetting phase
    static const int SwIdx = PVOffset + 1; //!< Saturation index of the wetting phase

    //! \cond \private
    //Set the types of the single models depending on the formulation
    static const int pressureType = pressureW;
    static const int saturationType = saturationW;

    static const int velocityDefault = velocityW;
    //! \endcond

    // indices of the equations
    static const int contiWEqIdx = PVOffset + 0; //!< Index of the continuity equation of the wetting phase
    static const int pressEqIdx = contiWEqIdx; //!< Index of the pressure equation (total mass balance)
    static const int contiNEqIdx = PVOffset + 1; //!< Index of the continuity equation of the non-wetting phase
    static const int satEqIdx = contiNEqIdx; //!< Index of the continuity equation of the non-wetting phase (saturation equation)
    static const int transportEqIdx = satEqIdx; //!< Index of the saturation transport equation
};

/*!
 * \brief The indices for the \f$p_n-S_n\f$ formulation of the
 *        isothermal two-phase model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int PVOffset>
struct DecoupledTwoPIndices<DecoupledTwoPCommonIndices::pnSn, PVOffset>
    : public DecoupledTwoPCommonIndices
{
    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!<  Index for the primary pressure vairable in a solution vector
    static const int saturationIdx = PVOffset + 1; //!< Index for the primary saturation variable in a solution vector

    // indices of the primary variables
    static const int pnIdx = PVOffset + 0; //!< Pressure index of the non-wetting phase
    static const int SnIdx = PVOffset + 1; //!< Saturation index of the non-wetting phase

    //! \cond \private
    //Set the types of the single models depending on the formulation
    static const int pressureType = pressureNW;
    static const int saturationType = saturationNW;

    static const int velocityDefault = velocityNW;
    //! \endcond

    // indices of the equations
    static const int contiNEqIdx = PVOffset + 0; //!< Index of the continuity equation of the non-wetting phase
    static const int pressEqIdx = contiNEqIdx; //!< Index of the pressure equation (total mass balance)
    static const int contiWEqIdx = PVOffset + 1; //!< Index of the continuity equation of the wetting phase
    static const int satEqIdx = contiWEqIdx; //!< Index of the continuity equation of the non-wetting phase (saturation equation)
    static const int transportEqIdx = satEqIdx; //!< Index of the saturation transport equation
};


/*!
 * \brief The indices for the \f$p_{global}-S_w\f$ formulation of the
 *        isothermal two-phase model.
 *
 * \tparam formulation Index of the formulation
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int PVOffset>
struct DecoupledTwoPIndices<DecoupledTwoPCommonIndices::pGlobalSw, PVOffset> : public DecoupledTwoPCommonIndices
{
    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!<  Index for the primary pressure vairable in a solution vector
    static const int saturationIdx = PVOffset + 1; //!< Index for the primary saturation variable in a solution vector

    // indices of the primary variables
    static const int pGlobalIdx = PVOffset + 0; //!< Global pressure index
    static const int SwIdx = PVOffset + 1; //!< Saturation index of the wetting phase

    //! \cond \private
    //Set the types of the single models depending on the formulation
    static const int pressureType = pressureGlobal;
    static const int saturationType = saturationW;

    static const int velocityDefault = velocityTotal;
    //! \endcond

    // indices of the equations
    static const int pressEqIdx = PVOffset + 0; //!< Index of the pressure equation (total mass balance)
    static const int satEqIdx = PVOffset + 1; //!< Index of the continuity equation of the non-wetting phase (saturation equation)
    static const int transportEqIdx = satEqIdx; //!< Index of the saturation transport equation
};

/*!
 * \brief The indices for the \f$p_{global}-S_n\f$ formulation of the
 *        isothermal two-phase model.
 *
 * \tparam PVOffset The first index in a primary variable vector.
 */
template <int PVOffset>
struct DecoupledTwoPIndices<DecoupledTwoPCommonIndices::pGlobalSn, PVOffset>
    : public DecoupledTwoPCommonIndices
{
    // Primary variable indices
    static const int pressureIdx = PVOffset + 0; //!<  Index for the primary pressure vairable in a solution vector
    static const int saturationIdx = PVOffset + 1; //!< Index for the primary saturation variable in a solution vector

    // indices of the primary variables
    static const int pGlobalIdx = PVOffset + 0; //!< Global pressure index
    static const int SnIdx = PVOffset + 1; //!< Saturation index of the non-wetting phase

    //! \cond \private
    //Set the types of the single models depending on the formulation
    static const int pressureType = pressureGlobal;
    static const int saturationType = saturationNW;

    static const int velocityDefault = velocityTotal;
    //! \endcond

    // indices of the equations
    static const int pressEqIdx = PVOffset + 0; //!< Index of the pressure equation (total mass balance)
    static const int satEqIdx = PVOffset + 1; //!< Index of the continuity equation of the non-wetting phase (saturation equation)
    static const int transportEqIdx = satEqIdx; //!< Index of the saturation transport equation
};

// \}
} // namespace Dumux


#endif
