// $Id: 1p2cfluidstate.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2010 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Calcultes the phase state from the primary variables in the
 *        1p2c model.
 */
#ifndef DUMUX_1P2C_PHASE_STATE_HH
#define DUMUX_1P2C_PHASE_STATE_HH

#include "1p2cproperties.hh"

#include <dumux/material/fluidstate.hh>

namespace Dumux
{
/*!
 * \brief Calcultes the phase state from the primary variables in the
 *        1p2c model.
 */
template <class TypeTag>
class OnePTwoCFluidState : public FluidState<typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)),
                                             OnePTwoCFluidState<TypeTag> >
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVarVector)) PrimaryVarVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;

    enum {
        konti = Indices::konti,
        transport = Indices::transport
    };

public:
    enum { numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)) };

    /*!
     * \brief Update the phase state from the primary variables.
     */
    void update(const PrimaryVarVector &primaryVars,
                Scalar temperature)
    {
        Valgrind::CheckDefined(primaryVars);

        temperature_ = temperature;

        phasePressure_[0] = primaryVars[konti];

        moleFraction_[0][konti] = 1.0 - primaryVars[transport];
        moleFraction_[0][transport] = primaryVars[transport];
    }

public:
    /*!
     * \brief Returns the molar fraction of a component in a fluid phase.
     */
    Scalar moleFrac(int phaseIdx, int compIdx) const
    {
        return moleFraction_[phaseIdx][compIdx];
    }

    /*!
     * \brief Returns the pressure of a fluid phase [Pa].
     */
    Scalar phasePressure(int phaseIdx) const
    { return phasePressure_[phaseIdx]; }

    /*!
     * \brief Returns the temperature of the fluids [K].
     *
     * Note that we assume thermodynamic equilibrium, so all fluids
     * and the rock matrix exhibit the same temperature.
     */
    Scalar temperature() const
    { return temperature_; };

public:
    Scalar moleFraction_[numPhases][numComponents];
    Scalar phasePressure_[numPhases];
    Scalar temperature_;
};

} // end namepace

#endif
