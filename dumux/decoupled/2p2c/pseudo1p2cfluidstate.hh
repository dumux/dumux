/*****************************************************************************
 *   Copyright (C) 2009 by Benjamin Faigle                                   *
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
 * \brief Calculates phase state for a single phase but two-component state.
 */
#ifndef DUMUX_PSEUDO1P2C_FLUID_STATE_HH
#define DUMUX_PSEUDO1P2C_FLUID_STATE_HH

#include <dumux/material/fluidstate.hh>
#include <dumux/decoupled/2p2c/2p2cproperties.hh>

namespace Dumux
{
/*!
 * \ingroup multiphysics
 * \brief Container for compositional variables in a 1p2c situation
 *
 *  This class represents a pseudo flash calculation when only single-phase situations
 *  occur. This is used in case of a multiphysics approach.
 *  \tparam TypeTag The property Type Tag
 */
template <class TypeTag>
class PseudoOnePTwoCFluidState : public FluidState<typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)),
                                           PseudoOnePTwoCFluidState<TypeTag> >
{
    typedef DecTwoPTwoCFluidState<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))      Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;


public:
    enum {     numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
            numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),
            wPhaseIdx = 0,
            nPhaseIdx = 1,};

public:
    /*!
     * \name flash calculation routines
     * Routine to determine the phase composition after the transport step.
     */
    //@{
    //! The simplest possible update routine for 1p2c "flash" calculations
    /*!
     * Routine goes as follows:
     * - Check if we are in single phase condition
     * - Assign total concentration to the present phase
     *
     * \param Z1 Feed mass fraction \f$\mathrm{[-]}\f$
     * \param press1p Pressure value for present phase \f$\mathrm{[Pa]}\f$
     * \param satW Saturation of the wetting phase \f$\mathrm{[-]}\f$
     */
    void update(Scalar Z1,Scalar press1p, Scalar satW)
    {
        Sw_ = satW;
        pressure1p_=press1p;

        if (satW == 1.)
        {
            massfracX1_[wPhaseIdx] = Z1;
            massfracX1_[nPhaseIdx] = 0.;
        }
        else if (satW == 0.)
        {
            massfracX1_[wPhaseIdx] = 0.;
            massfracX1_[nPhaseIdx] = Z1;
        }
        else
            Dune::dgrave << "Twophase conditions in simple-phase flash! Saturation is " << satW << std::endl;

        return;
    }
    //@}
    /*! \name Acess functions */
    //@{
    /*! \brief Returns the saturation of a phase.
     *  \param phaseIdx Index of the phase
     */
    Scalar saturation(int phaseIdx) const
    {
        if (phaseIdx == wPhaseIdx)
            return Sw_;
        else
            return Scalar(1.0) - Sw_;
    };

    /*! \brief Returns the pressure of a fluid phase \f$\mathrm{[Pa]}\f$.
     *  \param phaseIdx Index of the phase
     */
    Scalar phasePressure(int phaseIdx) const
    { return pressure1p_; }


    /*!
     * \brief Returns the mass fraction of a component in a phase.
     *  \param phaseIdx Index of the phase
     *  \param compIdx Index of the component
     */
    Scalar massFrac(int phaseIdx, int compIdx) const
    {
        if (compIdx == wPhaseIdx)
            return massfracX1_[phaseIdx];
        else
            return 1.-massfracX1_[phaseIdx];

    }

    /*!
     * \brief Returns the molar fraction of a component in a fluid phase.
     *  \param phaseIdx Index of the phase
     *  \param compIdx Index of the component
     */
    Scalar moleFrac(int phaseIdx, int compIdx) const
    {
        double moleFrac_(0.);
        // as the mass fractions are calculated, these are used to determine the mole fractions
        if (compIdx == wPhaseIdx)   //only interested in wComp => X1
        {
            moleFrac_ = ( massfracX1_[phaseIdx] / FluidSystem::molarMass(compIdx) );    // = moles of compIdx
            moleFrac_ /= ( massfracX1_[phaseIdx] / FluidSystem::molarMass(0)
                           + (1.-massfracX1_[phaseIdx]) / FluidSystem::molarMass(1) );    // /= total moles in phase
        }
        else    // interested in nComp => 1-X1
        {
            moleFrac_ = ( (1.-massfracX1_[phaseIdx]) / FluidSystem::molarMass(compIdx) );   // = moles of compIdx
            moleFrac_ /= ((1.- massfracX1_[phaseIdx] )/ FluidSystem::molarMass(0)
                           + massfracX1_[phaseIdx] / FluidSystem::molarMass(1) );    // /= total moles in phase
        }
        return moleFrac_;
    }
    //@}

public:
    Scalar massfracX1_[numPhases];
    Scalar Sw_;
    Scalar pressure1p_;
};

} // end namepace

#endif
