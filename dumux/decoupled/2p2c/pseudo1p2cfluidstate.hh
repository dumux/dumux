// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
class PseudoOnePTwoCFluidState
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)      Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

public:
    enum {     numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
            numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        wCompIdx = Indices::wPhaseIdx,
        nCompIdx = Indices::nPhaseIdx
    };

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
     * \param presentPhaseIdx Subdomain Index = Indication which phase is present
     * \param temperature Temperature \f$\mathrm{[K]}\f$
     */
    void update(const Scalar& Z1,const Scalar& press1p,const int presentPhaseIdx, const Scalar& temperature)
    {

        pressure1p_=press1p;
        temperature_ = temperature;

        if (presentPhaseIdx == wPhaseIdx)
        {
            Sw_ = 1.;
            massFractionWater_[wPhaseIdx] = Z1;
            massFractionWater_[nPhaseIdx] = 0.;

            moleFractionWater_[wPhaseIdx] = Z1 / FluidSystem::molarMass(0);
            moleFractionWater_[wPhaseIdx] /= ( Z1 / FluidSystem::molarMass(0)
                           + (1.-Z1) / FluidSystem::molarMass(1));
            moleFractionWater_[nPhaseIdx] = 0.;

            presentPhaseIdx_ = wPhaseIdx;
        }
        else if (presentPhaseIdx == nPhaseIdx)
        {
            Sw_ = 0.;
            massFractionWater_[wPhaseIdx] = 0.;
            massFractionWater_[nPhaseIdx] = Z1;

            // interested in nComp => 1-X1
            moleFractionWater_[nPhaseIdx] = ( Z1 / FluidSystem::molarMass(0) );   // = moles of compIdx
            moleFractionWater_[nPhaseIdx] /= (Z1/ FluidSystem::molarMass(0)
                           + (1.-Z1) / FluidSystem::molarMass(1) );    // /= total moles in phase
            moleFractionWater_[nPhaseIdx] = 0.;

            presentPhaseIdx_ = nPhaseIdx;
        }
        else
            Dune::dgrave << __FILE__ <<": Twophase conditions in single-phase flash!"
                << " Z1 is "<<Z1<< std::endl;

        density_ = FluidSystem::density(*this, presentPhaseIdx);

        aveMoMass_ =  moleFractionWater_[presentPhaseIdx]*FluidSystem::molarMass(wCompIdx)
        +   (1.-moleFractionWater_[presentPhaseIdx])*FluidSystem::molarMass(nCompIdx);

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

    int presentPhaseIdx() const
    {
        return presentPhaseIdx_;
    }

    /*! \brief Returns the pressure of a fluid phase \f$\mathrm{[Pa]}\f$.
     *  \param phaseIdx Index of the phase
     */
    Scalar pressure(int phaseIdx) const
    { return pressure1p_; }

    Scalar density(int phaseIdx) const
    {
        if(phaseIdx == presentPhaseIdx_)
            return density_;
        else
            return FluidSystem::density(*this, phaseIdx);
    }

    /*!
     * \brief Returns the mass fraction of a component in a phase.
     *  \param phaseIdx Index of the phase
     *  \param compIdx Index of the component
     */
    Scalar massFraction(int phaseIdx, int compIdx) const
    {
        if (compIdx == wPhaseIdx)
            return massFractionWater_[phaseIdx];
        else
            return 1.-massFractionWater_[phaseIdx];

    }

    /*!
     * \brief Returns the molar fraction of a component in a fluid phase.
     *  \param phaseIdx Index of the phase
     *  \param compIdx Index of the component
     */
    Scalar moleFraction(int phaseIdx, int compIdx) const
    {
        if (compIdx == wPhaseIdx)
            return moleFractionWater_[phaseIdx];
        else
            return 1.-moleFractionWater_[phaseIdx];
    }

    /*!
     * \brief Returns the viscosity of a phase TODO: \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param phaseIdx the index of the phase
     */
    Scalar viscosity(int phaseIdx) const
    {
        assert(phaseIdx == presentPhaseIdx_);
        return viscosity_;
    }
    void setViscosity(int phaseIdx, Scalar value)
    {
        assert(phaseIdx == presentPhaseIdx_);
        viscosity_ = value;
    }

    Scalar averageMolarMass(int phaseIdx) const
    {
        return aveMoMass_;
//        return moleFractionWater_[phaseIdx]*FluidSystem::molarMass(wCompIdx)
//                +   (1.-moleFractionWater_[phaseIdx])*FluidSystem::molarMass(nCompIdx);
    }

    /*!
     * \brief Returns the temperature of the fluids \f$\mathrm{[K]}\f$.
     *
     * Note that we assume thermodynamic equilibrium, so all fluids
     * and the rock matrix exhibit the same temperature.
     */
    Scalar temperature(int phaseIdx) const
    { return temperature_; };

    /*!
     * \brief Returns the total mass concentration of a component \f$\mathrm{[kg/m^3]}\f$.
     *
     * This is equivalent to the sum of the component concentrations for all
     * phases multiplied with the phase density.
     *
     * \param compIdx the index of the component
     */
    Scalar massConcentration(int compIdx) const
    {
        return massConcentration_[compIdx];
    };
    /*!
     * \brief Sets the total mass concentration of a component \f$\mathrm{[kg/m^3]}\f$.
     */
    void setMassConcentration(int compIdx, Scalar value)
    {
        massConcentration_[compIdx] = value;
    };
    //@}

public:
    Scalar aveMoMass_;
    Scalar massConcentration_[numComponents];
    Scalar massFractionWater_[numPhases];
    Scalar moleFractionWater_[numPhases];
    Scalar Sw_;
    Scalar pressure1p_;
    Scalar density_;
    Scalar viscosity_;
    Scalar temperature_;
    int presentPhaseIdx_;
};

} // end namepace

#endif
