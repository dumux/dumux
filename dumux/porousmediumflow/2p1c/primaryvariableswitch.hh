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
 * \ingroup TwoPOneCModel
 * \copydoc Dumux::TwoPOneCPrimaryVariableSwitch
 */
#ifndef DUMUX_2P1C_PRIMARY_VARIABLE_SWITCH_HH
#define DUMUX_2P1C_PRIMARY_VARIABLE_SWITCH_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/compositional/primaryvariableswitch.hh>

namespace Dumux {

/*!
 * \ingroup TwoPOneCModel
 * \brief The primary variable switch for the two-phase one-component model
 */
template<class TypeTag>
class TwoPOneCPrimaryVariableSwitch
: public PrimaryVariableSwitch<TwoPOneCPrimaryVariableSwitch<TypeTag>>
{
    using ParentType = PrimaryVariableSwitch<TwoPOneCPrimaryVariableSwitch<TypeTag>>;
    friend ParentType;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    enum {
        switchIdx = Indices::switchIdx,

        pressureIdx = Indices::pressureIdx,

        wPhaseIdx = FluidSystem::wPhaseIdx,
        nPhaseIdx = FluidSystem::nPhaseIdx,

        twoPhases = Indices::twoPhases,
        wPhaseOnly  = Indices::wPhaseOnly,
        nPhaseOnly  = Indices::nPhaseOnly
    };

public:
    using ParentType::ParentType;

protected:

    /*!
     * \brief Perform variable switch at a degree of freedom location.
     *
     * \param priVars The primary variables at the given degree of freedom (dof) location.
     * \param volVars The volume variables.
     * \param dofIdxGlobal The respective dof index.
     * \param globalPos The global position of the dof.
     */
    bool update_(PrimaryVariables& priVars,
                 const VolumeVariables& volVars,
                 IndexType dofIdxGlobal,
                 const GlobalPosition& globalPos)
    {
        // evaluate primary variable switch
        bool wouldSwitch = false;
        int phasePresence =  priVars.state();
        int newPhasePresence = phasePresence;

        // check if a primary var switch is necessary
        if (phasePresence == twoPhases)
        {
            Scalar Smin = 0;
            if (this->wasSwitched_[dofIdxGlobal])
                Smin = -0.01;

            if (volVars.saturation(nPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // gas phase disappears
                std::cout << "Gas phase disappears at vertex " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", sg: "
                          << volVars.saturation(nPhaseIdx) << std::endl;
                newPhasePresence = wPhaseOnly;

                priVars[switchIdx] = volVars.fluidState().temperature();
            }
            else if (volVars.saturation(wPhaseIdx) <= Smin)
            {
                wouldSwitch = true;
                // water phase disappears
                std::cout << "Wetting phase disappears at vertex " << dofIdxGlobal
                          << ", coordinates: " << globalPos << ", sw: "
                          << volVars.saturation(wPhaseIdx) << std::endl;
                newPhasePresence = nPhaseOnly;

                priVars[switchIdx] = volVars.fluidState().temperature();
            }

        }
        else if (phasePresence == wPhaseOnly)
        {
            const Scalar temp = volVars.fluidState().temperature();
            const Scalar tempVap = volVars.vaporTemperature();

            // if the the temperature would be larger than
            // the vapor temperature at the given pressure, gas phase appears
            if (temp >= tempVap)
            {
                wouldSwitch = true;
                // gas phase appears
                std::cout << "gas phase appears at vertex " << dofIdxGlobal
                          << ", coordinates: " << globalPos  << std::endl;

               newPhasePresence = twoPhases;
               priVars[switchIdx] = 0.9999; //wetting phase saturation
            }
        }


        else if (phasePresence == nPhaseOnly)
        {

            const Scalar temp = volVars.fluidState().temperature();
            const Scalar tempVap = volVars.vaporTemperature();

            if (temp < tempVap)
            {
                wouldSwitch = true;
                // wetting phase appears
                std::cout << "wetting phase appears at vertex " << dofIdxGlobal
                          << ", coordinates: " << globalPos  << std::endl;


               newPhasePresence = twoPhases;
               priVars[switchIdx] = 0.0001; //arbitrary small value
            }
    }
        priVars.setState(newPhasePresence);
        this->wasSwitched_[dofIdxGlobal] = wouldSwitch;
        return phasePresence != newPhasePresence;
    }
};

} // end namespace Dumux

#endif
