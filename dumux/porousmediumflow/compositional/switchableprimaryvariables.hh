// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup PorousmediumflowModels
 * \brief A primary variable vector with a state to allow variable switches.
 */

#ifndef DUMUX_SWITCHABLE_PRIMARY_VARIABLES_HH
#define DUMUX_SWITCHABLE_PRIMARY_VARIABLES_HH

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumflowModels
 * \brief A primary variable vector with a state to allow variable switches.
 */
template<class PrimaryVariables, class StateType>
class SwitchablePrimaryVariables : public PrimaryVariables
{
    using ParentType = PrimaryVariables;
public:
    //! Inherit the constructors
    using ParentType::ParentType;
    //! Use the assignment operators from the field vector
    using ParentType::operator=;

    //! Ask for the state of this primary variable object, e.g. the phase presence
    StateType state() const
    {
        if (!stateIsSet_)
            DUNE_THROW(Dune::InvalidStateException, "Model demands setting a primary variable state (like a phase presence)"
                                                 << " but none was set! Use its setState method to set the state.");

        return state_;
    }

    //! Set the state of this primary variable object, e.g. the phase presence.
    void setState(StateType state)
    {
        // NOTE: we use a copy instead of a reference in the signature to
        // avoid linker errors related to passing a static variable to this function
        state_ = std::move(state);
        stateIsSet_ = true;
    }

    //! Invalidate the state
    void invalidateState()
    {
        stateIsSet_ = false;
    }

private:
    StateType state_;
    bool stateIsSet_{false};
};

} // end namespace Dumux

#endif
