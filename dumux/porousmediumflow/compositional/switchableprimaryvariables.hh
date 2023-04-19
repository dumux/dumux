// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PorousmediumflowModels
 * \brief A primary variable vector with a state to allow variable switches.
 */

#ifndef DUMUX_SWITCHABLE_PRIMARY_VARIABLES_HH
#define DUMUX_SWITCHABLE_PRIMARY_VARIABLES_HH

#include <dune/common/ftraits.hh>
#include <dune/common/exceptions.hh>
#include <dumux/common/numeqvector.hh>

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

/*!
 * \ingroup PorousmediumflowModels
 * \brief The corresponding NumEqVectorTraits for the primary variables with switchable state
 */
template<class PrimaryVariables, class StateType>
struct NumEqVectorTraits<SwitchablePrimaryVariables<PrimaryVariables, StateType>>
{
    static constexpr std::size_t numEq = PrimaryVariables::size();
    using type = PrimaryVariables;
};

} // end namespace Dumux

// specialize field traits for this type
namespace Dune {

template <class PrimaryVariables, class StateType>
struct FieldTraits<Dumux::SwitchablePrimaryVariables<PrimaryVariables, StateType>>
: public FieldTraits<PrimaryVariables>
{};

} // end namespace Dune

#endif
