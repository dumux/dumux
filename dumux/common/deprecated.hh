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
 * \ingroup Common
 * \brief Helpers for deprecation
 */

#ifndef DUMUX_COMMON_DEPRECATED_HH
#define DUMUX_COMMON_DEPRECATED_HH

#include <utility>

#include <dune/common/ftraits.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/std/type_traits.hh>
#include <dumux/common/numeqvector.hh>

namespace Dumux {

#ifndef DOXYGEN // hide from doxygen
// Helper classes/functions for deprecation
// Each implementation has to state after which release
// it will be removed. Implementations in the Deprecated
// namespace will be removed without
// deprecation after their usage in the code expired,
// so most likely you don't want to use this in your code
namespace Deprecated {

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif // __clang__

template<class G>
using DetectFVGeometryHasSCVGeometry = decltype(std::declval<G>().geometry(std::declval<typename G::SubControlVolume>()));

template<class G>
using DetectFVGeometryHasSCVFGeometry = decltype(std::declval<G>().geometry(std::declval<typename G::SubControlVolumeFace>()));

template<class G>
constexpr inline bool hasSCVGeometryInterface()
{ return Dune::Std::is_detected<DetectFVGeometryHasSCVGeometry, G>::value; }

template<class G>
constexpr inline bool hasSCVFGeometryInterface()
{ return Dune::Std::is_detected<DetectFVGeometryHasSCVFGeometry, G>::value; }

#ifdef __clang__
#pragma clang diagnostic pop
#endif  // __clang__

template<bool enableWaterDiffusionInAir>
struct ExtendedRichardsHelper
{
    template<bool b>
    [[deprecated("Enabling the extended Richards model through a template parameter/properties is deprecated and will be removed after release (3.6). Use the new model ExtendedRichards instead.")]]
    static constexpr void extendedRichardsWithTemplateParameter() {}

    static constexpr bool isExtendedRichards()
    {
        if constexpr(enableWaterDiffusionInAir)
            extendedRichardsWithTemplateParameter<enableWaterDiffusionInAir>();
        return enableWaterDiffusionInAir;
    }
};

/*!
 * \ingroup RichardsModel
 * \brief A primary variable vector with a state to allow variable switches.
 */
template<class PrimaryVariables, class StateType>
class RichardsSwitchablePrimaryVariables : public PrimaryVariables
{
    using ParentType = PrimaryVariables;
public:
    //! Inherit the constructors
    using ParentType::ParentType;
    //! Use the assignment operators from the field vector
    using ParentType::operator=;

    [[deprecated("Will be removed after release (3.6). Normal Richards does not require setting a state. Use ExtendedRichards for extended model. "
                 "After removing all .setState()/.state() calls in the user files the remaining deprecation warnings can be removed by changing the PrimaryVariables property to Dune::FieldVector. "
                 "template<class TypeTag> "
                 "struct PrimaryVariables<TypeTag, TTag::Richards> "
                 "{ using type = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; }")]]
    StateType state() const
    {
        return state_;
    }

    //! Set the state of this primary variable object, e.g. the phase presence.
    [[deprecated("Will be removed after release (3.6). Normal Richards does not require setting a state. Use ExtendedRichards for extended model."
                 "After removing all .setState()/.state() calls in the user files the remaining deprecation warnings can be removed by changing the PrimaryVariables property to Dune::FieldVector. "
                 "template<class TypeTag> "
                 "struct PrimaryVariables<TypeTag, TTag::Richards> "
                 "{ using type = Dune::FieldVector<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::ModelTraits>::numEq()>; }")]]
    void setState(StateType state)
    {
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


} // end namespace Deprecated
#endif

/*!
 * \ingroup PorousmediumflowModels
 * \brief The corresponding NumEqVectorTraits for the primary variables with switchable state
 */
template<class PrimaryVariables, class StateType>
struct NumEqVectorTraits<Deprecated::RichardsSwitchablePrimaryVariables<PrimaryVariables, StateType>>
{
    static constexpr std::size_t numEq = PrimaryVariables::size();
    using type = PrimaryVariables;
};

} // end namespace Dumux

// specialize field traits for this type
namespace Dune {

template <class PrimaryVariables, class StateType>
struct FieldTraits<Dumux::Deprecated::RichardsSwitchablePrimaryVariables<PrimaryVariables, StateType>>
: public FieldTraits<PrimaryVariables>
{};

} // end namespace Dune
#endif
