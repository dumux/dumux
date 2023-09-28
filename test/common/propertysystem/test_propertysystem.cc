// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Common
 * \ingroup Tests
 * \brief Testing the Dumux property system
 */

#include <iostream>
#include <type_traits>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/properties/propertysystem.hh>

namespace Dumux::Properties {

// create some properties:
// the first type tag is the actual TypeTag for which the property will be obtained
// (can be used to make properties depend on other properties),
// the second type tag is for partial specialization
// the default property should be always undefined to produce a good error message
// if the user attempt to get an unset property
template<class TypeTag, class MyTypeTag>
struct Scalar { using type = UndefinedProperty;  };

template<class TypeTag, class MyTypeTag>
struct CoordinateType { using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct UseTpfaFlux { using type = UndefinedProperty; };

namespace TTag {
// create some type tags:
// the tuple is sorted by precedence, the first one overwriting the following
struct Base { };
struct Grid { };
struct CCTpfaDisc { using InheritsFrom = std::tuple<Grid, Base>; };
struct BoxDisc { using InheritsFrom = std::tuple<Grid, Base>; };
struct OnePModel { using InheritsFrom = std::tuple<Base>; };
struct OnePTest { using InheritsFrom = std::tuple<OnePModel, BoxDisc>; };

} // end namespace TTag

// set and overwrite some properties
template<class TypeTag>
struct Scalar<TypeTag, TTag::Base> { using type = float; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::OnePModel> { using type = double; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::OnePTest> { using type = int; };

template<class TypeTag>
struct CoordinateType<TypeTag, TTag::Grid> { using type = GetPropType<TypeTag, Scalar>; };

template<class TypeTag>
struct UseTpfaFlux<TypeTag, TTag::CCTpfaDisc> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

int main(int argc, char* argv[])
{
    using namespace Dumux;
    using namespace Properties;

    {
        using Scalar = GetPropType<TTag::Base, Scalar>;
        if (!std::is_same_v<Scalar, float>)
            DUNE_THROW(Dune::InvalidStateException, "Property Scalar in TTag::Base should be float but is " << Dune::className<Scalar>());
    }
    {
        using Scalar = GetPropType<TTag::OnePTest, Scalar>;
        if (!std::is_same_v<Scalar, int>)
            DUNE_THROW(Dune::InvalidStateException, "Property Scalar in TTag::OnePTest should be int but is " << Dune::className<Scalar>());
    }
    {
        using Scalar = GetPropType<TTag::OnePModel, Scalar>;
        if (!std::is_same_v<Scalar, double>)
            DUNE_THROW(Dune::InvalidStateException, "Property Scalar in TTag::OnePModel should be double but is " << Dune::className<Scalar>());
    }
    {
        static_assert(
            !hasDefinedType<TTag::Base, CoordinateType>(),
            "Property type should be undefined for TTag::Base"
        );
    }
    {
        using CoordinateType = GetPropTypeOr<TTag::Base, CoordinateType, double>;
        static_assert(
            std::is_same_v<CoordinateType, double>,
            "Property is expected to default to double"
        );
    }
    {
        using CoordinateType = GetPropType<TTag::OnePTest, CoordinateType>;
        if (!std::is_same_v<CoordinateType, int>)
            DUNE_THROW(Dune::InvalidStateException, "Property CoordinateType in TTag::OnePTest should be int but is " << Dune::className<CoordinateType>());
    }
    {
        using CoordinateType = GetPropType<TTag::CCTpfaDisc, CoordinateType>;
        if (!std::is_same_v<CoordinateType, float>)
            DUNE_THROW(Dune::InvalidStateException, "Property CoordinateType in TTag::CCTpfaDisc should be float but is " << Dune::className<CoordinateType>());
    }
    {
        constexpr bool useTpfaFlux = getPropValue<TTag::CCTpfaDisc, UseTpfaFlux>();
        if (!useTpfaFlux)
            DUNE_THROW(Dune::InvalidStateException, "Property UseTpfaFlux in TTag::CCTpfaDisc should be true!");
    }

    static_assert(inheritsFrom<TTag::Base, TTag::CCTpfaDisc>());
    static_assert(!inheritsFrom<TTag::CCTpfaDisc, TTag::Base>());
    static_assert(inheritsFrom<TTag::Base, TTag::OnePTest>());
    static_assert(!inheritsFrom<TTag::OnePTest, TTag::Base>());

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
