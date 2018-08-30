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
 * \ingroup Common
 * \ingroup Tests
 * \brief Testing the Dumux property system
 */

#include <iostream>
#include <type_traits>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/properties/propertysystem.hh>

namespace Dumux {
namespace Properties {

// create some Properties (equivalent to old macro NEW_PROP_TAG(...))
// the first type tag is the actual TypeTag for which the property will be obtained
// (can be used to make properties depend on other properties),
// the second type tag is for parital specialization (equivalent to old macro SET_PROP(...), see below)
// the default property should be always undefined to produce a good error message
// if the user attempt to get an unset property
template<class TypeTag, class MyTypeTag>
struct Scalar { using type = UndefinedProperty;  };

template<class TypeTag, class MyTypeTag>
struct CoordinateType { using type = UndefinedProperty; };

namespace TTag {
// create some TypeTags (equivalent to old macro NEW_TYPE_TAG(..., INHERITS_FROM(...)))
// the tuple is sorted by precedence, the first one overwriting the following
struct Base { };
struct Grid { };
struct CCTpfaDisc { using InheritsFrom = std::tuple<Grid, Base>; };
struct BoxDisc { using InheritsFrom = std::tuple<Grid, Base>; };
struct OnePModel { using InheritsFrom = std::tuple<Base>; };
struct OnePTestTypeTag { using InheritsFrom = std::tuple<OnePModel, BoxDisc>; };

} // end namespace TTag

// set and overwrite some properties (equivalent to old macro SET_PROP(...){};)
template<class TypeTag>
struct Scalar<TypeTag, TTag::Base> { using type = float; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::OnePModel> { using type = double; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::OnePTestTypeTag> { using type = int; };

template<class TypeTag>
struct CoordinateType<TypeTag, TTag::Grid> { using type = GetPropType<TypeTag, Scalar>; };

} // end namespace Properties
} // end namespace Dumux

//! the main function
int main(int argc, char* argv[]) try
{
    using namespace Dumux;
    using namespace Properties;

    {
        using Scalar = GetPropType<TTag::Base, Scalar>;
        if (!std::is_same<Scalar, float>::value)
            DUNE_THROW(Dune::InvalidStateException, "Property Scalar in TTag::Base should be float but is " << Dune::className<Scalar>());
    }
    {
        using Scalar = GetPropType<TTag::OnePTestTypeTag, Scalar>;
        if (!std::is_same<Scalar, int>::value)
            DUNE_THROW(Dune::InvalidStateException, "Property Scalar in TTag::OnePTestTypeTag should be int but is " << Dune::className<Scalar>());
    }
    {
        using Scalar = GetPropType<TTag::OnePModel, Scalar>;
        if (!std::is_same<Scalar, double>::value)
        DUNE_THROW(Dune::InvalidStateException, "Property Scalar in TTag::OnePModel should be double but is " << Dune::className<Scalar>());
    }
    {
        using CoordinateType = GetPropType<TTag::OnePTestTypeTag, CoordinateType>;
        if (!std::is_same<CoordinateType, int>::value)
            DUNE_THROW(Dune::InvalidStateException, "Property CoordinateType in TTag::OnePTestTypeTag should be int but is " << Dune::className<CoordinateType>());
    }
    {
        using CoordinateType = GetPropType<TTag::CCTpfaDisc, CoordinateType>;
        if (!std::is_same<CoordinateType, float>::value)
            DUNE_THROW(Dune::InvalidStateException, "Property CoordinateType in TTag::CCTpfaDisc should be float but is " << Dune::className<CoordinateType>());
    }

    std::cout << "All tests passed!" << std::endl;
    return 0;
}

// error handler
catch (const Dune::Exception& e)
{
    std::cerr << "Dune exception thrown: " << e << " --> Abort!" << std::endl;
    return 1;
}

catch (...)
{
    std::cerr << "Unknown exception thrown --> Abort!" << std::endl;
    return 2;
}
