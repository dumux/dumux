// $Id:$
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: andreas.lauser _at_ iws.uni-stuttgart.de                         *
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
 * \brief This file tests the properties system.
 *
 * We define a few type tags and property tags, then we attach values
 * to (TypeTag, PropertyTag) tuples and finally we use them in the
 * main function and print some diagnostic messages.
 */

#include <dumux/auxiliary/properties.hh>

namespace Dune {
namespace Properties {

///////////////////
// Define some hierarchy of type tags:
//
//  MyBase -- MyDerived --_
//              |          \ a
//              |            - MyDerived2 -_
//               \         /                \ a
//  MyBase2 --------------^                  \ a
//                 \__________________________ MoreDerived
//                                           / 
//                                          / 
//  MyBase3 -------------------------------^
///////////////////
NEW_TYPE_TAG(MyBase);
NEW_TYPE_TAG(MyBase2);
NEW_TYPE_TAG(MyBase3);

NEW_TYPE_TAG(MyDerived, INHERITS_FROM(MyBase));
NEW_TYPE_TAG(MyDerived2, INHERITS_FROM(MyBase2, MyDerived));

NEW_TYPE_TAG(MoreDerived, INHERITS_FROM(MyDerived2, MyDerived, MyBase3));

///////////////////
// Define the property tags: foo, bar, fooBar
///////////////////
NEW_PROP_TAG(foo);
NEW_PROP_TAG(bar);
NEW_PROP_TAG(fooBar);

///////////////////
// Define some values for the properties on the type tags:
//
// (MyBase, foo) = 0
// (MyBase3, foo) = 1
// (MyDerived, foo) = 2
// (MyDerived2, foo) = 3
//
// bar defaults to foo + 10, but is not defined for MyBase
// (MyBase2, bar) = 5
//
// fooBar defaults to bar + 1000
///////////////////

// some simple constant for the foo property
SET_INT_PROP(MyBase, foo, 0);
SET_INT_PROP(MyBase3, foo, 1);
SET_INT_PROP(MyDerived, foo, 2);
SET_INT_PROP(MyDerived2, foo, 3);

// Some default value for the bar property
SET_PROP_DEFAULT(bar)
{
    static const int value = GET_PROP_VALUE(TypeTag, PTAG(foo)) + 10;
};

// unset bar for MyBase
UNSET_PROP(MyBase, bar);

// use a different value for bar on MyBase2
SET_INT_PROP(MyBase2, bar, 5);

// A property which depends on a property only defined on a derived
// type tag
SET_PROP_DEFAULT(fooBar)
{
    static const int value = GET_PROP_VALUE(TypeTag, PTAG(bar)) + 1000;
};

}; // namespace Properties
}; // namespace Dune


int main()
{
    // print all properties for all type tags
    std::cout << "---------------------------------------\n";
    std::cout << "-- Property values\n";
    std::cout << "---------------------------------------\n";


    std::cout << "---------- Values for MyBase ----------\n";
    
    std::cout << "(MyBase, foo) = " << GET_PROP_VALUE(TTAG(MyBase), PTAG(foo)) << "\n";
    /* explcitly unset -> doesn't compile!
     */
    //std::cout << "(MyBase, bar) = " << GET_PROP_VALUE(TTAG(MyBase), PTAG(bar)) << "\n";
    /* fooBar is bar+10 by default, but bar was not set for MyBase ->
       doesn't compile
    */
    // std::cout << "(MyBase, fooBar) = " << GET_PROP_VALUE(TTAG(MyBase), PTAG(fooBar)) << "\n";

    std::cout << "---------- Values for MyBase2 ---------\n";
    /* foo was not set anywhere for MyBase2 -> doesn't compile
     */
    // std::cout << "(MyBase2, foo) = " << GET_PROP_VALUE(TTAG(MyBase2), PTAG(foo)) << "\n";
    std::cout << "(MyBase2, bar) = " << GET_PROP_VALUE(TTAG(MyBase2), PTAG(bar)) << "\n";
    std::cout << "(MyBase2, fooBar) = " << GET_PROP_VALUE(TTAG(MyBase2), PTAG(fooBar)) << "\n";

    std::cout << "--------- Values for MyBase3 ----------\n";
    std::cout << "(MyBase3, foo) = " << GET_PROP_VALUE(TTAG(MyBase3), PTAG(foo)) << "\n";
    std::cout << "(MyBase3, bar) = " << GET_PROP_VALUE(TTAG(MyBase3), PTAG(bar)) << "\n";
    std::cout << "(MyBase3, fooBar) = " << GET_PROP_VALUE(TTAG(MyBase3), PTAG(fooBar)) << "\n";

    std::cout << "--------- Values for MyDerived --------\n";
    std::cout << "(MyDerived, foo) = " << GET_PROP_VALUE(TTAG(MyDerived), PTAG(foo)) << "\n";
    /* 'bar' explicitly unset for only parent 'MyBase' -> doesn't
     * compile
     */
    //std::cout << "(MyDerived, bar) = " << GET_PROP_VALUE(TTAG(MyDerived), PTAG(bar)) << "\n";
    /* 'fooBar' requires 'bar' to be defined, but is unset for
     * 'MyDerived' -> doesn't compile
     */
    //std::cout << "(MyDerived, fooBar) = " << GET_PROP_VALUE(TTAG(MyDerived), PTAG(fooBar)) << "\n";
    
    std::cout << "-------- Values for MyDerived2 --------\n";
    std::cout << "(MyDerived2, foo) = " << GET_PROP_VALUE(TTAG(MyDerived2), PTAG(foo)) << "\n";
    std::cout << "(MyDerived2, bar) = " << GET_PROP_VALUE(TTAG(MyDerived2), PTAG(bar)) << "\n";
    std::cout << "(MyDerived2, fooBar) = " << GET_PROP_VALUE(TTAG(MyDerived2), PTAG(fooBar)) << "\n";
    
    std::cout << "-------- Values for MoreDerived -------\n";
    std::cout << "(MoreDerived, foo) = " << GET_PROP_VALUE(TTAG(MoreDerived), PTAG(foo)) << "\n";
    std::cout << "(MoreDerived, bar) = " << GET_PROP_VALUE(TTAG(MoreDerived), PTAG(bar)) << "\n";
    std::cout << "(MoreDerived, fooBar) = " << GET_PROP_VALUE(TTAG(MoreDerived), PTAG(fooBar)) << "\n";
    
    std::cout << "\n";
    std::cout << "---------------------------------------\n";
    std::cout << "-- Diagnostic messages\n";
    std::cout << "---------------------------------------\n";
    std::cout << "---- Message for (MoreDerived, foo) ---\n"
              << PROP_DIAGNOSTIC(TTAG(MoreDerived), PTAG(foo));
    std::cout << "---- Message for (MoreDerived, bar) ---\n"
              << PROP_DIAGNOSTIC(TTAG(MoreDerived), PTAG(bar));
    std::cout << "-- Message for (MoreDerived, fooBar) --\n"
              << PROP_DIAGNOSTIC(TTAG(MoreDerived), PTAG(fooBar));

 
    return 0;
};
