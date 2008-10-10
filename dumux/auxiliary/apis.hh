/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
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
 * \file Apis.hh
 *
 * This file provides tools to make sure a template class adheres to an
 * abstract API
 */
#ifndef APIS_HH
#define APIS_HH

namespace Dune
{
namespace Api
{
///////////////////////////////////////////////////////////////////////////////
// Low level stuff to make the API checks work.
//
// DON'T LOOK AT IT, IT'S _UGLY_!
///////////////////////////////////////////////////////////////////////////////
    class Interface { };

#define BEGIN_API_DEF(InterfaceName) \
    class InterfaceName : public Interface { }; \
    template <typename Implementation> \
    class CheckInterface<InterfaceName, Implementation>      \
    { \
    public: \
        CheckInterface(Implementation &impl, const Implementation &const_impl) {

#define END_API_DEF } };

    template <class Interface, class Implementation>
    class CheckInterface
    {
    public:
        CheckInterface(Implementation &impl,
                       const Implementation &const_impl)
        {
            // the requested API has not been declared in the
            // scope where you called require(). Look at the examples
            // below and declare your API the same way!
            //
            // To find out which API is missing, look at the
            // the source location which tried to instantiate
            // the Api::require template function. (Said location
            // can be deduced from the compiler error, at least
            // when using g++.)
            impl.ProgrammingError_RequestedApiUnknown__SeeCommentAbove();
        };
    };

///////////////////////////////////////////////////////////////////////////////
// End of the low level stuff
// (you may open your eyes again)
///////////////////////////////////////////////////////////////////////////////
    // make sure an object adheres to an interface definition
    template<class Interface, class Implementation>
    void require(const Implementation &const_impl)
    {
        while (0) {
            CheckInterface<Interface, Implementation> dummy((Implementation &) const_impl,
                                                            const_impl);
        }
    }

    // make sure a type adheres to an interface definition
    template<class Interface, class Implementation>
    void require()
    {
        while (0) {
            Implementation *dummyImpl;
            CheckInterface<Interface, Implementation> dummy(*dummyImpl,
                                                            *dummyImpl);
        }
    }

} // namespace Api
} // namespace Dune

#endif // APIS_HH
