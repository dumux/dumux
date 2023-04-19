// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Adds output fields to a given output module, this is the default if a
          model doesn't implement this functionality
 */
#ifndef DUMUX_IO_DEFAULT_IO_FIELDS_HH
#define DUMUX_IO_DEFAULT_IO_FIELDS_HH

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Adds output fields to a given output module
 */
class DefaultIOFields
{
public:
    template<class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        DUNE_THROW(Dune::NotImplemented, "This model doesn't implement default output fields!");
    }

    template <class FluidSystem = void, class SolidSystem = void>
    static std::string primaryVariableName(int pvIdx = 0, int state = 0)
    {
        DUNE_THROW(Dune::NotImplemented, "This model doesn't implement primaryVariableName!");
    }
};

} // end namespace Dumux

#endif
