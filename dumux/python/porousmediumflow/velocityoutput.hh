// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief TODO: docme!
 */

#ifndef DUMUX_PYTHON_POROUSMEDIUMFLOW_VELOCITYOUTPUT_HH
#define DUMUX_PYTHON_POROUSMEDIUMFLOW_VELOCITYOUTPUT_HH

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/common/typeregistry.hh>
#include <dumux/porousmediumflow/velocityoutput.hh>

namespace Dumux::Python {

// Python wrapper for the PorousMediumFlowVelocityOutput class
template<class GridVariables, class FluxVariables, class... options>
void registerPorousMediumFlowVelocityOutput(pybind11::handle scope,
                                            pybind11::class_<Dumux::PorousMediumFlowVelocityOutput<GridVariables, FluxVariables>, options...> cls)
{
    using pybind11::operator""_a;
    using namespace Dune::Python;
    using VelocityOutput = Dumux::PorousMediumFlowVelocityOutput<GridVariables, FluxVariables>;

    cls.def(pybind11::init([](const GridVariables& gridVariables){
        return std::make_shared<VelocityOutput>(gridVariables);
    }));
}

} // end namespace Dumux::Python

#endif
