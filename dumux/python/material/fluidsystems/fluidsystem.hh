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

#ifndef DUMUX_PYTHON_MATERIAL_FLUIDSYSTEMS_FLUIDSYSTEM_HH
#define DUMUX_PYTHON_MATERIAL_FLUIDSYSTEMS_FLUIDSYSTEM_HH

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

namespace Dumux::Python {

template <class FS, class... options>
void registerFluidSystem(pybind11::handle scope,
                         pybind11::class_<FS, options...> cls)
{
    cls.def(pybind11::init());
}

} // namespace Dumux::Python

#endif
