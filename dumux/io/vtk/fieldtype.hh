// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Vtk field types available in Dumux
 */
#ifndef DUMUX_IO_VTK_FIELD_TYPE_HH
#define DUMUX_IO_VTK_FIELD_TYPE_HH

namespace Dumux::Vtk {

/*!
 * \ingroup InputOutput
 * \brief Identifier for vtk field types.
 */
enum class FieldType : unsigned int
{
    element, vertex, automatic
};

} // end namespace Dumux::Vtk

#endif
