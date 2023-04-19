// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Convenience header that includes all grid manager specializations
 *
 * \todo add Petrel grids with dune-cornerpoint
 */
#ifndef DUMUX_IO_GRID_MANAGER_HH
#define DUMUX_IO_GRID_MANAGER_HH

#include <dumux/io/grid/gridmanager_base.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_oned.hh>
#include <dumux/io/grid/gridmanager_ug.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/io/grid/gridmanager_foam.hh>
#include <dumux/io/grid/gridmanager_sp.hh>
#include <dumux/io/grid/gridmanager_mmesh.hh>
#include <dumux/io/grid/gridmanager_sub.hh>

#endif
