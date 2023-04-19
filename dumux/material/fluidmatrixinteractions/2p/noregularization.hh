// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief   A tag to turn off regularization and it's overhead
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_TWOP_NO_REGULARIZATION_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_TWOP_NO_REGULARIZATION_HH

namespace Dumux::FluidMatrix {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A tag to turn off regularization and it's overhead
 */
struct NoRegularization
{
    //! Empty parameter structure
    template<class S> struct Params {};

    //! We are always equal to other instances of our kind
    bool operator== (const NoRegularization& o) const
    { return true; }
};

} // end namespace Dumux::FluidMatrix

#endif
