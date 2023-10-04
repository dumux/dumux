// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Adaptive
 * \brief Interface to be used by classes transferring grid data on adaptive grids
 */
#ifndef DUMUX_ADAPTIVE_GRIDDATATRANSFER_HH
#define DUMUX_ADAPTIVE_GRIDDATATRANSFER_HH

namespace Dumux {

/*!
 * \ingroup Adaptive
 * \brief Interface to be used by classes transferring grid data on adaptive grids
 */
template<class Grid>
class GridDataTransfer
{
public:
    //! pure virtual base class needs virtual destructor
    virtual ~GridDataTransfer() = default;

    //! store user data before grid adaption
    virtual void store(const Grid&) = 0;

    //! store user data after grid adaption
    virtual void reconstruct(const Grid&) = 0;
};
} // end namespace Dumux

#endif
