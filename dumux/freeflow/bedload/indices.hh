// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTransportModel
 * \copydoc Dumux::BedloadIndices
 */
#ifndef DUMUX_BEDLOAD_INDICES_HH
#define DUMUX_BEDLOAD_INDICES_HH

namespace Dumux {

/*!
 * \ingroup BedloadTransportModel
 * \brief The common indices for the bedload transport model.
 */
struct BedloadIndices
{
    static constexpr int dimXIdx = 0; //!< Index of the x-component of a vector of size dim
    static constexpr int dimYIdx = 1; //!< Index of the y-component of a vector of size dim

    // The indices for the sediment balance equations and the primary variables cannot be defined here,
    // because the number of grain classes is a runtime parameter. Rather than using predefined indices,
    // it is looped over the number of grainclasses where the indices are needed.
};

} // end namespace Dumux

#endif
