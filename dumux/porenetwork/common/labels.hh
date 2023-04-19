// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Defines labels for pores and throats.
 */
#ifndef DUMUX_PNM_LABELS_HH
#define DUMUX_PNM_LABELS_HH

namespace Dumux::PoreNetwork {
// \{

/*!
 * \ingroup PoreNetworkModels
 * \brief Labels for pores and throats.
 */
struct Labels
{
    // TODO revise this concept
    static constexpr int interior = -1; //!< Label for pores/throats not on a boundary
    static constexpr int noflow = 0; //!< Label for pores/throats with a now flow BC
    static constexpr int dirichlet = 1; //!< Label for pores/throats with Dirichlet BC
    static constexpr int inlet = 2; //!< Label for pores/throats on an inlet
    static constexpr int outlet = 3; //!< Label for pores/throats on an outlet
    static constexpr int source = 4; //!< Label for pores/throats with a sink/source term
};

// \}
} // end Dumux::PoreNetwork

#endif
