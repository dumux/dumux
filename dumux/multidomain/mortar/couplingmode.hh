// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MortarCoupling
 * \brief Mode in which mortar data enters the subdomain problems.
 */
#ifndef DUMUX_MULTIDOMAIN_MORTAR_COUPLING_MODE_HH
#define DUMUX_MULTIDOMAIN_MORTAR_COUPLING_MODE_HH

namespace Dumux::Mortar {

/*!
 * \ingroup MortarCoupling
 * \brief Mode in which mortar data enters the subdomain problems.
 *
 * In essential mode the mortar datum is imposed as a Dirichlet condition and the conjugate
 * (flux) trace is read back, so the mortar plays the role of an interface value (e.g. a
 * pressure); in natural mode it is imposed as a Neumann condition and the value trace is
 * read back, so the mortar plays the role of an interface flux, which is the flux-mortar
 * method \cite Boon2022 \cite Boon2023. An interface operator running in one mode has the
 * conjugate mode as its preconditioner.
 */
enum class CouplingMode { essential, natural };

} // end namespace Dumux::Mortar

#endif
