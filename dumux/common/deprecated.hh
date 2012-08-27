// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Provides macros for deprecating classes, functions and attributes.
 *
 * This file is a almost verbatim copy of Dune's dune/common/version.hh.
 *
 * \deprecated DUMUX_DEPRECATED_MSG is deprecated, use DUNE_DEPRECATED_MSG instead.
 */
#ifndef DUMUX_DEPRECATED_HH
#define DUMUX_DEPRECATED_HH

#ifndef DUMUX_DEPRECATED_MSG
#if defined(DOXYGEN) || !HAVE_ATTRIBUTE_DEPRECATED_MSG
//! Mark some entity as deprecated
#define DUMUX_DEPRECATED_MSG(text) DUNE_DEPRECATED
#else
#define DUMUX_DEPRECATED_MSG(text) __attribute__((deprecated(text)))
#endif
#endif

#endif
