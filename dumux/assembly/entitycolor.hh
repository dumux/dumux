// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \brief An enum class to define the colors of elements and vertices required
 *        for partial Jacobian reassembly.
 */
#ifndef DUMUX_ENTITY_COLOR_HH
#define DUMUX_ENTITY_COLOR_HH

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief The colors of elements and vertices required for partial
 *        Jacobian reassembly.
 */
enum class EntityColor {
    //! distance from last linearization is above the tolerance
    red,

    //! neighboring entity is red
    yellow,

    /*!
     * A yellow entity that has only non-green neighbor elements.
     *
     * This means that its relative error is below the tolerance,
     * but its defect can be linearized without any additional
     * cost.
     */
    orange,

    //! does not need to be reassembled
    green
};

} // end namespace Dumux

#endif
