// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCDiscretization
 * \brief Boundary types gathered on an element
 */
#ifndef DUMUX_CC_ELEMENT_BOUNDARY_TYPES_HH
#define DUMUX_CC_ELEMENT_BOUNDARY_TYPES_HH

namespace Dumux {

/*!
 * \ingroup CCDiscretization
 * \brief Boundary types gathered on an element
 * \note This class exists only for compatibility purposes with the
 *        box scheme. The cell-centered schemes and the box scheme use
 *        a common base local residual, which passes an ElementBoundaryTypes
 *        object to the implemented interfaces.
 */
class CCElementBoundaryTypes
{
public:
    /*!
     * \brief Update the boundary types for all vertices of an element.
     *
     * \param problem The problem object which needs to be simulated
     * \param element The DUNE Codim<0> entity for which the boundary
     *                types should be collected
     * \param fvGeometry The element finite volume geometry
     */
    template<class Problem, class Element, class FVElementGeometry>
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry)
    {}
};

} // namespace Dumux

#endif
