// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief The volumetric strain of a displacement field
 */
#ifndef DUMUX_TEST_MANDEL_VOLUMETRIC_STRAIN_HH
#define DUMUX_TEST_MANDEL_VOLUMETRIC_STRAIN_HH

#include <dumux/discretization/evalgradients.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief The volumetric strain \f$\nabla \cdot \mathbf{u}\f$ at a position in an element
 * \note The displacement components are the first entries of the primary variables.
 */
template<class Element, class GridGeometry, class ElementSolution, class GlobalPosition>
auto volumetricStrain(const Element& element,
                      const GridGeometry& gridGeometry,
                      const ElementSolution& elemSol,
                      const GlobalPosition& globalPos)
{
    const auto gradU = evalGradients(element, element.geometry(), gridGeometry, elemSol, globalPos);
    typename GlobalPosition::value_type divU = 0.0;
    for (int dir = 0; dir < GlobalPosition::dimension; ++dir)
        divU += gradU[dir][dir];
    return divU;
}

} // end namespace Dumux

#endif
