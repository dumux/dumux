// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedCoupling
 * \brief Global dof indices and basis functions of the local dofs of an element
 */

#ifndef DUMUX_MULTIDOMAIN_EMBEDDED_ELEMENTDOFS_HH
#define DUMUX_MULTIDOMAIN_EMBEDDED_ELEMENTDOFS_HH

#include <type_traits>
#include <utility>
#include <vector>

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux::EmbeddedCoupling {

/*!
 * \ingroup EmbeddedCoupling
 * \brief The global dof indices of the local dofs of the bound element of a control-volume
 *        finite element or finite element discretization, in the order of the basis functions
 * \param elemDisc the element discretization (local view) bound to the element
 */
template<class ElementDiscretization>
auto elementDofIndices(const ElementDiscretization& elemDisc)
{
    using LocalDof = std::decay_t<decltype(*localDofs(std::declval<const ElementDiscretization&>()).begin())>;
    std::vector<typename LocalDof::GridIndexType> dofIndices(Dumux::Detail::LocalDofs::numLocalDofs(elemDisc));
    for (const auto& localDof : localDofs(elemDisc))
        dofIndices[localDof.index()] = localDof.dofIndex();
    return dofIndices;
}

/*!
 * \ingroup EmbeddedCoupling
 * \brief The values of the basis functions of an element at a position
 * \param gridDiscretization the grid geometry or grid discretization of a control-volume
 *        finite element or finite element discretization
 * \param geometry the geometry of the element
 * \param globalPos the position
 * \param shapeValues the values, in the order of the local dof indices
 */
template<class GridDiscretization, class Geometry, class GlobalPosition, class ShapeValues>
void basisFunctionValues(const GridDiscretization& gridDiscretization,
                         const Geometry& geometry,
                         const GlobalPosition& globalPos,
                         ShapeValues& shapeValues)
{
    const auto ipLocal = geometry.local(globalPos);
    gridDiscretization.feCache().get(geometry.type()).localBasis().evaluateFunction(ipLocal, shapeValues);
}

} // end namespace Dumux::EmbeddedCoupling

#endif
