// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Core
 * \brief Class to specify the type of a boundary.
 */
#ifndef DUMUX_BOUNDARY_TYPES__HH
#define DUMUX_BOUNDARY_TYPES__HH

#include <algorithm>
#include <array>
#include <functional>

namespace Dumux::Experimental {

/*!
 * \ingroup Core
 * \brief Base class to specify the type of a boundary.
 *
 * The minimal interface of a boundary types class is to provide information if on a given face (scvf / boundaryFace)
 * boundary fluxes have be to evaluated. This corresponds to the `isFluxBoundary` function in this class.
 */
template <int numEq>
class BoundaryTypes
{
public:
    BoundaryTypes()
    { reset(); }

    //! we have a boundary condition for each equation
    static constexpr int size()
    { return numEq; }

    /*!
     * \brief Reset the boundary types for all equations.
     */
    void reset()
    { isFluxBoundary_.fill(false); }

    /*!
     * \brief Reset the boundary types for one equation.
     */
    void resetEq(int eqIdx)
    { isFluxBoundary_[eqIdx] = false; }

    /*!
     * \brief Set a flux boundary condition for a single equation.
     *
     * \param eqIdx The index of the equation
     */
    void setFluxBoundary(int eqIdx)
    {
        resetEq(eqIdx);
        isFluxBoundary_[eqIdx] = true;
    }

    /*!
     * \brief Set flux boundary conditions for all equations.
     */
    void setAllFluxBoundary()
    { isFluxBoundary_.fill(true); }

    /*!
     * \brief Returns true if an equation is used to specify a
     *        flux boundary condition.
     *
     * \param eqIdx The index of the equation
     */
    bool isFluxBoundary(unsigned eqIdx) const
    { return isFluxBoundary_[eqIdx]; }

    /*!
     * \brief Returns true if some equation is used to specify a
     *        flux boundary condition.
     */
    bool hasFluxBoundary() const
    { return std::ranges::any_of(isFluxBoundary_, std::identity{}); }

    /*!
     * \brief Returns true if all equations are used to specify a
     *        flux boundary condition.
     */
    bool hasOnlyFluxBoundary() const
    { return std::ranges::all_of(isFluxBoundary_, std::identity{}); }

protected:
    std::array<bool, numEq> isFluxBoundary_{};
};

} // end namespace Dumux::Experimental

#endif
