// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesIndices
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_INDICES_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_INDICES_HH

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The common indices for the isothermal Navier-Stokes model.
 *
 * \tparam dimension The dimension of the problem
 */
template <int dimension>
struct NavierStokesMomentumIndices // TODO specialize for staggered / diamond / etc
{
    static constexpr int dimXIdx = 0; //!< Index of the x-component of a vector of size dim
    static constexpr int dimYIdx = 1; //!< Index of the y-component of a vector of size dim
    static constexpr int dimZIdx = 2; //!< Index of the z-component of a vector of size dim

    static constexpr auto dim = dimension;

    static constexpr int momentumXBalanceIdx = 0; //!< Index of the momentum balance equation
    static constexpr int momentumYBalanceIdx = 1; //!< Index of the momentum balance equation
    static constexpr int momentumZBalanceIdx = 2; //!< Index of the momentum balance equation

    static constexpr int velocityXIdx = 0; //!< Index of the velocity in a solution vector
    static constexpr int velocityYIdx = 1; //!< Index of the velocity in a solution vector
    static constexpr int velocityZIdx = 2; //!< Index of the velocity in a solution vector

    /*!
     * \brief Index of the velocity in a solution vector given a certain direction.
     *
     * \param dirIdx The index of the direction.
     */
    static constexpr int velocity(int dirIdx)
    {
        return dirIdx;
    }

    /*!
     * \brief Index of the momentum balance equation given the direction
     *
     * \param dirIdx The index of the direction.
     */
    static constexpr int momentumBalanceIdx(int dirIdx)
    {
        return dirIdx;
    }
};

} // end namespace Dumux

#endif
