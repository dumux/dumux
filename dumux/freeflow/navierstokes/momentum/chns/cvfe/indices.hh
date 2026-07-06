// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMomentumCHNSCVFEIndices
 */
#ifndef DUMUX_NAVIERSTOKES_MOMENTUM_CHNS_CVFE_INDICES_HH
#define DUMUX_NAVIERSTOKES_MOMENTUM_CHNS_CVFE_INDICES_HH

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Indices for the co-located momentum + Cahn-Hilliard CVFE model.
 *
 * The primary-variable / equation layout is [ velocity_0 ... velocity_{dim-1}, c, mu ],
 * i.e. the dim momentum-balance components (Taylor-Hood P2 velocity) followed by the
 * phase field c and the chemical potential mu, both sharing the same P2 degree as the
 * velocity (Aland & Voigt P2-P1-P2-P2). numEq = dim + 2. The pressure/continuity lives
 * in a separate P1 (Box) subdomain, so there is NO conti0EqIdx here.
 *
 * \tparam dimension The dimension of the problem
 */
template <int dimension>
struct NavierStokesMomentumCHNSCVFEIndices
{
    static constexpr int dim = dimension;

    static constexpr int dimXIdx = 0; //!< Index of the x-component of a vector of size dim
    static constexpr int dimYIdx = 1; //!< Index of the y-component of a vector of size dim
    static constexpr int dimZIdx = 2; //!< Index of the z-component of a vector of size dim

    //! momentum balance equations (one per space dimension)
    static constexpr int momentumXBalanceIdx = 0;
    static constexpr int momentumYBalanceIdx = 1;
    static constexpr int momentumZBalanceIdx = 2;

    //! velocity primary-variable components (one per space dimension)
    static constexpr int velocityXIdx = 0;
    static constexpr int velocityYIdx = 1;
    static constexpr int velocityZIdx = 2;

    //! phase field c and chemical potential mu, appended after the dim velocity components
    static constexpr int phaseFieldIdx = dimension;          //!< primary variable index of c
    static constexpr int chemicalPotentialIdx = dimension+1; //!< primary variable index of mu
    static constexpr int phaseFieldEqIdx = dimension;         //!< equation index of the c-transport eq
    static constexpr int chemicalPotentialEqIdx = dimension+1;//!< equation index of the mu-definition eq

    //! Index of the velocity in a solution vector given a certain direction.
    static constexpr int velocity(int dirIdx)
    { return dirIdx; }

    //! Index of the momentum balance equation given the direction
    static constexpr int momentumBalanceIdx(int dirIdx)
    { return dirIdx; }
};

} // end namespace Dumux

#endif
