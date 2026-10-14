// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEFlux
 * \brief Darcy's law with the Boussinesq approximation for control-volume finite element schemes
 */
#ifndef DUMUX_FLUX_CVFE_BOUSSINESQ_DARCYS_LAW_HH
#define DUMUX_FLUX_CVFE_BOUSSINESQ_DARCYS_LAW_HH

#include <dumux/flux/cvfe/darcyslaw.hh>

namespace Dumux {

/*!
 * \ingroup CVFEFlux
 * \brief Density in the gravity term of Darcy's law with the Boussinesq approximation
 *
 * The density in the gravity term is \f$\rho + \Delta\rho\f$. The phase density \f$\rho\f$
 * of the fluid system is the constant reference density, which is also used in the storage
 * and advection terms. The problem returns the density deviation \f$\Delta\rho\f$ in
 * \f$\mathrm{kg/m^3}\f$ that causes buoyancy in phase `phaseIdx` with `problem.buoyantDensity(volVars, phaseIdx)`, for example
 * \f$\Delta\rho = -\rho \beta_T (T - T_0)\f$ for a thermal expansion coefficient \f$\beta_T\f$
 * or \f$\Delta\rho = \rho \beta_X X\f$ for a solute mass fraction \f$X\f$.
 */
struct BoussinesqGravityTermDensity
{
    template<class Problem, class VolumeVariables>
    static auto density(const Problem& problem, const VolumeVariables& volVars, int phaseIdx)
    { return volVars.density(phaseIdx) + problem.buoyantDensity(volVars, phaseIdx); }
};

/*!
 * \ingroup CVFEFlux
 * \brief Darcy's law with the Boussinesq approximation for control-volume finite element schemes
 *
 * Assembles \f$-|\sigma| \mathbf{n}^T \mathbf{K} \left( \nabla p - (\rho + \Delta\rho) \mathbf{g} \right)\f$,
 * see BoussinesqGravityTermDensity.
 */
template<class Scalar, class GridGeometry>
using CVFEBoussinesqDarcysLaw = CVFEDarcysLaw<Scalar, GridGeometry, BoussinesqGravityTermDensity>;

} // end namespace Dumux

#endif
