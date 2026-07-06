// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseFieldModels
 * \brief A smooth double-well free-energy density.
 */
#ifndef DUMUX_PHASEFIELD_FREEENERGY_DOUBLE_WELL_HH
#define DUMUX_PHASEFIELD_FREEENERGY_DOUBLE_WELL_HH

namespace Dumux::PhaseField::FreeEnergy {

/*!
 * \ingroup PhaseFieldModels
 * \brief Quartic double-well free-energy density \f$ f(c) = E c^2 (1-c)^2 \f$,
 *        with minima at \f$ c=0 \f$ and \f$ c=1 \f$ and a local maximum at \f$ c=1/2 \f$.
 */
template<class Scalar>
class DoubleWell
{
public:
    explicit DoubleWell(Scalar energyScale)
    : energyScale_(energyScale)
    {}

    //! the free-energy density f(c)
    Scalar value(Scalar c) const
    { return energyScale_*c*c*(1.0-c)*(1.0-c); }

    //! the derivative df/dc
    Scalar derivative(Scalar c) const
    { return energyScale_*2.0*c*(2.0*c*c - 3.0*c + 1.0); }

    //! the second derivative d^2f/dc^2
    Scalar secondDerivative(Scalar c) const
    { return energyScale_*(12.0*c*c - 12.0*c + 2.0); }

    /*!
     * \brief The derivative of the convex part of a convex-concave (Eyre)
     *        splitting \f$ f = f_{vex} - f_{cave} \f$, both convex.
     *
     * \f$ f''(c) = E(12c^2-12c+2) = 12E(c-1/2)^2 - E \f$ has global minimum
     * \f$ -E \f$ (at \f$ c=1/2 \f$), so the quadratic stabilizer
     * \f$ f_{vex}(c) = f(c) + (E/2)c^2 \f$, \f$ f_{cave}(c) = (E/2)c^2 \f$
     * makes \f$ f_{vex}'' = f''+E \geq 0 \f$ everywhere, using no parameter
     * beyond the potential's own \f$ E \f$.
     */
    Scalar convexDerivative(Scalar c) const
    { return derivative(c) + energyScale_*c; }

    //! the derivative of the concave part, see convexDerivative()
    Scalar concaveDerivative(Scalar c) const
    { return energyScale_*c; }

private:
    Scalar energyScale_;
};

} // end namespace Dumux::PhaseField::FreeEnergy

#endif
