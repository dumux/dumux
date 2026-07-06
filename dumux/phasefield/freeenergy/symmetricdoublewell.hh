// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseFieldModels
 * \brief A smooth double-well free-energy density on the symmetric [-1,1] range.
 */
#ifndef DUMUX_PHASEFIELD_FREEENERGY_SYMMETRIC_DOUBLE_WELL_HH
#define DUMUX_PHASEFIELD_FREEENERGY_SYMMETRIC_DOUBLE_WELL_HH

namespace Dumux::PhaseField::FreeEnergy {

/*!
 * \ingroup PhaseFieldModels
 * \brief Quartic double-well free-energy density \f$ f(\phi) = \frac{E}{4} (\phi^2-1)^2 \f$,
 *        with minima at \f$ \phi=-1 \f$ and \f$ \phi=1 \f$ and a local maximum at \f$ \phi=0 \f$.
 *
 * This is the symmetric ([-1,1]-ranged) counterpart of `DoubleWell` (which is
 * defined on [0,1]), the convention typically used for order-parameter-style
 * phase fields (e.g. the Hele-Shaw two-phase model) rather than
 * concentration-style fields.
 */
template<class Scalar>
class SymmetricDoubleWell
{
public:
    explicit SymmetricDoubleWell(Scalar energyScale)
    : energyScale_(energyScale)
    {}

    //! the free-energy density f(phi)
    Scalar value(Scalar phi) const
    { return 0.25*energyScale_*(phi*phi - 1.0)*(phi*phi - 1.0); }

    //! the derivative df/dphi
    Scalar derivative(Scalar phi) const
    { return energyScale_*phi*(phi*phi - 1.0); }

    //! the second derivative d^2f/dphi^2
    Scalar secondDerivative(Scalar phi) const
    { return energyScale_*(3.0*phi*phi - 1.0); }

    /*!
     * \brief The derivative of the convex part of a convex-concave (Eyre)
     *        splitting \f$ f = f_{vex} - f_{cave} \f$, both convex.
     *
     * \f$ f''(\phi) = E(3\phi^2-1) \f$ has global minimum \f$ -E \f$ (at
     * \f$ \phi=0 \f$), so the quadratic stabilizer
     * \f$ f_{vex}(\phi) = f(\phi) + (E/2)\phi^2 \f$,
     * \f$ f_{cave}(\phi) = (E/2)\phi^2 \f$ makes \f$ f_{vex}'' = f''+E \geq 0 \f$
     * everywhere, using no parameter beyond the potential's own \f$ E \f$.
     */
    Scalar convexDerivative(Scalar phi) const
    { return derivative(phi) + energyScale_*phi; }

    //! the derivative of the concave part, see convexDerivative()
    Scalar concaveDerivative(Scalar phi) const
    { return energyScale_*phi; }

private:
    Scalar energyScale_;
};

} // end namespace Dumux::PhaseField::FreeEnergy

#endif
