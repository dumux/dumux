// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Contact
 * \brief Contact pressure as a continuously differentiable function of the gap
 */
#ifndef DUMUX_SOLIDMECHANICS_CONTACT_PENALTY_HH
#define DUMUX_SOLIDMECHANICS_CONTACT_PENALTY_HH

namespace Dumux::Contact {

/*!
 * \ingroup Contact
 * \brief Penalty law with offset
 *
 * With the stiffness \f$ k \f$ and the offset \f$ \delta > 0 \f$, the contact pressure is
 * \f[
 *     p(g) = \begin{cases}
 *         0 & g > \delta, \\
 *         k (\delta - g)^2 / (2\delta) & 0 \leq g \leq \delta, \\
 *         k (\delta/2 - g) & g < 0,
 *     \end{cases}
 * \f]
 * so contact engages at the distance \f$ \delta \f$ and pressure and stiffness are continuous in
 * the gap \f$ g \f$. The energy per unit area is \f$ \psi(g) = \int_g^\delta p(s)\,\mathrm{d}s \f$.
 */
template<class Scalar>
class PenaltyLaw
{
public:
    PenaltyLaw(Scalar stiffness, Scalar offset)
    : k_(stiffness), delta_(offset)
    {}

    Scalar pressure(Scalar gap) const
    {
        if (gap >= delta_)
            return 0.0;
        else if (gap >= 0.0)
            return 0.5*k_*(delta_ - gap)*(delta_ - gap)/delta_;
        else
            return k_*(0.5*delta_ - gap);
    }

    //! The derivative of the pressure with respect to the gap
    Scalar pressureDerivative(Scalar gap) const
    {
        if (gap >= delta_)
            return 0.0;
        else if (gap >= 0.0)
            return -k_*(delta_ - gap)/delta_;
        else
            return -k_;
    }

    Scalar energyDensity(Scalar gap) const
    {
        if (gap >= delta_)
            return 0.0;
        else if (gap >= 0.0)
            return k_*(delta_ - gap)*(delta_ - gap)*(delta_ - gap)/(6.0*delta_);
        else
            return k_*(delta_*delta_/6.0 - 0.5*delta_*gap + 0.5*gap*gap);
    }

    Scalar stiffness() const
    { return k_; }

    Scalar offset() const
    { return delta_; }

private:
    Scalar k_, delta_;
};

} // end namespace Dumux::Contact

#endif
