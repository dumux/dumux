// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseFieldModels
 * \brief A degenerate mobility that vanishes in the pure phases.
 */
#ifndef DUMUX_PHASEFIELD_MOBILITY_DEGENERATE_HH
#define DUMUX_PHASEFIELD_MOBILITY_DEGENERATE_HH

#include <algorithm>

namespace Dumux::PhaseField::Mobility {

/*!
 * \ingroup PhaseFieldModels
 * \brief Degenerate mobility \f$ M(c) = M_0\, c (1-c) \f$, vanishing at
 *        \f$ c=0 \f$ and \f$ c=1 \f$ (the pure phases) and maximal at
 *        \f$ c=1/2 \f$.
 *
 * \f$ c \f$ is clamped to \f$ [0,1] \f$ before evaluation so the mobility
 * stays non-negative for solution values that slightly overshoot the
 * physical range (e.g. due to discretization).
 */
template<class Scalar>
class Degenerate
{
public:
    explicit Degenerate(Scalar mobilityScale)
    : mobilityScale_(mobilityScale)
    {}

    //! the mobility M(c)
    Scalar value(Scalar c) const
    {
        c = clamp_(c);
        return mobilityScale_*c*(1.0-c);
    }

    //! the derivative dM/dc
    Scalar derivative(Scalar c) const
    {
        c = clamp_(c);
        return mobilityScale_*(1.0-2.0*c);
    }

private:
    static Scalar clamp_(Scalar c)
    {
        using std::clamp;
        return clamp(c, Scalar(0.0), Scalar(1.0));
    }

    Scalar mobilityScale_;
};

} // end namespace Dumux::PhaseField::Mobility

#endif
