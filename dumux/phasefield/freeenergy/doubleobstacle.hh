// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseFieldModels
 * \brief A double-obstacle free-energy density.
 */
#ifndef DUMUX_PHASEFIELD_FREEENERGY_DOUBLE_OBSTACLE_HH
#define DUMUX_PHASEFIELD_FREEENERGY_DOUBLE_OBSTACLE_HH

namespace Dumux::PhaseField::FreeEnergy {

/*!
 * \ingroup PhaseFieldModels
 * \brief Double-obstacle free-energy density \f$ f(c) = \frac{E}{2} c (1-c) \f$
 *        (Blowey & Elliott), the smooth counterpart of the double well.
 *
 * The physical double-obstacle potential additionally restricts \f$ c \f$ to
 * \f$ [0,1] \f$ via a variational inequality; that constraint is not enforced
 * by this class and must be imposed by the model/solver if required.
 */
template<class Scalar>
class DoubleObstacle
{
public:
    explicit DoubleObstacle(Scalar energyScale)
    : energyScale_(energyScale)
    {}

    //! the free-energy density f(c)
    Scalar value(Scalar c) const
    { return 0.5*energyScale_*c*(1.0-c); }

    //! the derivative df/dc
    Scalar derivative(Scalar c) const
    { return 0.5*energyScale_*(1.0-2.0*c); }

    //! the second derivative d^2f/dc^2
    Scalar secondDerivative(Scalar c) const
    { return -energyScale_; }

private:
    Scalar energyScale_;
};

} // end namespace Dumux::PhaseField::FreeEnergy

#endif
