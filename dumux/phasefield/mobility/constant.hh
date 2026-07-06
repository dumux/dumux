// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseFieldModels
 * \brief A constant (concentration-independent) mobility.
 */
#ifndef DUMUX_PHASEFIELD_MOBILITY_CONSTANT_HH
#define DUMUX_PHASEFIELD_MOBILITY_CONSTANT_HH

namespace Dumux::PhaseField::Mobility {

/*!
 * \ingroup PhaseFieldModels
 * \brief Constant mobility \f$ M(c) = M_0 \f$.
 */
template<class Scalar>
class Constant
{
public:
    explicit Constant(Scalar mobilityScale)
    : mobilityScale_(mobilityScale)
    {}

    //! the mobility M(c)
    Scalar value(Scalar /*c*/) const
    { return mobilityScale_; }

    //! the derivative dM/dc
    Scalar derivative(Scalar /*c*/) const
    { return 0.0; }

private:
    Scalar mobilityScale_;
};

} // end namespace Dumux::PhaseField::Mobility

#endif
