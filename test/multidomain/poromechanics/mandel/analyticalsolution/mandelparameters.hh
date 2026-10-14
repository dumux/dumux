// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Material parameters of Mandel's problem
 */
#ifndef DUMUX_TEST_MANDEL_PARAMETERS_HH
#define DUMUX_TEST_MANDEL_PARAMETERS_HH

#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief Material parameters of Mandel's problem
 *
 * The parameters of linear poroelasticity are derived from the drained Young's modulus
 * \f$E\f$ and Poisson's ratio \f$\nu\f$ of the solid skeleton, the Biot coefficient
 * \f$\alpha\f$, the initial porosity \f$\phi_0\f$ and the fluid bulk modulus \f$K_f\f$.
 * The bulk modulus of the solid grains is \f$K_s = K/(1-\alpha)\f$, which gives the Biot modulus
 * \f$1/M = \phi_0/K_f + (\alpha - \phi_0)/K_s\f$ \cite Coussy2004.
 */
template<class Scalar>
class MandelParameters
{
public:
    MandelParameters()
    {
        const Scalar E = getParam<Scalar>("Problem.EModulus");
        nu_ = getParam<Scalar>("Problem.Nu");
        alpha_ = getParam<Scalar>("Problem.BiotCoefficient");
        phi0_ = getParam<Scalar>("Problem.InitialPorosity");
        Kf_ = getParam<Scalar>("Problem.FluidBulkModulus");

        G_ = E/(2.0*(1.0 + nu_));
        K_ = E/(3.0*(1.0 - 2.0*nu_));
        M_ = 1.0/(phi0_/Kf_ + (alpha_ - phi0_)*(1.0 - alpha_)/K_);
        Ku_ = K_ + alpha_*alpha_*M_;
        B_ = alpha_*M_/Ku_;
        nuu_ = (3.0*Ku_ - 2.0*G_)/(2.0*(3.0*Ku_ + G_));
    }

    //! Drained Poisson's ratio \f$\nu\f$
    Scalar poissonRatio() const { return nu_; }

    //! Undrained Poisson's ratio \f$\nu_u\f$
    Scalar undrainedPoissonRatio() const { return nuu_; }

    //! Shear modulus \f$G\f$
    Scalar shearModulus() const { return G_; }

    //! Drained bulk modulus \f$K\f$
    Scalar bulkModulus() const { return K_; }

    //! Undrained bulk modulus \f$K_u = K + \alpha^2 M\f$
    Scalar undrainedBulkModulus() const { return Ku_; }

    //! Biot coefficient \f$\alpha\f$
    Scalar biotCoefficient() const { return alpha_; }

    //! Biot modulus \f$M\f$
    Scalar biotModulus() const { return M_; }

    //! Skempton coefficient \f$B = \alpha M / K_u\f$
    Scalar skemptonCoefficient() const { return B_; }

    //! Porosity \f$\phi_0\f$ in the initial state
    Scalar initialPorosity() const { return phi0_; }

    //! Fluid bulk modulus \f$K_f\f$
    Scalar fluidBulkModulus() const { return Kf_; }

private:
    Scalar nu_, alpha_, phi0_, Kf_;
    Scalar G_, K_, M_, Ku_, B_, nuu_;
};

} // end namespace Dumux

#endif
