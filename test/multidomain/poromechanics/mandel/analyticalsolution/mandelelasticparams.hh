// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_TEST_MANDEL_ELASTIC_PARAMS_HH
#define DUMUX_TEST_MANDEL_ELASTIC_PARAMS_HH

#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup
 * \brief Parameters for the Mandel poroelastic problem and the analytical solution
 * \tparam Scalar
 */
template<class Scalar>
struct MandelPoroElasticParameters
{
    MandelPoroElasticParameters()
    {
        E_ = getParam<Scalar>("Problem.EModulus");
        nu_ = getParam<Scalar>("Problem.Nu");
        lambda_ = (E_ * nu_) / ((1.0 + nu_) * (1.0 - 2.0*nu_));
        G_ = E_ / 2.0 /(1.0 + nu_);
        Kb_ = E_ / 3.0 / (1.0 - 2.0*nu_);

        alpha_ = getParam<Scalar>("Problem.Alpha");
        M_ = getParam<Scalar>("Problem.MModulus");

        Kbu_ = Kb_ + alpha_*alpha_*M_;
        B_ = alpha_ * M_ / Kbu_;
        nuu_ = (3*Kbu_ - 2*G_) /2 /(3*Kbu_ + G_);
    };

    //! First Lamé parameter
    Scalar lambda() const { return lambda_; }

    //! Second Lamé parameter
    Scalar Nu() const { return nu_; }

    //! Young's modulus
    Scalar E() const { return E_; }

    //! Shear modulus
    Scalar G() const { return G_; }

    //! Bulk modulus
    Scalar Kb() const { return Kb_; }

    //! Biot's coefficient
    Scalar alpha() const { return alpha_; }

    //! Biot modulus
    Scalar M() const { return M_; }

    //! Skempton's coefficient
    Scalar B() const { return B_; }

    //! Undrained Poisson's ratio
    Scalar Nuu() const { return nuu_; }

    //! Undrained bulk modulus
    Scalar Kbu() const { return Kbu_; }

    //! Drained bulk modulus
    Scalar Kdr() const
    { return 2.0*lambda_; }

private:
    Scalar E_, nu_, G_, Kb_, lambda_, alpha_, M_, Kbu_, B_, nuu_;
};

} // end namespace Dumux

#endif
