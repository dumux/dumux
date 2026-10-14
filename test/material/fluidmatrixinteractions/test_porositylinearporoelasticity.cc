// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup MaterialTests
 * \brief Test for the porosity of a linear poroelastic medium.
 */

#include <config.h>

#include <cmath>
#include <iostream>

#include <dune/common/exceptions.hh>

#include <dumux/material/fluidmatrixinteractions/porositylinearporoelasticity.hh>

namespace Dumux {

void checkClose(double value, double reference, double relTol, const std::string& what)
{
    using std::abs;
    if (abs(value - reference) > relTol*abs(reference))
        DUNE_THROW(Dune::Exception, what << ": " << value << " != " << reference);
}

} // end namespace Dumux

int main()
{
    using namespace Dumux;

    const double phi0 = 0.2;
    const double K = 3.3e9;
    const double Kf = 2.2e9;
    const double rhoFluid0 = 1000.0;

    // incompressible grains
    {
        const PorosityLinearPoroElasticity<double> law(phi0, 1.0, K);
        checkClose(law.evaluatePorosity(0.0, 0.0), phi0, 1e-14, "reference porosity");
        checkClose(law.evaluatePorosity(1e-3, 0.0), phi0 + 1e-3, 1e-14, "alpha = 1, strain");
        checkClose(law.evaluatePorosity(1e-3, 1e6), phi0 + 1e-3, 1e-14, "alpha = 1, pressure");
        checkClose(law.biotModulus(Kf), Kf/phi0, 1e-14, "alpha = 1, Biot modulus");
    }

    // compressible grains
    for (const double alpha : {0.2, 0.7, 0.95})
    {
        const PorosityLinearPoroElasticity<double> law(phi0, alpha, K);
        const double Ks = K/(1.0 - alpha);

        checkClose(law.evaluatePorosity(0.0, 0.0), phi0, 1e-14, "reference porosity");
        checkClose(law.evaluatePorosity(1e-3, 0.0) - phi0, alpha*1e-3, 1e-12, "strain coefficient");
        if (alpha > phi0)
            checkClose(law.evaluatePorosity(0.0, 1e6) - phi0, (alpha - phi0)/Ks*1e6, 1e-12, "pressure coefficient");
        else
            checkClose(law.evaluatePorosity(0.0, 1e6), phi0, 1e-14, "pressure coefficient");

        // undrained loading conserves the fluid mass per unit volume of the undeformed medium
        const double M = law.biotModulus(Kf);
        const double strain = -1e-6;
        const double pressure = -alpha*M*strain;
        const double rhoFluid = rhoFluid0*(1.0 + pressure/Kf);
        const double massChange = rhoFluid*law.evaluatePorosity(strain, pressure) - rhoFluid0*phi0;
        using std::abs;
        if (abs(massChange) > 1e-5*rhoFluid0*alpha*abs(strain))
            DUNE_THROW(Dune::Exception, "undrained mass change " << massChange << " for alpha = " << alpha);
    }

    // the Biot coefficient has to be in [phi0, 1]
    for (const double alpha : {0.1, 1.1})
    {
        bool thrown = false;
        try { PorosityLinearPoroElasticity<double> law(phi0, alpha, K); }
        catch (const Dune::RangeError&) { thrown = true; }
        if (!thrown)
            DUNE_THROW(Dune::Exception, "No exception for alpha = " << alpha);
    }

    std::cout << "All checks passed" << std::endl;
    return 0;
}
