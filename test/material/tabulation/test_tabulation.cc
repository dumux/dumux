// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup MaterialTests
 * \brief This is a program to test the tabulation class for of
 *        individual components.
 *
 * It either prints "success" or "error", it does not do anything
 * else.
 */

#include <config.h>

#include <dune/common/float_cmp.hh>

#include <dumux/material/components/air.hh>
#include <dumux/material/components/benzene.hh>
#include <dumux/material/components/brine.hh>
#include <dumux/material/components/calcite.hh>
#include <dumux/material/components/calciumion.hh>
#include <dumux/material/components/cao.hh>
#include <dumux/material/components/cao2h2.hh>
#include <dumux/material/components/carbonateion.hh>
#include <dumux/material/components/ch4.hh>
#include <dumux/material/components/co2.hh>
#include <dumux/material/components/granite.hh>
#include <dumux/material/components/h2.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/heavyoil.hh>
#include <dumux/material/components/mesitylene.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/components/nacl.hh>
#include <dumux/material/components/o2.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/trichloroethene.hh>
#include <dumux/material/components/xylene.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/componenttraits.hh>

template<class Scalar>
void checkEquality(const std::string& name, Scalar tab, Scalar real, Scalar eps)
{
    if (!Dune::FloatCmp::eq(tab, real, eps))
        DUNE_THROW(Dune::InvalidStateException, "Tabulated and computed value for "
                             << name << " differs by more than " << eps*100 << "% "
                             << "(tabulated: " << tab << ", " << "actual value: " << real << ")");
}

int main(int argc, char *argv[])
{
    using namespace Dumux;
    using Scalar = double;

    // test IapwsH2O in detail
    {
        using IapwsH2O = Components::H2O<Scalar>;
        using TabulatedH2O = Components::TabulatedComponent<IapwsH2O>;

        // tabulation min/max
        const Scalar tempMin = 274.15;
        const Scalar tempMax = 622.15;
        const int nTemp = 800;

        const Scalar pMin = 1e4;
        const Scalar pMax = IapwsH2O::vaporPressure(tempMax*1.1);
        const int nPress = 200;

        std::cout << "Creating tabulation with " << nTemp*nPress << " entries per quantity\n";
        TabulatedH2O::init(tempMin, tempMax, nTemp, pMin, pMax, nPress);

        constexpr Scalar eps = 1e-1;
        std::cout << "Checking tabulation\n";
        const int m = nTemp*3;
        const int n = nPress*3;
        for (int i = 0; i < m; ++i)
        {
            const Scalar T = tempMin + (tempMax - tempMin)*Scalar(i)/Scalar(m-1);
            checkEquality("vaporPressure", TabulatedH2O::vaporPressure(T), IapwsH2O::vaporPressure(T), eps);

            for (int j = 0; j < n; ++j)
            {
                const Scalar p = pMin + (pMax - pMin)*Scalar(j)/Scalar(n-1);

                if (p < IapwsH2O::vaporPressure(T) * 1.001) {
                    Scalar rho = IapwsH2O::gasDensity(T,p);
                    checkEquality("Iapws::gasPressure", IapwsH2O::gasPressure(T,rho), p, eps);
                    checkEquality("gasPressure", TabulatedH2O::gasPressure(T,rho), p, eps);
                    checkEquality("gasEnthalpy", TabulatedH2O::gasEnthalpy(T,p), IapwsH2O::gasEnthalpy(T,p), eps);
                    checkEquality("gasInternalEnergy", TabulatedH2O::gasInternalEnergy(T,p), IapwsH2O::gasInternalEnergy(T,p), eps);
                    checkEquality("gasDensity", TabulatedH2O::gasDensity(T,p), rho, eps);
                    checkEquality("gasViscosity", TabulatedH2O::gasViscosity(T,p), IapwsH2O::gasViscosity(T,p), eps);
                }

                if (p > IapwsH2O::vaporPressure(T) / 1.001) {
                    Scalar rho = IapwsH2O::liquidDensity(T,p);
                    checkEquality("Iapws::liquidPressure", IapwsH2O::liquidPressure(T,rho), p, eps);
                    checkEquality("liquidPressure", TabulatedH2O::liquidPressure(T,rho), p, eps);
                    checkEquality("liquidEnthalpy", TabulatedH2O::liquidEnthalpy(T,p), IapwsH2O::liquidEnthalpy(T,p), eps);
                    checkEquality("liquidInternalEnergy", TabulatedH2O::liquidInternalEnergy(T,p), IapwsH2O::liquidInternalEnergy(T,p), eps);
                    checkEquality("liquidDensity", TabulatedH2O::liquidDensity(T,p), rho, eps);
                    checkEquality("liquidViscosity", TabulatedH2O::liquidViscosity(T,p), IapwsH2O::liquidViscosity(T,p), eps);
                }


            }

        }
    }

    // test if other components can be tabulated
    {
        Components::TabulatedComponent<Components::Air<Scalar>, false>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::Benzene<Scalar>, false>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::Calcite<Scalar>, false>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::CalciumIon<Scalar>, false>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::CaO<Scalar>, false>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::CaO2H2<Scalar>, false>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::CarbonateIon<Scalar>, false>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::CH4<Scalar>, false>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::Granite<Scalar>, false>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::H2O<Scalar>>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::HeavyOil<Scalar>>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::Mesitylene<Scalar>>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::N2<Scalar>>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::NaCl<Scalar>, false>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::O2<Scalar>>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::SimpleH2O<Scalar>>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::Trichloroethene<Scalar>>::init(273, 275, 3, 1e5, 1e6, 3);
        Components::TabulatedComponent<Components::Xylene<Scalar>>::init(273, 275, 3, 1e5, 1e6, 3);
    }

    return 0;
}
