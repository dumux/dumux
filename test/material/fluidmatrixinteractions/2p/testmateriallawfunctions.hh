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
 * \brief Test for the 2p van Genuchten material law
 */

#ifndef DUMUX_TEST_MATERIALLAW_FUNCTIONS_HH
#define DUMUX_TEST_MATERIALLAW_FUNCTIONS_HH

#include <iomanip>

#include <dune/common/float_cmp.hh>
#include <dumux/common/math.hh>
#include <dumux/common/numericdifferentiation.hh>

namespace Dumux::Test {

template<class F, class D>
void testDerivatives(std::string_view derivName,
                     const std::vector<double>& values,
                     const F& f, const D& deriv)
{
    static constexpr double eps = 1.0e-1;
    static constexpr double numEps = 1e-8;
    for (auto val : values)
    {
        double analyticDeriv = deriv(val);

        double numericDeriv = 0.0;
        Dumux::NumericDifferentiation::partialDerivative(f, val,
            numericDeriv, f(val), numEps, 0 /*central differences*/);

        if (Dune::FloatCmp::ne(analyticDeriv, numericDeriv, eps))
             DUNE_THROW(Dune::Exception, "Analytic derivative for " << derivName
                         << " doesn't match numerical derivative: "
                         << std::setprecision(10) << analyticDeriv << " != " << numericDeriv
                         << " evaluated at " << val);
    }
}

template<class F, class G>
void testValueEqualRange(std::string_view testName,
                         const std::vector<double>& values,
                         const F& f, const G& g)
{
    static constexpr double eps = 1e-7;
    for (auto val : values)
    {
        const auto a = f(val);
        const auto b = g(val);

        if (Dune::FloatCmp::ne(a, b, eps))
            DUNE_THROW(Dune::Exception, "Test: " << testName << ": Function values do not match: "
                       << a << " != " << b << " evaluated at " << val << "\n");

    }
}


template<class Law, class RegLaw>
void runMaterialLawTest(const std::string& name,
                        const Law& law,
                        const RegLaw& regLaw,
                        const std::vector<typename Law::Scalar>& sw,
                        const std::vector<typename Law::Scalar>& swNonReg)
{
    const auto pc = [&](){ auto pc = sw;
        for (int i = 0; i < sw.size(); ++i)
            pc[i] = regLaw.pc(sw[i]);
        return pc;
    }();

    testDerivatives("dpc_dsw", sw, [&](auto sw){ return regLaw.pc(sw); }, [&](auto sw){ return regLaw.dpc_dsw(sw); });
    testDerivatives("dkrw_dsw", sw, [&](auto sw){ return regLaw.krw(sw); }, [&](auto sw){ return regLaw.dkrw_dsw(sw); });
    testDerivatives("dkrn_dsw", sw, [&](auto sw){ return regLaw.krn(sw); }, [&](auto sw){ return regLaw.dkrn_dsw(sw); });
    testDerivatives("dsw_dpc", pc, [&](auto pc){ return regLaw.sw(pc); }, [&](auto pc){ return regLaw.dsw_dpc(pc); });
    testValueEqualRange("Checking sw == sw(pc(sw))", sw, [](auto sw){ return sw; }, [&](auto sw) { return regLaw.sw(regLaw.pc(sw)); });
    testValueEqualRange("Checking 1.0 == dsw_dpc*dpc_dsw^-1", sw, [](auto sw){ return 1.0; }, [&](auto sw) { return regLaw.dpc_dsw(sw)*regLaw.dsw_dpc(regLaw.pc(sw)); });

    // check that regularized and unregularized are the same in the region without regularization
    testValueEqualRange("Checking NoReg::pc == Reg::pc", swNonReg, [&](auto sw){ return law.pc(sw); }, [&](auto sw) { return regLaw.pc(sw); });

    // test pc-sw curve against some precomputed values
    writeContainerToFile(pc, "test_pcsw_" + name + ".dat", 100);
}


template<class FluidMatrix>
void runEffToAbsTest(const std::string& name, const FluidMatrix& fmLaw,
                     const std::vector<typename FluidMatrix::Scalar>& sw)

{
    using EffToAbs = typename FluidMatrix::EffToAbs;
    const auto params = fmLaw.effToAbsParams();

    testValueEqualRange("Checking 1.0 == Abs::swToSwe(1-snr)", sw, [](auto sw){ return 1.0; }, [&](auto sw) { return EffToAbs::swToSwe(1-params.snr(), params); });
    testValueEqualRange("Checking 0.0 == Abs::swToSwe(snr)", sw, [](auto sw){ return 0.0; }, [&](auto sw) { return EffToAbs::swToSwe(params.snr(), params); });
    testValueEqualRange("Checking sw == sweToSw(swToSwe(sw))", sw, [](auto sw){ return sw; }, [&](auto sw) { return EffToAbs::sweToSw(EffToAbs::swToSwe(sw, params), params); });
}

} // end namespace Dumux

#endif
