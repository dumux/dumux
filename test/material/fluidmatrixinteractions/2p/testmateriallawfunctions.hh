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

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/numericdifferentiation.hh>

namespace Dumux::Test {

template<class F, class D>
void testDerivatives(std::string_view derivName,
                     const std::vector<double>& values,
                     const F& f, const D& deriv)
{
    static constexpr double eps = 1e-7;
    static constexpr double numEps = 1e-8;
    for (auto val : values)
    {
        const double analyticDeriv = deriv(val);

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
void runTest(const std::string& name, const typename RegLaw::Params& params,
             const std::vector<typename Law::Scalar>& sw,
             const std::vector<typename Law::Scalar>& swNonReg)
{
    const auto pc = [&](){ auto pc = sw;
        for (int i = 0; i < sw.size(); ++i)
            pc[i] = RegLaw::pc(params, sw[i]);
        return pc;
    }();

    testDerivatives("dpc_dsw", sw, [&](auto sw){ return RegLaw::pc(params, sw); }, [&](auto sw){ return RegLaw::dpc_dsw(params, sw); });
    //testDerivatives("dkrw_dsw", sw, [&](auto sw){ return RegLaw::krw(params, sw); }, [&](auto sw){ return RegLaw::dkrw_dsw(params, sw); });
    //testDerivatives([&](auto sw){ return RegLaw::krn(params, sw); }, [&](auto sw){ return RegLaw::dkrn_dsw(params, sw); }, sw, "dkrn_dsw");
    //testDerivatives("dsw_dpc", [&](auto sw){ return RegLaw::sw(params, sw); }, [&](auto sw){ return RegLaw::dsw_dpc(params, sw); }, pc);
    testValueEqualRange("Checking sw == sw(pc(sw))", sw, [](auto sw){ return sw; }, [&](auto sw) { return RegLaw::sw(params, RegLaw::pc(params, sw)); });
    testValueEqualRange("Checking 1.0 == dsw_dpc*dpc_dsw^-1", sw, [](auto sw){ return 1.0; }, [&](auto sw) { return RegLaw::dpc_dsw(params, sw)*RegLaw::dsw_dpc(params, RegLaw::pc(params, sw)); });

    // check that regularized and unregularized are the same in the region without regularization
    testValueEqualRange("Checking NoReg::pc == Reg::pc", swNonReg, [&](auto sw){ return RegLaw::pc(params, sw); }, [&](auto sw) { return Law::pc(params, sw); });

    // test pc-sw curve against some precomputed values
    writeContainerToFile(pc, "test_pcsw_" + name + ".dat", 100);
}


template<class EffLaw, class AbsLaw>
void runEffToAbsTest(const std::string& name, const typename EffLaw::Params& effParams, const typename AbsLaw::Params& absParams,
             const std::vector<typename EffLaw::Scalar>& sw)

{
    testValueEqualRange("Checking 1.0 == Abs::swToSwe(1-snr)", sw, [](auto sw){ return 1.0; }, [&](auto sw) { return AbsLaw::swToSwe(absParams, 1-absParams.snr()); });
    testValueEqualRange("Checking 0.0 == Abs::swToSwe(snr)", sw, [](auto sw){ return 0.0; }, [&](auto sw) { return AbsLaw::swToSwe(absParams, absParams.snr()); });
    testValueEqualRange("Checking 1.0 == Abs::snToSne(1-swr)", sw, [](auto sw){ return 1.0; }, [&](auto sw) { return AbsLaw::snToSne(absParams, 1-absParams.swr()); });
    testValueEqualRange("Checking 0.0 == Abs::snToSne(swr)", sw, [](auto sw){ return 0.0; }, [&](auto sw) { return AbsLaw::snToSne(absParams, absParams.swr()); });

    testValueEqualRange("Checking Abs::pc(1-snr) == Eff::pc(1.0)", sw, [&](auto pc){ return AbsLaw::pc(absParams, 1-absParams.snr()); }, [&](auto pc) { return EffLaw::pc(effParams, 1.0); });
    testValueEqualRange("Checking Abs::endPointPc == Eff::pc(1.0)", sw, [&](auto pc){ return AbsLaw::endPointPc(absParams); }, [&](auto pc) { return EffLaw::pc(absParams, 1.0); });
}

} // end namespace Dumux

#endif
