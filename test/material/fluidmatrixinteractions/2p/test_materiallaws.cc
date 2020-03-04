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
 * \brief Test for the material laws Brooks-Corey and van-Genuchten.
 */

 #include <config.h>
 #include <dumux/common/math.hh>

 #include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
 #include <dumux/material/fluidmatrixinteractions/2p/brookscoreyparams.hh>
 #include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
 #include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscoreyparams.hh>

 #include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
 #include <dumux/material/fluidmatrixinteractions/2p/vangenuchtenparams.hh>
 #include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
 #include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchtenparams.hh>

 #include <dumux/common/numericdifferentiation.hh>
 #include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

 #include <dumux/io/plotmateriallaw.hh>
 #include <dumux/io/gnuplotinterface.hh>
 #include <dumux/io/container.hh>

template <class Scalar, class MaterialLaw, class MaterialLawRegularized, class MaterialLawParams, class MaterialLawParamsRegularized>
void testMaterialLawCommon(MaterialLaw &materialLaw,
                          MaterialLawRegularized &materialLawRegularized,
                          MaterialLawParams &materialLawParams,
                          MaterialLawParamsRegularized &materialLawParamsRegularized,
                          Scalar sweValue,
                          Scalar pcValue,
                          Scalar eps)
{
    // test if derivatives match numerical derivatives for dpc_dsw
    Scalar analyticDerivativePcSw = materialLaw.dpc_dswe(materialLawParams, sweValue);
    Scalar numericDerivativePcSw = 0.0;

    auto materialLawPcNumericPcSw = [&](Scalar sw)
    {
        return materialLaw.pc(materialLawParams, sw);
    };

    Dumux::NumericDifferentiation::partialDerivative(materialLawPcNumericPcSw, sweValue,
                                                     numericDerivativePcSw,
                                                     materialLaw.pc(materialLawParams, sweValue),
                                                     eps, 0 /*central*/);

    if (Dune::FloatCmp::eq(analyticDerivativePcSw, numericDerivativePcSw, eps*analyticDerivativePcSw))
        DUNE_THROW(Dune::Exception, "analytic derivative dpc_dswe doesn't match numerical derivative dpc_dswe: "
        << analyticDerivativePcSw << " != " << numericDerivativePcSw << "\n");

    // test if derivatives match numerical derivatives for dsw_dpc
    Scalar analyticDerivativeSwPc = materialLaw.dswe_dpc(materialLawParams, pcValue);
    Scalar numericDerivativeSwPc = 0.0;

    auto materialLawPcNumericSwPc = [&](Scalar pc)
    {
       return materialLaw.sw(materialLawParams, pc);
    };

    Dumux::NumericDifferentiation::partialDerivative(materialLawPcNumericSwPc, pcValue,
                                                     numericDerivativeSwPc,
                                                     materialLaw.sw(materialLawParams, pcValue),
                                                     eps, 0 /*central*/);

    if (Dune::FloatCmp::eq(analyticDerivativeSwPc, numericDerivativeSwPc, eps*analyticDerivativeSwPc))
        DUNE_THROW(Dune::Exception, "analytic derivative dswe_dpc doesn't match numerical derivative dsw_dpc: "
        << analyticDerivativeSwPc << " != " << numericDerivativeSwPc << "\n");

    // test if the unregularized pc is the same as regularized pc in the non-regularized range
    Scalar pC = materialLaw.pc(materialLawParams, sweValue);
    Scalar pCReg = materialLawRegularized.pc(materialLawParamsRegularized, sweValue);

    if (Dune::FloatCmp::ne(pC, pCReg, eps*pC))
        DUNE_THROW(Dune::Exception, "regulized pc doesn't match unregularized pc: "
        << pCReg << " != " << pC << "\n");

    // test if endPointPc() is the same as evaluation at Sw=1 and entryPressure
    // for the un-regularized material law
    Scalar pCSat = materialLaw.pc(materialLawParams, 1.0);
    Scalar endPointPc = materialLaw.endPointPc(materialLawParams);

    if (Dune::FloatCmp::ne(pCSat, endPointPc, eps*pCSat))
        DUNE_THROW(Dune::Exception, "pc(Sw=1) doesn't match endPointPc: " << pCSat << " != "
        << endPointPc << "\n");

    // for the regularized material law
    Scalar pCSatReg = materialLawRegularized.pc(materialLawParamsRegularized, 1.0);
    Scalar endPointPcReg = materialLawRegularized.endPointPc(materialLawParamsRegularized);

    if (Dune::FloatCmp::ne(pCSatReg, endPointPcReg, eps*pCSatReg))
        DUNE_THROW(Dune::Exception, "regularized pc(Sw=1) doesn't match regularized endPointPc: "
        << pCSatReg << " != " << endPointPcReg << "\n");
}

 int main(int argc, char** argv)
 {
    using namespace Dumux;

    using Scalar = double;

    using BrooksCorey = BrooksCorey<Scalar>;
    using BrooksCoreyParams = typename BrooksCorey::Params;
    using BrooksCoreyRegularized = RegularizedBrooksCorey<Scalar>;
    using BrooksCoreyParamsRegularized = typename BrooksCoreyRegularized::Params;

    BrooksCorey brooksCorey;
    BrooksCoreyParams brooksCoreyParams;
    BrooksCoreyRegularized brooksCoreyRegularized;
    BrooksCoreyParamsRegularized brooksCoreyParamsRegularized;

    using VanGenuchten = VanGenuchten<Scalar>;
    using VanGenuchtenParams = typename VanGenuchten::Params;
    using VanGenuchtenRegularized = RegularizedVanGenuchten<Scalar>;
    using VanGenuchtenParamsRegularized = typename VanGenuchtenRegularized::Params;

    // set Params Brooks-Corey
    brooksCoreyParams.setPe(1e4);
    brooksCoreyParams.setLambda(2.0);
    brooksCoreyParamsRegularized.setPe(1e4);
    brooksCoreyParamsRegularized.setLambda(2.0);

    VanGenuchten vanGenuchten;
    VanGenuchtenParams vanGenuchtenParams;
    VanGenuchtenRegularized vanGenuchtenRegularized;
    VanGenuchtenParamsRegularized vanGenuchtenParamsRegularized;

    // set Params van-Genuchten
    vanGenuchtenParams.setVgAlpha(6.66e-5);
    vanGenuchtenParams.setVgn(3.652);
    vanGenuchtenParamsRegularized.setVgAlpha(6.66e-5);
    vanGenuchtenParamsRegularized.setVgn(3.652);

    // test common functions of both material laws
    // define ranges for sw and pc for testing the material laws
    const int n = 10; //size of the following vectors
    const auto sweValues = Dumux::linspace(0.005, 1.0, n); // range of tested sw values
    const auto pcValues = Dumux::linspace(1.0, 10000.0, n); // range of tested pc values
    Scalar eps = 1.0e-6; // threshhold value

    for( int i = 0; i < n; i++)
    {
    // test Brooks-Corey
        testMaterialLawCommon<Scalar, BrooksCorey, BrooksCoreyRegularized, BrooksCoreyParams, BrooksCoreyParamsRegularized>
        (brooksCorey, brooksCoreyRegularized, brooksCoreyParams, brooksCoreyParamsRegularized, sweValues[i], pcValues[i], eps);
    // test van-Genuchten
        testMaterialLawCommon<Scalar, VanGenuchten, VanGenuchtenRegularized, VanGenuchtenParams, VanGenuchtenParamsRegularized>
        (vanGenuchten, vanGenuchtenRegularized, vanGenuchtenParams, vanGenuchtenParamsRegularized, sweValues[i], pcValues[i], eps);
    }

    // test Brooks-Corey specific functions
    // test if endPointPc() is the same as entryPressure
    Scalar pCSat = brooksCorey.pc(brooksCoreyParams, 1.0);
    Scalar pCSatReg = brooksCoreyRegularized.pc(brooksCoreyParamsRegularized, 1.0);
    if (Dune::FloatCmp::ne(pCSat, brooksCoreyParams.pe(), eps*pCSat))
        DUNE_THROW(Dune::Exception, "pc(Sw=1) doesn't match entryPressure: " << pCSat << " != " << brooksCoreyParams.pe() << "\n");
    if (Dune::FloatCmp::ne(pCSatReg, brooksCoreyParamsRegularized.pe(), eps*pCSatReg))
        DUNE_THROW(Dune::Exception, "regularized pc(Sw=1) doesn't match regularized entryPressure: " << pCSatReg << " != " << brooksCoreyParamsRegularized.pe() << "\n");

// 4) test against some precomputed reference values (a good regression test)
    GnuplotInterface<double> gnuplot;
    gnuplot.setOpenPlotWindow(false);
    using PlotMaterialLaw = PlotMaterialLaw<Scalar, BrooksCorey>;
    PlotMaterialLaw plotMaterialLaw;
    const std::string fileName = "pcswcurve.dat";

    plotMaterialLaw.addpcswcurve(gnuplot, brooksCoreyParams, 0.0 + eps, 1.0, fileName);

// 5) //TODO test eff to abs law
//     using MaterialLaw = EffToAbsLaw<BrooksCorey>;
//     using MaterialLawRegularized = EffToAbsLaw<EffectiveLawRegularized>;


    return 0;
 }