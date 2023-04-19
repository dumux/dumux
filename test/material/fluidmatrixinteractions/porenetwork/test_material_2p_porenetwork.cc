//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <dune/common/float_cmp.hh>

#include <dumux/common/parameters.hh>
#include <dumux/io/container.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/localrulesforplatonicbody.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/multishapelocalrules.hh>

#include "testmateriallawfunctions.hh"

int main(int argc, char** argv) try
{
    using namespace Dumux;

    const double poreRadius = 1e-6;
    const double surfaceTension = 0.0725;
    const auto sw = Dumux::linspace(0.0, 1.0, 100);
    const auto swForDerivatives = Dumux::linspace(1e-8, 0.7, 100); // derivatives can't be numerically reproduced for sw > 0.7. TODO: find better numeric epsilon?
    const auto swNonReg = Dumux::linspace(1e-2, 1.0-1e-2, 100);

    GnuplotInterface<double> gnuplotPc;

    using LocalRulesCubeNoReg = PoreNetwork::FluidMatrix::TwoPLocalRulesPlatonicBodyNoReg<PoreNetwork::Pore::Shape::cube>;
    using LocalRulesCube = PoreNetwork::FluidMatrix::TwoPLocalRulesPlatonicBodyDefault<PoreNetwork::Pore::Shape::cube>;
    using Method = LocalRulesCube::RegularizationParams::HighSwRegularizationMethod;

    const auto params = LocalRulesCubeNoReg::BasicParams().setPoreShape(PoreNetwork::Pore::Shape::cube).setPoreInscribedRadius(poreRadius).setSurfaceTension(surfaceTension);
    LocalRulesCubeNoReg cubeLawNoReg(params);

    for (const auto method : std::array{Method::linear, Method::powerLaw, Method::spline})
    {
        LocalRulesCube::RegularizationParams regularizationParams;
        regularizationParams.setHighSwRegularizationMethod(method);

        const std::string name = [&]
        {
            if (method == Method::linear) return "linear";
            else if (method == Method::powerLaw) return "powerlaw";
            else return "spline";
        }();

        LocalRulesCube cubeLaw(params, regularizationParams);
        Dumux::Test::runMaterialLawTest("cube_with_" + name + "_regularization", cubeLaw, sw, swForDerivatives);

        // check that regularized and unregularized are the same in the region without regularization
        Dumux::Test::testValueEqualRange("Checking NoReg::pc == Reg::pc", swNonReg, [&](auto sw){ return cubeLawNoReg.pc(sw); }, [&](auto sw) { return cubeLaw.pc(sw); });

        // plot regularization part for high Sw
        if (argc > 1 && std::atoi(argv[1]) == 1)
        {
            const auto swPlot = Dumux::linspace(0.989, 1.01, 100);
            auto pc = swPlot;
            std::transform(swPlot.begin(), swPlot.end(), pc.begin(), [&](auto s){ return cubeLaw.pc(s); });

            // extend the x range by 10% on each side
            gnuplotPc.setXRange(0.98, 1.01);
            // gnuplotPc.setYRange(-1e-1, cubeLaw.pc(0.96)*1.01);
            gnuplotPc.setXlabel("wetting phase saturation [-]");
            gnuplotPc.setYlabel("capillary pressure [Pa]");
            gnuplotPc.addDataSetToPlot(swPlot, pc, name);
            gnuplotPc.plot("pc-Sw");
        }
    }

    std::cout << "\n\nChecking multi shape" << std::endl;
    using LocalRulesMultiShape = PoreNetwork::FluidMatrix::MultiShapeTwoPLocalRules<double>;
    using Shape = Dumux::PoreNetwork::Pore::Shape;

    for (const auto shape : std::array{Shape::tetrahedron, Shape::cube, Shape::octahedron, Shape::dodecahedron, Shape::icosahedron})
    {
        std::cout << "\nChecking " << Dumux::PoreNetwork::Pore::shapeToString(shape) << std::endl;
        LocalRulesMultiShape::BasicParams multiParams;
        const auto params = LocalRulesCubeNoReg::BasicParams().setPoreShape(shape).setPoreInscribedRadius(poreRadius).setSurfaceTension(surfaceTension);
        multiParams.setParams(params);
        LocalRulesMultiShape multiShapeLaw(multiParams);
        auto name = Dumux::PoreNetwork::Pore::shapeToString(shape);
        name[0] = std::tolower(name[0]);
        Dumux::Test::runMaterialLawTest(name, multiShapeLaw, sw, swForDerivatives);
        std::cout << "done\n" << std::endl;

    }

    return 0;
}
// error handler
catch (const Dune::Exception& e)
{
    std::cerr << "Test failed with exception: " << e << std::endl;
    return 1;
}
