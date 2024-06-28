#include <config.h>

#include <dune/common/float_cmp.hh>

#include <dumux/common/parameters.hh>

#include <dumux/material/fluidmatrixinteractions/2p/tabulated.hh>

#include <dumux/io/container.hh>
#include "testmateriallawfunctions.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using Scalar = double;
    Parameters::init(argc, argv);

//    using VGReg = FluidMatrix::VanGenuchtenDefault<double>;
//    using VG = FluidMatrix::VanGenuchtenNoReg<double>;
//
//    // set some parameters
//    const double alpha = 6.66e-5;
//    const double n = 3.652;
//    const double l = 0.5;
//    VGReg::BasicParams params(alpha, n, l);
//
//    VGReg::EffToAbsParams eaParams;
//    eaParams.setSwr(0.1);
//    eaParams.setSnr(0.1);
//
//    VGReg::RegularizationParams regParams;
//    regParams.setPcLowSwe(0.01);
//    regParams.setPcHighSwe(0.99);
//    regParams.setKrnLowSwe(0.1);
//    regParams.setKrwHighSwe(0.9);
//
//    VGReg vgRegLaw(params, eaParams, regParams);
//    VG vgLaw(params, eaParams);
//
//    Test::checkEndPointPc(vgRegLaw);
//    Test::checkEndPointPc(vgLaw);
//
//    const auto sw = linspace(0.0, 1.0, 100);
//    const auto swNonReg = linspace(VGReg::EffToAbs::sweToSw(regParams.pcLowSwe(), eaParams), VGReg::EffToAbs::sweToSw(regParams.pcHighSwe(), eaParams), 100);
//
//    Test::runMaterialLawTest("vangenuchten", vgLaw, vgRegLaw, sw, swNonReg);
//    Test::runEffToAbsTest("vangenuchten-efftoabs", vgLaw, sw);
//    Test::runEffToAbsTest("vangenuchten-reg-efftoabs", vgRegLaw, sw);

    using Tab = FluidMatrix::TabulatedPropertiesDefault<61, Scalar>;
    // set some parameters
    const double alpha = 6.66e-5;
    const double n = 3.652;
    const double l = 0.5;
    Tab::BasicParams params("X-dir",alpha, n, l);
    Tab::EffToAbsParams eaParams;
    eaParams.setSwr(0.0);
    eaParams.setSnr(0.0);
    Tab table(params, eaParams);
    int testSteps = 200;
    for (int i = 0; i <= testSteps; ++i)
    {
        Scalar sw = Scalar(i)/testSteps;
        Scalar pc = table.pc(sw);
        Scalar sw_pc = table.sw(pc);
        Scalar krw = table.krw(sw);
        Scalar krn = table.krn(sw);
        std::cout << sw << ": " << pc << ", " << sw_pc
            << ", " << krw << ", " << krn << std::endl;
    }
    std::cout << "endPointPc = " << table.endPointPc() << std::endl;
    Tab tableParams("");
    for (int i = 0; i <= testSteps; ++i)
    {
        Scalar sw = Scalar(i)/testSteps;
        Scalar pc = tableParams.pc(sw);
        Scalar sw_pc = tableParams.sw(pc);
        Scalar krw = tableParams.krw(sw);
        Scalar krn = tableParams.krn(sw);
        std::cout << sw << ": " << pc << ", " << sw_pc
            << ", " << krw << ", " << krn << std::endl;
    }
    std::cout << "endPointPc = " << tableParams.endPointPc() << std::endl;

    return 0;
}
