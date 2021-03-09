#include <config.h>

#include <dune/common/float_cmp.hh>

#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

#include <dumux/io/container.hh>
#include "testmateriallawfunctions.hh"

namespace Dumux::Test {

// test if endPointPc() is the same as evaluation at sw=1
template<class Law>
void checkEndPointPc(const Law& law)
{
    const auto pcSat = law.pc(Law::EffToAbs::sweToSw(1.0, law.effToAbsParams()));
    const auto endPointPc = law.endPointPc();
    static constexpr double eps = 1e-7;

    if (Dune::FloatCmp::lt(eps, std::abs(pcSat-endPointPc)))
        DUNE_THROW(Dune::Exception, "pc(sw=1) != endPointPc(): " << pcSat << " != " << endPointPc);
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    using VGReg = FluidMatrix::VanGenuchtenDefault<double>;
    using VG = FluidMatrix::VanGenuchtenNoReg<double>;

    // set some parameters
    const double alpha = 6.66e-5;
    const double n = 3.652;
    const double l = 0.5;
    VGReg::BasicParams params(alpha, n, l);

    VGReg::EffToAbsParams eaParams;
    eaParams.setSwr(0.1);
    eaParams.setSnr(0.1);

    VGReg::RegularizationParams regParams;
    regParams.setPcLowSwe(0.01);
    regParams.setPcHighSwe(0.99);
    regParams.setKrnLowSwe(0.1);
    regParams.setKrwHighSwe(0.9);

    VGReg vgRegLaw(params, eaParams, regParams);
    VG vgLaw(params, eaParams);

    Test::checkEndPointPc(vgRegLaw);
    Test::checkEndPointPc(vgLaw);

    const auto sw = linspace(0.0, 1.0, 100);
    const auto swNonReg = linspace(VGReg::EffToAbs::sweToSw(regParams.pcLowSwe(), eaParams), VGReg::EffToAbs::sweToSw(regParams.pcHighSwe(), eaParams), 100);

    Test::runMaterialLawTest("vangenuchten", vgLaw, vgRegLaw, sw, swNonReg);
    Test::runEffToAbsTest("vangenuchten-efftoabs", vgLaw, sw);
    Test::runEffToAbsTest("vangenuchten-reg-efftoabs", vgRegLaw, sw);

    return 0;
}
