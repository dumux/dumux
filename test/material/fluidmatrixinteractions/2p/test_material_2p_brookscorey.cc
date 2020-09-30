#include <config.h>

#include <dune/common/float_cmp.hh>

#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>

#include <dumux/io/container.hh>
#include "testmateriallawfunctions.hh"

namespace Dumux::Test {

// test if endPointPc() is the same as evaluation at sw=1
template<class Law>
void checkEndPointPc(const Law& law, const double entryPressure)
{
    const auto pcSat = law.pc(Law::EffToAbs::sweToSw(1.0, law.effToAbsParams()));
    const auto endPointPc = law.endPointPc();
    static constexpr double eps = 1e-10;

    if (Dune::FloatCmp::ne(pcSat, endPointPc, eps))
        DUNE_THROW(Dune::Exception, "pc(sw=1) != endPointPc(): " << pcSat << " != " << endPointPc);
    if (Dune::FloatCmp::ne(pcSat, entryPressure, eps))
        DUNE_THROW(Dune::Exception, "pc(sw=1) != entryPressure: " << pcSat << " != " << entryPressure);
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    using BCReg = FluidMatrix::BrooksCoreyDefault<double>;
    using BC = FluidMatrix::BrooksCoreyNoReg<double>;

    // set some parameters
    BCReg::BasicParams params;
    const double entryPressure = 1e4;
    params.pe = entryPressure;
    params.lambda = 2.0;

    BCReg::EffToAbsParams eaParams;
    eaParams.swr = 0.1;
    eaParams.snr = 0.1;

    BCReg::RegularizationParams regParams;
    const double thresholdSw = 0.01;
    regParams.pcLowSwe = thresholdSw;

    BCReg bcRegLaw(params, eaParams, regParams);
    BC bcLaw(params, eaParams);

    Test::checkEndPointPc(bcRegLaw, entryPressure);
    Test::checkEndPointPc(bcLaw, entryPressure);

    const auto sw = linspace(0.0, 1.0, 100);
    const auto swNonReg = linspace(BCReg::EffToAbs::sweToSw(thresholdSw, eaParams), BCReg::EffToAbs::sweToSw(1.0, eaParams), 100);

    Test::runMaterialLawTest("brookscorey", bcLaw, bcRegLaw, sw, swNonReg);
    Test::runEffToAbsTest("brookscorey-efftoabs", bcLaw, sw);
    Test::runEffToAbsTest("brookscorey-reg-efftoabs", bcRegLaw, sw);

    return 0;
}
