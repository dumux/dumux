#include <config.h>

#include <dune/common/float_cmp.hh>

#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscoreyparams.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscoreyparams.hh>

#include <dumux/io/container.hh>
#include "testmateriallawfunctions.hh"

namespace Dumux::Test {

// test if endPointPc() is the same as evaluation at sw=1
template<class Law>
void checkEndPointPc(const typename Law::Params& params)
{
    const auto pcSat = Law::pc(params, Law::sweToSw_(params, 1.0));
    const auto endPointPc = Law::endPointPc(params);
    const auto entryPressure = params.pe();
    static constexpr double eps = 1e-10;

    if (Dune::FloatCmp::ne(pcSat, endPointPc, eps))
        DUNE_THROW(Dune::Exception, "pc(sw=1) != endPointPc(): " << pcSat << " != " << endPointPc);
    if (Dune::FloatCmp::ne(pcSat, entryPressure, eps))
        DUNE_THROW(Dune::Exception, "pc(sw=1) != entryPressure: " << pcSat << " != " << entryPressure);
}

} // end namespace Dumux

int main(int argc, char** argv) try
{
    using namespace Dumux;

    using BCReg = EffToAbsLaw<RegularizedBrooksCorey<double>>;
    using BC = EffToAbsLaw<BrooksCorey<double>, BCReg::Params>;

    // set some parameters
    BCReg::Params params;
    params.setPe(1e4);
    params.setLambda(2.0);
    params.setSwr(0.1);
    params.setSnr(0.1);
    params.setThresholdSw(0.01);

    Test::checkEndPointPc<BC>(params);
    Test::checkEndPointPc<BCReg>(params);

    const auto sw = Dumux::linspace(0.0, 1.0, 100);
    const auto swNonReg = Dumux::linspace(BCReg::sweToSw_(params, params.thresholdSw()), BCReg::sweToSw_(params, 1.0), 100);

    Test::runTest<BC, BCReg>("brookscorey", params, sw, swNonReg);

    return 0;
}
// error handler
catch (const Dune::Exception& e)
{
    std::cerr << "Test failed with exception: " << e << std::endl;
    return 1;
}
