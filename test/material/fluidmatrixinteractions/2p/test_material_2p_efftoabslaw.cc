#include <config.h>

#include <dune/common/float_cmp.hh>

#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchtenparams.hh>

#include <dumux/io/container.hh>
#include "testmateriallawfunctions.hh"


int main(int argc, char** argv) try
{
    using namespace Dumux;

    using VGEff = VanGenuchten<double>;
    using VGAbs = EffToAbsLaw<VanGenuchten<double>>;


    // set some parameters
    VGEff::Params effParams;
    VGAbs::Params absParams;
    absParams.setVgAlpha(6.66e-5);
    absParams.setVgn(3.652);
    absParams.setVgl(0.5);
    absParams.setSwr(0.1);
    absParams.setSnr(0.1);

    effParams.setVgAlpha(6.66e-5);
    effParams.setVgn(3.652);
    effParams.setVgl(0.5);

    const auto sw = Dumux::linspace(0.0, 1.0, 100);

    Test::runEffToAbsTest<VGEff, VGAbs>("efftoabs", effParams, absParams, sw);


    return 0;
}
// error handler
catch (const Dune::Exception& e)
{
    std::cerr << "Test failed with exception: " << e << std::endl;
    return 1;
}
