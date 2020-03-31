#include <config.h>

#include <dune/common/float_cmp.hh>

#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchtenparams.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchtenparams.hh>

#include <dumux/io/container.hh>
#include "testmateriallawfunctions.hh"

namespace Dumux::Test {

// test if endPointPc() is the same as evaluation at sw=1
template<class Law>
void checkEndPointPc(const typename Law::Params& params)
{
    const auto pcSat = Law::pc(params, Law::sweToSw_(params, 1.0));
    const auto endPointPc = Law::endPointPc(params);
    static constexpr double eps = 1e-10;

    if (Dune::FloatCmp::ne(pcSat, endPointPc, eps))
        DUNE_THROW(Dune::Exception, "pc(sw=1) != endPointPc(): " << pcSat << " != " << endPointPc);
}

} // end namespace Dumux

int main(int argc, char** argv) try
{
    using namespace Dumux;

    using VGRegEff = RegularizedVanGenuchten<double>;
    using VGEff = VanGenuchten<double>;
    using VGReg = EffToAbsLaw<VGRegEff>;
    using VG = EffToAbsLaw<VGEff, VGReg::Params>;

    // set some parameters
    VGReg::Params params;
    params.setVgAlpha(6.66e-5);
    params.setVgn(3.652);
    params.setVgl(0.5);
    params.setSwr(0.1);
    params.setSnr(0.1);
    params.setPcLowSw(0.01);
    params.setPcHighSw(0.99);
    params.setKrnLowSw(0.1);
    params.setKrwHighSw(0.9);

    Test::checkEndPointPc<VG>(params);
    Test::checkEndPointPc<VGReg>(params);

    const auto sw = Dumux::linspace(0.0, 1.0, 100);
    const auto swNonReg = Dumux::linspace(VGReg::sweToSw_(params, params.pcLowSw()), VGReg::sweToSw_(params, params.pcHighSw()), 100);

    Test::runMaterialLawTest<VG, VGReg>("vangenuchten", params, sw, swNonReg);
    Test::runEffToAbsTest<VGRegEff, VGReg>("vangenuchten-efftoabs", params, sw);

    return 0;
}
// error handler
catch (const Dune::Exception& e)
{
    std::cerr << "Test failed with exception: " << e << std::endl;
    return 1;
}
