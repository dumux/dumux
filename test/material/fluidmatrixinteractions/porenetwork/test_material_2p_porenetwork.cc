#include <config.h>

#include <dune/common/float_cmp.hh>

#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchtenparams.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchtenparams.hh>

#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/regularizedlocalrules.hh>



#include <dumux/io/container.hh>
#include "testmateriallawfunctions.hh"

namespace Dumux::Test {

// test if endPointPc() is the same as evaluation at sw=1
template<class Law>
void checkEndPointPc(const typename Law::Params& params)
{
    const auto pcSat = Law::pc(params, Law::sweToSw(params, 1.0));
    const auto endPointPc = Law::endPointPc(params);
    static constexpr double eps = 1e-10;

    if (Dune::FloatCmp::ne(pcSat, endPointPc, eps))
        DUNE_THROW(Dune::Exception, "pc(sw=1) != endPointPc(): " << pcSat << " != " << endPointPc);
}

} // end namespace Dumux

int main(int argc, char** argv) try
{
    using namespace Dumux;


    const double poreRadius = 1e-5;
    const double contactAngle = 0.0;
    const double surfaceTension = 0.0725;
    const auto shape = Pore::Shape::cube; // todo more shapes
    const auto sw = Dumux::linspace(0.0, 1.0, 100);


    using LocalRules = TwoPLocalRules<double>;
    using RegLocalRules = RegularizedTwoPLocalRules<double>;
    const auto params = RegLocalRules::makeParams(poreRadius, contactAngle, surfaceTension, shape);
    const auto swNonReg = Dumux::linspace(params.lowSw, params.highSw, 100);

    // set some parameters


    // Test::checkEndPointPc<VG>(params);
    // Test::checkEndPointPc<VGReg>(params);


    Dumux::Test::testValueEqualRange("Checking sw == sw(pc(sw))", sw, [](auto sw){ return sw; }, [&](auto sw) { return LocalRules::sw(params, LocalRules::pc(params, sw)); });
    // Dumux::Test::testValueEqualRange("Checking sw == sw(pc(sw))", sw, [](auto sw){ return sw; }, [&](auto sw) { return RegLocalRules::sw(params, RegLocalRules::pc(params, sw)); });
    // check that regularized and unregularized are the same in the region without regularization
    Dumux::Test::testValueEqualRange("Checking NoReg::pc == Reg::pc", swNonReg, [&](auto sw){ return RegLocalRules::pc(params, sw); }, [&](auto sw) { return LocalRules::pc(params, sw); });

    // Test::runMaterialLawTest<VG, RegLocalRules>("vangenuchten", params, sw, swNonReg);

    return 0;
}
// error handler
catch (const Dune::Exception& e)
{
    std::cerr << "Test failed with exception: " << e << std::endl;
    return 1;
}
