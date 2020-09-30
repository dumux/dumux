#include <config.h>

#include <iostream>
#include <iomanip>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dumux/nonlinear/findscalarroot.hh>

int main(int argc, char* argv[])
{
    using namespace Dumux;

    const auto func = [](const double x){ return x*x - 5.0; };
    const auto deriv = [](const double x){ return 2*x; };
    const auto absTol = 1e-15;
    const auto relTol = 1e-15;

    const auto rootNewton = findScalarRootNewton(5.0, func, deriv, absTol);
    const auto rootNewtonNumDiff = findScalarRootNewton(5.0, func, absTol);
    const auto rootBrent = findScalarRootBrent(0.0, 5.0, func, relTol);
    const auto rootExact = std::sqrt(5.0);

    if (Dune::FloatCmp::ne(rootExact, rootNewton, absTol))
        DUNE_THROW(Dune::Exception, "Wrong root for Newton algorithm (exact derivative)! Relative error: " << std::abs(rootExact-rootNewton)/rootExact);
    if (Dune::FloatCmp::ne(rootExact, rootNewtonNumDiff, absTol))
        DUNE_THROW(Dune::Exception, "Wrong root for Newton algorithm (numeric differentiation)! Relative error: " << std::abs(rootExact-rootNewtonNumDiff)/rootExact);
    if (Dune::FloatCmp::ne(rootExact, rootBrent, absTol))
        DUNE_THROW(Dune::Exception, "Wrong root for Brent algorithm! Relative error: " << std::abs(rootExact-rootBrent)/rootExact);

    std::cout << "Positive root of f(x) = (x*x - 5) is\n" << std::setprecision(15)
              << "-- Exact: " << rootExact << "\n"
              << "-- Newton (exact derivative): " << rootNewton << "\n"
              << "-- Newton (numeric derivative): " << rootNewtonNumDiff << "\n"
              << "-- Brent: " << rootBrent << std::endl;

    return 0;

}
