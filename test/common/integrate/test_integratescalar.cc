#include <config.h>

#include <iostream>
#include <iomanip>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/math.hh>
#include <dumux/common/integrate.hh>
#include <dumux/nonlinear/findscalarroot.hh>

namespace Dumux {

class VanGenuchtenKrw
{
public:
    VanGenuchtenKrw(double alpha, double n, double l = 0.5)
    : alpha_(alpha), n_(n), m_(1.0 - 1.0/n), l_(l)
    {}

    double operator()(double pw) const
    {
        using std::pow;
        const auto pc = 1.0e5-pw;
        const auto swe = pow(pow(alpha_*pc, n_) + 1.0, -m_);
        const auto r = 1.0 - pow(1.0 - pow(swe, 1.0/m_), m_);
        return pow(swe, l_)*r*r;
    }

private:
    double alpha_, n_, m_, l_;
};

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;

    {
        // test some integral where we know the analytical expression
        const auto func = [](const double x){ return std::exp(2.0*x); };
        const auto exactIntegral = 0.5*(std::exp(10.0)-std::exp(-2.0));
        const auto numIntegral = integrateScalarFunction(func, -1.0, 5.0, 1e-8);
        if (Dune::FloatCmp::ne<double, Dune::FloatCmp::CmpStyle::absolute>(exactIntegral, numIntegral, 1e-8))
            DUNE_THROW(Dune::Exception, "Wrong integral: " << numIntegral << ", should be " << exactIntegral);
    }

    {
        // test some integral where we know the analytical expression
        const auto func = [](const double x){ return Dune::power(1.0 - x, 5)*pow(x, -1.0/3.0); };
        const auto exactIntegral = 2187.0/5236.0;
        const auto numIntegral = integrateScalarFunction(func, 0.0, 1.0, 1e-8);
        if (Dune::FloatCmp::ne<double, Dune::FloatCmp::CmpStyle::absolute>(exactIntegral, numIntegral, 1e-8))
            DUNE_THROW(Dune::Exception, "Wrong integral: " << numIntegral << ", should be " << exactIntegral);
    }

    {
        // relative permeability function of water pressure
        const auto krFunc = VanGenuchtenKrw(/*alpha=*/1e-4, /*n=*/3.0);

        // introduce a transformed variable u(p) = int_{-1.5e6}^{p} kr(q)dq
        const auto kirchhoffTransform = [&](const double p){
            return integrateScalarFunction(krFunc, -1.5e6, p, 1e-16);
        };

        // compute the inverse of the transformation numerically
        const auto inverseTransform = [&](const auto u){
            const auto residual = [&, u](const auto& p){ return u - kirchhoffTransform(p); };
            return findScalarRootBrent(-1.6e6, 1.0e5, residual, 1e-10, 200000);
        };

        // compute transformation and inverse for several pressures
        for (const auto p : {0.9e5, 0.1e5, 0.0, -5.0e5, -1.5e6})
        {
            const auto u = kirchhoffTransform(p);
            const auto pNum = inverseTransform(u);
            if (Dune::FloatCmp::ne<double, Dune::FloatCmp::CmpStyle::absolute>(p, pNum, 0.001))
                DUNE_THROW(Dune::Exception, "Inverse did not recover original pressure value: "
                           << std::setprecision(15) << pNum << ", should be " << p);
        }
    }

    return 0;
}
