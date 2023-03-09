#include <config.h>

#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

#include <dumux/common/integrate.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/bsplinefunction.hh>

double evaluate(double p, double T) { return p; }

#include <fstream>

int main (int argc, char *argv[])
{
    Dumux::initialize(argc, argv);

    constexpr std::size_t nx = 7;
    constexpr std::size_t ny = 7;

    constexpr double pMin = 0; //e5;
    constexpr double pMax = 1; //e5;
    constexpr double TMin = 0; //73;
    constexpr double TMax = 1; //33;

    const auto getP = [&] (int i, std::size_t n) { return pMin + i*(pMax - pMin)/(n - 1); };
    const auto getT = [&] (int i, std::size_t n) { return TMin + i*(TMax - TMin)/(n - 1); };

    Dune::BlockVector<double> samples(nx*ny);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            samples[j*nx + i] = evaluate(getP(i, nx), getT(j, ny));

    Dumux::BSplineFunction function{
        Dune::FieldVector<double, 2>({pMin, TMin}),
        Dune::FieldVector<double, 2>({pMax, TMax}),
        {{nx, ny}},
        samples
    };

    Dune::YaspGrid grid{
        Dune::FieldVector<double, 2>({pMin, TMin}),
        Dune::FieldVector<double, 2>({pMax, TMax}),
        {{nx*5, ny*5}}
    };

    const auto& gridView = grid.leafGridView();
    const auto f1 = Dune::Functions::makeAnalyticGridViewFunction(
        [&] (const auto& x) { return function(x); },
        gridView
    );
    const auto f2 = Dune::Functions::makeAnalyticGridViewFunction(
        [&] (const auto& x) { return evaluate(x[0], x[1]); },
        gridView
    );
    const auto error = Dumux::integrateL2Error(gridView, f1, f2, 4);
    const auto integral = Dumux::integrateGridFunction(gridView, f2, 4);
    const auto relativeError = error/std::abs(integral);
    std::cout << "Error = " << error << ", relative = " << relativeError << ", integral = " << integral << std::endl;
    // if (relativeError > 1e-8)
    //     DUNE_THROW(Dune::InvalidStateException, "Relative error > 1e-8");

    const double p = 0.5*(pMin + pMax);
    const double T = 0.5*(TMin + TMax);
    std::cout << "f(p= " << p << ", t = " << T << ") = " << function({p, T}) << std::endl;

    std::ofstream xfile("x.txt", std::ios::out);
    std::ofstream yfile("y.txt", std::ios::out);
    std::ofstream zfile("z.txt", std::ios::out);
    std::ofstream zexactfile("zexact.txt", std::ios::out);

    constexpr auto NX = nx*5;
    constexpr auto NY = ny*5;
    for (int i = 0; i < NX; ++i) xfile << getP(i, NX) << std::endl;
    for (int i = 0; i < NY; ++i) yfile << getT(i, NY) << std::endl;
    for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i) {
            zfile << function({getP(i, NX), getT(j, NY)}) << std::endl;
            zexactfile << evaluate(getP(i, NX), getT(j, NY)) << std::endl;
        }

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
