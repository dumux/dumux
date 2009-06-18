#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>

#include "dumux/material/matrixproperties.hh"
#include "dumux/diffusion/fv/fvwettingvelocity2p.hh"
#include "dumux/diffusion/mpfa/mpfaodiffusionvelocity.hh"
#include "dumux/diffusion/mimetic/mimeticpressure2p.hh"
#include "dumux/diffusion/fe/fepressure2p.hh"
#include "dumux/fractionalflow/variableclass2p.hh"

#include "fvca5test5problem.hh"
#include "benchmarkresult.hh"

int main(int argc, char** argv)
{
    try{
        // define the problem dimensions
        const int dim=2;

        typedef double Scalar;
#if HAVE_UG
        typedef Dune::UGGrid<dim> Grid;
#else
        typedef Dune::SGrid<dim,dim> Grid;
#endif

#if HAVE_PARDISO
        std::string solver = "Loop";
        std::string preconditioner = "SeqPardiso";
#else
        std::string solver = "CG";
        std::string preconditioner = "SeqILU0";
#endif

        if (argc < 2 || argc > 3) {
            std::cout << "\nUsage: test_diffusion refinementsteps [delta]\n" << std::endl;
            std::cout << "- refinementsteps: number of uniform refinements of the initial grid." << std::endl;
            std::cout << "   The initial grid is grids/fvca5test5.dgf.\n" << std::endl;
            std::cout << "- delta: parameter in (0, 1] for the permeability tensor K." << std::endl;
            std::cout << "   Default is 1e-3. K becomes singular for delta = 0, and unity for delta = 1.\n" << std::endl;
            return (1);
        }
        std::stringstream dgfFileName;
        dgfFileName << "./grids/fvca5test5.dgf";
        Dune::GridPtr<Grid> gridPtr( dgfFileName.str() );
        Grid& grid = *gridPtr;

        int refinementSteps;
        std::string arg2(argv[1]);
        std::istringstream is2(arg2);
        is2 >> refinementSteps;

        double delta = 1e-3;
        if (argc > 2) {
            std::string arg3(argv[2]);
            std::istringstream is3(arg3);
            is3 >> delta;
        }

        if (refinementSteps)
            grid.globalRefine(refinementSteps);

        typedef Grid::LevelGridView GridView;
        GridView gridView(grid.levelView(grid.maxLevel()));

        Dune::FVariableClassA5Test5Soil<Grid, Scalar> soil(delta);

        typedef Dune::VariableClass<GridView, Scalar> VariableClass;

        Dune::Timer timer;

        VariableClass fvVariables(gridView);
        Dune::FVariableClassA5Test5Problem<GridView, Scalar, VariableClass> fvProblem(fvVariables, soil, delta);
        timer.reset();
        Dune::FVWettingPhaseVelocity2P<GridView, Scalar, VariableClass> fvDiffusion(gridView, fvProblem, "pw","Sw", solver, preconditioner);
        fvDiffusion.pressure(true, 0, false);
        fvDiffusion.calculateVelocity();
        double fvTime = timer.elapsed();
        fvVariables.vtkout("fv-2pfa", 0);
        Dune::P0Function<GridView, Scalar, 1> fvPressureP0(gridView);
        *fvPressureP0 = fvVariables.pressure();
        Dune::ResultEvaluation fvResult;
        fvResult.evaluate(gridView, fvProblem, fvPressureP0, fvVariables.velocity());

#if HAVE_UG
        VariableClass mpfaVariables(gridView);
        Dune::FVariableClassA5Test5Problem<GridView, Scalar, VariableClass> mpfaProblem(mpfaVariables, soil, delta);
        timer.reset();
        Dune::MPFAODiffusionVelocity<GridView, Scalar, VariableClass> mpfaDiffusion(gridView, mpfaProblem, solver, preconditioner);
        mpfaDiffusion.pressure();
        mpfaDiffusion.calculateVelocity();
        double mpfaTime = timer.elapsed();
        mpfaVariables.vtkout("fv-mpfa", 0);
        Dune::P0Function<GridView, Scalar, 1> mpfaPressureP0(gridView);
        *mpfaPressureP0 = mpfaVariables.pressure();
        Dune::ResultEvaluation mpfaResult;
        mpfaResult.evaluate(gridView, mpfaProblem, mpfaPressureP0, mpfaVariables.velocity());
#endif

        VariableClass mimeticVariables(gridView);
        Dune::FVariableClassA5Test5Problem<GridView, Scalar, VariableClass> mimeticProblem(mimeticVariables, soil, delta);
        timer.reset();
        Dune::MimeticPressure2P<GridView, Scalar, VariableClass> mimeticDiffusion(gridView, mimeticProblem, solver, preconditioner);
        mimeticDiffusion.pressure();
        mimeticDiffusion.calculateVelocity();
        double mimeticTime = timer.elapsed();
        mimeticVariables.vtkout("mimetic", 0);
        Dune::BenchmarkResult mimeticResult;
        mimeticResult.evaluate(grid, mimeticProblem, mimeticDiffusion);

        VariableClass feVariables(gridView, dim);
        Dune::FVariableClassA5Test5Problem<GridView, Scalar, VariableClass> feProblem(feVariables, soil, delta);
        timer.reset();
        Dune::FEPressure2P<GridView, Scalar, VariableClass> feDiffusion(gridView, feProblem, solver, preconditioner);
        feDiffusion.pressure();
        feDiffusion.calculateVelocity();
        double feTime = timer.elapsed();
        feVariables.vtkout("fe", 0);
        Dune::ResultEvaluation feResult;
        feResult.evaluate(gridView, feProblem, feDiffusion.pressP1, feVariables.velocity());

        std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
        std::cout.precision(2);
        std::cout << "\t error press \t error grad\t sumflux\t erflm\t\t uMin\t\t uMax\t\t time" << std::endl;
        std::cout << "fv-2pfa\t " << fvResult.relativeL2Error << "\t " << fvResult.ergrad << "\t " << fvResult.sumflux
			<< "\t " << fvResult.erflm << "\t " << fvResult.uMin << "\t " << fvResult.uMax << "\t " << fvTime << std::endl;
#if HAVE_UG
        std::cout << "fv-mpfa\t " << mpfaResult.relativeL2Error	<< "\t " << mpfaResult.ergrad << "\t " << mpfaResult.sumflux
			<< "\t " << mpfaResult.erflm << "\t " << mpfaResult.uMin << "\t " << mpfaResult.uMax << "\t " << mpfaTime << std::endl;
#else
        std::cout << "fv-mpfa\t please install UGGrid, MPFAODiffusion currently does not work with SGrid" << std::endl;
#endif
        std::cout << "mimetic\t " << mimeticResult.relativeL2Error << "\t " << mimeticResult.ergrad << "\t " << mimeticResult.sumflux
			<< "\t " << mimeticResult.erflm << "\t " << mimeticResult.uMin << "\t " << mimeticResult.uMax << "\t " << mimeticTime << std::endl;
        std::cout << "fe\t " << feResult.relativeL2Error << "\t " << feResult.ergrad << "\t " << feResult.sumflux
			<< "\t " << feResult.erflm << "\t " << feResult.uMin << "\t " << feResult.uMax << "\t " << feTime << std::endl;

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}

