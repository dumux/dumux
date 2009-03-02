// commented lines 1454, 1464-1467 in istl/communicator.hh
#include <config.h>
#include <iostream>
#undef DUMMY
#ifdef DUMMY
#include<mpi.h>
#include <dune/grid/alugrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/common/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dumux/timedisc/timeloop.hh>
#include "parallelboxdiffusion.hh"

namespace Dune
{
template<int dim>
struct P1Layout
{
    bool contains (Dune::GeometryType gt)
    {
        return gt.dim() == 0;
    }
};
template<class Grid, class Solution, class Problem>
double discreteError(const Grid& grid, const Solution& solution, const Problem& problem)
{
    enum{dim=Grid::dimension};
    typedef typename Grid::LeafGridView GV;
    typedef typename GV::IndexSet IS;
    typedef typename GV::template Codim<dim>::Iterator VertexIterator;
    typedef MultipleCodimMultipleGeomTypeMapper<Grid,IS,P1Layout> VM;

    VM vertexMapper(grid, grid.leafIndexSet());
    double error = 0.0;

    VertexIterator endIt = grid.leafView().template end<dim>();
    VertexIterator it = grid.leafView().template begin<dim>();
    for (; it != endIt; ++it)
        {
            const PartitionType& partitionType = (*it).partitionType();
            //std::cout << "type = " << partitionType << std::endl;

            if (partitionType != GhostEntity)
                {
                    // get exact solution at vertex
                    FieldVector<double,dim> globalCoord = (*it).geometry().corner(0);
                    double exact = problem.exact(globalCoord);

                    // get approximate solution at vertex
                    int globalId = vertexMapper.map(*it);
                    double approximate = (*solution)[globalId];

                    //          std::cout << grid.comm().rank() << ": id = " << globalId << ", coord = " << globalCoord << ", type = " << partitionType << std::endl;

                    error += (exact - approximate)*(exact - approximate);
                }
        }

    return sqrt(error)/(*solution).two_norm();
}
}

int main(int argc, char** argv)
{
    try{

        Dune::MPIHelper::instance(argc, argv);

        // define the problem dimensions
        const int dim=3;
        typedef double NumberType;
        if (argc != 2 && argc != 3) {
            std::cout << "usage: test_parallel dgffilename/basefilename [refinementsteps]" << std::endl;
            return 1;
        }
        int refinementSteps = 0;
        if (argc == 3) {
            std::string arg2(argv[2]);
            std::istringstream is2(arg2);
            is2 >> refinementSteps;
        }

        // instantiate a distributed grid with overlap
        //    Dune::FieldVector<double,dim> length(8.0);
        //    Dune::FieldVector<int,dim> size(refinementSteps);
        //    Dune::FieldVector<bool,dim> periodic(false);
        //    int overlap = 0;
        //    typedef Dune::YaspGrid<dim,dim> GridType;
        //    GridType grid(MPI_COMM_WORLD, length, size, periodic, overlap);

        // create a grid object
        typedef Dune::ALUSimplexGrid<dim,dim> GridType;
        //typedef Dune::ALUCubeGrid<dim,dim> GridType;

        // create grid pointer
        Dune::GridPtr<GridType> gridPtr( argv[1], Dune::MPIHelper::getCommunicator() );
        // grid reference
        GridType& grid = *gridPtr;

        grid.loadBalance();

        if (refinementSteps)
            grid.globalRefine(refinementSteps);

        Dune::gridinfo(grid);

        Dune::Timer timer;
        timer.reset();

        DiffusionParameters<GridType,NumberType> problem;

        typedef Dune::LeafP1ParallelBoxDiffusion<GridType, NumberType> Diffusion;
        Diffusion diffusion(grid, problem);
        //    discreteError(grid, *diffusion, problem);

        Dune::TimeLoop<GridType, Diffusion> timeloop(0, 1, 1, "test_parallel", 1);

        timeloop.execute(diffusion);

        double discreteErr = discreteError(grid, *diffusion, problem);
        grid.comm().sum(&discreteErr, 1);
        double elapsedTime = timer.elapsed();
        grid.comm().max(&elapsedTime, 1);

        if (grid.comm().rank() == 0) {
            std::cout << "discrete error = " << discreteErr << std::endl;
            std::cout << "Calculation took " << elapsedTime << " seconds." << std::endl;
        }
        //    char buffer[128];
        //    sprintf(buffer, "rank %d :", grid.comm().rank());
        //    printvector(std::cout, *(*diffusion), "solution", buffer, 200, 1, 3);

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
#else

int main (int argc , char **argv) try
    {
        std::cout << "This test is not finished yet." << std::endl;
        //  std::cout << "Please install MPI." << std::endl;

        return 1;
    }
 catch (...)
     {
         std::cerr << "Generic exception!" << std::endl;
         return 2;
     }
#endif

