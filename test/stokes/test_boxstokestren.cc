#include "config.h"
#include <iostream>
#define DUMMY
#ifdef DUMMY
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfalberta.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/timedisc/timeloop.hh"
#include "boxstokes.hh"
#include "yxproblem.hh"
#include "sinproblem.hh"
#include "curlproblem.hh"
#include "dumux/material/phaseproperties/phaseproperties_waterair.hh"
#include "dumux/material/multicomponentrelations.hh"

#include "boxstokestren.hh"
#include "FFproblem.hh"
#include "steproblem.hh"

template<int dim>
struct VertexLayout
{
    bool contains (Dune::GeometryType gt)
    {
        return (gt.dim() == 0);
    }
};

template<int dim>
struct ElementLayout
{
    bool contains (Dune::GeometryType gt)
    {
        return (gt.dim() == dim);
    }
};

template<class Vector, class Grid, class Problem>
void calculateError(const Grid& grid, const Problem& problem, Vector& solution)
{
    typedef typename Grid::ctype Scalar;
    enum {dim=Grid::dimension};
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::template Codim<0>::Iterator ElementIterator;
    typedef typename Grid::LeafGridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename Grid::LeafGridView::IndexSet IS;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,IS,ElementLayout> ElementMapper;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,IS,VertexLayout> VertexMapper;

    VertexMapper vertexMapper(grid, grid.leafView().indexSet());
    ElementMapper elementMapper(grid, grid.leafView().indexSet());

    Scalar errPressure = 0;
    Scalar errVelocity = 0;
    Scalar constant = 0;
    Scalar errPartialDensity = 0;
    Scalar errTemperature = 0;

    ElementIterator endEIt = grid.template leafend<0>();
    for (ElementIterator eIt = grid.template leafbegin<0>(); eIt != endEIt; ++eIt)
    {
        const Element& element = *eIt;

        Dune::GeometryType geomType = element.geometry().type();

        const Dune::FieldVector<Scalar,dim>& local = Dune::ReferenceElements<Scalar,dim>::general(geomType).position(0, 0);
        Dune::FieldVector<Scalar,dim> global = element.geometry().global(local);

        Scalar volume = element.geometry().integrationElement(local)
                *Dune::ReferenceElements<Scalar,dim>::general(geomType).volume();

        int eIdx = elementMapper.map(element);

        Scalar approxPressure = solution.evallocal (4, element, local);
        Scalar exactPressure = problem.pressure(global);
        if (eIdx == 0)
            constant = exactPressure - approxPressure;
        approxPressure += constant;
        errPressure += volume*(approxPressure - exactPressure)*(approxPressure - exactPressure);

        //transport
        Scalar approxPartialDensity = solution.evallocal (2, element, local);
        Scalar exactPartialDensity = problem.partialdensity(global);
        errPartialDensity += volume*(approxPartialDensity - exactPartialDensity)*(approxPartialDensity - exactPartialDensity);

        //energy
        Scalar approxTemperature = solution.evallocal (3, element, local);
        Scalar exactTemperature = problem.temperature(global);
        errTemperature += volume*(approxTemperature - exactTemperature)*(approxTemperature - exactTemperature);

        //velocity
        Scalar approxXV = solution.evallocal (0, element, local);
        Scalar exactXV = problem.velocity(global)[0];
        errVelocity += volume*(approxXV - exactXV)*(approxXV - exactXV);

        Scalar approxYV = solution.evallocal (1, element, local);
        Scalar exactYV = problem.velocity(global)[1];
        errVelocity += volume*(approxYV - exactYV)*(approxYV - exactYV);
    }

    errPressure = sqrt(errPressure);
    errVelocity = sqrt(errVelocity);
    errPartialDensity = sqrt(errPartialDensity);
    errTemperature = sqrt(errTemperature);

    std::cout << "Error in discrete L2 norm:\nPressure: " << errPressure << "\nVelocity: " << errVelocity << "\nPdensity: " << errPartialDensity << "\nTemperature: " << errTemperature <<std::endl;
}

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim = 2;
    double tStart = 0;
    double tEnd = 1;
    double dt = 0.1;
    int modulo = 1; // print every modulo'th time step

    // create a grid object
    typedef double NumberType;
    typedef Dune::SGrid<dim,dim> GridType;
    //typedef Dune::AlbertaGrid<dim,dim> GridType;

    if (argc != 2 && argc != 3) {
        std::cout << "Usage: test_boxstokestren dgffilename [refinementsteps]" << std::endl;
        return (1);
    }
    int refinementSteps = 0;
    if (argc == 3) {
        std::string arg2(argv[2]);
        std::istringstream is2(arg2);
        is2 >> refinementSteps;
    }

    Dune::GridPtr<GridType> gridPtr( argv[1] );
    GridType& grid = *gridPtr;

    if (refinementSteps)
        grid.globalRefine(refinementSteps);

    Dune::Liq_WaterAir dummyPhase;
    Dune::Gas_WaterAir gasPhase;

    Dune::CWaterAir multicomp(dummyPhase, gasPhase);

    Dune::SteProblem<GridType, double> problem(gasPhase, multicomp);
    typedef Dune::LeafP1BoxStokesTrEn<GridType, NumberType, dim> BoxStokesTrEn;
    BoxStokesTrEn boxStokesTrEn(grid, problem);

    Dune::TimeLoop<GridType, BoxStokesTrEn> timeloop(tStart, tEnd, dt, "test_boxstokestren", modulo);

    timeloop.execute(boxStokesTrEn);

    //printvector(std::cout, *(boxStokes.u), "solution", "row", 200, 1, 3);
    calculateError(grid, problem, boxStokesTrEn.u);

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

  return 1;
}
catch (...)
{
    std::cerr << "Generic exception!" << std::endl;
    return 2;
}
#endif
