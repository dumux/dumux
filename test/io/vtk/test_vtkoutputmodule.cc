#include <config.h>

#include <array>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/parameters.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dune::MPIHelper::instance(argc, argv);

    Parameters::init([](Dune::ParameterTree& params)
    {
        params["Double.Vtk.Precision"] = "Float64";
        params["Single.Vtk.Precision"] = "Float32";
        params["DoubleCoord.Vtk.CoordPrecision"] = "Float64";
    });

    using Grid = Dune::YaspGrid<2>;
    Grid grid({1.0,1.0}, {2,2});

    using GridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView>;
    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());

    std::vector<int> integers({0, 1, 2, 3});
    std::vector<float> floats({0.1, 1.2, 2.3, 3.4});
    std::vector<double> doubles({0.123456789101112131415, 1.3, 2.5, 3.9});

    // single precision
    {
        VtkOutputModuleBase<GridGeometry> vtkWriter(*gridGeometry, "test_vtkoutputmodule_float", "Single");

        vtkWriter.addField(integers, "integer", Vtk::Precision::int32);
        vtkWriter.addField(floats, "float", Vtk::Precision::float32);
        vtkWriter.addField(doubles, "double", Vtk::Precision::float64);
        vtkWriter.addField(doubles, "default");
        vtkWriter.write(0.0);
    }
    // double precision
    {
        VtkOutputModuleBase<GridGeometry> vtkWriter(*gridGeometry, "test_vtkoutputmodule_double", "Double");

        vtkWriter.addField(integers, "integer", Vtk::Precision::int32);
        vtkWriter.addField(floats, "float", Vtk::Precision::float32);
        vtkWriter.addField(doubles, "double", Vtk::Precision::float64);
        vtkWriter.addField(doubles, "default");
        vtkWriter.write(0.0);
    }
    // single precision, coordinates double
    {
        VtkOutputModuleBase<GridGeometry> vtkWriter(*gridGeometry, "test_vtkoutputmodule_doublecoord", "DoubleCoord");

        vtkWriter.addField(integers, "integer", Vtk::Precision::int32);
        vtkWriter.addField(floats, "float", Vtk::Precision::float32);
        vtkWriter.addField(doubles, "double", Vtk::Precision::float64);
        vtkWriter.addField(doubles, "default");
        vtkWriter.write(0.0);
    }
    // single precision, all float32
    {
        VtkOutputModuleBase<GridGeometry> vtkWriter(*gridGeometry, "test_vtkoutputmodule_allfloat", "Single");

        vtkWriter.addField(integers, "integer", Vtk::Precision::float32);
        vtkWriter.addField(floats, "float", Vtk::Precision::float32);
        vtkWriter.addField(doubles, "double", Vtk::Precision::float32);
        vtkWriter.addField(doubles, "default");
        vtkWriter.write(0.0);
    }

    Parameters::print();

    return 0;
} // end main
