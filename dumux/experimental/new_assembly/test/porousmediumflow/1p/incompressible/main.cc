#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>

#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/experimental/new_assembly/dumux/common/variables.hh>


template<int numEq,
         typename ElementGeometry,
         typename ElementVariables,
         typename StorageOperator,
         typename... FluxOperator>
class CCLocalOperator
{
    using Scalar = Variables::ScalarType<typename ElementVariables::GridVariables>;

public:
    using LocalResidual = Dune::FieldVector<Scalar, numEq>;

    CCLocalOperator(const ElementGeometry& elemGeom,
                    const ElementVariables& elemVars)
    : elemGeom_(elemGeom)
    , elemVars_(elemVars)
    {}

    void evalStorage(LocalResidual& res) const
    {
        res = 0.0;
        for (const auto& scv : scvs(elemGeom_))
            res += storageOperator_();
        res *= scv.volume()*exFac;
    }

    void evalFluxes(LocalResidual& res) const
    {
        res = 0.0;
        evalFluxes_(res);
    }

    void evalSources(LocalResidual& res) const
    {
        res = 0.0;
        for (const auto& scv : scvs(elemGeom_))
            res += problem.source(...)*scv.volume()*exFac;
        res *= scv.volume()*exFac;
    }

private:
    const ElementGeometry& elemGeom_;
    const ElementVariables& elemVars_;
};


int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    Dumux::Parameters::init(argc, argv);

    static constexpr int dim = 2;
    using Scalar = double;

    // make grid
    using Grid = Dune::YaspGrid<dim>;
    const Grid grid{
        Dumux::getParam<Dune::FieldVector<Scalar, dim>>("Grid.UpperRight"),
        Dumux::getParam<std::array<int, dim>>("Grid.Cells")
    };

    // make grid geometry of the scheme to be used
    using GridGeometry = Dumux::CCTpfaFVGridGeometry<typename Grid::LeafGridView>;
    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());

    std::cout << "Test passed" << std::endl;
    return 0;
}
