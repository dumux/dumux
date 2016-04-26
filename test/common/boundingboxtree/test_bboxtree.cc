#include <config.h>
#include <iostream>
#include <unordered_map>

#include <dune/common/parametertreeparser.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/io/gridcreator.hh>
#include <dumux/common/basicproperties.hh>
#include <dumux/common/boundingboxtree.hh>

namespace Dumux {

namespace Properties
{
    NEW_TYPE_TAG(BBoxTreeTest, INHERITS_FROM(NumericModel));
#if HAVE_UG
    SET_TYPE_PROP(BBoxTreeTest, Grid, Dune::UGGrid<3>);
#else
    SET_TYPE_PROP(BBoxTreeTest, Grid, Dune::YaspGrid<3>);
#endif
}

template<class TypeTag>
class BBoxTreeTests
{
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridView = typename Grid::LeafGridView;
    using Scalar = typename Grid::ctype;
    enum { dimWorld = Grid::dimensionworld };
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    int construct(const GridView& gv)
    {
        std::cout << std::endl
                  << "Construct bounding box tree from grid view:" << std::endl
                  << "*******************************************"
                  << std::endl;

        // construct a bounding box tree
        tree_ = std::make_shared<Dumux::BoundingBoxTree<GridView>>();
        tree_->build(gv);

        return 0;
    }

    int intersectPoint(const GlobalPosition& p, std::size_t expectedCollisions)
    {
        std::cout << std::endl
                  << "Intersect with point ("<< p <<"):" << std::endl
                  << "************************************************"
                  << std::endl;

        auto entities = tree_->computeEntityCollisions(p);

        std::cout << entities.size() << " intersection(s) found ("
                  << expectedCollisions << " expected)." << std::endl;
        if (entities.size() != expectedCollisions)
        {
            std::cerr << "Point intersection failed: Expected "
                      << expectedCollisions << " and got "
                      << entities.size() << "!" <<std::endl;
            return 1;
        }
        return 0;
    }

    template <class OtherGridView>
    int intersectTree(const Dumux::BoundingBoxTree<OtherGridView>& otherTree,
                      const OtherGridView& otherGridView,
                      std::size_t expectedIntersections)
    {
        Dune::Timer timer;
        auto intersections = tree_->computeEntityCollisions(otherTree);
        std::cout << "Computed tree intersections in " << timer.elapsed() << std::endl;
        timer.reset();

        std::vector<std::vector<std::vector<GlobalPosition>>> map;
        map.resize(otherGridView.size(0));
        std::vector<std::vector<GlobalPosition>> uniqueIntersections;
        for (const auto& is : intersections)
        {
            bool add = true;
            for (const auto& i : map[is.second()])
            {
                if (intersectionsEqual(is.corners(), i))
                {
                    add = false;
                    break;
                }
            }
            if(add)
            {
                map[is.second()].push_back(is.corners());
                uniqueIntersections.push_back(is.corners());
            }
        }

        std::cout << "Found " << uniqueIntersections.size() << " unique intersections "
                  << "in " << timer.elapsed() << std::endl;

        if (uniqueIntersections.size() != expectedIntersections)
        {
            std::cerr << "BoundingBoxTree intersection failed: Expected "
                      << expectedIntersections << " and got "
                      << uniqueIntersections.size() << "!" <<std::endl;
            return 1;
        }
        return 0;
    }

private:
    template<class Intersection>
    bool intersectionsEqual(const Intersection& is1, const Intersection& is2)
    {
        auto eps = 1e-7*std::max((is1[0] - is1[1]).two_norm(), (is2[0] - is2[1]).two_norm());
        return (is1[0] - is2[0]).two_norm() < eps && (is1[1] - is2[1]).two_norm() < eps;
    }
    std::shared_ptr<Dumux::BoundingBoxTree<GridView>> tree_;
};

} // end namespace Dumux

int main (int argc, char *argv[]) try
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // Some aliases two type tags for tests using two grids
    using TypeTag = TTAG(BBoxTreeTest);
    using Grid = GET_PROP_TYPE(TypeTag, Grid);
    using Scalar = Grid::ctype;
    enum { dimWorld = Grid::dimensionworld };
    enum { dim = Grid::dimension };
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    // collect returns to determine exit code
    std::vector<int> returns;
    Dumux::BBoxTreeTests<TypeTag> test;

    for (const auto scaling : {1e10, 1.0, 1e-3, 1e-10})
    {
        std::cout << std::endl
                      << "Testing with scaling = " << scaling << std::endl
                      << "***************************************"
                      << std::endl;

        // create a cube grid
        const GlobalPosition lowerLeft(0.0);
        const GlobalPosition upperRight(1.0*scaling);
        std::array<unsigned int, dim> elems; elems.fill(3);
        auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft, upperRight, elems);

        Dune::VTKWriter<Grid::LeafGridView> vtkWriter(grid->leafGridView());
        vtkWriter.write("grid", Dune::VTK::ascii);

        // bboxtree tests using one bboxtree
        returns.push_back(test.construct(grid->leafGridView()));
        returns.push_back(test.intersectPoint(GlobalPosition(0.0), 1));
        returns.push_back(test.intersectPoint(GlobalPosition(1e-3*scaling), 1));
        returns.push_back(test.intersectPoint(GlobalPosition(1.0/3.0*scaling), 8));

#if HAVE_DUNE_FOAMGRID
        using NetworkGrid = Dune::FoamGrid<1, 3>;
        using NetworkGridView = NetworkGrid::LeafGridView;

        std::cout << std::endl
                      << "Intersect with other bounding box tree:" << std::endl
                      << "***************************************"
                      << std::endl;

        auto networkGrid = std::shared_ptr<NetworkGrid>(Dune::GmshReader<NetworkGrid>::read("network1d.msh", false, false));

        // scaling
        for (const auto& vertex : vertices(networkGrid->leafGridView()))
        {
            auto newPos = vertex.geometry().corner(0);
            newPos *= scaling;
            networkGrid->setPosition(vertex, newPos);
        }

        Dune::VTKWriter<NetworkGridView> lowDimVtkWriter(networkGrid->leafGridView());
        lowDimVtkWriter.write("network", Dune::VTK::ascii);

        std::cout << "Constructed " << networkGrid->leafGridView().size(0) <<
                                  " element 1d network grid." << std::endl;

        Dumux::BoundingBoxTree<NetworkGridView> networkTree;
        networkTree.build(networkGrid->leafGridView());
        returns.push_back(test.intersectTree(networkTree, networkGrid->leafGridView(), 9));
#endif
    }

    std::cout << std::endl;

    // determine the exit code
    if (std::any_of(returns.begin(), returns.end(), [](int i){ return i==1; }))
        return 1;

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dumux::ParameterException &e) {
    std::cerr << e << ". Abort!\n";
    return 1;
}
catch (const Dune::Exception& e) {
    std::cout << e << std::endl;
    return 1;
}
