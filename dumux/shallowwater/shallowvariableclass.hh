// $Id$

#ifndef DUNE_SHALLOWVARIABLECLASS_HH
#define DUNE_SHALLOWVARIABLECLASS_HH

/**
 * @file
 * @brief  class including the variables
 * @author Markus Wolff
 */

namespace Dune
{

template<class Grid, class Scalar> class ShallowVariableClass
{
    template<int dim> struct ElementLayout
    {
        bool contains(Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    enum
        {   dim=Grid::dimension};

    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,dim> GlobalPosition;
    typedef typename Grid::LeafGridView GridView;
    typedef typename GridView::IndexSet IndexSet;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,IndexSet,ElementLayout>
    ElementMapper;

public:
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,1> > ScalarType;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,dim> > VelType;
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,dim+1> > SolutionType;

    ElementMapper elementMapper;

    ScalarType wDepth;
    VelType velocity;
    SolutionType globalSolution;

    Grid& grid;

    //Constructor: defines the initital values for the variables and the methods includes in the class

    ShallowVariableClass(Grid& grid) :
        elementMapper(grid, grid.leafIndexSet()), grid(grid),
        size(elementMapper.size())
    {
        sizeInitWDepth(size);
        sizeInitVel(size);
        sizeInitGlobalSolution(size);
    }

    //Definition of the methods to resize the vectors and fill with initial data. !!!Initial data is not the boundary condition, it´s just to initialize sth

    void sizeInitWDepth(int size)
    {
        wDepth.resize(size);
        wDepth=0;
        return;
    }

    void sizeInitVel(int size)
    {
        velocity.resize(size);
        velocity=0;
        return;
    }

    void sizeInitGlobalSolution(int size)
    {
        globalSolution.resize(size);
        globalSolution=0;
        return;
    }

    //Declaration of methods to return variable values

    ScalarType& returnWDepth()
    {
        return wDepth;
    }

    VelType& returnVel()
    {
        return velocity;
    }

    SolutionType& returnGlobalSol()
    {
        return globalSolution;
    }

    //Definition of the return methods

    const Scalar& returnWDepth(const GlobalPosition& globalPos,
                               const Element& element, const LocalPosition& localPos)
    {
        return wDepth[elementMapper.map(element)];
    }

    const VelType& returnVel(const GlobalPosition& globalPos,
                             const Element& element, const LocalPosition& localPos)
    {
        return velocity[elementMapper.map(element)];
    }

    const SolutionType& returnGlobalSol(const GlobalPosition& globalPos,
                                        const Element& element, const LocalPosition& localPos)
    {
        return globalSolution[elementMapper.map(element)];
    }

    void vtkout(const char* name, int k)
    {
        VTKWriter<typename Grid::LeafGridView> vtkwriter(grid.leafView());
        Dune::BlockVector<Dune::FieldVector<Scalar, 1> > vX(size);
        Dune::BlockVector<Dune::FieldVector<Scalar, 1> > vY(size);
        vX = 0;
        vY = 0;

        for (int i = 0; i < size; i++)
        {
            std::cout<<"wDepth = "<<wDepth[i]<<std::endl;
            std::cout<<"velX = "<<velocity[i][0]<<std::endl;
            std::cout<<"velY = "<<velocity[i][1]<<std::endl;
            vX[i] = velocity[i][0];
            vY[i] = velocity[i][1];
        }
        char fname[128];
        sprintf(fname, "%s-%05d", name, k);
        vtkwriter.addCellData(vX, "Velocity_X");
        vtkwriter.addCellData(vY, "Velocity_Y");
        vtkwriter.addCellData(wDepth, "wDepth");
        vtkwriter.write(fname, VTKOptions::ascii);

        return;
    }

    int size;
private:
    int i;
    double eps;
};
}
#endif
