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
    ScalarType waterDepth;
    ScalarType bottomElevation;
    ScalarType waterLevel;
    VelType velocity;
    SolutionType globalSolution;

    Grid& grid;

    ShallowVariableClass(Grid& grid, Scalar& initialWDepth = *(new Scalar(0)), Scalar& initialBottomElevation = *(new Scalar(0)), Scalar& initialWaterLevel = *(new Scalar(0)), Dune::FieldVector<Scalar, dim>& initialVel = *(new Dune::FieldVector<Scalar,dim>(0))) :
        elementMapper(grid, grid.leafIndexSet()), grid(grid),
                size(elementMapper.size())
    {
        sizeInitWDepth(initialWDepth, size);
        sizeInitBottomElevation(initialBottomElevation, size);
        sizeInitWaterLevel(initialWaterLevel, size);
        sizeInitVel(initialVel, size);
        sizeInitGlobalSolution(size);
    }

    //Definition of the methods to resize the vectors and fill with initial data. !!!Initial data is not the boundary condition, it´s just to initialize sth

    void sizeInitWDepth(Scalar& initialWDepth, int size)
    {
        waterDepth.resize(size);
        waterDepth=initialWDepth;
        return;
    }
    void sizeInitBottomElevation(Scalar& initialBottomElevation, int size)
    {
        bottomElevation.resize(size);
        bottomElevation=initialBottomElevation;
        return;
    }
    void sizeInitWaterLevel(Scalar& initialWaterLevel, int size)
    {
        waterLevel.resize(size);
        waterLevel=initialWaterLevel;
        return;
    }

    void sizeInitVel(Dune::FieldVector<Scalar, dim>& initialVel, int size)
    {
        velocity.resize(size);
        velocity[0]=initialVel[0];
        velocity[1]=initialVel[1];
        return;
    }

    void sizeInitGlobalSolution(int size)
    {
        globalSolution.resize(size);
        globalSolution=0;
        return;
    }

    //Declaration of return methods 

    ScalarType& returnWDepth()
    {
        return waterDepth;
    }
    ScalarType& returnBottomElevation()
    {
        return bottomElevation;
    }
    ScalarType& returnWaterLevel()
    {
        return waterLevel;
    }

    VelType& returnVel()
    {
        return velocity;
    }

    SolutionType& returnGlobalSol()
    {
        return globalSolution;
    }

    //Definition of return methods

    const ScalarType& returnWDepth(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos)
    {
        return waterDepth[elementMapper.map(element)];

    }

    const ScalarType& returnBottomElevation(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos)
    {
        return bottomElevation[elementMapper.map(element)];

    }
    const ScalarType& returnWaterLevel(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos)
    {
        waterLevel = bottomElevation[elementMapper.map(element)]+ waterDepth[elementMapper.map(element)];
        return waterLevel;
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
        //Dune::BlockVector<Dune::FieldVector<Scalar, 1> > vY(size);
        vX = 0;
        //  vY = 0;

        for (int i = 0; i < size; i++)
        {
            //std::cout<<"waterDepth = "<<waterDepth[i]<<std::endl;
            //std::cout<<"velX = "<<velocity[i][0]<<std::endl;
            //std::cout<<"velY = "<<velocity[i][1]<<std::endl;
            vX[i] = velocity[i][0];
            // vY[i] = velocity[i][1];
        }
        char fname[128];
        sprintf(fname, "%s-%05d", name, k);
        vtkwriter.addCellData(vX, "Velocity_X");
        //vtkwriter.addCellData(vY, "Velocity_Y");
        vtkwriter.addCellData(waterDepth, "waterDepth");
        vtkwriter.addCellData(bottomElevation, "bottomElevation");
        vtkwriter.addCellData(waterLevel, "waterLevel");
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
