// $Id$
#ifndef DUNE_VARIABLECLASSSUBPROBS_HH
#define DUNE_VARIABLECLASSSUBPROBS_HH

/**
 * @file
 * @brief  class including the variables
 * @author Markus Wolff
 */

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar> class VariableClassSubProbs
{
    template<int dim> struct ElementLayout
    {
        bool contains(Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };
    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::Traits::LevelIndexSet IndexSet;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,IndexSet,ElementLayout>
            ElementMapper;

public:
    typedef Dune::BlockVector< Dune::FieldVector<Scalar,1> > ScalarVectorType;
     typedef Dune::BlockVector< FieldVector<FieldVector<Scalar, dim>, 2*dim> > VelType;


    Grid& grid;
    ElementMapper diffMapper;
    ElementMapper transMapper;
    int diffSize;
    int transSize;
    int diffLevel;
    int transLevel;

    ScalarVectorType saturation;
    ScalarVectorType pressure;
    VelType velocity;

    VariableClassSubProbs(Grid& grid, int lev)
    : grid(grid), diffMapper(grid, grid.levelIndexSet(lev)), transMapper(grid, grid.levelIndexSet(lev)),
      diffSize(diffMapper.size()),transSize(transMapper.size()),
      diffLevel(lev),transLevel(lev)
    {
        initsat(transSize);
        initpress(diffSize);
        initvel(transSize);
    }

    void initsat(int size)
    {
        saturation.resize(size);
        saturation=0;
        return;
    }

    void initpress(int size)
    {
        pressure.resize(size);
        pressure=0;
        return;
    }
    void initvel(int size)
    {
        Dune::FieldVector<Scalar, dim> vel(0);
        velocity.resize(size);
        velocity = vel;

        return;
    }

    ScalarVectorType& sat() const
    {
        return saturation;
    }
    ScalarVectorType& press() const
    {
        return pressure;
    }
    VelType& vel() const
    {
        return velocity;
    }

    const Dune::FieldVector<Scalar,1>& sat(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        return saturation[transMapper.map(element)];;
    }

    const Dune::FieldVector<Scalar,1>& press(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        return pressure[diffMapper.map(element)];
    }

    const Dune::FieldVector<Scalar,dim>& vTotal(const Element& element,
            const int numberInSelf) const
    {
        int elemId = transMapper.map(element);

        return (velocity[elemId][numberInSelf]);
    }
};
}
#endif
