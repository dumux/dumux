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

template<class G, class RT> class VariableClassSubProbs
{
    template<int dim> struct ElementLayout
    {
        bool contains(Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };
public:
    enum
    {    n=G::dimension};
    typedef typename G::ctype DT;
    typedef Dune::BlockVector< Dune::FieldVector<RT,1> > ScalarType;
    typedef Dune::BlockVector< Dune::FieldVector<RT,n> > SlopeType;
    typedef Dune::BlockVector< Dune::FieldVector<Dune::FieldVector<double, n>, 2*n> >
            VelType;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename G::Traits::LevelIndexSet IS;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout>
            ElementMapper;

    G& grid;
    ElementMapper diffMapper;
    ElementMapper transMapper;
    int diffSize;
    int transSize;
    int diffLevel;
    int transLevel;

    ScalarType saturation;
    ScalarType pressure;
    VelType velocity;
    SlopeType slope;

    VariableClassSubProbs(G& g, int lev)
    : grid(g), diffMapper(g, g.levelIndexSet(lev)), transMapper(g, g.levelIndexSet(lev)),
      diffSize(diffMapper.size()),transSize(transMapper.size()),
      diffLevel(lev),transLevel(lev)
    {
        initsat(transSize);
        initpress(diffSize);
        initvel(transSize);
        initslopes(transSize);
    }

    void initsat(int size)
    {
        saturation.resize(size);
        saturation=0;
        return;
    }

    void initslopes(int size)
        {
            slope.resize(size);
            slope=0;
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
        Dune::FieldVector<double, n> vel(0);
        velocity.resize(size);
        velocity = vel;

        return;
    }

    ScalarType& sat() const
    {
        return saturation;
    }
    ScalarType& press() const
    {
        return pressure;
    }
    VelType& vel() const
    {
        return velocity;
    }

    const Dune::FieldVector<RT,1>& sat(const Dune::FieldVector<DT,n>& x,
            const Entity& e, const Dune::FieldVector<DT,n>& xi) const
    {
        return saturation[transMapper.map(e)];;
    }

    const Dune::FieldVector<RT,1>& press(const Dune::FieldVector<DT,n>& x,
            const Entity& e, const Dune::FieldVector<DT,n>& xi) const
    {
        return pressure[diffMapper.map(e)];
    }

    const Dune::FieldVector<DT,n>& vTotal(const Entity& e,
            const int numberInSelf) const
    {
        int elemId = transMapper.map(e);

        return (velocity[elemId][numberInSelf]);
    }
};
}
#endif
