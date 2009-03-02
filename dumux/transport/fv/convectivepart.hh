// $Id$

#ifndef DUNE_CONVECTIVEPART_HH
#define DUNE_CONVECTIVEPART_HH

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class ConvectivePart
{
private:
    enum{dim = Grid::dimension, dimWorld = Grid::dimensionworld};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:

    virtual Scalar operator() (const Element& element, const Scalar sat, const GlobalPosition& faceGlobal) const
    {
        Scalar trivial(0);
        return trivial;
    }

    virtual ~ConvectivePart()
    { }
};
}

#endif
