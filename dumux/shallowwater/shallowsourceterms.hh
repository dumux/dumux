// $Id$

#ifndef DUNE_SHALLOW_SOURCE_TERMS_HH
#define DUNE_SHALLOW_SOURCE_TERMS_HH

namespace Dune
{

template<class Grid, class Scalar> class ShallowSourceTerms
{
public:

    enum
    {   dim=Grid::dimension};
    enum
    {   dimWorld = Grid::dimensionworld};
    typedef Dune::FieldVector<Scalar, dim> VelType;
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;

    virtual VelType operator()(Scalar waterLevelI, Scalar waterLevelJ, VelType nVec,
            Scalar gravity, Scalar bottomElevationI, Scalar bottomElevationJ)
    {
        VelType trivial(0);
        return trivial;
    }

    virtual VelType operator()(VelType velocityI, VelType velocityJ, Scalar waterDepthL,
            Scalar waterDepthR, VelType nVec, Scalar gravity)
    {
        VelType trivial(0);
        return trivial;
    }

    virtual ~ShallowSourceTerms()
    {
    }
};
}
#endif

