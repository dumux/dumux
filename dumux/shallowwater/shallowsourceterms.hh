// $Id$

#ifndef DUNE_SHALLOW_SOURCE_TERMS_HH
#define DUNE_SHALLOW_SOURCE_TERMS_HH

namespace Dune
{

template<class GridView, class Scalar> class ShallowSourceTerms
{
public:

    enum
    {   dim=GridView::dimension};
    enum
    {   dimWorld = GridView::dimensionworld};
    typedef Dune::FieldVector<Scalar, dim> VelType;
    typedef typename GridView::Traits::template Codim<0>::Entity Entity;

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

