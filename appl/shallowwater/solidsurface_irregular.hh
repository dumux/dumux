// $Id: solidsurfaceplain.hh 1501 2009-03-26 15:52:22Z anneb $

#ifndef SOLIDSURFACE_IRREGULAR_HH
#define SOLIDSURFACE_IRREGULAR_HH

#include <dumux/shallowwater/solidsurfacebase.hh>

namespace Dune
{

template<class Grid, class Scalar> class SolidSurfacePlain :
    public SolidSurfaceBase<Grid,Scalar>
{
    enum
    {   dim=Grid::dimension};

    enum
    {   dimWorld = Grid::dimensionworld};

    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,dim> GlobalPosition;

public:

    Scalar evalBottomElevation(const GlobalPosition& globalPos)
    {

        if (globalPos[0]<8 ||globalPos[0]> 12)
        {
            bottomElevation_ = 0*globalPos[0];
            return bottomElevation_;
        }
        if (8 <= globalPos[0]<= 12)
        {
            bottomElevation_= 0.2 - 0.05*((globalPos[0]-10)*(globalPos[0]-10));
            return bottomElevation_;
        }
        return bottomElevation_;

    }
    FieldVector<Scalar, dim> evalBottomSlopes()
    {
        switch (dim)
        {
        case 1:
            bottomSlope_=crossSlope_;
            break;
        case 2:
            bottomSlope_[0]=crossSlope_;
            bottomSlope_[1]=longSlope_;
            break;
        }

        return bottomSlope_;
    }

    SolidSurfacePlain() :
        SolidSurfaceBase<Grid, Scalar>(), bottomElevation_(0),
                crossSlope_(-0.0), longSlope_(0.0), eps_(1e-12)
    {

    }

private:
    //Definition der Variablen für Quer- und Längsneigung
    Scalar bottomElevation_;
    Scalar crossSlope_;
    Scalar longSlope_;
    FieldVector<Scalar, dim> bottomSlope_;
    Scalar eps_;

};

} // end namespace
#endif
