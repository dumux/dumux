// $Id$

#ifndef SOLIDSURFACEPLAIN_HH
#define SOLIDSURFACEPLAIN_HH

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

    Scalar evalBottomElevation(const GlobalPosition& globalPos,
                               const Element& element, const LocalPosition& localPos)
    {
        bottomElevation_ = 100 + crossSlope_*globalPos[0]+ longSlope_
            *globalPos[1];

        return bottomElevation_;

    }
    FieldVector<Scalar, dim> calcBottomSlopes(const GlobalPosition& globalPos,
                                              const Element& element, const LocalPosition& localPos)
    {

        bottomSlope_[0]=crossSlope_;
        bottomSlope_[1]=longSlope_;

        return bottomSlope_;
    }

    SolidSurfacePlain() :
        SolidSurfaceBase<Grid, Scalar>(), crossSlope_(-0.02), longSlope_(0.0),
        eps_(1e-12)
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
