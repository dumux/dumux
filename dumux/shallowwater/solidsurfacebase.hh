// $Id$

#ifndef SOLIDSURFACEBASE_HH
#define SOLIDSURFACEBASE_HH

#include <dune/common/fvector.hh>
#include <vector>

namespace Dune
{

template<class Grid, class Scalar> class SolidSurfaceBase
{
public:
    enum
    {   dim=Grid::dimension};

    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef FieldVector<Scalar,dim> LocalPosition;
    typedef FieldVector<Scalar,dim> GlobalPosition;

    virtual Scalar frictionRelationType(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos)
    {
        DUNE_THROW(NotImplemented, "friction term not implemented!");
    }

    virtual Scalar friction(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos)
    {
        DUNE_THROW(NotImplemented, "friction term not implemented!");
    }

    virtual Scalar evalBottomElevation(const GlobalPosition& globalPos) = 0;
    
    virtual FieldVector<Scalar, dim> evalBottomSlopes() = 0;

    virtual FieldVector<Scalar, dim> calcBottomSlopes(
            const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos)
            {
        DUNE_THROW(NotImplemented, "calculation of bottom slopes not implemented!");
            }

    virtual ~SolidSurfaceBase()
    {
    }
};

} // end namespace
#endif /*PROPERTY_BASECLASSES*/
