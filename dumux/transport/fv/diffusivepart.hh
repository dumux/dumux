// $Id$

#ifndef DUNE_DIFFUSIVEPART_HH
#define DUNE_DIFFUSIVEPART_HH

//! \ingroup transport
//! \defgroup diffPart Diffusive transport
/**
 * @file
 * @brief  Base class for defining the diffusive part of an advection-diffusion equation
 * @author Bernd Flemisch
 */
namespace Dune
{
/*!\ingroup diffPart
 * @brief  Base class for defining the diffusive part of an advection-diffusion equation
 */
template<class GridView, class Scalar>
class DiffusivePart
{
private:
    enum{dim = GridView::dimension};
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> FieldVector;

public:
    virtual FieldVector operator() (const Element& element, const int numberInSelf, Scalar satI, Scalar satJ, const FieldVector& satGradient) const
    {
        FieldVector trivial(0);
        return trivial;
    }

    virtual FieldVector operator() (const Element& element, const int numberInSelf,
                                    const Scalar satIntersection, const FieldVector& satGradient, const Scalar time) const
    {
        FieldVector trivial(0);
        return trivial;
    }

    virtual FieldVector operator() (const Element& element, const int numberInSelf,
                                    const Scalar satIntersection, const FieldVector& satGradient, const Scalar time,
                                    Scalar satI, Scalar satJ) const
    {
        FieldVector trivial(0);
        return trivial;
    }

    virtual ~DiffusivePart()
    { }
};
}

#endif
