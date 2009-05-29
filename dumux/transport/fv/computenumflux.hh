// $Id$

#ifndef DUNE_COMPUTENUMFLUX_HH
#define DUNE_COMPUTENUMFLUX_HH


//! \ingroup transport
//! \defgroup numerical flux transport
/**
 * @file
 * @brief  Base class for defining the numerical flux of an advection-diffusion equation
 * @author Yufei Cao
 */

namespace Dune
{
//! \ingroup transport
//! Compute the numerical flux of the transport equation
template<class GridView, class Scalar>
class ComputeNumFlux {
public:
    /*! \brief Realizes the numerical flux function.
     *
     *  \param element cell I
     *  \param faceNumber index of faceIJ between cell I and cell J
     *  \param satI usually the saturation value of cell I
     *  \param satJ usually the saturation value of cell J
     */
private:
    enum{dim = GridView::dimension};
    typedef typename GridView::Traits::template Codim<0>::Entity Element;

public:
    virtual Scalar operator() (const Element& element, const int faceNumber, const Scalar satI, const Scalar satJ) const
    {
        Scalar trivial(0);
        return trivial;
    }

    virtual ~ComputeNumFlux()
    { }
};
}

#endif
