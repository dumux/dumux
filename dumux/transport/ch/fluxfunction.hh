// $Id$

#ifndef FLUXFUNCTION_HH
#define FLUXFUNCTION_HH

//! \ingroup transport
//! \defgroup flux function transport
/**
 * @file
 * @brief  Base class for defining the flux function of an advection-diffusion equation
 * @author Yufei Cao
 */

namespace Dune
{
//! \ingroup transport
//! Compute the flux function of the transport equation
template<class Grid, class Scalar>
class FluxFunction {
public:
    /*! \brief Realizes the numerical flux function.
     *
     *  \param element cell I
     *  \param faceNumber index of faceIJ between cell I and cell J
     *  \param satI usually the saturation value of cell I
     *  \param satJ usually the saturation value of cell J
     */
private:
    enum{dim = Grid::dimension, dimWorld = Grid::dimensionworld};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    virtual Scalar operator() (Scalar saturationW, const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, Scalar T=283.15, Scalar p=1e5) const = 0;

    virtual ~FluxFunction()
    { }
};
}

#endif
