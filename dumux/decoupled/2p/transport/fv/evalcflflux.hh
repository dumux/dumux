// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                 *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_EVALCFLFLUX_HH
#define DUMUX_EVALCFLFLUX_HH

/**
 * @file
 * @brief  Base class for implementations of different kinds of fluxes to evaluate a CFL-Condition
 * @author Markus Wolff
 */
namespace Dumux
{
/*!\ingroup Transport2p
 * @brief  Base class for implementations of different kinds of fluxes to evaluate a CFL-Condition
 *
 *  Base class for implementations of different kinds of fluxes to evaluate a CFL-Condition of the form
 *
 *  \f[
 *    \frac{F_i \Delta t}{V_p_i} < 1
 *  \f]
 *  where \f(V_p_i\f) is the pore volume of cell i and
 *
 *  \f[
 *    F_i = \sum f_{ij}
 *  \f]
 *
 * with \f(f_{ij}\f) being the CFL-flux over edge \f(ij\f).
 *
 * Template parameters are:

 - TypeTag PropertyTag of the problem implementation
 */
template<class TypeTag>
class EvalCflFlux
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:

    void addFlux(Scalar& lambdaW, Scalar& lambdaNW, Scalar& viscosityW, Scalar& viscosityNW, Scalar flux, const Intersection& intersection, int phaseIdx)
    {}

    void addFlux(Scalar& lambdaW, Scalar& lambdaNW, Scalar& viscosityW, Scalar& viscosityNW, Scalar flux, const Element& element, int phaseIdx)
    {}

    Scalar getCFLFluxFunction(int phaseIdx = 0)
    {
        return 0.0;
    }

    Scalar getCFLFluxIn(int phaseIdx = 0)
    {
        return 0.0;
    }

    Scalar getCFLFluxOut(int phaseIdx = 0)
    {
        return 0.0;
    }

    Scalar getCFLFluxFunction(const GlobalPosition& globalPos, const Element& element)
    {
        return 0.0;
    }

    /*! @brief Constructs a EvalCflFlux instance */
    EvalCflFlux ()
    {}
private:

};
}

#endif
