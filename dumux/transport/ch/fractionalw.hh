// $Id:$
/*****************************************************************************
 *   Copyright (C) <YEARS> by <ADD_AUTHOR_HERE>                              *
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

#ifndef DUNE_FRACTIONALW_HH
#define DUNE_FRACTIONALW_HH

#include "dumux/transport/ch/fluxfunction.hh"
#include "dumux/transport/transportproblem.hh"

//! \ingroup transport
//! \defgroup flux function transport
/**
 * @file
 * @brief  Base class for defining the flux function of an advection-diffusion equation
 * @author Yufei Cao
 */

namespace Dune
{
/*!\ingroup transport
 * @brief  Base class for defining the flux function of an advection-diffusion equation
 */
template<class GridView, class Scalar, class VC, class Problem = TransportProblem<GridView, Scalar, VC> >
class FractionalW : public FluxFunction<GridView,Scalar> {
    enum{dim = GridView::dimension,dimWorld = GridView::dimensionworld};
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     *  \param saturationW the saturation of the wetting phase
     *  \param T temperature
     *  \param p pressure
     *  \return the fractional flow function of the wetting phase
     */
    virtual Scalar operator() (Scalar saturationW, const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, Scalar T=283.15, Scalar p=1e5) const
    {
        Scalar result = problem_.materialLaw().fractionalW(saturationW, globalPos, element,localPos, T, p);
        return result;
    }

    FractionalW (Problem& problem)
        : problem_(problem)
    { }

private:
    Problem& problem_;
};
}

#endif
