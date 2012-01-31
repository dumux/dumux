// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_CONVECTIVEPART_HH
#define DUMUX_CONVECTIVEPART_HH

#include <dumux/decoupled/2p/2pproperties.hh>

/**
 * @file
 * @brief  Base class for defining a convective part of an advection-diffusion equation
 * @author Markus Wolff
 */

namespace Dumux
{

/*!\ingroup Saturation2p
 * @brief  Base class for defining the convective part of an advection-diffusion equation
 *
 * @tparam TypeTag The Type Tag
 */

template<class TypeTag>
class ConvectivePart
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    enum{dim = GridView::dimension, dimWorld = GridView::dimensionworld};
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> FieldVector;

public:
    //! Returns convective term
    /*! Returns convective term for current element face
     *  @param[in] element        entity of codim 0
     *  @param[in] indexInInside  face index in reference element
     *  @param[in] sat           saturation of current element
     *  \return     convective term of an advection-diffusion equation
     */
    Scalar getFlux(const Intersection& intersection, const Scalar sat) const
    {
        return 0.0;
    }
    //! Returns convective term
    /*! Returns convective term for current element face
     *  @param[in] element        entity of codim 0
     *  @param[in] indexInInside  face index in reference element
     *  @param[in] satI           saturation of current element
     *  @param[in] satJ           saturation of neighbor element
     *  \return     convective term of an advection-diffusion equation
     */
    void getFlux(FieldVector& flux, const Intersection& intersection, const Scalar satI, const Scalar satJ) const
    {}

    //! The constructor
    /*
     *  \param problem object including the problem definition
     */
    ConvectivePart(Problem& problem)
    {}

    //! always define virtual destructor in abstract base class
    virtual ~ConvectivePart()
    { }
};
}

#endif
