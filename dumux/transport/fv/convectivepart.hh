// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
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
#ifndef DUNE_CONVECTIVEPART_HH
#define DUNE_CONVECTIVEPART_HH

//! \ingroup transport
//! \defgroup convPart Convective transport
/**
 * @file
 * @brief  Base class for defining a convective part of an advection-diffusion equation
 * @author Markus Wolff
 */

namespace Dune
{

/*!\ingroup convPart
 * @brief  Base class for defining the convective part of an advection-diffusion equation
 *
 * Template parameters are:

 - GridView      a DUNE gridview type
 - Scalar        type used for scalar quantities
 */

template<class GridView, class Scalar>
class ConvectivePart
{
private:
    enum{dim = GridView::dimension, dimWorld = GridView::dimensionworld};
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    //! Returns convective term
    /*! Returns convective term for current element face
     *  @param[in] element        entity of codim 0
     *  @param[in] sat           saturation of current element
     *  @param[in] indexInInside  face index in reference element
     *  \return     convective term of an advection-diffusion equation
     */
    virtual Scalar operator() (const Element& element, const Scalar sat, const int indexInInside) const
    {
        Scalar trivial(0);
        return trivial;
    }
    //! Returns convective term
    /*! Returns convective term for current element face
     *  @param[in] element        entity of codim 0
     *  @param[in] satI           saturation of current element
     *  @param[in] satJ           saturation of neighbor element
     *  @param[in] indexInInside  face index in reference element
     *  \return     convective term of an advection-diffusion equation
     */
    virtual FieldVector<Scalar, dimWorld> operator() (const Element& element, const Scalar satI, const Scalar satJ, const int indexInInside) const
    {
        FieldVector<Scalar, dimWorld> trivial(0);
        return trivial;
    }

    //! always define virtual destructor in abstract base class
    virtual ~ConvectivePart()
    { }
};
}

#endif
