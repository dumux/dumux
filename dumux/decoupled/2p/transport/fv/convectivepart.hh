// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
 *   Institute of Hydraulic Engineering                                      *
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
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    enum{dim = GridView::dimension, dimWorld = GridView::dimensionworld};
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    //! Returns convective term
    /*! Returns convective term for current element face
     *  @param[in] element        entity of codim 0
     *  @param[in] indexInInside  face index in reference element
     *  @param[in] sat           saturation of current element
     *  \return     convective term of an advection-diffusion equation
     */
    Scalar operator() (const Element& element, const int indexInInside, const Scalar sat) const
    {
        Scalar trivial(0);
        return trivial;
    }
    //! Returns convective term
    /*! Returns convective term for current element face
     *  @param[in] element        entity of codim 0
     *  @param[in] indexInInside  face index in reference element
     *  @param[in] satI           saturation of current element
     *  @param[in] satJ           saturation of neighbor element
     *  \return     convective term of an advection-diffusion equation
     */
    Dune::FieldVector<Scalar, dimWorld> operator() (const Element& element, const int indexInInside, const Scalar satI, const Scalar satJ) const
    {
        Dune::FieldVector<Scalar, dimWorld> trivial(0);
        return trivial;
    }

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
