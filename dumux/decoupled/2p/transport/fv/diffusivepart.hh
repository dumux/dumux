// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
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
#ifndef DUMUX_DIFFUSIVEPART_HH
#define DUMUX_DIFFUSIVEPART_HH

/**
 * @file
 * @brief  Base class for defining the diffusive part of an advection-diffusion equation
 * @author Bernd Flemisch, Markus Wolff
 */
namespace Dumux
{
/*!\ingroup Saturation2p
 * @brief  Base class for defining the diffusive part of an advection-diffusion equation
 *
 * @tparam TypeTag The Type Tag
 */
template<class TypeTag>
class DiffusivePart
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

    enum{dim = GridView::dimension};
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> FieldVector;

public:
    //! Returns diffusive term
    /*! Returns diffusive term for current element face
     *  @param[in] element        entity of codim 0
     *  @param[in] indexInInside  face index in reference element
     *  @param[in] satI           saturation of current element
     *  @param[in] satJ           saturation of neighbor element
     *  @param[in] pcGradient     gradient of capillary pressure between element I and J
     *  \return     diffusive term of an advection-diffusion equation
     */
    FieldVector operator() (const Element& element, const int indexInInside, Scalar satI, Scalar satJ, const FieldVector& pcGradient) const
    {
        FieldVector trivial(0);
        return trivial;
    }

    //! Returns diffusive term
    /*! Returns diffusive term for current element face
     *  @param[in] element          entity of codim 0
     *  @param[in] indexInInside    face index in reference element
     *  @param[in] satIntersection  saturation at the face between element I and J
     *  @param[in] satGradient       gradient of saturation between element I and J
     *  @param[in] time             time
     *  \return     diffusive term of an advection-diffusion equation
     */
    FieldVector operator() (const Element& element, const int indexInInside,
                                    const Scalar satIntersection, const FieldVector& satGradient, const Scalar time) const
    {
        FieldVector trivial(0);
        return trivial;
    }

    //! Returns diffusive term
    /*! Returns diffusive term for current element face
     *  @param[in] element          entity of codim 0
     *  @param[in] indexInInside    face index in reference element
     *  @param[in] satIntersection  saturation at the face between element I and J
     *  @param[in] satGradient       gradient of saturation between element I and J
     *  @param[in] time             time
     *  @param[in] satI             saturation of current element
     *  @param[in] satJ             saturation of neighbor element
     *  \return     diffusive term of an advection-diffusion equation
     */
    FieldVector operator() (const Element& element, const int indexInInside,
                                    const Scalar satIntersection, const FieldVector& satGradient, const Scalar time,
                                    Scalar satI, Scalar satJ) const
    {
        FieldVector trivial(0);
        return trivial;
    }


    //! The constructor
    /*
     *  \param problem object including the problem definition
     */
    DiffusivePart(Problem& problem)
    {}
    //! always define virtual destructor in abstract base class
    ~DiffusivePart()
    { }
};
}

#endif
