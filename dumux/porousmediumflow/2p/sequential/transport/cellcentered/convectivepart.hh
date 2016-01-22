// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_CONVECTIVEPART_HH
#define DUMUX_CONVECTIVEPART_HH

#include <dumux/porousmediumflow/2p/sequential/properties.hh>

/**
 * \file
 * \brief  Base class for defining a convective part of the saturation transport equation
 */

namespace Dumux
{

/*!\ingroup FVSaturation2p
 * \brief  Base class for defining a convective part of the saturation transport equation
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class ConvectivePart
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    enum{dimWorld = GridView::dimensionworld};
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;

public:
    //! For initialization
    void initialize()
    {}

    /*! \brief Returns convective term for current element face
     *  \param intersection  Intersection of two grid elements/global boundary
     *  \param sat           Saturation of current element
     *  \return     Convective flux
     */
    Scalar getFlux(const Intersection& intersection, const Scalar sat) const
    {
        return 0.0;
    }

    /*! \brief Returns convective term for current intersection
     *
     *  \param flux        Flux vector (gets the flux from the function)
     *  \param intersection  Intersection of two grid elements/global boundary
     *  \param satI           Saturation of current element
     *  \param satJ           Saturation of neighbor element
     */
    void getFlux(DimVector& flux, const Intersection& intersection, const Scalar satI, const Scalar satJ) const
    {}

    /*! Constructs a ConvectivePart object
     *
     *  \param problem A problem class object
     */
    ConvectivePart(Problem& problem)
    {}

};
}

#endif
