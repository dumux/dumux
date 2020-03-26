// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
/*!
 * \file
 * \ingroup SequentialTwoPModel
 * \brief  Base class for defining a convective part of the saturation transport equation.
 */
#ifndef DUMUX_CONVECTIVEPART_HH
#define DUMUX_CONVECTIVEPART_HH

#include <dumux/porousmediumflow/2p/sequential/properties.hh>

namespace Dumux {

/*!
 * \ingroup SequentialTwoPModel
 * \brief  Base class for defining a convective part of the saturation transport equation.
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class ConvectivePart
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
      using Scalar = GetPropType<TypeTag, Properties::Scalar>;
      using Problem = GetPropType<TypeTag, Properties::Problem>;

    enum{dimWorld = GridView::dimensionworld};
    using Intersection = typename GridView::Intersection;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;

public:
    //! For initialization
    void initialize()
    {}

    /*!
     * \brief Returns convective term for current element face
     *  \param intersection  Intersection of two grid elements/global boundary
     *  \param sat           Saturation of current element
     *  \return     Convective flux
     */
    Scalar getFlux(const Intersection& intersection, const Scalar sat) const
    {
        return 0.0;
    }

    /*!
     * \brief Returns convective term for current intersection
     *
     *  \param flux        Flux vector (gets the flux from the function)
     *  \param intersection  Intersection of two grid elements/global boundary
     *  \param satI           Saturation of current element
     *  \param satJ           Saturation of neighbor element
     */
    void getFlux(DimVector& flux, const Intersection& intersection, const Scalar satI, const Scalar satJ) const
    {}

    /*!
     * \brief Constructs a ConvectivePart object
     *
     *  \param problem A problem class object
     */
    ConvectivePart(Problem& problem)
    {}

};
} // end namespace Dumux

#endif
