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
 * \brief  Base class for defining a diffusive part of the saturation transport equation
 */
#ifndef DUMUX_DIFFUSIVEPART_HH
#define DUMUX_DIFFUSIVEPART_HH

#include <dumux/porousmediumflow/2p/sequential/transport/properties.hh>

namespace Dumux {
/*!
 * \ingroup SequentialTwoPModel
 * \brief Base class for defining the diffusive part of the saturation transport equation
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class DiffusivePart
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
      using Scalar = GetPropType<TypeTag, Properties::Scalar>;
      using Problem = GetPropType<TypeTag, Properties::Problem>;

    enum{dim = GridView::dimension};
    using Intersection = typename GridView::Intersection;
    using DimVector = Dune::FieldVector<Scalar, dim>;

public:

    //! For initialization
    void initialize()
    {}

    /*!
     * \brief Returns diffusive term for current element face
     *
     *  \param flux        Flux vector (gets the flux from the function)
     *  \param intersection  Intersection of two grid elements/global boundary
     *  \param satI           saturation of current element
     *  \param satJ           saturation of neighbor element
     *  \param pcGradient     gradient of capillary pressure between element I and J
     */
    void getFlux(DimVector& flux, const Intersection& intersection, Scalar satI, Scalar satJ, const DimVector& pcGradient) const
    {}

    /*!
     * \brief Returns diffusive term for current element face
     *
     *  \param flux        Flux vector (gets the flux from the function)
     *  \param intersection  Intersection of two grid elements/global boundary
     *  \param satIntersection  saturation at the face between element I and J
     *  \param satGradient       gradient of saturation between element I and J
     *  \param time             time
     */
    void getFlux(DimVector& flux, const Intersection& intersection,
                                    const Scalar satIntersection, const DimVector& satGradient, const Scalar time) const
    {}

    /*!
     * \brief Returns diffusive term for current element face
     *
     *  \param flux        Flux vector (gets the flux from the function)
     *  \param intersection  Intersection of two grid elements/global boundary
     *  \param satIntersection  saturation at the face between element I and J
     *  \param satGradient       gradient of saturation between element I and J
     *  \param time             time
     *  \param satI             saturation of current element
     *  \param satJ             saturation of neighbor element
     */
    void getFlux(DimVector& flux, const Intersection& intersection,
                                    const Scalar satIntersection, const DimVector& satGradient, const Scalar time,
                                    Scalar satI, Scalar satJ) const
    {}

    /*!
     * \brief Constructs a DiffusivePart object
     *
     *  \param problem problem class object
     */
    DiffusivePart(Problem& problem)
    {}

};
} // end namespace Dumux

#endif
