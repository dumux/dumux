// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
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
#ifndef DUMUX_DIFFUSIVEPART_HH
#define DUMUX_DIFFUSIVEPART_HH

#include <dumux/decoupled/2p/transport/transportproperties2p.hh>

/**
 * \file
 * \brief  Base class for defining a diffusive part of the saturation transport equation
 * \author Bernd Flemisch, Markus Wolff
 */
namespace Dumux
{
/*!\ingroup FVSaturation2p
 * \brief  Base class for defining the diffusive part of the saturation transport equation
 *
 * \tparam TypeTag The problem TypeTag
 */
template<class TypeTag>
class DiffusivePart
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    enum{dim = GridView::dimension};
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dim> FieldVector;

public:
    /*! \brief Returns diffusive term for current element face
     *
     *  \param flux        Flux vector (gets the flux from the function)
     *  \param intersection  Intersection of two grid elements/global boundary
     *  \param satI           saturation of current element
     *  \param satJ           saturation of neighbor element
     *  \param pcGradient     gradient of capillary pressure between element I and J
     */
    void getFlux(FieldVector& flux, const Intersection& intersection, Scalar satI, Scalar satJ, const FieldVector& pcGradient) const
    {}

    /*! \brief Returns diffusive term for current element face
     *
     *  \param flux        Flux vector (gets the flux from the function)
     *  \param intersection  Intersection of two grid elements/global boundary
     *  \param satIntersection  saturation at the face between element I and J
     *  \param satGradient       gradient of saturation between element I and J
     *  \param time             time
     */
    void getFlux(FieldVector& flux, const Intersection& intersection,
                                    const Scalar satIntersection, const FieldVector& satGradient, const Scalar time) const
    {}

    /*! \brief Returns diffusive term for current element face
     *
     *  \param flux        Flux vector (gets the flux from the function)
     *  \param intersection  Intersection of two grid elements/global boundary
     *  \param satIntersection  saturation at the face between element I and J
     *  \param satGradient       gradient of saturation between element I and J
     *  \param time             time
     *  \param satI             saturation of current element
     *  \param satJ             saturation of neighbor element
     */
    void getFlux(FieldVector& flux, const Intersection& intersection,
                                    const Scalar satIntersection, const FieldVector& satGradient, const Scalar time,
                                    Scalar satI, Scalar satJ) const
    {}


    //! Constructs a DiffusivePart object
    /*
     *  \param A problem class object
     */
    DiffusivePart(Problem& problem)
    {}

};
}

#endif
