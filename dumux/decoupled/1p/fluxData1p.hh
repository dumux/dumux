/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
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
#ifndef DUMUX_FLUXDATA1P_HH
#define DUMUX_FLUXDATA1P_HH

#include "1pproperties.hh"

/**
 * @file
 * @brief  Class including the variables and data of discretized data of the constitutive relations
 * @author Markus Wolff
 */

namespace Dumux
{
/*!
 * \ingroup IMPES
 */
//! Class including the variables and data of discretized data of the constitutive relations.
/*! The variables of two-phase flow, which are one pressure and one saturation are stored in this class.
 * Additionally, a velocity needed in the transport part of the decoupled two-phase flow is stored, as well as discretized data of constitutive relationships like
 * mobilities, fractional flow functions and capillary pressure. Thus, they have to be callculated just once in every time step or every iteration step.
 *
 * @tparam TypeTag The Type Tag
 1*/
template<class TypeTag>
class FluxData1P
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dim> FieldVector;
    typedef Dune::FieldVector<FieldVector, 2 * dim> VelocityVector;

    VelocityVector velocity_;
    Scalar potential_[2 * dim];
    bool velocityMarker_[2 * dim];

public:

    //! Constructs a FluxData object
    /**
     *
     */

    FluxData1P()
    {
        for (int face = 0;  face < 2*dim; face++)
        {
            velocity_[face] = FieldVector(0.0);
            potential_[face] = 0.0;
            velocityMarker_[face] = false;
        }
    }

    ////////////////////////////////////////////////////////////
    // functions returning the vectors of the primary variables
    ////////////////////////////////////////////////////////////

    //! Return the velocity
    const FieldVector& velocity(int indexInInside)
    {
        return velocity_[indexInInside];
    }

    const FieldVector& velocity(int indexInInside) const
    {
        return velocity_[indexInInside];
    }

    void setVelocity(int indexInInside, FieldVector& velocity)
    {
        velocity_[indexInInside] = velocity;
    }

    void resetVelocity()
    {
        for (int i = 0; i < 2 * dim; i++)
        {
            velocity_[i] = 0.;
            velocityMarker_[i] = false;
        }
    }

    void setVelocityMarker(int indexInInside)
    {
        velocityMarker_[indexInInside] = true;
    }

    bool haveVelocity(int indexInInside)
    {
        return velocityMarker_[indexInInside];
    }

    void resetVelocityMarker()
    {
        for (int i = 0; i < 2*dim; i++)
            velocityMarker_[i] = false;
    }

    bool isUpwindCell(int indexInInside)
    {
        return (potential_[indexInInside] >= 0.);
    }

    bool isUpwindCell(int indexInInside) const
    {
        return (potential_[indexInInside] >= 0.);
    }

    Scalar potential(int indexInInside)
    {
        return potential_[indexInInside];
    }

    Scalar potential(int indexInInside) const
    {
        return potential_[indexInInside];
    }

    void setPotential(int indexInInside, Scalar pot)
    {
        potential_[indexInInside] = pot;
    }

};
}
#endif
