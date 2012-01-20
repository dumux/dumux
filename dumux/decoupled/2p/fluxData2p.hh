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
#ifndef DUMUX_FLUXDATA2P_HH
#define DUMUX_FLUXDATA2P_HH

#include "2pproperties.hh"

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
class FluxData2P
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

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef Dune::FieldVector<Scalar, dim> FieldVector;
    typedef Dune::FieldVector<FieldVector, 2 * dim> VelocityVector;

    VelocityVector velocity_[numPhases];
    Scalar potential_[2 * dim][numPhases];
    bool velocityMarker_[2 * dim];

public:

    //! Constructs a VariableClass object
    /**
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */

    FluxData2P()
    {
        for (int face = 0;  face < 2*dim; face++)
        {
            for (int phase = 0; phase < numPhases; phase++)
            {
                velocity_[phase][face] = FieldVector(0.0);

                potential_[face][phase] = 0.0;
            }
            velocityMarker_[face] = false;
        }
    }

    ////////////////////////////////////////////////////////////
    // functions returning the vectors of the primary variables
    ////////////////////////////////////////////////////////////

    //! Return the velocity
    const FieldVector& velocity(int phaseIdx, int indexInInside)
    {
        return velocity_[phaseIdx][indexInInside];
    }

    const FieldVector& velocity(int phaseIdx, int indexInInside) const
    {
        return velocity_[phaseIdx][indexInInside];
    }

    void setVelocity(int phaseIdx, int indexInInside, FieldVector& velocity)
    {
        velocity_[phaseIdx][indexInInside] = velocity;
    }

    void resetVelocity()
    {
        for (int i = 0; i < 2 * dim; i++)
        {
            for (int j = 0; j < numPhases; j++)
            {
                velocity_[j][i] = 0.;
            }
            velocityMarker_[i] = false;
        }
    }

    //! Return the velocity
    FieldVector velocityTotal(int indexInInside)
    {
        return velocity_[wPhaseIdx][indexInInside]
                + velocity_[nPhaseIdx][indexInInside];
    }

    FieldVector velocityTotal(int indexInInside) const
    {
        return velocity_[wPhaseIdx][indexInInside]
                + velocity_[nPhaseIdx][indexInInside];
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

    bool isUpwindCell(int phaseIdx, int indexInInside)
    {
        return (potential_[indexInInside][phaseIdx] >= 0.);
    }

    bool isUpwindCell(int phaseIdx, int indexInInside) const
    {
        return (potential_[indexInInside][phaseIdx] >= 0.);
    }

    Scalar potential(int phaseIdx, int indexInInside)
    {
        return potential_[indexInInside][phaseIdx];
    }

    Scalar potential(int phaseIdx, int indexInInside) const
    {
        return potential_[indexInInside][phaseIdx];
    }

    void setPotential(int phaseIdx, int indexInInside, Scalar pot)
    {
        potential_[indexInInside][phaseIdx] = pot;
    }

};
}
#endif
