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
#ifndef DUMUX_FVPRESSUREVELOCITY2P_HH
#define DUMUX_FVPRESSUREVELOCITY2P_HH


// dumux environment
#include "pressure.hh"
#include <dumux/porousmediumflow/2p/sequential/properties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/velocity.hh>

/**
 * \file
 * \brief  Two Phase Finite Volume Model
 */

namespace Dumux
{

//! \ingroup FV2p
//! \brief Two Phase Finite Volume Model
/*!
 *
 * \tparam TypeTag The Type Tag
 *
 */
template<class TypeTag> class FVPressureVelocity2P: public FVPressure2P<TypeTag>
{
    typedef FVPressure2P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
public:
    /*! \brief Initializes the pressure model
     *
     * \copydetails ParentType::initialize()
     */
    void initialize()
    {
        ParentType::initialize();
        velocity_.initialize();
        velocity_.calculateVelocity();
    }

    /*! \brief Pressure update
     *
     * \copydetails ParentType::update()
     *
     */
    void update()
    {
        ParentType::update();
        velocity_.calculateVelocity();
    }

    /*! \brief Adds velocity output to the output file
     *
     * Adds the velocities to the output.
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     *
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        ParentType::addOutputVtkFields(writer);
        velocity_.addOutputVtkFields(writer);
    }

    //! Constructs a FVPressure2P object
    /**
     * \param problem A problem class object
     */
    FVPressureVelocity2P(Problem& problem) :
        ParentType(problem), velocity_(problem)
    {}

private:
    Dumux::FVVelocity<TypeTag, typename GET_PROP_TYPE(TypeTag, Velocity) > velocity_;
};

}
#endif
