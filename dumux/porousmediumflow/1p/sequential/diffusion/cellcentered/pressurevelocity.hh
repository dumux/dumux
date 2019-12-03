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
 * \ingroup SequentialOnePModel
 * \brief  Single Phase Finite Volume Model
 */

#ifndef DUMUX_FVPRESSUREVELOCITY1P_HH
#define DUMUX_FVPRESSUREVELOCITY1P_HH

// dumux environment
#include "pressure.hh"
#include <dumux/porousmediumflow/1p/sequential/properties.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/velocity.hh>


namespace Dumux {

/*!
 * \ingroup SequentialOnePModel
 * \brief Single Phase Finite Volume Model
 *
 * This model solves equations of the form
 * \f[
 *  \textbf{div}\, \boldsymbol v = q.
 * \f]
 * The velocity \f$ \boldsymbol v \f$ is the single phase Darcy velocity:
 * \f[
 *  \boldsymbol v = -\frac{1}{\mu} \boldsymbol K \left(\textbf{grad}\, p + \rho \, g  \, \textbf{grad}\, z\right),
 * \f]
 * where \f$ p \f$ is the pressure, \f$ \boldsymbol K \f$ the absolute permeability,
 *       \f$ \mu \f$ the viscosity, \f$ \rho \f$ the density, and \f$ g \f$ the gravity constant,
 * and \f$ q \f$ is the source term.
 * At the boundary, \f$ p = p_D \f$ on \f$ \Gamma_{Dirichlet} \f$, and \f$ \boldsymbol v \cdot \boldsymbol n = q_N\f$
 * on \f$ \Gamma_{Neumann} \f$.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FVPressureVelocity1P: public FVPressure1P<TypeTag>
{
    using ParentType = FVPressure1P<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
public:
    /*!
     * \brief Initializes the pressure model
     *
     * \copydetails FVPressure::initialize()
     */
    void initialize()
    {
        ParentType::initialize(false);
        velocity_.calculateVelocity();
    }

    /*!
     * \brief Pressure update
     *
     * \copydetails FVPressure::update()
     */
    void update()
    {
        ParentType::update();
        velocity_.calculateVelocity();
    }

    /*!
     * \brief Adds velocity output to the output file
     *
     * Adds the velocities to the output.
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        ParentType::addOutputVtkFields(writer);
        velocity_.addOutputVtkFields(writer);
    }

    /*!
     * Constructs a FVPressure1P object
     *
     * \param problem A problem class object
     */
    FVPressureVelocity1P(Problem& problem) :
        ParentType(problem), velocity_(problem)
    {}

private:
    FVVelocity<TypeTag, GetPropType<TypeTag, Properties::Velocity> > velocity_;
};

} // end namespace Dumux
#endif
