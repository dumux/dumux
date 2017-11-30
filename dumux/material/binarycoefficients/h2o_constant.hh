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
/*!
 * \file
 *
 * \brief Binary coefficients for water and a "constant" component.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_CONSTANT_HH
#define DUMUX_BINARY_COEFF_H2O_CONSTANT_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/basicproperties.hh>

#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/constant.hh>

namespace Dumux
{

namespace BinaryCoeff
{

/*!
 * \brief Binary coefficients for water and another component.
 * \todo All other binary coefficient could be generalized like this
 */
template<class TypeTag, class Component>
class H2O_Component
{
    H2O_Component()
    {
        DUNE_THROW(Dune::NotImplemented, "The binary coefficients for H2O and your "
                     << "component are not implemented! Please implement the needed specialization.");
    }
};

/*!
 * \brief Binary coefficients for water and a constant component
 */
template<class TypeTag, int id>
class H2O_Component<TypeTag, Components::Constant<id, typename GET_PROP_TYPE(TypeTag, Scalar)>>
{
public:
    /*!
     * \brief Henry coefficent \f$N/m^2\f$  for heavy oil in liquid water.
     *
     * See:
     *
     */

    template <class Scalar>
    static Scalar henryCompInWater(Scalar temperature)
    {
        static const Scalar h = getParamFromGroup<Scalar>(std::to_string(id), "Component.HenryComponentInWater", 1.0);
        return h;
    }

    /*!
     * \brief Henry coefficent \f$N/m^2\f$  for water in liquid heavy oil.
     *
     * See:
     *
     */

    template <class Scalar>
    static Scalar henryWaterInComp(Scalar temperature)
    {
        static const Scalar h = getParamFromGroup<Scalar>(std::to_string(id), "Component.HenryWaterInComponent", 1.0);
        return h;
    }


    /*!
     * \brief Binary diffusion coefficent \f$m^2/s\f$ for molecular water and heavy oil.
     *
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        static const Scalar D = getParamFromGroup<Scalar>(std::to_string(id), "Component.GasDiffusionCoefficient", 1.0);
        return D;
    }

    /*!
     * \brief Diffusion coefficent \f$m^2/s\f$ for tce in liquid water.
     *
     * \todo
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        static const Scalar D = getParamFromGroup<Scalar>(std::to_string(id), "Component.LiquidDiffusionCoefficient", 1.0);
        return D;
    }
};

} // end namespace BinaryCoeff

} // end namespace Dumux

#endif
