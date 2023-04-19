// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and a "constant" component.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_CONSTANT_HH
#define DUMUX_BINARY_COEFF_H2O_CONSTANT_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/constant.hh>

namespace Dumux {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and another component.
 * \todo All other binary coefficient could be generalized like this
 */
template<class Scalar, class Component>
class H2O_Component
{
    H2O_Component()
    {
        DUNE_THROW(Dune::NotImplemented, "The binary coefficients for H2O and your "
                     << "component are not implemented! Please implement the needed specialization.");
    }
};

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and a constant component
 */
template<class Scalar, int id>
class H2O_Component<Scalar, Components::Constant<id, Scalar>>
{
public:
    /*!
     * \brief Henry coefficient \f$N/m^2\f$  for the constant component in liquid water.
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     */
    static Scalar henryCompInWater(Scalar temperature)
    {
        static const Scalar h = getParamFromGroup<Scalar>(std::to_string(id), "Component.HenryComponentInWater", 1.0);
        return h;
    }

    /*!
     * \brief Henry coefficient \f$N/m^2\f$  for water in the constant component.
     *
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     */
    static Scalar henryWaterInComp(Scalar temperature)
    {
        static const Scalar h = getParamFromGroup<Scalar>(std::to_string(id), "Component.HenryWaterInComponent", 1.0);
        return h;
    }


    /*!
     * \brief Binary diffusion coefficient \f$m^2/s\f$ for molecular water and the constant component.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        static const Scalar D = getParamFromGroup<Scalar>(std::to_string(id), "Component.GasDiffusionCoefficient", 1.0);
        return D;
    }

    /*!
     * \brief Diffusion coefficient \f$m^2/s\f$ for the constant component in liquid water.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        static const Scalar D = getParamFromGroup<Scalar>(std::to_string(id), "Component.LiquidDiffusionCoefficient", 1.0);
        return D;
    }
};

} // end namespace BinaryCoeff
} // end namespace Dumux

#endif
