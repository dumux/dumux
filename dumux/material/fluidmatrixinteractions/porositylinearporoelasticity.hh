// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Fluidmatrixinteractions
 * \brief The porosity of a linear poroelastic medium.
 */
#ifndef DUMUX_MATERIAL_POROSITY_LINEAR_POROELASTICITY_HH
#define DUMUX_MATERIAL_POROSITY_LINEAR_POROELASTICITY_HH

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief The porosity of a linear poroelastic medium per unit volume of the undeformed medium.
 *
 * In linear poroelasticity \cite Coussy2004, the porosity relative to a reference state
 * with porosity \f$\phi_0\f$ is
 * \f[
 *     \phi = \phi_0 + \alpha \Delta\varepsilon_v + \frac{\Delta p}{N},
 *     \qquad
 *     \frac{1}{N} = \frac{\alpha - \phi_0}{K_s} = \frac{(\alpha - \phi_0)(1 - \alpha)}{K},
 * \f]
 * where \f$\Delta\varepsilon_v\f$ and \f$\Delta p\f$ are the changes of the volumetric strain
 * \f$\varepsilon_v = \nabla \cdot \mathbf{u}\f$ and of the pore pressure with respect to the
 * reference state, \f$\alpha = 1 - K/K_s\f$ is the Biot coefficient, \f$K\f$ is the drained bulk
 * modulus of the solid skeleton and \f$K_s\f$ is the bulk modulus of the solid grains.
 *
 * This porosity (Lagrangian porosity) is the pore volume per unit volume of the undeformed medium,
 * which is the porosity entering the storage term of a fluid mass balance on the undeformed grid.
 * With a pore fluid of bulk modulus \f$K_f\f$, this storage term corresponds to the Biot modulus
 * \f$M\f$ with \f$1/M = \phi_0/K_f + 1/N\f$.
 * For incompressible grains (\f$\alpha = 1\f$), the porosity is \f$\phi = \phi_0 + \Delta\varepsilon_v\f$.
 *
 * \tparam Scalar The type used for scalar values
 */
template<class Scalar>
class PorosityLinearPoroElasticity
{
public:
    /*!
     * \brief Constructor
     * \param refPorosity The porosity \f$\phi_0\f$ in the reference state
     * \param biotCoefficient The Biot coefficient \f$\alpha \in [\phi_0, 1]\f$
     * \param drainedBulkModulus The drained bulk modulus \f$K\f$ of the solid skeleton
     */
    PorosityLinearPoroElasticity(Scalar refPorosity, Scalar biotCoefficient, Scalar drainedBulkModulus)
    : refPorosity_(refPorosity)
    , biotCoefficient_(biotCoefficient)
    , inverseN_((biotCoefficient - refPorosity)*(1.0 - biotCoefficient)/drainedBulkModulus)
    {
        if (biotCoefficient < refPorosity || biotCoefficient > 1.0)
            DUNE_THROW(Dune::RangeError, "The Biot coefficient " << biotCoefficient
                        << " is not in [" << refPorosity << ", 1]");
    }

    /*!
     * \brief The porosity
     * \param volumetricStrainChange The change \f$\Delta\varepsilon_v\f$ of the volumetric strain with respect to the reference state
     * \param pressureChange The change \f$\Delta p\f$ of the pore pressure with respect to the reference state
     */
    Scalar evaluatePorosity(Scalar volumetricStrainChange, Scalar pressureChange) const
    { return refPorosity_ + biotCoefficient_*volumetricStrainChange + inverseN_*pressureChange; }

    /*!
     * \brief The Biot modulus \f$M\f$
     * \param fluidBulkModulus The bulk modulus \f$K_f\f$ of the pore fluid
     */
    Scalar biotModulus(Scalar fluidBulkModulus) const
    { return 1.0/(refPorosity_/fluidBulkModulus + inverseN_); }

private:
    Scalar refPorosity_;
    Scalar biotCoefficient_;
    Scalar inverseN_;
};

} // end namespace Dumux

#endif
