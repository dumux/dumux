// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTransport
 * \copydoc Dumux::BedloadFormulaGrass
 */
#ifndef DUMUX_MATERIAL_BEDLOADFORMULA_GRASS_HH
#define DUMUX_MATERIAL_BEDLOADFORMULA_GRASS_HH

#include "bedloadformula.hh"
#include <dumux/flux/bedload/extendvector.hh>

namespace Dumux {
/*!
 * \ingroup BedloadFormulas
 * \brief Implementation of the Grass bedload transport formula.
 *
 * This is one of the simplest formulas to calculate the bedload transport rate \f$q_b\f$.
 * It just depends on the velocity (\f$u\f$ and \f$v\f$) and a constant coefficient \f$A_g\f$ which
 * depends on the properties of the sediment. In x-direction:
 *
 * \f[
 * q_b,x = A_g * u * (u^2 + v^2)
 * \f]
 *
 * Analogous in y-direction. This is a simplified version of the general
 * Grass formula \f$ \mathbf{q_b} = A_g*\mathbf{u}*|\mathbf{u}|^{(m_g-1)}\f$ with \f$m_g = 3\f$.
 */

template <typename VolumeVariables>
class BedloadFormulaGrass : public BedloadFormula<VolumeVariables>
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
public:
    /*!
     * \brief Constructor
     *
     * \param grassAlpha Grass bedload transport coefficient [-]
     * \param hidingExposureCoefficientKarimHollyYang Hiding and exposure coefficient of Karim, Holly and Yang [-]
     */
    BedloadFormulaGrass(const Scalar grassAlpha, const Scalar hidingExposureCoefficientKarimHollyYang=1.0)
        : grassAlpha_(grassAlpha), hidingExposureCoefficientKarimHollyYang_(hidingExposureCoefficientKarimHollyYang) {}

    /*!
     * \brief Compute the bedload transport rate.
     *
     * \param volVars Volume Variables.
     *
     * Compute the bedload transport rate in x- and y-direction.
     *
     * \return bedload transport rate \f$[m^3/s/m]\f$. First entry is the x-component, the second the y-component.
     */
    Dune::FieldVector<Scalar, 2> bedloadDischarge(const VolumeVariables& volVars) const final
    {
        Dune::FieldVector<Scalar, 2> bedloadDischarge(0.0);

        // Consider the influence of the secondary currents on the bedload discharge. Actually, the secondary currents
        // influence the bed shear stress, but since the Grass formula does not use the bed shear stress we modify the
        // velocity by the secondary currents angle.
        Dune::FieldVector<Scalar, 2> velocity = Sediment::getExtendedVector(volVars.velocity(), volVars.tanAngleSecondaryCurrents());

        for (int i=0; i<2; i++)
        {
            bedloadDischarge[i] = grassAlpha_ * velocity[i] * (velocity[0] * velocity[0] + velocity[1] * velocity[1]);
        }

        // Consider the influence of hiding & exposure using the approach of Karim, Holly and Yang (1987)
        bedloadDischarge *= hidingExposureCoefficientKarimHollyYang_;

        return bedloadDischarge;
    }

private:
    Scalar grassAlpha_;
    Scalar hidingExposureCoefficientKarimHollyYang_;
};

} // end namespace Dumux

#endif // DUMUX_MATERIAL_BEDLOADFORMULA_GRASS_HH
