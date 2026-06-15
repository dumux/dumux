// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTransport
 * \copydoc Dumux::BedloadFormulaMeyerPeterMueller
 */
#ifndef DUMUX_MATERIAL_BEDLOADFORMULA_MEYERPETERMUELLER_HH
#define DUMUX_MATERIAL_BEDLOADFORMULA_MEYERPETERMUELLER_HH

#include "bedloadformula.hh"
#include <dumux/flux/bedload/extendvector.hh>

namespace Dumux {
/*!
 * \ingroup BedloadFormulas
 * \brief Implementation of Meyer-Peter and Muellers bedload transport formula.
 *
 * This is a common formula to calculate the bedload transport rate \f$q_b\f$.
 *
 * \f[
 * q_b = 8 * \sqrt{g * (\rho_s / \rho_w - 1) * d_{50}^3} * (max(0, \theta - \theta_{c}))^{1,5}
 * \f]
 *
 * with:
 *
 * \f$g\f$: gravity \f$[m/s^2]\f$\n
 * \f$\rho_s\f$: bedload grain density \f$[kg/m^3]\f$\n
 * \f$\rho_w\f$: water density \f$[kg/m^3]\f$\n
 * \f$d_{50}\f$: a representative grain diameter \f$[m]\f$\n
 * \f$\theta\f$: dimensionless bed shear stress \f$[-]\f$\n
 * \f$\theta_{c}\f$: critical Shields stress \f$[-]\f$
 */

template <typename VolumeVariables>
class BedloadFormulaMeyerPeterMueller : public BedloadFormula<VolumeVariables>
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
public:
    /*!
     * \brief Constructor
     *
     * \param criticalShieldsStress critical Shields stress, also known as critical non-dimensional bed shear stress.
     *        Beneath this value no bedload transport occurs \f$[-]\f$.
     * \param grainDensity Grain density of the transported material \f$[kg/m^3]\f$
     * \param gravity Gravity \f$[m/s^2]\f$
     * \param mpmCoefficient Meyer-Peter Mueller coefficient \f$[-]\f$
     * \param representativeGrainDiameter Representative grain diameter \f$[m]\f$
     * \param waterDensity Water density \f$[kg/m^3]\f$
     * \param hidingExposureCoefficientKarimHollyYang Hiding and exposure coefficient of Karim, Holly and Yang \f$[-]\f$
     */
    BedloadFormulaMeyerPeterMueller(const Scalar criticalShieldsStress,
                                    const Scalar grainDensity,
                                    const Scalar gravity,
                                    const Scalar mpmCoefficient,
                                    const Scalar representativeGrainDiameter,
                                    const Scalar waterDensity,
                                    const Scalar hidingExposureCoefficientKarimHollyYang=1.0)
        : criticalShieldsStress_(criticalShieldsStress), grainDensity_(grainDensity), gravity_(gravity),
          hidingExposureCoefficientKarimHollyYang_(hidingExposureCoefficientKarimHollyYang),
          meyerPeterMuellerCoefficient_(mpmCoefficient), representativeGrainDiameter_(representativeGrainDiameter),
          waterDensity_(waterDensity)
    {
        specificDensity_ = grainDensity_  / waterDensity_;
    }

    /*!
     * \brief Compute the dimensionless bottom shear stress.
     *
     * \param volVars Volume Variables.
     *
     * Compute the dimensionless bottom shear stress in x- and y-direction.
     *
     * \return dimensionless bottom shear stress \f$[-]\f$. First entry is the x-component, the second the y-component.
     */
    Dune::FieldVector<Scalar, 2> dimensionlessBottomShearStress(const VolumeVariables& volVars) const
    {

        // Consider the influence of the secondary currents on the bottom shear stress.
        Dune::FieldVector<Scalar, 2> bottomShearStress = Sediment::getExtendedVector(volVars.bottomShearStress(), volVars.tanAngleSecondaryCurrents());

        return bottomShearStress / (gravity_ * representativeGrainDiameter_ * (grainDensity_ - waterDensity_));
    }

    /*!
     * \brief Compute the bedload transport rate.
     *
     * \param volVars Volume Variables.
     *
     * Compute the bedload transport rate in x- and y-direction.
     *
     * \return Bedload transport rate \f$[m^3/s/m]\f$. First entry is the x-component, the second the y-component.
     */
    Dune::FieldVector<Scalar, 2> bedloadDischarge(const VolumeVariables& volVars) const final
    {
        using std::sqrt;
        using std::pow;
        using std::max;

        Dune::FieldVector<Scalar, 2> bedloadDischarge(0.0);

        Dune::FieldVector<Scalar, 2> dimensionlessBottomShearStress = this->dimensionlessBottomShearStress(volVars);
        Scalar scalarDimensionlessBottomShearStress = dimensionlessBottomShearStress.two_norm();

        Scalar scalarBedloadDischarge = meyerPeterMuellerCoefficient_ * sqrt(gravity_ * (specificDensity_ - 1.0) * pow(representativeGrainDiameter_, 3))
                                        * pow(max(scalarDimensionlessBottomShearStress - criticalShieldsStress_, 0.0), 1.5);

        if (scalarDimensionlessBottomShearStress < eps_) { // check for small bottom shear stress to avoid divison by zero
            return bedloadDischarge;          // both elements of bedloadDischarge are zero
        }
        else {
            bedloadDischarge[0] = scalarBedloadDischarge * dimensionlessBottomShearStress[0]/scalarDimensionlessBottomShearStress;
            bedloadDischarge[1] = scalarBedloadDischarge * dimensionlessBottomShearStress[1]/scalarDimensionlessBottomShearStress;
        }

        // Consider the influence of hiding & exposure using the approach of Karim, Holly and Yang (1987)
        bedloadDischarge *= hidingExposureCoefficientKarimHollyYang_;

        return bedloadDischarge;
    }

private:
    Scalar criticalShieldsStress_;
    static constexpr Scalar eps_ = 1e-6;
    Scalar grainDensity_;
    Scalar gravity_;
    Scalar hidingExposureCoefficientKarimHollyYang_;
    Scalar meyerPeterMuellerCoefficient_;
    Scalar representativeGrainDiameter_;
    Scalar specificDensity_;
    Scalar waterDensity_;
};

} // end namespace Dumux

#endif // DUMUX_MATERIAL_FLUIDMATRIX_BEDLOADFORMULA_MEYERPETERMUELLER_HH
