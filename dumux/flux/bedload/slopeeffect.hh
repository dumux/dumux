// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadFlux
 * \brief Formulas concidering the slope effect
 *
 */
#ifndef DUMUX_FLUX_SLOPE_EFFECT_HH
#define DUMUX_FLUX_SLOPE_EFFECT_HH

namespace Dumux {

/*!
 * \ingroup BedloadFlux
 * \brief Calculate the deviation of the bedload discharge from the main flow due to the bed slope.
 *
 *  The deviation \phi_b of the bedload discharge from the main flow velocity is described by the following formula
 *  (see e.g. Parker et al. (1983) "Bedload and Size Distribution in Paved Gravel-Bed Streams" or
 *  Sekine and Kikkawa (1992) "Mechanics of saltating grains")
 *
 * \f[
 * \phi_b = \arctan( - f(\theta)\mathbf{s} \mathbf{n_q})
 * \f]
 *
 * with:
 *
 * \f$\theta\f$: dimensionless bed shear stress \f$[-]\f$\n
 * \f$\mathbf{s}\f$: bed slope vector \f$[-]\f$\n
 * \f$\mathbf{n_q}\f$: flow unit normal \f$[-]\f$\n
 *
 *  Talmon et al. (1995) "Laboratory measurements of the direction of sediment transport on transverse alluvial-bed slopes"
 *  defines the term  \f$f(\theta)\f$ as:
 *
 * \f[
 * f(\theta) = \frac{1}{a \sqrt{\theta}}
 * \f]
 *
 * For laboratory conditions he found good results for \f$a = 1.7 \f$,
 * while for natural rivers fair results were obtained with \f$a = 0.85 \f$.
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \param bedSlope bed slope vector at the edge
 * \param bedShearStress bed shear stress vector
 * \param velocity flow velocity vector
 *
 * \return The angle by which the bedload discharge deviates from the main flow velocity.
 */
template<class Scalar>
const Scalar calculateDirectionCorrection(const Dune::FieldVector<Scalar, 2>& bedSlope, const Dune::FieldVector<Scalar, 2>& bedShearStress, const Dune::FieldVector<Scalar, 2>& velocity, const Scalar representativeGrainDiameter)
{
    static const std::string formula = getParam<std::string>("Sediment.SlopeEffectDirectionCorrection", "None");
    // Use the approach of Talmon et al. (1995)
    if (formula == "Talmon")
    {
        Scalar eps = 1e-12;
        static const Scalar gravity = getParam<Scalar>("Problem.Gravity", 9.81);
        static const Scalar waterDensity = getParam<Scalar>("Sediment.WaterDensity", 1000);
        static const Scalar grainDensity = getParam<Scalar>("Sediment.GrainDensity", 2650.0);
        static const Scalar talmonParameter = getParam<Scalar>("Sediment.TalmonParameter");

        // The dimensionless shear stress is also called Shields parameter
        Scalar dimensionlessShearStress = bedShearStress.two_norm() / (gravity * representativeGrainDiameter * (grainDensity - waterDensity));

        // unit vector perpendicular to the flow direction
        Scalar scalarVelocity = velocity.two_norm();
        Dune::FieldVector<Scalar, 2> flowUnitNormal {-velocity[1]/scalarVelocity, velocity[0]/scalarVelocity};

        // check for small shear stresses
        if (dimensionlessShearStress < eps) {
            return 0.0;
        }
        else {
            return atan(-1/(talmonParameter*sqrt(dimensionlessShearStress)) * bedSlope * flowUnitNormal);
        }
    }
    else if (formula == "None") { return 0.0; }
    else
    {
        DUNE_THROW(Dune::InvalidStateException, "The parameter 'Sediment.SlopeEffectDirectionCorrection' is set to '"<<formula
                   <<"', which is not a valid value. Valid values 'Talmon' and 'None'!");
    }
}
/*!
 * \ingroup BedloadFlux
 * \brief Correct the magnitude of the bedload discharge due to the bed slope.
 *
 * Koch (1980) "bed level computations for axisymmetric curved channels" proposed an approach to modify the
 * magnitude of bedload discharge, which was rewritten by Struiksma et al. (1985) "Bed deformation in curved
 * alluvial channels" using a coefficient \f$c_{grav}\f$ by which the bedload discharge have to be multiplied to
 * consider the influence of the gravitational transport.
 *
 * \f[
 * c_{grav} = 1 - \zeta \frac{\partial{z}}{\partial{s}}
 * \f]
 *
 * with:
 *
 * \f$\zeta\f$: empirical factor (default=1.3) \f$[-]\f$\n
 * \f$\partial{z}/\partial{s}\f$: bed slope in flow direction \f$[-]\f$\n
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \param bedloadDischarge bedload discharge vector
 * \param bedSlope bed slope vector at the edge
 *
 */
template<class Scalar>
const void correctBedloadMagnitude(Dune::FieldVector<Scalar, 2>& bedloadDischarge, const Dune::FieldVector<Scalar, 2>& bedSlope)
{
    static const std::string formula = getParam<std::string>("Sediment.SlopeEffectMagnitudeCorrection", "None");
    Scalar eps = 1e-12;
    // Use the approach of Koch (1980)
    if (formula == "Koch")
    {
        if (bedloadDischarge.two_norm()>eps) {
        static const Scalar beta = getParam<Scalar>("Sediment.KochParameter",1.3);
        // Use for the correction of the magnitude only the part of bed slope, which is in direction of the bedload discharge.
        // Therefore the scalar product of the bed slope and the normalized bedload discharge is used for the correction.
        bedloadDischarge *= 1 - beta * bedSlope * bedloadDischarge / bedloadDischarge.two_norm();
        }
    }
    else if (formula == "None") { }
    else
    {
        DUNE_THROW(Dune::InvalidStateException, "The parameter 'Sediment.SlopeEffectMagnitudeCorrection' is set to '"<<formula
                   <<"', which is not a valid value. Valid values are 'Koch' and 'None'!");
    }
}
} // end namespace Dumux

#endif
