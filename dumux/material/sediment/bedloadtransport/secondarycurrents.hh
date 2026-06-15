// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTransport
 * \brief Formulas concidering the influence of secondary currents
 *
 */
#ifndef DUMUX_MATERIAL_SECONDARY_CURRENTS_HH
#define DUMUX_MATERIAL_SECONDARY_CURRENTS_HH

namespace Dumux {

/*!
 * \ingroup BedloadFlux
 * \brief Calculate the deviation of the bed shear stress from the main flow due to secondary currents.
 *
 * We use the approach proposed by Engelund (1974) in <em> "Flow and bed topography in channel bend" </em> to quantify this influence.
 * This approach yields the tanget of the angle \f$\delta\f$ by which the bed shear stress deviates from the main flow velocity.
 *
 *  \f[
 * tan(\delta) = -N * h / R
 * \f]
 *
 * With
 *
 * \f$N\f$: curvature factor \f$[-]\f$ \n
 * \f$h\f$: water depth \f$[m]\f$ \n
 * \f$R\f$: curvature radius \f$[m]\f$
 *
 * Engelund proposes a fixed curvature factor of 7.
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \param waterDepth water depth
 * \param curvatureRadius curvature radius. A positive radius means a counterclockwise curvature, a negative radius a clockwise curvature.
 *
 * \return The tangent of the angle by which the bed shear stress deviates from the main flow velocity. The angle is defined from the main flow velocity to the bed shear shress.
 */
template<class Scalar>
const Scalar calculateTanAngleSecondaryCurrents(const Scalar waterDepth, const Scalar curvatureRadius)
{
    const Scalar eps = 1e-8;
    static const std::string formula = getParam<std::string>("Sediment.SecondaryCurrents", "None");
    if (formula=="Engelund") {
        if (abs(curvatureRadius) > eps) {
            static const Scalar curvatureFactor = getParam<Scalar>("Sediment.CurvatureFactor");
            return curvatureFactor * waterDepth / curvatureRadius;
        }
        else {
            return 0.0;
        }
    }
    else if (formula == "None") {
        return 0.0;
    }
    else {
        DUNE_THROW(Dune::InvalidStateException, "The parameter 'Sediment.SecondaryCurrents' is set to '"<<formula
                   <<"', which is not a valid value. Valid values 'Engelund' and 'None'!");
    }
}
} // end namespace Dumux

#endif
