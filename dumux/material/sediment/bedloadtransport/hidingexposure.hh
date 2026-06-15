// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTransport
 * \brief Formulas concidering the influence of hiding and exposure
 *
 */
#ifndef DUMUX_MATERIAL_HIDING_EXPOSURE_HH
#define DUMUX_MATERIAL_HIDING_EXPOSURE_HH

#include <vector>

namespace Dumux {

/*!
 * \ingroup BedloadFlux
 * \brief Calculate the hiding and exposure coefficients (Karim, Holly Yang).
 *
 * We use the approach proposed by Karim, Holly and Yang (1987) in <em> "Ialluvial Numerical Simulation of Mobile- Bed Rivers:
 * Part I. Theoretical and Numerical Principles." </em>.
 * This approach yields a coefficient \f$c_{hide}\f$ for each grain class, by which the bedload transport rate of the corresponding grain class has to be multiplied.
 *
 *  \f[
 *      c_{hide} = c_1 \left(\frac{d}{d_{50}}\right)^{c_2}
 * \f]
 *
 * With
 *
 * \f$c_1\f$: bed-material dependend parameter\f$[-]\f$ \n
 * \f$c_2\f$: bed-material dependend parameter\f$[-]\f$ \n
 * \f$d\f$: representative grain diameter \f$[m]\f$ \n
 * \f$d_{50}\f$: median grain grain diameter \f$[m]\f$
 *
 * Karim et al. propose values of \f$c_1 = 1.0\f$ and \f$c_1\f$ = 0.8 for the typical Missouri River bed-material they used.
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \param representativeGrainDiameters representative grain diameters of all grain classes
 * \param medianGrainDiameters median grain diameter
 *
 * \return hiding and exposure coefficients
 */
template<typename Scalar>
const std::vector<Scalar> calculateHidingExposureCoefficientsKarimHollyYang(const std::vector<Scalar>& representativeGrainDiameters, const Scalar& medianGrainDiameter)
{
    using std::pow;
    static const Scalar c1 = getParam<Scalar>("Sediment.KarimHollyYangParameter1",1.0);
    static const Scalar c2 = getParam<Scalar>("Sediment.KarimHollyYangParameter2",0.8);
    std::vector<Scalar> hidingexposureCoefficients {};
    for (Scalar diameter : representativeGrainDiameters) {
        hidingexposureCoefficients.push_back(c1 * pow(diameter / medianGrainDiameter, c2));
    }
    return hidingexposureCoefficients;
}
} // end namespace Dumux

#endif
