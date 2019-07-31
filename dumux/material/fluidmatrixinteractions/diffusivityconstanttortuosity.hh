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
 * \ingroup Fluidmatrixinteractions
 * \brief   Relation for the saturation-dependent effective diffusion coefficient
 */
#ifndef DUMUX_MATERIAL_DIFFUSIVITY_CONSTANT_TORTUOSITY_HH
#define DUMUX_MATERIAL_DIFFUSIVITY_CONSTANT_TORTUOSITY_HH

#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Relation for the saturation-dependent effective diffusion coefficient
 *
 * The material law is:
 * \f[
 *  D_\text{eff,pm} = \phi * S_w * \tau * D
 * \f]
 *
 * with a constant tau.
 *
 * The default value is 0.5, empirically obtained in Carman 1937:
 * <i>Fluid flow through granular beds</i> \cite carman1937
 * Additionally, Bear 1972 \cite bear1972 mentions values 0.4 and in the
 * range of 0.56 to 0.8.
 */
template<class Scalar>
class DiffusivityConstantTortuosity
{
public:
    /*!
     * \brief Returns the effective diffusion coefficient \f$\mathrm{[m^2/s]}\f$ based
     *        on a constant tortuosity value
     * \param porosity The porosity
     * \param saturation The saturation of the wetting phase
     * \param diffCoeff The diffusion coefficient of the phase in \f$\mathrm{[m^2/s]}\f$
     */
    [[deprecated("Signature deprecated. Use signature with volume variables!")]]
    static Scalar effectiveDiffusivity(const Scalar porosity,
                                       const Scalar saturation,
                                       const Scalar diffCoeff)
    {
        static const Scalar tau = getParam<Scalar>("SpatialParams.Tortuosity", 0.5);

        return porosity * saturation * tau * diffCoeff;
    }

    /*!
     * \brief Returns the effective diffusion coefficient \f$\mathrm{[m^2/s]}\f$ based
     *        on a constant tortuosity value
     * \param volVars The Volume Variables
     * \param phaseIdx the index of the phase
     * \param compIdx the component index
     */
    template<class VolumeVariables>
    static Scalar effectiveDiffusivity(const VolumeVariables& volVars,
                                       const Scalar diffCoeff,
                                       const int phaseIdx)
    {
        static const Scalar tau = getParam<Scalar>("SpatialParams.Tortuosity", 0.5);

        return volVars.porosity() * volVars.saturation(phaseIdx) * tau * diffCoeff;
    }
};

} // end namespace Dumux

#endif // DIFFUSIVITY_CONSTANT_TORTUOSITY_HH
