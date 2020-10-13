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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesFluxVariablesImpl
 */
#ifndef DUMUX_NAVIERSTOKES_MASS_1PNC_FLUXVARIABLES_HH
#define DUMUX_NAVIERSTOKES_MASS_1PNC_FLUXVARIABLES_HH

#include <dumux/common/properties.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/flux/upwindscheme.hh>
#include <dumux/freeflow/navierstokes/scalarvolumevariables.hh>
#include <dumux/freeflow/navierstokes/scalarfluxvariables.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief The flux variables class for the Navier-Stokes model using the staggered grid discretization.
 */
template<class TypeTag, class UpwindScheme = UpwindScheme<GetPropType<TypeTag, Properties::GridGeometry>>>
class NavierStokesMassOnePNCFluxVariables
: public NavierStokesScalarConservationModelFluxVariables<TypeTag, UpwindScheme>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    static constexpr bool enableMolecularDiffusion = ModelTraits::enableMolecularDiffusion();

public:

    static constexpr auto numComponents = ModelTraits::numFluidComponents();
    static constexpr bool useMoles = ModelTraits::useMoles();
    using MolecularDiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>;

    /*!
     * \brief Returns the diffusive fluxes computed by the respective law.
     */
    Dune::FieldVector<Scalar, numComponents> molecularDiffusionFlux(int phaseIdx = 0) const
    {
        if constexpr (enableMolecularDiffusion)
            return MolecularDiffusionType::flux(this->problem(),
                                                this->element(),
                                                this->fvGeometry(),
                                                this->elemVolVars(),
                                                this->scvFace(),
                                                phaseIdx,
                                                this->elemFluxVarsCache());
        else
            return {0.0};
    }

};

} // end namespace Dumux

#endif // DUMUX_NAVIERSTOKES_STAGGERED_FLUXVARIABLES_HH
