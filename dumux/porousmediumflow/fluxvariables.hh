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
 * \ingroup PorousmediumflowModels
 * \brief Base class for the flux variables in porous medium models
 */

#ifndef DUMUX_POROUSMEDIUMFLOW_FLUXVARIABLES_HH
#define DUMUX_POROUSMEDIUMFLOW_FLUXVARIABLES_HH

#include <bitset>
#include <array>

#include <dumux/common/properties.hh>
#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/flux/upwindscheme.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumflowModels
 * \brief The porous medium flux variables class that computes advective / convective,
 *        molecular diffusive and heat conduction fluxes.
 *
 * \param TypeTag The type tag for access to type traits
 * \param UpwindScheme The upwind scheme to be applied to advective fluxes
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag,
         class UpwindScheme = UpwindScheme<GetPropType<TypeTag, Properties::GridGeometry>> >
class PorousMediumFluxVariables
: public FluxVariablesBase<GetPropType<TypeTag, Properties::Problem>,
                           typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView,
                           typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    enum
    {
        numPhases = ModelTraits::numFluidPhases(),
        numComponents = ModelTraits::numFluidComponents()
    };

public:
    using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
    using MolecularDiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using HeatConductionType = GetPropType<TypeTag, Properties::HeatConductionType>;

    static constexpr bool enableAdvection = ModelTraits::enableAdvection();
    static constexpr bool enableMolecularDiffusion = ModelTraits::enableMolecularDiffusion();
    static constexpr bool enableEnergyBalance = ModelTraits::enableEnergyBalance();
    static constexpr bool enableThermalNonEquilibrium = getPropValue<TypeTag, Properties::EnableThermalNonEquilibrium>();

    //! The constructor
    PorousMediumFluxVariables()
    {
        advFluxIsCached_.reset();
        advFluxBeforeUpwinding_.fill(0.0);
    }

    /*!
     * \brief Returns the advective flux computed by the respective law.
     */
    template<typename FunctionType>
    Scalar advectiveFlux([[maybe_unused]] const int phaseIdx, [[maybe_unused]] const FunctionType& upwindTerm) const
    {
        if constexpr (enableAdvection)
        {
            if (!advFluxIsCached_[phaseIdx])
            {

                advFluxBeforeUpwinding_[phaseIdx] = AdvectionType::flux(this->problem(),
                                                                        this->element(),
                                                                        this->fvGeometry(),
                                                                        this->elemVolVars(),
                                                                        this->scvFace(),
                                                                        phaseIdx,
                                                                        this->elemFluxVarsCache());
                advFluxIsCached_.set(phaseIdx, true);
            }

            //! Give the upwind scheme access to the cached variables
            return UpwindScheme::apply(*this, upwindTerm, advFluxBeforeUpwinding_[phaseIdx], phaseIdx);
        }
        else
            return 0.0;
    }

    /*!
     * \brief Returns the diffusive fluxes computed by the respective law.
     */
    Dune::FieldVector<Scalar, numComponents> molecularDiffusionFlux([[maybe_unused]] const int phaseIdx) const
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
            return Dune::FieldVector<Scalar, numComponents>(0.0);
    }

    /*!
     * \brief Returns the conductive flux computed by the respective law.
     * \note This overload is used in models considering local thermal equilibrium
     */
    Scalar heatConductionFlux() const
    {
        static_assert(!enableThermalNonEquilibrium, "This only works for thermal equilibrium");
        if constexpr (enableEnergyBalance)
            return HeatConductionType::flux(this->problem(),
                                            this->element(),
                                            this->fvGeometry(),
                                            this->elemVolVars(),
                                            this->scvFace(),
                                            this->elemFluxVarsCache());
        else
            return 0.0;
    }

    /*!
     * \brief Returns the conductive flux computed by the respective law.
     * \note This overload is used in models considering local thermal nonequilibrium
     */
    Scalar heatConductionFlux([[maybe_unused]] const int phaseIdx) const
    {
        if constexpr (enableEnergyBalance)
            return HeatConductionType::flux(this->problem(),
                                            this->element(),
                                            this->fvGeometry(),
                                            this->elemVolVars(),
                                            this->scvFace(),
                                            phaseIdx,
                                            this->elemFluxVarsCache());
        else
            return 0.0;
    }

private:
    //! simple caching if advection flux is used twice with different upwind function
    mutable std::bitset<numPhases> advFluxIsCached_;
    mutable std::array<Scalar, numPhases> advFluxBeforeUpwinding_;
};

} // end namespace Dumux

#endif
