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
 * \brief Helper class to compute fluxes for porous medium flow models
 */

#ifndef DUMUX_POROUSMEDIUMFLOW_FLUX_HH
#define DUMUX_POROUSMEDIUMFLOW_FLUX_HH

#include <bitset>
#include <array>

#include <dumux/common/properties.hh>
#include <dumux/flux/upwindscheme.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumflowModels
 * \brief The porous medium flux helper that computes advective / convective,
 *        molecular diffusive and heat conduction fluxes.
 *
 * \param TypeTag The type tag for access to type traits
 */
template<class TypeTag>
class PorousMediumFlux
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr std::size_t numPhases = ModelTraits::numFluidPhases();
    static constexpr std::size_t numComponents = ModelTraits::numFluidComponents();

public:
    using UpwindScheme = UpwindScheme<GetPropType<TypeTag, Properties::GridGeometry>;
    using AdvectionType = GetPropType<TypeTag, Properties::AdvectionType>;
    using DispersionFluxType = GetPropType<TypeTag, Properties::DispersionFluxType>;
    using MolecularDiffusionType = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    using HeatConductionType = GetPropType<TypeTag, Properties::HeatConductionType>;

    static constexpr bool enableAdvection = ModelTraits::enableAdvection();
    static constexpr bool enableMolecularDiffusion = ModelTraits::enableMolecularDiffusion();
    static constexpr bool enableCompositionalDispersion = ModelTraits::enableCompositionalDispersion();
    static constexpr bool enableThermalDispersion = ModelTraits::enableThermalDispersion();
    static constexpr bool enableEnergyBalance = ModelTraits::enableEnergyBalance();
    static constexpr bool enableThermalNonEquilibrium = getPropValue<TypeTag, Properties::EnableThermalNonEquilibrium>();

    //! The constructor resets the cache
    PorousMediumFlux()
    {
        advFluxIsCached_.reset();
        advFluxBeforeUpwinding_.fill(0.0);
    }

    /*!
     * \brief Returns the advective flux computed by the respective law.
     */
    template<typename FunctionType, typename Context>
    Scalar advectiveFlux([[maybe_unused]] const int phaseIdx,
                         [[maybe_unused]] const FunctionType& upwindTerm,
                         [[maybe_unused]] const Context& context) const
    {
        if constexpr (enableAdvection)
        {
            if (!advFluxIsCached_[phaseIdx])
            {

                advFluxBeforeUpwinding_[phaseIdx] = AdvectionType::flux(
                    context.problem(), context.element(), context.fvGeometry(),
                    context.elemVolVars(), context.scvFace(), phaseIdx,
                    context.elemFluxVarsCache()
                );
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
    template<typename Context>
    Dune::FieldVector<Scalar, numComponents> molecularDiffusionFlux([[maybe_unused]] const int phaseIdx,
                                                                    [[maybe_unused]] const Context& context) const
    {
        if constexpr (enableMolecularDiffusion)
            return MolecularDiffusionType::flux(
                context.problem(), context.element(), context.fvGeometry(),
                context.elemVolVars(), context.scvFace(), phaseIdx,
                context.elemFluxVarsCache()
            );
        else
            return Dune::FieldVector<Scalar, numComponents>(0.0);
    }

    /*!
     * \brief Returns the compositional dispersion flux computed by the respective law.
     */
    template<typename Context>
    Dune::FieldVector<Scalar, numComponents> compositionalDispersionFlux([[maybe_unused]] const int phaseIdx,
                                                                         [[maybe_unused]] const Context& context) const
    {
        if constexpr (enableCompositionalDispersion)
        {
            return DispersionFluxType::compositionalDispersionFlux(
                context.problem(), context.element(), context.fvGeometry(),
                context.elemVolVars(), context.scvFace(), phaseIdx,
                context.elemFluxVarsCache()
            );
        }
        else
            return Dune::FieldVector<Scalar, numComponents>(0.0);
    }

    /*!
     * \brief Returns the thermal dispersion flux computed by the respective law.
     */
    template<typename Context>
    Dune::FieldVector<Scalar, 1> thermalDispersionFlux([[maybe_unused]] const int phaseIdx = 0,
                                                       [[maybe_unused]] const Context& context) const
    {
        if constexpr (enableThermalDispersion)
        {
            return DispersionFluxType::thermalDispersionFlux(
                context.problem(), context.element(), context.fvGeometry(),
                context.elemVolVars(), context.scvFace(), phaseIdx,
                context.elemFluxVarsCache()
            );
        }
        else
            return Dune::FieldVector<Scalar, 1>(0.0);
    }

    /*!
     * \brief Returns the conductive flux computed by the respective law.
     * \note This overload is used in models considering local thermal equilibrium
     */
    template<typename Context>
    Scalar heatConductionFlux([[maybe_unused]] const Context& context) const
    {
        static_assert(!enableThermalNonEquilibrium, "This only works for thermal equilibrium");
        if constexpr (enableEnergyBalance)
            return HeatConductionType::flux(
                context.problem(), context.element(), context.fvGeometry(),
                context.elemVolVars(), context.scvFace(), context.elemFluxVarsCache()
            );
        else
            return 0.0;
    }

    /*!
     * \brief Returns the conductive flux computed by the respective law.
     * \note This overload is used in models considering local thermal nonequilibrium
     */
    template<typename Context>
    Scalar heatConductionFlux([[maybe_unused]] const int phaseIdx, [[maybe_unused]] const Context& context) const
    {
        if constexpr (enableEnergyBalance)
            return HeatConductionType::flux(
                context.problem(), context.element(), context.fvGeometry(),
                context.elemVolVars(), context.scvFace(), phaseIdx,
                context.elemFluxVarsCache()
            );
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
