// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
 * \param UpScheme The upwind scheme to be applied to advective fluxes
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag,
         class UpScheme = UpwindScheme<GetPropType<TypeTag, Properties::GridGeometry>> >
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
    using UpwindScheme = UpScheme;
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
     * \brief Returns the compositional dispersion flux computed by the respective law.
     */
    Dune::FieldVector<Scalar, numComponents> compositionalDispersionFlux([[maybe_unused]] const int phaseIdx) const
    {
        if constexpr (enableCompositionalDispersion)
        {
            return DispersionFluxType::compositionalDispersionFlux(this->problem(),
                                                                   this->element(),
                                                                   this->fvGeometry(),
                                                                   this->elemVolVars(),
                                                                   this->scvFace(),
                                                                   phaseIdx,
                                                                   this->elemFluxVarsCache());
        }
        else
            return Dune::FieldVector<Scalar, numComponents>(0.0);
    }

    /*!
     * \brief Returns the thermal dispersion flux computed by the respective law.
     */
    Dune::FieldVector<Scalar, 1> thermalDispersionFlux([[maybe_unused]] const int phaseIdx = 0) const
    {
        if constexpr (enableThermalDispersion)
        {
            return DispersionFluxType::thermalDispersionFlux(this->problem(),
                                                             this->element(),
                                                             this->fvGeometry(),
                                                             this->elemVolVars(),
                                                             this->scvFace(),
                                                             phaseIdx,
                                                             this->elemFluxVarsCache());
        }
        else
            return Dune::FieldVector<Scalar, 1>(0.0);
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
