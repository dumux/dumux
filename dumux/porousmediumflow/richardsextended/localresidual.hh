// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ExtendedRichardsModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the extended Richards fully implicit models.
 */

#ifndef DUMUX_RICHARDSEXTENDED_LOCAL_RESIDUAL_HH
#define DUMUX_RICHARDSEXTENDED_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/porousmediumflow/richards/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup ExtendedRichardsModel
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the extended Richards fully implicit models.
 */
template<class TypeTag>
class ExtendedRichardsLocalResidual : public RichardsLocalResidual<TypeTag>
{
    using Implementation = GetPropType<TypeTag, Properties::LocalResidual>;

    using ParentType = RichardsLocalResidual<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using EnergyLocalResidual = GetPropType<TypeTag, Properties::EnergyLocalResidual>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    // first index for the mass balance
    enum { conti0EqIdx = Indices::conti0EqIdx };

    // phase & component indices
    static constexpr auto liquidPhaseIdx = FluidSystem::phase0Idx;
    static constexpr auto gasPhaseIdx = FluidSystem::phase1Idx;
    static constexpr auto liquidCompIdx = FluidSystem::comp0Idx;

public:
    using ParentType::ParentType;

    /*!
     * \brief Evaluates the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub-control
     *        volume of a finite volume element for the immiscible models.
     * \param problem The problem
     * \param scv The sub control volume
     * \param volVars The current or previous volVars
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage = ParentType::computeStorage(problem, scv, volVars);
        // for extended Richards we consider water in air
        storage[conti0EqIdx] += volVars.porosity()
                                * volVars.molarDensity(gasPhaseIdx)
                                * volVars.moleFraction(gasPhaseIdx, liquidCompIdx)
                                * FluidSystem::molarMass(liquidCompIdx)
                                * volVars.saturation(gasPhaseIdx);

        return storage;
    }


    /*!
     * \brief Evaluates the mass flux over a face of a sub control volume.
     *
     * \param problem The problem
     * \param element The current element.
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache The cache related to flux computation
     */
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        NumEqVector flux = ParentType::computeFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        // for extended Richards we consider water vapor diffusion in air
        //check for the reference system and adapt units of the diffusive flux accordingly.
        if (FluxVariables::MolecularDiffusionType::referenceSystemFormulation() == ReferenceSystemFormulation::massAveraged)
            flux[conti0EqIdx] += fluxVars.molecularDiffusionFlux(gasPhaseIdx)[liquidCompIdx];
        else
            flux[conti0EqIdx] += fluxVars.molecularDiffusionFlux(gasPhaseIdx)[liquidCompIdx]*FluidSystem::molarMass(liquidCompIdx);

        return flux;
    }
};

} // end namespace Dumux

#endif
