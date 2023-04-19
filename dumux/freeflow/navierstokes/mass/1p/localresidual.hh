// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesResidualImpl
 */
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_MASS_1P_LOCAL_RESIDUAL_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_MASS_1P_LOCAL_RESIDUAL_HH

#include <type_traits>

#include <dumux/common/numeqvector.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief Traits class to be specialized for problems to add auxiliary fluxes
 * \note This can be used, for example, to implement flux stabilization terms
 */
template<class Problem>
struct ImplementsAuxiliaryFluxNavierStokesMassOneP
: public std::false_type
{};

/*!
 * \ingroup NavierStokesModel
 * \brief Element-wise calculation of the Navier-Stokes residual for single-phase flow.
 */
template<class TypeTag>
class NavierStokesMassOnePLocalResidual : public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using ElementFluxVariablesCache = typename GridFluxVariablesCache::LocalView;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

public:
    //! Use the parent type's constructor
    using ParentType::ParentType;

    /*!
     * \brief Calculate the storage term of the equation
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);
        storage[ModelTraits::Indices::conti0EqIdx] = volVars.density();

        // consider energy storage for non-isothermal models
        if constexpr (ModelTraits::enableEnergyBalance())
            storage[ModelTraits::Indices::energyEqIdx] = volVars.density() * volVars.internalEnergy();

        return storage;
    }

    /*!
     * \brief Evaluate the mass flux over a face of a sub control volume.
     *
     * \param problem The problem
     * \param element The element
     * \param fvGeometry The finite volume geometry context
     * \param elemVolVars The volume variables for all flux stencil elements
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
        auto flux = fluxVars.flux(0);

        // the auxiliary flux is enabled if the trait is specialized for the problem
        // this can be used, for example, to implement flux stabilization terms
        if constexpr (ImplementsAuxiliaryFluxNavierStokesMassOneP<Problem>::value)
            flux += problem.auxiliaryFlux(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);

        return flux;
    }
};

} // end namespace Dumux

#endif
