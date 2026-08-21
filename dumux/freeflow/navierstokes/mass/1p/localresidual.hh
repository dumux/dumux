// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesMassOnePLocalResidual
 */
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_MASS_1P_LOCAL_RESIDUAL_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_MASS_1P_LOCAL_RESIDUAL_HH

#include <type_traits>

#include <dumux/common/numeqvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/concepts/variables_.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

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
class NavierStokesMassOnePLocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVariablesCache = Concept::GridVariablesCache_t<GridVariables>;
    using ElementVariables = typename GridVariablesCache::LocalView;
    using Variables = Concept::Variables_t<GridVariables>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using SubControlVolume = typename ElementDiscretization::SubControlVolume;
    using SubControlVolumeFace = typename ElementDiscretization::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    using Extrusion = Extrusion_t<GridGeometry>;

public:
    //! Use the parent type's constructor
    using ParentType::ParentType;

    /*!
     * \brief Calculate the storage term of the equation
     */
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const Variables& vars) const
    {
        NumEqVector storage(0.0);
        storage[ModelTraits::Indices::conti0EqIdx] = vars.density();

        // consider energy storage for non-isothermal models
        if constexpr (ModelTraits::enableEnergyBalance())
            storage[ModelTraits::Indices::energyEqIdx] = vars.density() * vars.internalEnergy();

        return storage;
    }

    /*!
     * \brief Calculate the storage term of the equation
     *
     * \param elemDisc The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param scv The sub control volume
     * \param isPreviousTimeLevel If set to true, the storage term is evaluated on the previous time level.
     *
     */
    NumEqVector storageIntegral(const ElementDiscretization& elemDisc,
                                const ElementVariables& elemVars,
                                const SubControlVolume& scv,
                                bool isPreviousTimeLevel) const
    {
        const auto& vars = elemVars[scv];
        // We apply mass lumping here
        NumEqVector storage(0.0);
        storage[ModelTraits::Indices::conti0EqIdx] = vars.density();

        // consider energy storage for non-isothermal models
        if constexpr (ModelTraits::enableEnergyBalance())
            storage[ModelTraits::Indices::energyEqIdx] = vars.density() * vars.internalEnergy();

        storage *= Extrusion::volume(elemDisc, scv) * vars.extrusionFactor();

        return storage;
    }

    /*!
     * \brief Evaluate the mass flux over a face of a sub control volume.
     *
     * \param problem The problem
     * \param element The element
     * \param elemDisc The finite volume geometry context
     * \param elemVars The variables for all local dofs of the element
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache The cache related to flux computation
     */
    template<class ElementFluxVariablesCache>
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const ElementDiscretization& elemDisc,
                            const ElementVariables& elemVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        FluxVariables fluxVars;
        fluxVars.init(problem, element, elemDisc, elemVars, scvf, elemFluxVarsCache);
        auto flux = fluxVars.flux(0);

        // the auxiliary flux is enabled if the trait is specialized for the problem
        // this can be used, for example, to implement flux stabilization terms
        if constexpr (ImplementsAuxiliaryFluxNavierStokesMassOneP<Problem>::value)
            flux += problem.auxiliaryFlux(element, elemDisc, elemVars, elemFluxVarsCache, scvf);

        return flux;
    }

    /*!
     * \brief Calculates the flux integral over a face of a sub control volume.
     *
     * \param elemDisc The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param scvf The sub control volume face
     *
     */
    NumEqVector fluxIntegral(const ElementDiscretization& elemDisc,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        const auto& problem = this->asImp().problem();
        NumEqVector flux(0.0);
        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scvf))
        {
            const auto& faceIpData = qpData.ipData();
            flux += qpData.weight() * (problem.velocity(elemDisc, faceIpData) * faceIpData.unitOuterNormal());
        }

        static const auto upwindWeight
            = getParamFromGroup<Scalar>(this->problem().paramGroup(), "Flux.UpwindWeight", 1.0);

        const auto& insideVars = elemVars[elemDisc.scv(scvf.insideScvIdx())];
        const auto& outsideVars = elemVars[elemDisc.scv(scvf.outsideScvIdx())];

        flux *= (upwindWeight * insideVars.density() + (1.0 - upwindWeight) * outsideVars.density());

        // the auxiliary flux is enabled if the trait is specialized for the problem
        // this can be used, for example, to implement flux stabilization terms
        if constexpr (ImplementsAuxiliaryFluxNavierStokesMassOneP<Problem>::value)
            flux += problem.auxiliaryFlux(elemDisc.element(), elemDisc, elemVars, scvf);

        return flux * insideVars.extrusionFactor();
    }
};

} // end namespace Dumux

#endif
