// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CahnHilliardModel
 * \brief Element-wise calculation of the local residual for the Cahn-Hilliard model
 *        using hybrid control-volume/finite-element (CVFE) discretizations.
 */
#ifndef DUMUX_PHASEFIELD_CAHNHILLIARD_LOCAL_RESIDUAL_HH
#define DUMUX_PHASEFIELD_CAHNHILLIARD_LOCAL_RESIDUAL_HH

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

#include "felocalresidual.hh"
#include "flux.hh"

namespace Dumux {

/*!
 * \ingroup CahnHilliardModel
 * \brief Element-wise calculation of the Cahn-Hilliard local residual for hybrid CVFE discretizations.
 *
 * Local dofs with an associated sub-control volume (e.g. the vertex dofs of PQ2/PQ3) get
 * classic control-volume storage/flux contributions (`storageIntegral`/`fluxIntegral`).
 * Local dofs without an associated sub-control volume (e.g. the edge dofs of PQ2/PQ3) get
 * finite-element (weak form) contributions added on top, via
 * `CahnHilliardFELocalResidualTerms` (see `addToElementStorageResidual`/`addToElementFluxAndSourceResidual`).
 */
template<class TypeTag>
class CahnHilliardLocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVariablesCache = typename GridVariables::GridVariablesCache;
    using ElementVariables = typename GridVariablesCache::LocalView;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    using Extrusion = Extrusion_t<GridGeometry>;
    using FeResidual = CahnHilliardFELocalResidualTerms<NumEqVector, Indices>;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    /*!
     * \brief Calculate the storage integral for a control-volume local dof.
     *
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param scv The sub-control volume over which we integrate the storage term
     * \param isPreviousTimeLevel If set to true, the storage term is evaluated on the previous time level.
     */
    NumEqVector storageIntegral(const FVElementGeometry& fvGeometry,
                                const ElementVariables& elemVars,
                                const SubControlVolume& scv,
                                bool isPreviousTimeLevel) const
    {
        const auto& vars = elemVars[scv];

        NumEqVector storage(0.0);
        // the chemical potential equation has no storage term
        storage[Indices::massBalanceEqIdx] = vars.concentration();
        storage *= Extrusion::volume(fvGeometry, scv) * vars.extrusionFactor();

        return storage;
    }

    /*!
     * \brief Calculate the source integral for a control-volume local dof.
     *
     * In addition to the (possibly solution-dependent) source provided by the
     * problem (e.g. the derivative of a double-well free energy), the chemical
     * potential equation has a model-intrinsic algebraic identity contribution
     * \f$ \mu \f$, closing the system \f$ \mu = \partial f/\partial c - \gamma \nabla^2 c \f$.
     *
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param scv The sub-control volume over which we integrate the source term
     */
    NumEqVector sourceIntegral(const FVElementGeometry& fvGeometry,
                               const ElementVariables& elemVars,
                               const SubControlVolume& scv) const
    {
        const auto& problem = this->asImp().problem();

        // Note: a control-volume dof's source is evaluated "nodally" (a
        // single representative value for the whole scv, via
        // Concept::LocalDofIpData), so unlike the element (FE-dof) reaction
        // integral below, bumping this quadrature's order would not add any
        // accuracy -- the integrand is constant over the scv regardless of
        // how many points it's sampled at.
        NumEqVector source(0.0);
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scv))
        {
            NumEqVector sourceAtIp = problem.source(fvGeometry, elemVars, qpData.ipData());
            sourceAtIp[Indices::chemicalPotentialEqIdx] += elemVars[scv].chemicalPotential();
            source += qpData.weight()*sourceAtIp;
        }

        source *= elemVars[scv].extrusionFactor();
        return source;
    }

    /*!
     * \brief Calculate the flux integral over a sub-control volume face.
     *
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables for all local dofs of the element
     * \param scvf The sub-control volume face
     */
    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        const auto& problem = this->asImp().problem();

        NumEqVector flux(0.0);
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
        {
            const auto& ipData = qpData.ipData();
            const auto& ipCache = cache(elemVars, ipData);
            const CahnHilliardFluxFunctionContext context(fvGeometry, elemVars, ipCache);

            NumEqVector fluxAtIp(0.0);
            fluxAtIp[Indices::massBalanceEqIdx] = -1.0*vtmv(ipData.unitOuterNormal(), problem.mobility(), context.gradChemicalPotential());
            fluxAtIp[Indices::chemicalPotentialEqIdx] = -1.0*vtmv(ipData.unitOuterNormal(), problem.surfaceTension(), context.gradConcentration());
            flux += qpData.weight()*fluxAtIp;
        }

        flux *= elemVars[fvGeometry.scv(scvf.insideScvIdx())].extrusionFactor();
        return flux;
    }

    void addToElementStorageResidual(ElementResidualVector& residual,
                                     const Problem& problem,
                                     const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVariables& prevElemVars,
                                     const ElementVariables& curElemVars) const
    {
        FeResidual::addStorageTerms(
            residual, problem, fvGeometry, prevElemVars, curElemVars, this->timeLoop().timeStepSize()
        );
    }

    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVariables& elemVars) const
    {
        FeResidual::addFluxAndSourceTerms(
            residual, problem, fvGeometry, elemVars
        );
    }
};

} // end namespace Dumux

#endif
