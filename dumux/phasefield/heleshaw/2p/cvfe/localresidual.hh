// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup HeleShawModel
 * \brief Element-wise calculation of the local residual for the hybrid CVFE
 *        Hele-Shaw two-phase model.
 */
#ifndef DUMUX_PHASEFIELD_HELESHAW_2P_CVFE_LOCAL_RESIDUAL_HH
#define DUMUX_PHASEFIELD_HELESHAW_2P_CVFE_LOCAL_RESIDUAL_HH

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
 * \ingroup HeleShawModel
 * \brief Element-wise calculation of the hybrid CVFE Hele-Shaw local residual.
 *
 * Local dofs with an associated sub-control volume (e.g. the vertex dofs of PQ2/PQ3) get
 * classic control-volume storage/flux contributions (`storageIntegral`/`fluxIntegral`),
 * with the Darcy-advective phase-field term upwinded exactly as in the classic Box
 * model. Local dofs without an associated sub-control volume (e.g. the edge dofs of
 * PQ2/PQ3) get finite-element (weak form) contributions added on top, via
 * `HeleShawTwoPCVFEFELocalResidualTerms` (see
 * `addToElementStorageResidual`/`addToElementFluxAndSourceResidual`); there the
 * advective term is a plain (non-upwinded) Galerkin discretization.
 */
template<class TypeTag>
class HeleShawTwoPCVFELocalResidual
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
    using FeResidual = HeleShawTwoPCVFEFELocalResidualTerms<NumEqVector, Indices>;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    /*!
     * \brief Calculate the storage integral for a control-volume local dof.
     *        Only the phase field is transported (the pressure equation is
     *        elliptic, the chemical potential equation is algebraic).
     */
    NumEqVector storageIntegral(const FVElementGeometry& fvGeometry,
                                const ElementVariables& elemVars,
                                const SubControlVolume& scv,
                                bool isPreviousTimeLevel) const
    {
        const auto& vars = elemVars[scv];

        NumEqVector storage(0.0);
        storage[Indices::phaseFieldEqIdx] = vars.phaseField();
        storage *= Extrusion::volume(fvGeometry, scv) * vars.extrusionFactor();

        return storage;
    }

    /*!
     * \brief Calculate the source integral for a control-volume local dof.
     *
     * In addition to the (possibly solution-dependent) source provided by the
     * problem (e.g. the derivative of the double-well free energy), the chemical
     * potential equation has a model-intrinsic algebraic identity contribution
     * \f$ \mu \f$, closing the system \f$ \mu = \partial f/\partial \phi - \gamma \nabla^2 \phi \f$.
     */
    NumEqVector sourceIntegral(const FVElementGeometry& fvGeometry,
                               const ElementVariables& elemVars,
                               const SubControlVolume& scv) const
    {
        const auto& problem = this->asImp().problem();

        NumEqVector source(0.0);
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scv))
        {
            NumEqVector sourceAtIp = problem.source(fvGeometry, elemVars, qpData.ipData());
            sourceAtIp[Indices::chemPotEqIdx] += elemVars[scv].chemicalPotential();
            source += qpData.weight()*sourceAtIp;
        }

        source *= elemVars[scv].extrusionFactor();
        return source;
    }

    /*!
     * \brief Calculate the flux integral over a sub-control volume face.
     *
     * Mirrors the classic Box model's `computeFlux` term-for-term (Darcy
     * continuity, upwinded phase-field advection plus Cahn-Hilliard
     * diffusion, and the chemical-potential Laplacian term), generalized
     * from a single evaluation point to a proper quadrature rule (needed
     * since PQ2 gradients vary within a face, unlike Box's piecewise-constant
     * gradients).
     */
    NumEqVector fluxIntegral(const FVElementGeometry& fvGeometry,
                             const ElementVariables& elemVars,
                             const SubControlVolumeFace& scvf) const
    {
        const auto& problem = this->asImp().problem();
        const auto& insideVars = elemVars[fvGeometry.scv(scvf.insideScvIdx())];
        const auto& outsideVars = elemVars[fvGeometry.scv(scvf.outsideScvIdx())];

        NumEqVector flux(0.0);
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
        {
            const auto& ipData = qpData.ipData();
            const auto& ipCache = cache(elemVars, ipData);
            const HeleShawTwoPCVFEFluxFunctionContext context(fvGeometry, elemVars, ipCache);

            const auto phi = context.phaseField();
            const auto mu = context.chemicalPotential();
            const auto rhoMix = problem.mixtureDensity(phi);
            const auto etaMix = problem.mixtureViscosity(phi);
            const auto lambda = problem.permeability()/etaMix;
            const auto& n = ipData.unitOuterNormal();
            const auto& g = problem.gravity();

            // Darcy volume flux per unit area: v.n = -lambda*(gradP - mu*gradPhi - rhoMix*g).n
            const auto darcyFluxPerArea =
                -lambda*(context.gradPressure()*n - mu*(context.gradPhaseField()*n) - rhoMix*(g*n));

            NumEqVector fluxAtIp(0.0);
            fluxAtIp[Indices::continuityEqIdx] = darcyFluxPerArea;

            // upwind the phase field based on the Darcy flux direction
            const auto phiUpwind = (darcyFluxPerArea >= 0.0) ? insideVars.phaseField() : outsideVars.phaseField();
            fluxAtIp[Indices::phaseFieldEqIdx] = phiUpwind*darcyFluxPerArea
                                                - problem.chMobility()*(context.gradChemicalPotential()*n);

            fluxAtIp[Indices::chemPotEqIdx] = -problem.surfaceTensionCoeff()*(context.gradPhaseField()*n);

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
