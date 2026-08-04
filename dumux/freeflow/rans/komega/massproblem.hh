// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::KOmegaMassProblem
 */
#ifndef DUMUX_RANS_KOMEGA_MASS_PROBLEM_HH
#define DUMUX_RANS_KOMEGA_MASS_PROBLEM_HH

#include <cmath>
#include <memory>
#include <vector>
#include <bitset>

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Mixin that adds the mass-domain side of the two-equation k-omega RANS turbulence
 *        closure: the production/destruction/cross-diffusion source term (through the
 *        standard problem.source() extension point), the once-per-time-step (lagged) update
 *        of the k/omega gradients needed for the cross-diffusion term, and the wall-adjacent
 *        analytical omega Dirichlet condition ω_wall = 6ν/(β_ω·y²), enforced as an internal
 *        (cell, not boundary-face) Dirichlet constraint via the standard
 *        problem.hasInternalDirichletConstraint()/internalDirichlet() extension point.
 *
 * \tparam TypeTag The mass sub-model's TypeTag.
 * \tparam MomentumProblem The concrete momentum problem type (a Dumux::KOmegaMomentumProblem),
 *         explicit template parameter, exactly as for Dumux::RANSMassOneEqProblem.
 *
 * stressTensorScalarProduct() is *not* recomputed here - it is a frozen (lagged) value
 * computed once per time step by the momentum domain's Dumux::RANSMomentumProblem and
 * forwarded, read-only, through the stored momentum-problem pointer - the same lagging
 * strategy Dumux::RANSMassOneEqProblem already established.
 */
template<class TypeTag, class MomentumProblem>
class KOmegaMassProblem : public NavierStokesMassProblem<TypeTag>
{
    using ParentType = NavierStokesMassProblem<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    static constexpr int dim = GridView::dimension;
    using DimVector = Dune::FieldVector<Scalar, dim>;

public:
    template<class... Args>
    KOmegaMassProblem(Args&&... args)
    : ParentType(std::forward<Args>(args)...)
    // k/omega gradients are read (KOmegaMassVolumeVariables::update()) as soon as the grid
    // variables are first initialized in main.cc, before updateDynamicWallProperties() is ever
    // called (once per time step, inside the time loop) - so this must be pre-sized here, to
    // zero, exactly as for Dumux::RANSMassOneEqProblem's viscosityTildeGradient_.
    , crossDiffusionGradientProduct_(this->gridGeometry().elementMapper().size(), 0.0)
    {}

    //! Must be called once, after both sub-problems have been constructed.
    void setMomentumProblem(std::shared_ptr<const MomentumProblem> momentumProblem)
    { momentumProblem_ = momentumProblem; }

    Scalar wallDistance(std::size_t eIdx) const
    { return momentumProblem_->wallDistance(eIdx); }

    Scalar stressTensorScalarProduct(std::size_t eIdx) const
    { return momentumProblem_->stressTensorScalarProduct(eIdx); }

    Scalar storedCrossDiffusionGradientProduct(std::size_t eIdx) const
    { return crossDiffusionGradientProduct_[eIdx]; }

    //! The von Karman constant and the beta_omega closure constant, needed for the wall
    //! omega value below (kept here, not just in volumevariables.hh, since this method needs
    //! it without a volVars object at hand).
    Scalar betaOmega() const
    { return 0.0708; }

    /*!
     * \brief Updates the (lagged) k/omega gradients used by the cross-diffusion source term, via
     *        a central difference over each element's two neighbors along every axis - reusing
     *        the momentum domain's neighborIndex()/cellCenter() bookkeeping - reproducing
     *        releases/3.10:dumux/freeflow/rans/twoeq/komega/problem.hh's own gradient
     *        computation verbatim (a Green-Gauss reconstruction over CCTpfa's face/neighbor
     *        connectivity was used here previously; it gives a shorter effective distance at
     *        boundary-adjacent cells, see Dumux::SSTMassProblem for where that was found to
     *        matter). Called once per time step from main.cc, before the next assembly.
     */
    template<class SolutionVector>
    void updateDynamicWallProperties(const SolutionVector& sol)
    {
        const auto& gridGeometry = this->gridGeometry();
        const auto numElements = gridGeometry.elementMapper().size();
        crossDiffusionGradientProduct_.assign(numElements, 0.0);

        std::vector<Scalar> k(numElements);
        std::vector<Scalar> omega(numElements);
        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            fvGeometry.bind(element);
            const auto& insideScv = *scvs(fvGeometry).begin();
            const auto eIdx = insideScv.elementIndex();
            k[eIdx] = sol[insideScv.dofIndex()][Indices::turbulentKineticEnergyIdx];
            omega[eIdx] = sol[insideScv.dofIndex()][Indices::dissipationIdx];
        }

        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);

            DimVector kGradient(0.0);
            DimVector omegaGradient(0.0);
            for (int axisIdx = 0; axisIdx < dim; ++axisIdx)
            {
                const auto neighborIdx0 = momentumProblem_->neighborIndex(eIdx, axisIdx, 0);
                const auto neighborIdx1 = momentumProblem_->neighborIndex(eIdx, axisIdx, 1);
                const auto distance = momentumProblem_->cellCenter(neighborIdx1)[axisIdx]
                                     - momentumProblem_->cellCenter(neighborIdx0)[axisIdx];
                if (std::abs(distance) < 1e-8)
                    continue;

                kGradient[axisIdx] = (k[neighborIdx1] - k[neighborIdx0])/distance;
                omegaGradient[axisIdx] = (omega[neighborIdx1] - omega[neighborIdx0])/distance;
            }

            crossDiffusionGradientProduct_[eIdx] = kGradient*omegaGradient;
        }
    }

    /*!
     * \brief The k-omega production/destruction/cross-diffusion source term, ported from
     *        releases/3.10:dumux/freeflow/rans/twoeq/komega/staggered/localresidual.hh's
     *        computeSourceForCellCenter(). Added through the standard problem.source(...)
     *        extension point.
     */
    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        using std::max;

        NumEqVector source(0.0);
        const auto& volVars = elemVolVars[scv];

        const Scalar productionTerm = 2.0*volVars.dynamicEddyViscosity()*volVars.stressTensorScalarProduct();
        source[Indices::turbulentKineticEnergyEqIdx] += productionTerm;
        source[Indices::dissipationEqIdx] += volVars.alpha()*volVars.dissipation()/volVars.turbulentKineticEnergy()*productionTerm;

        source[Indices::turbulentKineticEnergyEqIdx] -= volVars.betaK()*volVars.density()*volVars.turbulentKineticEnergy()*volVars.dissipation();
        source[Indices::dissipationEqIdx] -= volVars.betaOmega()*volVars.density()*volVars.dissipation()*volVars.dissipation();

        const Scalar gradientProduct = volVars.storedCrossDiffusionGradientProduct();
        if (gradientProduct > 0.0)
            source[Indices::dissipationEqIdx] += 0.125*volVars.density()/volVars.dissipation()*gradientProduct;

        return source;
    }

    //! \name Wall-adjacent-cell analytical omega Dirichlet condition (an internal, not
    //! boundary-face, constraint - see releases/3.10's isDirichletCell_/
    //! dirichletTurbulentTwoEq_ komega branches).
    // \{

    static constexpr bool enableInternalDirichletConstraints()
    { return true; }

    std::bitset<ModelTraits::numEq()> hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<ModelTraits::numEq()> constraints;

        auto fvGeometry = localView(this->gridGeometry());
        fvGeometry.bindElement(element);
        for (const auto& scvf : scvfs(fvGeometry))
            if (scvf.boundary() && this->asImp_().isOnWallAtPos(scvf.center()))
                constraints.set(Indices::dissipationEqIdx);

        return constraints;
    }

    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);
        const auto eIdx = scv.elementIndex();
        const auto y = wallDistance(eIdx);
        values[Indices::dissipationEqIdx] = 6.0*momentumProblem_->kinematicViscosity(eIdx)/(betaOmega()*y*y);
        return values;
    }

    // \}

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    std::shared_ptr<const MomentumProblem> momentumProblem_;
    std::vector<Scalar> crossDiffusionGradientProduct_;
};

} // end namespace Dumux

#endif
