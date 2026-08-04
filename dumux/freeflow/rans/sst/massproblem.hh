// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::SSTMassProblem
 */
#ifndef DUMUX_RANS_SST_MASS_PROBLEM_HH
#define DUMUX_RANS_SST_MASS_PROBLEM_HH

#include <cmath>
#include <memory>
#include <vector>
#include <bitset>

#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/freeflow/turbulencemodel.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Mixin that adds the mass-domain side of the two-equation SST RANS turbulence closure:
 *        the production/destruction/cross-diffusion source term (through the standard
 *        problem.source() extension point, branching on the runtime RANS.SSTModelVersion), the
 *        once-per-time-step (lagged) update of the k/omega cross-diffusion gradient product
 *        needed by F1() and the cross-diffusion source term, and the wall-adjacent analytical
 *        omega Dirichlet condition omega_wall = 6*nu/(beta_omega*y^2), enforced as an internal
 *        (cell, not boundary-face) Dirichlet constraint - identical in every respect to
 *        Dumux::KOmegaMassProblem's own wall condition (SST does not override betaOmega(); the
 *        resulting 0.0708-vs-0.075 discrepancy is a deliberately preserved quirk of the ported
 *        formula, not a bug).
 *
 * \tparam TypeTag The mass sub-model's TypeTag.
 * \tparam MomentumProblem The concrete momentum problem type (a Dumux::SSTMomentumProblem),
 *         explicit template parameter, exactly as for Dumux::KOmegaMassProblem.
 *
 * Ported from the deleted releases/3.10:dumux/freeflow/rans/twoeq/sst/problem.hh
 * (RANSProblemImpl<TypeTag, TurbulenceModel::sst>) and
 * releases/3.10:test/freeflow/rans/problem.hh's isDirichletCell_/dirichletTurbulentTwoEq_
 * komega/sst-shared branches.
 */
template<class TypeTag, class MomentumProblem>
class SSTMassProblem : public NavierStokesMassProblem<TypeTag>
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
    SSTMassProblem(Args&&... args)
    : ParentType(std::forward<Args>(args)...)
    // The k/omega cross-diffusion gradient product is read (SSTMassVolumeVariables::update())
    // as soon as the grid variables are first initialized in main.cc, before
    // updateDynamicWallProperties() is ever called (once per time step, inside the time loop) -
    // so this must be pre-sized here, to zero, exactly as for Dumux::KOmegaMassProblem's
    // crossDiffusionGradientProduct_.
    , crossDiffusionGradientProduct_(this->gridGeometry().elementMapper().size(), 0.0)
    , storedTurbulentKineticEnergy_(this->gridGeometry().elementMapper().size(), 0.0)
    , storedDissipation_(this->gridGeometry().elementMapper().size(), 0.0)
    , sstModelVersion_(sstModelFromString(getParamFromGroup<std::string>(this->paramGroup(), "RANS.SSTModelVersion", "SST")))
    {}

    //! Must be called once, after both sub-problems have been constructed.
    void setMomentumProblem(std::shared_ptr<const MomentumProblem> momentumProblem)
    { momentumProblem_ = momentumProblem; }

    Scalar wallDistance(std::size_t eIdx) const
    { return momentumProblem_->wallDistance(eIdx); }

    //! \note Guarded the same way as Dumux::LowReKEpsilonMassProblem::kinematicViscosity() (see
    //! its comment): RANSMomentumProblem's stored molecular density/viscosity are only
    //! zero-*initialized* by updateStaticWallProperties(), so this quotient would otherwise
    //! evaluate to NaN in the narrow window before the first updateDynamicWallProperties() call.
    Scalar kinematicViscosity(std::size_t eIdx) const
    {
        const auto density = momentumProblem_->storedDensity(eIdx);
        return density > 0.0 ? momentumProblem_->kinematicViscosity(eIdx) : 0.0;
    }

    Scalar stressTensorScalarProduct(std::size_t eIdx) const
    { return momentumProblem_->stressTensorScalarProduct(eIdx); }

    Scalar vorticityTensorScalarProduct(std::size_t eIdx) const
    { return momentumProblem_->vorticityTensorScalarProduct(eIdx); }

    Scalar storedTurbulentKineticEnergy(std::size_t eIdx) const
    { return storedTurbulentKineticEnergy_[eIdx]; }

    Scalar storedDissipation(std::size_t eIdx) const
    { return storedDissipation_[eIdx]; }

    Scalar storedCrossDiffusionGradientProduct(std::size_t eIdx) const
    { return crossDiffusionGradientProduct_[eIdx]; }

    //! The runtime-selected SST model version (RANS.SSTModelVersion = "SST" (default) or "BSL").
    SSTModel sstModelVersion() const
    { return sstModelVersion_; }

    //! The von Karman constant and the beta_omega closure constant, needed for the wall omega
    //! value below (kept here, not just in volumevariables.hh, since this method needs it
    //! without a volVars object at hand) - identical, unblended constant to Dumux::KOmegaMassProblem's.
    Scalar betaOmega() const
    { return 0.0708; }

    /*!
     * \brief Updates the (lagged) stored k/omega and their cross-diffusion gradient product, via
     *        a central difference over each element's two neighbors along every axis - reusing
     *        the momentum domain's neighborIndex()/cellCenter() bookkeeping - reproducing
     *        releases/3.10:dumux/freeflow/rans/twoeq/sst/problem.hh's own gradient computation
     *        verbatim (a Green-Gauss reconstruction over CCTpfa's face/neighbor connectivity was
     *        used here previously; it gives a shorter effective distance at boundary-adjacent
     *        cells that was found to needlessly amplify SST's already gradient-sensitive
     *        near-wall blending). Called once per time step from main.cc, before the next
     *        assembly.
     */
    template<class SolutionVector>
    void updateDynamicWallProperties(const SolutionVector& sol)
    {
        const auto& gridGeometry = this->gridGeometry();

        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            fvGeometry.bind(element);
            const auto& insideScv = *scvs(fvGeometry).begin();
            const auto eIdx = insideScv.elementIndex();
            storedTurbulentKineticEnergy_[eIdx] = sol[insideScv.dofIndex()][Indices::turbulentKineticEnergyIdx];
            storedDissipation_[eIdx] = sol[insideScv.dofIndex()][Indices::dissipationIdx];
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

                kGradient[axisIdx] = (storedTurbulentKineticEnergy_[neighborIdx1] - storedTurbulentKineticEnergy_[neighborIdx0])/distance;
                omegaGradient[axisIdx] = (storedDissipation_[neighborIdx1] - storedDissipation_[neighborIdx0])/distance;
            }

            crossDiffusionGradientProduct_[eIdx] = kGradient*omegaGradient;
        }
    }

    /*!
     * \brief The SST production/destruction/cross-diffusion source term, ported from
     *        releases/3.10:dumux/freeflow/rans/twoeq/sst/staggered/localresidual.hh's
     *        computeSourceForCellCenter(), branching on the runtime SST model version. Added
     *        through the standard problem.source(...) extension point.
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

        const Scalar dynamicEddyViscosity = volVars.dynamicEddyViscosity(*this);
        const Scalar kinematicEddyViscosity = dynamicEddyViscosity/volVars.density();
        const Scalar productionTerm = 2.0*dynamicEddyViscosity*volVars.stressTensorScalarProduct();
        source[Indices::turbulentKineticEnergyEqIdx] += productionTerm;

        if (sstModelVersion() == SSTModel::BSL)
        {
            source[Indices::dissipationEqIdx] += volVars.gammaBSL()/kinematicEddyViscosity*productionTerm;
            source[Indices::turbulentKineticEnergyEqIdx] -= volVars.betaStarBSL()*volVars.density()*volVars.turbulentKineticEnergy()*volVars.dissipation();
            source[Indices::dissipationEqIdx] -= volVars.betaBSL()*volVars.density()*volVars.dissipation()*volVars.dissipation();
        }
        else
        {
            source[Indices::dissipationEqIdx] += volVars.gammaSST()/kinematicEddyViscosity*productionTerm;
            source[Indices::turbulentKineticEnergyEqIdx] -= volVars.betaStarSST()*volVars.density()*volVars.turbulentKineticEnergy()*volVars.dissipation();
            source[Indices::dissipationEqIdx] -= volVars.betaSST()*volVars.density()*volVars.dissipation()*volVars.dissipation();
        }

        // cross-diffusion term - identical in both branches, always weighted by (1-F1)
        const Scalar gradientProduct = volVars.storedCrossDiffusionGradientProduct();
        source[Indices::dissipationEqIdx] += 2.0*volVars.density()*(1.0 - volVars.F1())*volVars.sigmaOmega2()/volVars.dissipation()*gradientProduct;

        return source;
    }

    //! \name Wall-adjacent-cell analytical omega Dirichlet condition - identical mechanism and
    //! formula to Dumux::KOmegaMassProblem's (see class docs above for why SST reuses it as-is).
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
    std::vector<Scalar> storedTurbulentKineticEnergy_;
    std::vector<Scalar> storedDissipation_;
    SSTModel sstModelVersion_;
};

} // end namespace Dumux

#endif
