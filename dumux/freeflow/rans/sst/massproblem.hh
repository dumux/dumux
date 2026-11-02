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
 *        (cell, not boundary-face) Dirichlet constraint.
 *
 * \tparam TypeTag The mass sub-model's TypeTag.
 * \tparam MomentumProblem The concrete momentum problem type (a Dumux::SSTMomentumProblem)
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
    // The stored fields are read by SSTMassVolumeVariables::update() when the grid
    // variables are initialized, i.e. before updateDynamicWallProperties() is first called,
    // so they are sized and zero-initialized here.
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

    //! \note Returns zero as long as the stored density of the momentum problem is zero (i.e.
    //! before the first updateDynamicWallProperties() call), to avoid a division by zero.
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

    //! The (unblended) beta_omega closure constant, needed for the wall value of omega
    //! where no volume variables are available.
    Scalar betaOmega() const
    { return 0.0708; }

    /*!
     * \brief Updates the (lagged) stored k and omega and the product of their gradients.
     *        The gradients are computed by central differences over the two neighboring cells
     *        along each axis, using the neighborIndex() and cellCenter() of the momentum
     *        problem. To be called once per time step before the next assembly.
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
     * \brief The SST production/destruction/cross-diffusion source term, using the
     *        coefficients of the selected SST model version.
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

        const Scalar dynamicEddyViscosity = volVars.dynamicEddyViscosity();
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

        // cross-diffusion term, the same for both model versions, weighted by (1-F1)
        const Scalar gradientProduct = volVars.storedCrossDiffusionGradientProduct();
        source[Indices::dissipationEqIdx] += 2.0*volVars.density()*(1.0 - volVars.F1())*volVars.sigmaOmega2()/volVars.dissipation()*gradientProduct;

        return source;
    }

    //! \name Analytical omega Dirichlet condition in the wall-adjacent cells (an internal
    //! constraint on the cell, not a boundary face condition).
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
