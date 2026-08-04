// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::KEpsilonMassProblem
 */
#ifndef DUMUX_RANS_KEPSILON_MASS_PROBLEM_HH
#define DUMUX_RANS_KEPSILON_MASS_PROBLEM_HH

#include <memory>
#include <vector>
#include <bitset>
#include <cmath>
#include <iostream>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Mixin that adds the mass-domain side of the two-equation, high-Reynolds-number
 *        wall-function k-epsilon RANS turbulence closure: the per-time-step (lagged)
 *        matching-point/near-wall-region classification, the blended two-layer eddy viscosity,
 *        the production/destruction source term (with the matching-point production/
 *        destruction skip for k), and the near-wall-region/matching-point internal (whole-cell)
 *        Dirichlet constraints for k and epsilon.
 *
 * \tparam TypeTag The mass sub-model's TypeTag.
 * \tparam MomentumProblem The concrete momentum problem type (a Dumux::KEpsilonMomentumProblem),
 *         explicit template parameter, exactly as for Dumux::RANSMassOneEqProblem/
 *         Dumux::KOmegaMassProblem.
 *
 * The "matching point" is the cell whose wall distance is closest to the target y+ at which
 * the log-law and viscous sublayer profiles are matched; the "near-wall region" comprises the
 * cells between the wall and the matching point, where k/epsilon are not solved for but
 * instead fixed via internal Dirichlet constraints from the blended two-layer closure.
 *
 * wallDistance()/stressTensorScalarProduct()/wallElementIndex()/neighborIndex()/
 * velocityGradient()/kinematicViscosity()/karmanConstant() are *not* recomputed here - they are
 * frozen (lagged) values computed once per time step by the momentum domain's
 * Dumux::RANSMomentumProblem and forwarded, read-only, through the stored momentum-problem
 * pointer - the same lagging strategy already established for the one-equation and k-omega
 * models.
 */
template<class TypeTag, class MomentumProblem>
class KEpsilonMassProblem : public NavierStokesMassProblem<TypeTag>
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

public:
    template<class... Args>
    KEpsilonMassProblem(Args&&... args)
    : ParentType(std::forward<Args>(args)...)
    // Sized here (to zero/false), not lazily inside updateDynamicWallProperties(), since
    // KEpsilonMassVolumeVariables::update() is read as soon as the grid variables are first
    // initialized in main.cc - before updateDynamicWallProperties() is ever called (only
    // happens once per time step, inside the time loop). Same reasoning as
    // Dumux::RANSMassOneEqProblem's viscosityTildeGradient_/Dumux::KOmegaMassProblem's
    // crossDiffusionGradientProduct_.
    , matchingPointIdx_(this->gridGeometry().elementMapper().size(), 0)
    , storedTurbulentKineticEnergy_(this->gridGeometry().elementMapper().size(), 0.0)
    , storedDissipation_(this->gridGeometry().elementMapper().size(), 0.0)
    , storedDynamicEddyViscosity_(this->gridGeometry().elementMapper().size(), 0.0)
    , zeroEqDynamicEddyViscosity_(this->gridGeometry().elementMapper().size(), 0.0)
    {}

    //! Must be called once, after both sub-problems have been constructed.
    void setMomentumProblem(std::shared_ptr<const MomentumProblem> momentumProblem)
    { momentumProblem_ = momentumProblem; }

    Scalar wallDistance(std::size_t eIdx) const
    { return momentumProblem_->wallDistance(eIdx); }

    Scalar stressTensorScalarProduct(std::size_t eIdx) const
    { return momentumProblem_->stressTensorScalarProduct(eIdx); }

    Scalar kinematicViscosity(std::size_t eIdx) const
    { return momentumProblem_->kinematicViscosity(eIdx); }

    Scalar karmanConstant() const
    { return momentumProblem_->karmanConstant(); }

    Scalar cMu() const
    { return 0.09; }

    Scalar yPlusThreshold() const
    {
        static const Scalar threshold = getParamFromGroup<Scalar>(this->paramGroup(), "KEpsilon.YPlusThreshold", 30.0);
        return threshold;
    }

    Scalar storedTurbulentKineticEnergy(std::size_t eIdx) const
    { return storedTurbulentKineticEnergy_[eIdx]; }

    Scalar storedDissipation(std::size_t eIdx) const
    { return storedDissipation_[eIdx]; }

    Scalar storedDynamicEddyViscosity(std::size_t eIdx) const
    { return storedDynamicEddyViscosity_[eIdx]; }

    Scalar zeroEqDynamicEddyViscosity(std::size_t eIdx) const
    { return zeroEqDynamicEddyViscosity_[eIdx]; }

    std::size_t matchingPointIdx(std::size_t wallElementIdx) const
    { return matchingPointIdx_[wallElementIdx]; }

    //! \brief Returns if an element is the matching point of its wall-normal column.
    bool isMatchingPoint(std::size_t eIdx) const
    { return matchingPointIdx(momentumProblem_->wallElementIndex(eIdx)) == eIdx; }

    //! \brief Returns if an element is located in the near-wall (algebraic wall-function) region.
    bool inNearWallRegion(std::size_t eIdx) const
    {
        const auto wallElementIdx = momentumProblem_->wallElementIndex(eIdx);
        const auto matchingPointIndex = matchingPointIdx(wallElementIdx);
        return (wallElementIdx == matchingPointIndex) ? yPlusNominal(eIdx) < yPlusThreshold()
                                                      : yPlus(eIdx) < yPlusThreshold();
    }

    Scalar yPlus(std::size_t eIdx) const
    { return wallDistance(eIdx)*uStar(eIdx)/kinematicViscosity(eIdx); }

    Scalar yPlusNominal(std::size_t eIdx) const
    { return wallDistance(eIdx)*uStarNominal(eIdx)/kinematicViscosity(eIdx); }

    //! \brief The wall shear-stress velocity, from the actual velocity gradient at the wall.
    Scalar uStar(std::size_t eIdx) const
    {
        using std::abs;
        using std::sqrt;
        const auto wallElementIdx = momentumProblem_->wallElementIndex(eIdx);
        const auto wallNormalAxis = momentumProblem_->wallNormalAxis();
        const auto flowDirectionAxis = momentumProblem_->flowDirectionAxis();
        return sqrt(kinematicViscosity(wallElementIdx)
                    * abs(momentumProblem_->velocityGradient(wallElementIdx, flowDirectionAxis, wallNormalAxis)));
    }

    //! \brief The nominal wall shear-stress velocity, from the matching point's stored k
    //!        (accounts for a poor mesh approximation of the viscous sublayer).
    Scalar uStarNominal(std::size_t eIdx) const
    {
        using std::pow;
        using std::sqrt;
        const auto matchingPointIndex = matchingPointIdx(momentumProblem_->wallElementIndex(eIdx));
        return pow(cMu(), 0.25)*sqrt(storedTurbulentKineticEnergy(matchingPointIndex));
    }

    //! \brief The dissipation wall-function value (log-law).
    Scalar dissipationWallFunction(std::size_t eIdx) const
    {
        const auto uStarN = uStarNominal(eIdx);
        return uStarN*uStarN*uStarN/karmanConstant()/wallDistance(eIdx);
    }

    //! \brief The turbulent-kinetic-energy wall-function value: pinned to the matching point's
    //!        own stored k for every cell in the same column (near-wall region and the matching
    //!        point itself).
    Scalar turbulentKineticEnergyWallFunction(std::size_t eIdx) const
    { return storedTurbulentKineticEnergy(matchingPointIdx(momentumProblem_->wallElementIndex(eIdx))); }

    //! \brief The Van Driest (zero-equation) kinematic eddy viscosity, used in the near-wall
    //!        region, blended/scaled at the matching point ("two-layer model").
    Scalar zeroEqEddyViscosityModel(std::size_t eIdx) const
    {
        using std::abs;
        using std::exp;
        using std::sqrt;

        const Scalar yPlusValue = yPlus(eIdx);
        Scalar mixingLength = 0.0;
        if (yPlusValue > 0.0)
        {
            mixingLength = karmanConstant()*wallDistance(eIdx)
                           *(1.0 - exp(-yPlusValue/26.0))/sqrt(1.0 - exp(-0.26*yPlusValue));
        }

        const auto wallNormalAxis = momentumProblem_->wallNormalAxis();
        const auto flowDirectionAxis = momentumProblem_->flowDirectionAxis();
        const auto velocityGradient = momentumProblem_->velocityGradient(eIdx, flowDirectionAxis, wallNormalAxis);
        return mixingLength*mixingLength*abs(velocityGradient)*momentumProblem_->storedDensity(eIdx);
    }

    /*!
     * \brief Updates the (lagged) matching-point/near-wall-region classification and the
     *        blended two-layer eddy viscosity, from the given (typically: the last converged)
     *        mass solution. Called once per time step from main.cc, before the next assembly -
     *        not re-evaluated within Newton, matching how releases/3.10 treated all of this
     *        (see class documentation above).
     */
    template<class SolutionVector>
    void updateDynamicWallProperties(const SolutionVector& sol)
    {
        const auto& gridGeometry = this->gridGeometry();

        // update stored k/epsilon and the standard (non-blended) eddy viscosity
        auto fvGeometry = localView(gridGeometry);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            fvGeometry.bind(element);
            const auto& scv = *scvs(fvGeometry).begin();
            const auto eIdx = scv.elementIndex();
            storedTurbulentKineticEnergy_[eIdx] = sol[scv.dofIndex()][Indices::turbulentKineticEnergyIdx];
            storedDissipation_[eIdx] = sol[scv.dofIndex()][Indices::dissipationIdx];
            storedDynamicEddyViscosity_[eIdx] = cMu()*storedTurbulentKineticEnergy_[eIdx]*storedTurbulentKineticEnergy_[eIdx]
                                                / storedDissipation_[eIdx]*momentumProblem_->storedDensity(eIdx);
        }

        // find the matching point for each wall-normal column - ported verbatim from
        // releases/3.10:dumux/freeflow/rans/twoeq/kepsilon/problem.hh::updateDynamicWallProperties()
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);
            const auto wallNormalAxis = momentumProblem_->wallNormalAxis();
            const auto neighborIndex0 = momentumProblem_->neighborIndex(eIdx, wallNormalAxis, 0);
            const auto neighborIndex1 = momentumProblem_->neighborIndex(eIdx, wallNormalAxis, 1);

            if ((!inNearWallRegion(eIdx) && (inNearWallRegion(neighborIndex0) || inNearWallRegion(neighborIndex1)))
                || (!inNearWallRegion(eIdx) && eIdx == momentumProblem_->wallElementIndex(eIdx))
                || (inNearWallRegion(eIdx) && (momentumProblem_->wallElementIndex(neighborIndex0) != momentumProblem_->wallElementIndex(neighborIndex1))))
            {
                matchingPointIdx_[momentumProblem_->wallElementIndex(eIdx)] = eIdx;
            }
        }

        // compute the potential zero-eq (Van Driest) eddy viscosities for the two-layer model
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);
            zeroEqDynamicEddyViscosity_[eIdx] = zeroEqEddyViscosityModel(eIdx);
        }

        // scale the zero-eq eddy viscosity to match the two-equation model's value at the
        // matching point, then pin the matching point itself to the exact two-equation value
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);
            const auto matchingPointIndex = matchingPointIdx(momentumProblem_->wallElementIndex(eIdx));
            const Scalar scalingFactor = storedDynamicEddyViscosity(matchingPointIndex)/zeroEqDynamicEddyViscosity_[matchingPointIndex];
            if (!isMatchingPoint(eIdx) && !std::isnan(scalingFactor) && !std::isinf(scalingFactor))
                zeroEqDynamicEddyViscosity_[eIdx] *= scalingFactor;
        }
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);
            const auto matchingPointIndex = matchingPointIdx(momentumProblem_->wallElementIndex(eIdx));
            if (isMatchingPoint(eIdx))
                zeroEqDynamicEddyViscosity_[matchingPointIndex] = storedDynamicEddyViscosity(matchingPointIndex);
        }
    }

    /*!
     * \brief The k-epsilon production/destruction source term, ported from
     *        releases/3.10:dumux/freeflow/rans/twoeq/kepsilon/staggered/localresidual.hh's
     *        computeSourceForCellCenter(). k's production/destruction are skipped exactly at
     *        the matching point (local equilibrium hypothesis: production = dissipation there).
     */
    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);
        const auto& volVars = elemVolVars[scv];
        const Scalar productionTerm = 2.0*volVars.dynamicEddyViscosity()*volVars.stressTensorScalarProduct();

        if (!volVars.isMatchingPoint())
        {
            source[Indices::turbulentKineticEnergyEqIdx] += productionTerm;
            source[Indices::turbulentKineticEnergyEqIdx] -= volVars.dissipation()*volVars.density();
        }

        source[Indices::dissipationEqIdx] += volVars.cOneEpsilon()*volVars.dissipation()/volVars.turbulentKineticEnergy()*productionTerm;
        source[Indices::dissipationEqIdx] -= volVars.cTwoEpsilon()*volVars.dissipation()*volVars.dissipation()/volVars.turbulentKineticEnergy()*volVars.density();

        return source;
    }

    //! \name Near-wall-region/matching-point internal (whole-cell) Dirichlet constraints for
    //! k and epsilon - see releases/3.10's isDirichletCell_/dirichletTurbulentTwoEq_ kepsilon
    //! branches.
    // \{

    static constexpr bool enableInternalDirichletConstraints()
    { return true; }

    std::bitset<ModelTraits::numEq()> hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<ModelTraits::numEq()> constraints;
        const auto eIdx = scv.elementIndex();

        if (inNearWallRegion(eIdx))
        {
            constraints.set(Indices::turbulentKineticEnergyEqIdx);
            constraints.set(Indices::dissipationEqIdx);
        }
        else if (isMatchingPoint(eIdx))
            constraints.set(Indices::dissipationEqIdx);

        return constraints;
    }

    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);
        const auto eIdx = scv.elementIndex();
        values[Indices::turbulentKineticEnergyEqIdx] = turbulentKineticEnergyWallFunction(eIdx);
        values[Indices::dissipationEqIdx] = dissipationWallFunction(eIdx);
        return values;
    }

    // \}

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    std::shared_ptr<const MomentumProblem> momentumProblem_;
    std::vector<std::size_t> matchingPointIdx_;
    std::vector<Scalar> storedTurbulentKineticEnergy_;
    std::vector<Scalar> storedDissipation_;
    std::vector<Scalar> storedDynamicEddyViscosity_;
    std::vector<Scalar> zeroEqDynamicEddyViscosity_;
};

} // end namespace Dumux

#endif
