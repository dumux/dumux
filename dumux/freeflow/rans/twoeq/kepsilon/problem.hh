// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup KEpsilonModel
 * \copydoc Dumux::KEpsilonProblem
 */
#ifndef DUMUX_KEPSILON_PROBLEM_HH
#define DUMUX_KEPSILON_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/staggeredfvproblem.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/freeflow/rans/problem.hh>

#include "model.hh"

namespace Dumux
{

/*!
 * \ingroup KEpsilonModel
 * \brief K-epsilon turbulence problem base class.
 *
 * This implements some base functionality for k-epsilon models.
 */
template<class TypeTag>
class KEpsilonProblem : public RANSProblem<TypeTag>
{
    using ParentType = RANSProblem<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, Problem);

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;

    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using Indices = typename ModelTraits::Indices;

    static constexpr bool enableEnergyBalance = ModelTraits::enableEnergyBalance();
    static constexpr bool isCompositional = ModelTraits::numComponents() > 1;

    // account for the offset of the cell center privars within the PrimaryVariables container
    static constexpr auto cellCenterOffset = ModelTraits::numEq() - CellCenterPrimaryVariables::dimension;
    static_assert(cellCenterOffset == ModelTraits::dim(), "cellCenterOffset must equal dim for staggered NavierStokes");

public:

    //! The constructor sets the gravity, if desired by the user.
    KEpsilonProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, paramGroup)
    {
        yPlusThreshold_ = getParamFromGroup<Scalar>(this->paramGroup(), "KEpsilon.YPlusThreshold", 30);
        useStoredEddyViscosity_ = getParamFromGroup<bool>(this->paramGroup(), "RANS.UseStoredEddyViscosity", false);
    }

    /*!
     * \brief Correct size of the static (solution independent) wall variables
     */
    void updateStaticWallProperties()
    {
        ParentType::updateStaticWallProperties();

        // update size and initial values of the global vectors
        matchingPointIdx_.resize(this->fvGridGeometry().elementMapper().size(), 0);
        storedDensity_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        storedDissipation_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        storedTurbulentKineticEnergy_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        storedDynamicEddyViscosity_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        zeroEqDynamicEddyViscosity_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
    }

    /*!
     * \brief Update the dynamic (solution dependent) relations to the walls
     *
     * \param curSol The solution vector.
     */
    void updateDynamicWallProperties(const SolutionVector& curSol)
    {
        ParentType::updateDynamicWallProperties(curSol);

        // update the stored eddy viscosities
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);

            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                const int dofIdx = scv.dofIndex();
                const auto& cellCenterPriVars = curSol[FVGridGeometry::cellCenterIdx()][dofIdx];
                PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);
                auto elemSol = elementSolution<typename FVGridGeometry::LocalView>(std::move(priVars));
                // NOTE: first update the turbulence quantities
                storedDissipation_[elementIdx] = elemSol[0][Indices::dissipationEqIdx];
                storedTurbulentKineticEnergy_[elementIdx] = elemSol[0][Indices::turbulentKineticEnergyEqIdx];
                // NOTE: then update the volVars
                VolumeVariables volVars;
                volVars.update(elemSol, asImp_(), element, scv);
                storedDynamicEddyViscosity_[elementIdx] = volVars.calculateEddyViscosity();
                storedDensity_[elementIdx] = volVars.density();
            }
        }

        // get matching point for k-epsilon wall function
        unsigned int numElementsInNearWallRegion = 0;
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
            unsigned int wallNormalAxis = asImp_().wallNormalAxis_[elementIdx];
            unsigned int neighborIdx0 = asImp_().neighborIdx_[elementIdx][wallNormalAxis][0];
            unsigned int neighborIdx1 = asImp_().neighborIdx_[elementIdx][wallNormalAxis][1];
            numElementsInNearWallRegion = inNearWallRegion(elementIdx)
                                          ? numElementsInNearWallRegion + 1
                                          : numElementsInNearWallRegion + 0;
            if ((!inNearWallRegion(elementIdx) && (inNearWallRegion(neighborIdx0) || inNearWallRegion(neighborIdx1)))
                || (!inNearWallRegion(elementIdx) && elementIdx == asImp_().wallElementIdx_[elementIdx])
                || (inNearWallRegion(elementIdx) && (asImp_().wallElementIdx_[neighborIdx0] != asImp_().wallElementIdx_[neighborIdx1])))
            {
                matchingPointIdx_[asImp_().wallElementIdx_[elementIdx]] = elementIdx;
            }
        }
        std::cout << "numElementsInNearWallRegion: " << numElementsInNearWallRegion << std::endl;

        // calculate the potential zeroeq eddy viscosities for two-layer model
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
            zeroEqDynamicEddyViscosity_[elementIdx] = zeroEqEddyViscosityModel(elementIdx);
        }

        // then make them match at the matching point
        static const auto enableZeroEqScaling
            = getParamFromGroup<bool>(this->paramGroup(), "KEpsilon.EnableZeroEqScaling", true);
        if (enableZeroEqScaling)
        {
            for (const auto& element : elements(this->fvGridGeometry().gridView()))
            {
                unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
                unsigned int matchingPointIdx = matchingPointIdx_[asImp_().wallElementIdx_[elementIdx]];

                Scalar scalingFactor = storedDynamicEddyViscosity_[matchingPointIdx]
                                       / zeroEqDynamicEddyViscosity_[matchingPointIdx];
                if (!isMatchingPoint(elementIdx)
                    && !std::isnan(scalingFactor) && !std::isinf(scalingFactor))
                {
                    zeroEqDynamicEddyViscosity_[elementIdx] *= scalingFactor;
                }
            }
            for (const auto& element : elements(this->fvGridGeometry().gridView()))
            {
                unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
                unsigned int matchingPointIdx = matchingPointIdx_[asImp_().wallElementIdx_[elementIdx]];
                if (isMatchingPoint(elementIdx))
                {
                    zeroEqDynamicEddyViscosity_[matchingPointIdx] = storedDynamicEddyViscosity_[matchingPointIdx];
                }
            }
        }
    }

    /*!
     * \brief Returns if an element is located in the near-wall region
     */
    const bool inNearWallRegion(unsigned int elementIdx) const
    {
        unsigned int wallElementIdx = asImp_().wallElementIdx_[elementIdx];
        unsigned int matchingPointIdx = matchingPointIdx_[wallElementIdx];
        return wallElementIdx == matchingPointIdx ? yPlusNominal(elementIdx) < yPlusThreshold_
                                                : yPlus(elementIdx) < yPlusThreshold_;
    }

    /*!
     * \brief Returns if an element is the matching point
     */
    const bool isMatchingPoint(unsigned int elementIdx) const
    { return matchingPointIdx_[asImp_().wallElementIdx_[elementIdx]] == elementIdx; }

    /*!
     * \brief Returns the \f$ y^+ \f$ value at an element center
     */
    const Scalar yPlus(unsigned int elementIdx) const
    {
        return asImp_().wallDistance_[elementIdx] * uStar(elementIdx)
               / asImp_().kinematicViscosity_[elementIdx];
    }
    /*!
     * \brief Returns the nominal \f$ y^+ \f$ value at an element center
     */
    const Scalar yPlusNominal(unsigned int elementIdx) const
    {
        return asImp_().wallDistance_[elementIdx] * uStarNominal(elementIdx)
               / asImp_().kinematicViscosity_[elementIdx];
    }

    /*!
     * \brief Returns the kinematic eddy viscosity of a 0-Eq. model
     */
    const Scalar zeroEqEddyViscosityModel(unsigned int elementIdx) const
    {
        using std::abs;
        using std::exp;
        using std::sqrt;

        // use VanDriest's model
        Scalar yPlusValue = yPlus(elementIdx);
        Scalar mixingLength = 0.0;
        if (yPlusValue > 0.0)
        {
            mixingLength = asImp_().karmanConstant() * asImp_().wallDistance_[elementIdx]
                           * (1.0 - exp(-yPlusValue / 26.0 ))
                           / sqrt(1.0 - exp(-0.26 * yPlusValue));
        }

        unsigned int wallNormalAxis = asImp_().wallNormalAxis_[elementIdx];
        unsigned int flowNormalAxis = asImp_().flowNormalAxis_[elementIdx];
        Scalar velocityGradient = asImp_().velocityGradients_[elementIdx][flowNormalAxis][wallNormalAxis];
        return mixingLength * mixingLength * abs(velocityGradient) * storedDensity_[elementIdx];
    }

    //! \brief Returns the wall shear stress velocity
    const Scalar uStar(unsigned int elementIdx) const
    {
        using std::abs;
        using std::sqrt;
        unsigned int wallElementIdx = asImp_().wallElementIdx_[elementIdx];
        unsigned int wallNormalAxis = asImp_().wallNormalAxis_[elementIdx];
        unsigned int flowNormalAxis = asImp_().flowNormalAxis_[elementIdx];
        return sqrt(asImp_().kinematicViscosity_[wallElementIdx]
                    * abs(asImp_().velocityGradients_[wallElementIdx][flowNormalAxis][wallNormalAxis]));
    }

    //! \brief Returns the nominal wall shear stress velocity (accounts for poor approximation of viscous sublayer)
    const Scalar uStarNominal(unsigned int elementIdx) const
    {
        using std::pow;
        using std::sqrt;
        unsigned int matchingPointIdx = matchingPointIdx_[asImp_().wallElementIdx_[elementIdx]];
        return pow(cMu(), 0.25)
               * sqrt(storedTurbulentKineticEnergy_[matchingPointIdx]);
    }

    /*!
     * \brief Returns the dissipation calculated from the wall function consideration
     */
    const Scalar dissipationWallFunction(unsigned int elementIdx) const
    {
        return uStarNominal(elementIdx) * uStarNominal(elementIdx) * uStarNominal(elementIdx)
               / asImp_().karmanConstant() / asImp_().wallDistance_[elementIdx];
    }

    /*!
     * \brief Returns the turbulentKineticEnergy calculated from the wall function consideration
     */
    const Scalar turbulentKineticEnergyWallFunction(unsigned int elementIdx) const
    {
        unsigned int wallElementIdx = asImp_().wallElementIdx_[elementIdx];
        unsigned int matchingPointIdx = matchingPointIdx_[wallElementIdx];
        return storedTurbulentKineticEnergy_[matchingPointIdx];
    }

    //! \brief Returns the nominal wall shear stress (accounts for poor approximation of viscous sublayer)
    const Scalar tangentialMomentumWallFunction(unsigned int elementIdx, Scalar velocity) const
    {
        using std::log;
        Scalar velocityNominal = uStarNominal(elementIdx) * (1.0 / asImp_().karmanConstant() * log(yPlusNominal(elementIdx)) + 5.0);
        return uStarNominal(elementIdx) * uStarNominal(elementIdx)
               * velocity / velocityNominal;
    }

    //! \brief Checks whether a wall function should be used
    bool useWallFunction(const Element& element,
                         const SubControlVolumeFace& localSubFace,
                         const int& eqIdx) const
    {
        unsigned int elementIdx = asImp_().fvGridGeometry().elementMapper().index(element);
        auto bcTypes = asImp_().boundaryTypes(element, localSubFace);
        return asImp_().isOnWall(localSubFace)
               && bcTypes.isDirichlet(eqIdx)
               && isMatchingPoint(elementIdx);
    }

    //! \brief Returns an additional wall function momentum flux (only needed for RANS models)
    FacePrimaryVariables wallFunction(const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVolumeVariables& elemVolVars,
                                      const ElementFaceVariables& elemFaceVars,
                                      const SubControlVolumeFace& scvf,
                                      const SubControlVolumeFace& localSubFace) const
    {
        unsigned int elementIdx = asImp_().fvGridGeometry().elementMapper().index(element);
        return FacePrimaryVariables(asImp_().tangentialMomentumWallFunction(elementIdx, elemFaceVars[scvf].velocitySelf())
                                    * elemVolVars[scvf.insideScvIdx()].density());
    }

    //! \brief Returns the flux for non-isothermal and compositional RANS models
    template<bool eB = enableEnergyBalance, bool compositional = isCompositional,
             typename std::enable_if_t<eB && compositional, int> = 0>
    CellCenterPrimaryVariables wallFunction(const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const ElementFaceVariables& elemFaceVars,
                                            const SubControlVolumeFace& scvf) const
    {
        return wallFunctionComponent(element, fvGeometry, elemVolVars, elemFaceVars, scvf)
               + wallFunctionEnergy(element, fvGeometry, elemVolVars, elemFaceVars, scvf);
    }

    //! \brief Returns the flux for isothermal and compositional RANS models
    template<bool eB = enableEnergyBalance, bool compositional = isCompositional,
             typename std::enable_if_t<!eB && compositional, int> = 0>
    CellCenterPrimaryVariables wallFunction(const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const ElementFaceVariables& elemFaceVars,
                                            const SubControlVolumeFace& scvf) const
    { return wallFunctionComponent(element, fvGeometry, elemVolVars, elemFaceVars, scvf); }

    //! \brief Returns the flux for non-isothermal RANS models
    template<bool eB = enableEnergyBalance, bool compositional = isCompositional,
             typename std::enable_if_t<eB && !compositional, int> = 0>
    CellCenterPrimaryVariables wallFunction(const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const ElementFaceVariables& elemFaceVars,
                                            const SubControlVolumeFace& scvf) const
    { return wallFunctionEnergy(element, fvGeometry, elemVolVars, elemFaceVars, scvf); }

    //! \brief Returns the flux for isothermal RANS models
    template<bool eB = enableEnergyBalance, bool compositional = isCompositional,
             typename std::enable_if_t<!eB && !compositional, int> = 0>
    CellCenterPrimaryVariables wallFunction(const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const ElementFaceVariables& elemFaceVars,
                                            const SubControlVolumeFace& scvf) const
    { return CellCenterPrimaryVariables(0.0); }

    //! \brief Returns the component wall-function flux
    CellCenterPrimaryVariables wallFunctionComponent(const Element& element,
                                                     const FVElementGeometry& fvGeometry,
                                                     const ElementVolumeVariables& elemVolVars,
                                                     const ElementFaceVariables& elemFaceVars,
                                                     const SubControlVolumeFace& scvf) const
    {
        using std::log;
        auto wallFunctionFlux = CellCenterPrimaryVariables(0.0);
        unsigned int elementIdx = asImp_().fvGridGeometry().elementMapper().index(element);

        // component mass fluxes
        for (int compIdx = 0; compIdx < ModelTraits::numComponents(); ++compIdx)
        {
            if (ModelTraits::replaceCompEqIdx() == compIdx)
                continue;

            Scalar schmidtNumber = elemVolVars[scvf.insideScvIdx()].kinematicViscosity()
                                   / elemVolVars[scvf.insideScvIdx()].diffusionCoefficient(0, compIdx);
            Scalar moleToMassConversionFactor = ModelTraits::useMoles()
                                                ? 1.0 : FluidSystem::molarMass(compIdx);
            wallFunctionFlux[compIdx] +=
                -1.0 * (asImp_().dirichlet(element, scvf)[Indices::conti0EqIdx + compIdx]
                        - elemVolVars[scvf.insideScvIdx()].moleFraction(compIdx))
                * elemVolVars[scvf.insideScvIdx()].molarDensity()
                * moleToMassConversionFactor
                * uStarNominal(elementIdx)
                / asImp_().turbulentSchmidtNumber()
                / (1. / asImp_().karmanConstant() * log(yPlusNominal(elementIdx) * 9.793)
                    + pFunction(schmidtNumber, asImp_().turbulentSchmidtNumber()));
        }

        if (ModelTraits::replaceCompEqIdx() < ModelTraits::numComponents())
        {
            wallFunctionFlux[ModelTraits::replaceCompEqIdx()] =
                -std::accumulate(wallFunctionFlux.begin(), wallFunctionFlux.end(), 0.0);
        }

        return wallFunctionFlux;
    }

    //! \brief Returns the energy wall-function flux
    CellCenterPrimaryVariables wallFunctionEnergy(const Element& element,
                                                  const FVElementGeometry& fvGeometry,
                                                  const ElementVolumeVariables& elemVolVars,
                                                  const ElementFaceVariables& elemFaceVars,
                                                  const SubControlVolumeFace& scvf) const
    {
        using std::log;
        auto wallFunctionFlux = CellCenterPrimaryVariables(0.0);
        unsigned int elementIdx = asImp_().fvGridGeometry().elementMapper().index(element);
        // energy fluxes
        Scalar prandtlNumber = elemVolVars[scvf.insideScvIdx()].kinematicViscosity()
                               * elemVolVars[scvf.insideScvIdx()].density()
                               * elemVolVars[scvf.insideScvIdx()].heatCapacity()
                               / elemVolVars[scvf.insideScvIdx()].thermalConductivity();
        wallFunctionFlux[Indices::energyEqIdx - cellCenterOffset] +=
            -1.0 * (asImp_().dirichlet(element, scvf)[Indices::temperatureIdx]
                    - elemVolVars[scvf.insideScvIdx()].temperature())
            * elemVolVars[scvf.insideScvIdx()].density()
            * elemVolVars[scvf.insideScvIdx()].heatCapacity()
            * uStarNominal(elementIdx)
            / asImp_().turbulentPrandtlNumber()
            / (1. / asImp_().karmanConstant() * log(yPlusNominal(elementIdx) * 9.793)
                + pFunction(prandtlNumber, asImp_().turbulentPrandtlNumber()));

        return wallFunctionFlux;
    }

    //! \brief Returns the value of the P-function after Jayatilleke \cite Versteeg2009a
    const Scalar pFunction(Scalar molecularNumber, Scalar turbulentNumber) const
    {
        using std::pow;
        using std::exp;
        return 9.24
               * (pow(molecularNumber / turbulentNumber, 0.75) - 1.0)
               * (1.0 + 0.28 * exp(-0.007 * molecularNumber / turbulentNumber));
    }

    //! \brief Returns the \$f C_{\mu} \$f constant
    const Scalar cMu() const
    {
        return 0.09;
    }

public:
    std::vector<unsigned int> matchingPointIdx_;
    std::vector<Scalar> storedDensity_;
    std::vector<Scalar> storedDissipation_;
    std::vector<Scalar> storedTurbulentKineticEnergy_;
    std::vector<Scalar> storedDynamicEddyViscosity_;
    std::vector<Scalar> zeroEqDynamicEddyViscosity_;
    bool useStoredEddyViscosity_;
    Scalar yPlusThreshold_;

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

}

#endif
