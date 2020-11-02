// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \brief K-epsilon turbulence problem base class.
 */
#ifndef DUMUX_KEPSILON_PROBLEM_HH
#define DUMUX_KEPSILON_PROBLEM_HH

#include <numeric>

#include <dumux/common/properties.hh>
#include <dumux/common/staggeredfvproblem.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/rans/problem.hh>
#include <dumux/freeflow/turbulencemodel.hh>

#include "model.hh"

namespace Dumux {

/*!
 * \ingroup KEpsilonModel
 * \brief K-epsilon turbulence problem base class.
 *
 * This implements some base functionality for k-epsilon models.
 */
template<class TypeTag>
class RANSProblemImpl<TypeTag, TurbulenceModel::kepsilon> : public RANSProblemBase<TypeTag>
{
    using ParentType = RANSProblemBase<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;

    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FacePrimaryVariables = GetPropType<TypeTag, Properties::FacePrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    static constexpr bool enableEnergyBalance = ModelTraits::enableEnergyBalance();
    static constexpr bool isCompositional = ModelTraits::numFluidComponents() > 1;

    // account for the offset of the cell center privars within the PrimaryVariables container
    static constexpr auto cellCenterOffset = ModelTraits::numEq() - CellCenterPrimaryVariables::dimension;
    static_assert(cellCenterOffset == ModelTraits::dim(), "cellCenterOffset must equal dim for staggered NavierStokes");

public:

    //! The constructor sets the gravity, if desired by the user.
    RANSProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    { }

    /*!
     * \brief Correct size of the static (solution independent) wall variables
     */
    void updateStaticWallProperties()
    {
        if (!ParentType::isFlatWallBounded())
        {
            DUNE_THROW(Dune::NotImplemented, "\n Due to grid/geometric concerns, k-epsilon models should only be used for "
                                          << " wall bounded flows with flat channel geometries. "
                                          << "\n If your geometry is a flat channel, please set the runtime parameter RANS.IsFlatWallBounded to true. \n");
        }

        ParentType::updateStaticWallProperties();
        // update size and initial values of the global vectors
        matchingPointIdx_.resize(this->gridGeometry().elementMapper().size(), 0);
        storedDissipation_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        storedTurbulentKineticEnergy_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        storedDynamicEddyViscosity_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        zeroEqDynamicEddyViscosity_.resize(this->gridGeometry().elementMapper().size(), 0.0);
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
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);

            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                const int dofIdx = scv.dofIndex();
                const auto& cellCenterPriVars = curSol[GridGeometry::cellCenterIdx()][dofIdx];
                PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);
                auto elemSol = elementSolution<typename GridGeometry::LocalView>(std::move(priVars));
                // NOTE: first update the turbulence quantities
                storedDissipation_[elementIdx] = elemSol[0][Indices::dissipationEqIdx];
                storedTurbulentKineticEnergy_[elementIdx] = elemSol[0][Indices::turbulentKineticEnergyEqIdx];
                // NOTE: then update the volVars
                VolumeVariables volVars;
                volVars.update(elemSol, asImp_(), element, scv);
                storedDynamicEddyViscosity_[elementIdx] = volVars.calculateEddyViscosity();
            }
        }

        // get matching point for k-epsilon wall function
        unsigned int numElementsInNearWallRegion = 0;
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            unsigned int wallNormalAxis = asImp_().wallNormalAxis(elementIdx);
            unsigned int neighborIndex0 = asImp_().neighborIndex(elementIdx, wallNormalAxis, 0);
            unsigned int neighborIndex1 = asImp_().neighborIndex(elementIdx, wallNormalAxis, 1);
            numElementsInNearWallRegion = inNearWallRegion(elementIdx)
                                          ? numElementsInNearWallRegion + 1
                                          : numElementsInNearWallRegion + 0;
            if ((!inNearWallRegion(elementIdx) && (inNearWallRegion(neighborIndex0) || inNearWallRegion(neighborIndex1)))
                || (!inNearWallRegion(elementIdx) && elementIdx == asImp_().wallElementIndex(elementIdx))
                || (inNearWallRegion(elementIdx) && (asImp_().wallElementIndex(neighborIndex0) != asImp_().wallElementIndex(neighborIndex1))))
            {
                matchingPointIdx_[asImp_().wallElementIndex(elementIdx)] = elementIdx;
            }
        }
        std::cout << "numElementsInNearWallRegion: " << numElementsInNearWallRegion << std::endl;

        // calculate the potential zeroeq eddy viscosities for two-layer model
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            zeroEqDynamicEddyViscosity_[elementIdx] = zeroEqEddyViscosityModel(elementIdx);
        }

        // then make them match at the matching point
        static const auto enableZeroEqScaling
            = getParamFromGroup<bool>(this->paramGroup(), "KEpsilon.EnableZeroEqScaling", true);
        if (enableZeroEqScaling)
        {
            for (const auto& element : elements(this->gridGeometry().gridView()))
            {
                unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
                unsigned int matchingPointIndex = matchingPointIdx(asImp_().wallElementIndex(elementIdx));

                Scalar scalingFactor = storedDynamicEddyViscosity(matchingPointIndex)
                                       / zeroEqDynamicEddyViscosity_[matchingPointIndex];
                if (!isMatchingPoint(elementIdx)
                    && !std::isnan(scalingFactor) && !std::isinf(scalingFactor))
                {
                    zeroEqDynamicEddyViscosity_[elementIdx] *= scalingFactor;
                }
            }
            for (const auto& element : elements(this->gridGeometry().gridView()))
            {
                unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
                unsigned int matchingPointIndex = matchingPointIdx(asImp_().wallElementIndex(elementIdx));
                if (isMatchingPoint(elementIdx))
                {
                    zeroEqDynamicEddyViscosity_[matchingPointIndex] = storedDynamicEddyViscosity(matchingPointIndex);
                }
            }
        }
    }

    /*!
     * \brief Returns if an element is located in the near-wall region
     */
    bool inNearWallRegion(unsigned int elementIdx) const
    {
        unsigned int wallElementIdx = asImp_().wallElementIndex(elementIdx);
        unsigned int matchingPointIndex = matchingPointIdx(wallElementIdx);
        return (wallElementIdx == matchingPointIndex) ? yPlusNominal(elementIdx) < yPlusThreshold()
                                                      : yPlus(elementIdx) < yPlusThreshold();
    }

    /*!
     * \brief Returns if an element is the matching point
     */
    bool isMatchingPoint(unsigned int elementIdx) const
    { return matchingPointIdx(asImp_().wallElementIndex(elementIdx)) == elementIdx; }

    /*!
     * \brief Returns the \f$ y^+ \f$ value at an element center
     */
    const Scalar yPlus(unsigned int elementIdx) const
    {
        return asImp_().wallDistance(elementIdx) * uStar(elementIdx)
             / asImp_().kinematicViscosity(elementIdx);
    }
    /*!
     * \brief Returns the nominal \f$ y^+ \f$ value at an element center
     */
    const Scalar yPlusNominal(unsigned int elementIdx) const
    {
        return asImp_().wallDistance(elementIdx) * uStarNominal(elementIdx)
             / asImp_().kinematicViscosity(elementIdx);
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
            mixingLength = asImp_().karmanConstant() * asImp_().wallDistance(elementIdx)
                           * (1.0 - exp(-yPlusValue / 26.0 ))
                           / sqrt(1.0 - exp(-0.26 * yPlusValue));
        }

        unsigned int wallNormalAxis = asImp_().wallNormalAxis(elementIdx);
        unsigned int flowDirectionAxis = asImp_().flowDirectionAxis(elementIdx);
        Scalar velocityGradient = asImp_().velocityGradient(elementIdx, flowDirectionAxis, wallNormalAxis);
        return mixingLength * mixingLength * abs(velocityGradient) * asImp_().storedDensity(elementIdx);
    }

    //! \brief Returns the wall shear stress velocity
    const Scalar uStar(unsigned int elementIdx) const
    {
        using std::abs;
        using std::sqrt;
        unsigned int wallElementIdx = asImp_().wallElementIndex(elementIdx);
        unsigned int wallNormalAxis = asImp_().wallNormalAxis(elementIdx);
        unsigned int flowDirectionAxis = asImp_().flowDirectionAxis(elementIdx);
        return sqrt(asImp_().kinematicViscosity(wallElementIdx)
                    * abs(asImp_().velocityGradient(wallElementIdx, flowDirectionAxis, wallNormalAxis)));
    }

    //! \brief Returns the nominal wall shear stress velocity (accounts for poor approximation of viscous sublayer)
    const Scalar uStarNominal(unsigned int elementIdx) const
    {
        using std::pow;
        using std::sqrt;
        unsigned int matchingPointIndex = matchingPointIdx(asImp_().wallElementIndex(elementIdx));
        return pow(cMu(), 0.25) * sqrt(storedTurbulentKineticEnergy(matchingPointIndex));
    }

    /*!
     * \brief Returns the dissipation calculated from the wall function consideration
     */
    const Scalar dissipationWallFunction(unsigned int elementIdx) const
    {
        return uStarNominal(elementIdx) * uStarNominal(elementIdx) * uStarNominal(elementIdx)
               / asImp_().karmanConstant() / asImp_().wallDistance(elementIdx);
    }

    /*!
     * \brief Returns the turbulentKineticEnergy calculated from the wall function consideration
     */
    const Scalar turbulentKineticEnergyWallFunction(unsigned int elementIdx) const
    {
        unsigned int matchingPointIndex = matchingPointIdx(asImp_().wallElementIndex(elementIdx));
        return storedTurbulentKineticEnergy(matchingPointIndex);
    }

    //! \brief Returns the nominal wall shear stress (accounts for poor approximation of viscous sublayer)
    const Scalar tangentialMomentumWallFunction(unsigned int elementIdx, Scalar velocity) const
    {
        using std::log;
        Scalar velocityNominal = uStarNominal(elementIdx) * (1.0 / asImp_().karmanConstant() * log(yPlusNominal(elementIdx)) + 5.0);
        return uStarNominal(elementIdx) * uStarNominal(elementIdx) * velocity / velocityNominal;
    }

    //! \brief Checks whether a wall function should be used
    bool useWallFunction(const Element& element,
                         const SubControlVolumeFace& localSubFace,
                         const int& eqIdx) const
    {
        unsigned int elementIdx = asImp_().gridGeometry().elementMapper().index(element);
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
        unsigned int elementIdx = asImp_().gridGeometry().elementMapper().index(element);
        return FacePrimaryVariables(asImp_().tangentialMomentumWallFunction(elementIdx, elemFaceVars[scvf].velocitySelf())
                                    * asImp_().storedDensity(elementIdx) );
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
        unsigned int elementIdx = asImp_().gridGeometry().elementMapper().index(element);

        // component mass fluxes
        for (int compIdx = 0; compIdx < ModelTraits::numFluidComponents(); ++compIdx)
        {
            if (ModelTraits::replaceCompEqIdx() == compIdx)
                continue;

            Scalar schmidtNumber = elemVolVars[scvf.insideScvIdx()].kinematicViscosity()
                                   / elemVolVars[scvf.insideScvIdx()].diffusionCoefficient(0, 0, compIdx);
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

        if (ModelTraits::replaceCompEqIdx() < ModelTraits::numFluidComponents())
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
        unsigned int elementIdx = asImp_().gridGeometry().elementMapper().index(element);
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

    //! \brief Returns the \f$ C_{\mu} \f$ constant
    const Scalar cMu() const
    { return 0.09; }

    Scalar yPlusThreshold() const
    {
        static const Scalar yPlusThreshold = getParamFromGroup<Scalar>(this->paramGroup(), "KEpsilon.YPlusThreshold", 30);
        return yPlusThreshold;
    }

    bool useStoredEddyViscosity() const
    {
        static const bool useStoredEddyViscosity = getParamFromGroup<bool>(this->paramGroup(), "RANS.UseStoredEddyViscosity", false);
        return useStoredEddyViscosity;
    }

    Scalar storedDissipation(const int elementIdx) const
    { return storedDissipation_[elementIdx]; }

    Scalar storedTurbulentKineticEnergy(const int elementIdx) const
    { return storedTurbulentKineticEnergy_[elementIdx]; }

    Scalar storedDynamicEddyViscosity(const int elementIdx) const
    { return storedDynamicEddyViscosity_[elementIdx]; }

    Scalar zeroEqDynamicEddyViscosity(const int elementIdx) const
    { return zeroEqDynamicEddyViscosity_[elementIdx]; }

    unsigned int matchingPointIdx(const int elementIdx) const
    { return matchingPointIdx_[elementIdx]; }

private:
    std::vector<unsigned int> matchingPointIdx_;
    std::vector<Scalar> storedDissipation_;
    std::vector<Scalar> storedTurbulentKineticEnergy_;
    std::vector<Scalar> storedDynamicEddyViscosity_;
    std::vector<Scalar> zeroEqDynamicEddyViscosity_;

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace Dumux

#endif
