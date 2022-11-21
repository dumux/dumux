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
 * \ingroup PoromechanicsTests
 * \brief Definition of the spatial parameters for the poro-elastic problem.
 */
#ifndef DUMUX_POROELASTIC_SPATIAL_PARAMS_HH
#define DUMUX_POROELASTIC_SPATIAL_PARAMS_HH

#include <algorithm>

#include <dumux/geomechanics/lameparams.hh>
#include <dumux/geomechanics/poroelastic/fvspatialparams.hh>
#include <dumux/geomechanics/stressstate/stressdroplawparams.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>
#include <dumux/geomechanics/stressstate/constantcementmodel.hh>

#include <dumux/geomechanics/stressstate/stressdroplawparams.hh>
#include <dumux/geomechanics/stressstate/math/mohrspace.hh>
#include "regionindex.hh"
namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief Definition of the spatial parameters for the poro-elastic
 *        sub-problem in the coupled poro-mechanical el2p problem.
 */
template<class Scalar, class GridGeometry, class CouplingManager, class FluidSystem>
class PoroElasticSpatialParams : public FVPoroElasticSpatialParams< GridGeometry,
                                                                    Scalar,
                                                                    PoroElasticSpatialParams<Scalar, GridGeometry, CouplingManager,
                                                                    FluidSystem> >
{
    using ThisType = PoroElasticSpatialParams<Scalar, GridGeometry, CouplingManager, FluidSystem>;
    using ParentType = FVPoroElasticSpatialParams<GridGeometry, Scalar, ThisType>;

    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ConstantCementModel = ConstantCementModel<Scalar>;

    using StressDropLawParams = StressDropLawParams<Scalar,Line<Scalar>>;

public:
    //! Export the type of the lame parameters
    using LameParams = Dumux::LameParams<Scalar>;

    PoroElasticSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                             std::shared_ptr<CouplingManager> couplingManagerPtr)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManagerPtr)
    , rockModel_{makeModel("Rock")}
    , faultModel_{makeModel("Fault")}
    , stressDropLawParam_{makeStressDropLawParams()}
    {
        for (const auto& zoneType : AllZoneTypes)
        {
            initPorosity_[zoneType] = getParam<Scalar>(ZoneName.at(zoneType)+".Porosity");
        }

        hasFailure_.resize(gridGeometry->gridView().size(0));
        std::fill(hasFailure_.begin(), hasFailure_.end(),false);
    }

    //! Defines the Lame parameters.
    template<class ElemVolVars, class FluxVarsCache>
    const LameParams lameParams(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElemVolVars& elemVolVars,
                                 const FluxVarsCache& fluxVarsCach) const
    {
        const GlobalPosition& center = element.geometry().center();

        const auto zoneType = zoneTypeAtPos(center);
        const Scalar& initPorosity = initPorosity_.at(zoneType);

        const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
        const Scalar porosity = PorosityDeformation<Scalar>::evaluatePorosity(this->gridGeometry(), element, center, elemSol, initPorosity);

        if( zoneType == ZoneType::LeftFault || zoneType == ZoneType::RightFault)
            return faultModel_.effectiveLameModuli(porosity);
        else
            return rockModel_.effectiveLameModuli(porosity);

    }

    //! Returns the porosity of the porous medium.
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const

    {
        const Scalar initPorosity = initPorosity_.at(zoneTypeAtPos(scv.center()));
        return PorosityDeformation<Scalar>::evaluatePorosity(this->gridGeometry(), element, scv, elemSol, initPorosity);
    }

    /*!
     * \brief Returns the effective fluid density.
     */
    Scalar effectiveFluidDensity(const Element& element, const SubControlVolume& scv) const
    {
        // get porous medium flow volume variables from coupling manager
        const auto& pmFlowVolVars = couplingManager().getPMFlowVolVars(element);

        Scalar wPhaseDensity = pmFlowVolVars.density(FluidSystem::phase0Idx);
        Scalar nPhaseDensity = pmFlowVolVars.density(FluidSystem::phase1Idx);
        Scalar Sw = pmFlowVolVars.saturation(FluidSystem::phase0Idx);
        Scalar Sn = pmFlowVolVars.saturation(FluidSystem::phase1Idx);
        return (wPhaseDensity * Sw + nPhaseDensity * Sn);
    }

    /*!
     * \brief Returns the effective pore pressure.
     */
    template<class ElementVolumeVariables, class FluxVarsCache >
    Scalar effectivePorePressure(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const FluxVarsCache& fluxVarsCache) const
    {
        if (!initialized_)
        {
            return -9.81*1e3*element.geometry().center()[1];
        }
        else{
            // get porous medium flow volume variables from coupling manager
            const auto& pmFlowVolVars = couplingManager().getPMFlowVolVars(element);

            Scalar pw = pmFlowVolVars.pressure(FluidSystem::phase0Idx);
            Scalar pn = pmFlowVolVars.pressure(FluidSystem::phase1Idx);
            Scalar Sw = pmFlowVolVars.saturation(FluidSystem::phase0Idx);
            Scalar Sn = pmFlowVolVars.saturation(FluidSystem::phase1Idx);
            return (pw * Sw + pn * Sn);
        }
    }

    //! Returns the Biot coefficient of the porous medium.
    Scalar biotCoefficientAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    //! Returns the temperature in the domain at the given position
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 273.15 + 22.5 - 50.0/2000*(globalPos[1] + 500); }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    StressDropLawParams stressDropLawParam(const Element& element) const
    { return stressDropLawParam_; }

    void setFailure(const Element& element) const
    { std::cout << "shear failure detected." << std::endl;
      const auto eIdx = this->gridGeometry().elementMapper().index(element);
      hasFailure_[eIdx] = true;
    }

    const auto& getFailureState()
    { return hasFailure_; }

    void resetFailureState()
    {
        std::fill(hasFailure_.begin(), hasFailure_.end(),false);
    }

    bool hasFailure() const
    {
        return std::any_of(hasFailure_.begin(), hasFailure_.end(),
                           [&](const auto& val){return val;});
    }

    void isInitialized(const bool& state)
    {
        initialized_ = state;
    }

private:
    std::shared_ptr<const CouplingManager> couplingManager_;


    ConstantCementModel makeModel(const std::string& groupname)
    {
        Scalar phiCrit{getParam<Scalar>(groupname+".CriticalPorosity")};
        Scalar phib{getParam<Scalar>(groupname+".WellSortedPorosity")};

        Scalar Ks{getParam<Scalar>("MaterialParameters.GrainBulkModulus")};
        Scalar Gs{getParam<Scalar>("MaterialParameters.GrainShearModulus")};
        Scalar Kc{getParam<Scalar>("MaterialParameters.CementBulkModulus")};
        Scalar Gc{getParam<Scalar>("MaterialParameters.CementShearModulus")};

        return ConstantCementModel(phiCrit, phib, Ks, Gs, Kc, Gc);
    }

    StressDropLawParams makeStressDropLawParams()
    {
        Scalar phi = getParam<Scalar>("StressDropLaw.InternalFrictionAngle");
        Scalar sigma = getParam<Scalar>("StressDropLaw.Cohesion");
        Scalar stressDrop = getParam<Scalar>("StressDropLaw.StressDrop");
        return StressDropLawParams(phi, sigma, stressDrop);
    }

    ConstantCementModel rockModel_;
    ConstantCementModel faultModel_;
    std::vector<bool> isFracture_;
    std::map<ZoneType,Scalar> initPorosity_;
    StressDropLawParams stressDropLawParam_;
    mutable std::vector<bool> hasFailure_;
    bool initialized_ = false;
    Scalar eps_ = 3e-6;
};

} // end namespace Dumux

#endif
