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

#include <dumux/geomechanics/lameparams.hh>
#include <dumux/geomechanics/poroelastic/fvspatialparams.hh>
#include <dumux/geomechanics/stressdroplawparams.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>
#include <dumux/material/fluidmatrixinteractions/constantcementmodel.hh>
namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief Definition of the spatial parameters for the poro-elastic
 *        sub-problem in the coupled poro-mechanical el2p problem.
 */
template<class Scalar, class GridGeometry>
class PoroElasticSpatialParams : public FVPoroElasticSpatialParams< GridGeometry,
                                                                    Scalar,
                                                                    PoroElasticSpatialParams<Scalar, GridGeometry> >
{
    using ThisType = PoroElasticSpatialParams<Scalar, GridGeometry>;
    using ParentType = FVPoroElasticSpatialParams<GridGeometry, Scalar, ThisType>;

    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ConstantCementModel = ConstantCementModel<Scalar>;
    using StressDropLawParams = StressDropLawParams<Scalar>;
    static constexpr int regionNum = 3;
    enum Region{
        Fault, Caprock, Reservoir
    };
    static constexpr Region AllRegions[] = {Fault, Caprock, Reservoir};

public:
    //! Export the type of the lame parameters
    using LameParams = Dumux::LameParams<Scalar>;

    PoroElasticSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                             std::shared_ptr<CouplingManager> couplingManagerPtr)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManagerPtr)
    , rockModel_{makeModel("Rock")}
    , fractureModel_{makeModel("Fracture")}
    , stressDropLawParams_{makeStressDropLawParams()}
    {
        for (const auto& region : AllRegions)
        {
            initPorosity_[region] = getParam<Scalar>("SpatialParams."+regionName[region]+".Porosity");
        }
    }

    //! Defines the Lame parameters.
    template<class ElemVolVars, class FluxVarsCache>
    const LameParams& lameParams(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElemVolVars& elemVolVars,
                                 const FluxVarsCache& fluxVarsCach) const
    {
        const auto region = regionAtPos(element.geometry().center());
        const Scalar initPorosity = initPorosity_[region];

        const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
        const GlobalPosition center = element.geometry().center();
        const Scalar porosity = PorosityDeformation<Scalar>::evaluatePorosity(this->gridGeometry(), element, center, elemSol, initPorosity);

        if( region == Region::Fault)
            return fractureModel_.effectiveLameModuli(porosity);
        else
            return rockModel_.effectiveLameModuli(porosity);

    }

    //! Returns the porosity of the porous medium.
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const

    {
        return 0.0; // no porosity needed
        // const Scalar initPorosity = initPorosity_[regionAtPos(scv.center())];
        // return PorosityDeformation<Scalar>::evaluatePorosity(this->gridGeometry(), element, scv, elemSol, initPorosity);
    }

    /*!
     * \brief Returns the effective fluid density.
     */
    Scalar effectiveFluidDensity(const Element& element, const SubControlVolume& scv) const
    {
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
        // get porous medium flow volume variables from coupling manager
        const auto& pmFlowVolVars = couplingManager().getPMFlowVolVars(element);

        Scalar pw = pmFlowVolVars.pressure(FluidSystem::phase0Idx);
        Scalar pn = pmFlowVolVars.pressure(FluidSystem::phase1Idx);
        Scalar Sw = pmFlowVolVars.saturation(FluidSystem::phase0Idx);
        Scalar Sn = pmFlowVolVars.saturation(FluidSystem::phase1Idx);
        return (pw * Sw + pn * Sn);
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

    const StressDropLawParams& stressDropLawParams() const
    { return stressDropLawParams_; }

    const StressDropLawParams makeStressDropLawParams()
    {
        Scalar phi = getParam<Scalar>("StressDropLaw.InternalFrictionAngle");
        Scalar sigma = getParam<Scalar>("StressDropLaw.Cohesion");
        Scalar stressDrop = getParam<Scalar>("StressDropLaw.StressDrop");
        return StressDropLawParams(phi, sigma, stressDrop);
    }
private:
    std::shared_ptr<const CouplingManager> couplingManager_;

    const std::array<std::string, regionNum> regionName = {"Fault", "Caprock", "Reservoir"};

    ConstantCementModel makeModel(const std::string& groupname)
    {
        Scalar phiCrit{getParam<Scalar>("MaterialParameters."+groupname+".CriticalPorosity")};
        Scalar phib{getParam<Scalar>("MaterialParameters."+groupname+".WellSortedPorosity")};

        Scalar Ks{getParam<Scalar>("MaterialParameters.GrainBulkModulus")};
        Scalar Gs{getParam<Scalar>("MaterialParameters.GrainShearModulus")};
        Scalar Kc{getParam<Scalar>("MaterialParameters.CementBulkModulus")};
        Scalar Gc{getParam<Scalar>("MaterialParameters.CementShearModulus")};

        return ConstantCementModel(phiCrit, phib, Ks, Gs, Kc, Gc);
    }

    Region regionAtPos(const GlobalPosition& globalPos) const
    {
        using std::tan;
        using std::abs;

        // check on Fault
        // Fault 1000m long
        // (x - (-1500))
        // each fault has width of 10 m
        static const Scalar width = 10;
        static const GlobalPosition leftFaultCenter = {500,-1500};
        static const GlobalPosition rightFaultCenter = {1500,-1500};

        if (abs(globalPos[1] - leftFaultCenter[1]) < 500 + eps_)
        {
            static const Scalar tan_80 = tan(80.0/180*M_PI);
            //std::cout << M_PI << std::endl;
            //std::cout << tan_80 << std::endl;
            Scalar localLeftFractureCenter = (globalPos[1] - leftFaultCenter[1]) /tan_80 + leftFaultCenter[0];

            // the left fault
            if (std::abs(localLeftFractureCenter - globalPos[0]) < width/2 + eps_)
                return Region::Fault;

            // the right fault
            if (std::abs(localLeftFractureCenter + rightFaultCenter[0] - leftFaultCenter[0] - globalPos[0]) < width/2 + eps_)
                return Region::Fault;
        }

        // the upper caprock
        if (globalPos[1] > -1450 -eps_ && globalPos[1] < -1300 + eps_ )
            return Region::Caprock;

        // the lower caprock
        if (globalPos[1] > -1700 -eps_ && globalPos[1] < -1550 + eps_ )
            return Region::Caprock;

        // else
        return Region::Reservoir;
    }
    ConstantCementModel rockModel_;
    ConstantCementModel fractureModel_;
    std::array<Scalar,regionNum> initPorosity_;
    LameParams lameParams_;
    StressDropLawParams stressDropLawParams_;
    Scalar eps_ = 3e-6;
};

} // end namespace Dumux

#endif
