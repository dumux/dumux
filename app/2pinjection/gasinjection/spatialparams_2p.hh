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
 * \brief The spatial parameters class for the two-phase sub problem in the el2p test problem.
 */

#ifndef DUMUX_2P_TEST_SPATIALPARAMS_HH
#define DUMUX_2P_TEST_SPATIALPARAMS_HH

#include <array>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include <dumux/material/gstatrandomfield.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>
#include <dumux/material/fluidmatrixinteractions/permeabilitykozenycarman.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief The spatial parameters class for the two-phase sub problem in the el2p test problem.
 */
template<class GridGeometry, class Scalar, class CouplingManager>
class TwoPSpatialParams : public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                                               TwoPSpatialParams<GridGeometry, Scalar, CouplingManager>>
{
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using ThisType = TwoPSpatialParams<GridGeometry, Scalar, CouplingManager>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;

    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;
    static constexpr int regionNum = 3;
    static constexpr int curveNum = 2;
    enum Region{
        Fault, Caprock, Reservoir
    };
    enum Curve{
        Shale, OtherLayers
    };
    static constexpr Region AllRegions[] = {Fault, Caprock, Reservoir};

public:
    // export permeability type
    using PermeabilityType = Scalar;

    TwoPSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::shared_ptr<CouplingManager> couplingManagerPtr)
    : ParentType(gridGeometry)
    , couplingManagerPtr_(couplingManagerPtr)
    {
        for (const auto& region : AllRegions)
        {
            initPermeability_[region] = getParam<Scalar>("SpatialParams."+regionName[region]+".Permeability");
            initPorosity_[region] = getParam<Scalar>("SpatialParams."+regionName[region]+".Porosity");
        }

        //pcKrSwCurve_.resize(curveNum);
        for(int curve = 0; curve < curveNum; ++curve) {
            pcKrSwCurve_.push_back(PcKrSwCurve("SpatialParams."+curveName[curve]));
        }

        onlyForOutput_.resize(this->gridGeometry().numDofs());
        for(const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bind(element);

            for( const auto& scv : scvs(fvGeometry))
            {
                onlyForOutput_[scv.dofIndex()] = regionAtPos(scv.center());
            }
        }

    }

    const auto& getOutPut()
    {return onlyForOutput_;}
    //! Returns the porosity for a sub-control volume.
    template< class ElementSolution >
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        static constexpr auto poroMechId = CouplingManager::poroMechId;

        const auto& poroMechGridGeom = couplingManagerPtr_->problem(poroMechId).gridGeometry();
        const auto& poroMechSolution = couplingManagerPtr_->getSolution(poroMechId);
        const auto poroMechElemSol = elementSolution(element, poroMechSolution, poroMechGridGeom);

        const auto region = regionAtPos(scv.center());
        const Scalar initPorosity = initPorosity_[region];
        // evaluate the deformation-dependent porosity at the scv center
        return PorosityDeformation<Scalar>::evaluatePorosity(poroMechGridGeom, element, scv.center(), poroMechElemSol, initPorosity);
    }

    //! Functions for defining the (intrinsic) permeability \f$[m^2]\f$.
     template< class ElementSolution >
     PermeabilityType permeability(const Element& element,
                                   const SubControlVolume& scv,
                                   const ElementSolution& elemSol) const
     {
         PermeabilityKozenyCarman<PermeabilityType> permLaw;
         const auto region = regionAtPos(scv.center());
         const Scalar initPermeability = initPermeability_[region];
         const Scalar initPorosity = initPorosity_[region];
         return permLaw.evaluatePermeability(initPermeability, initPorosity, porosity(element, scv, elemSol));
     }

     /*!
      * \brief Returns the parameters for the material law at a given location
      *
      * \param globalPos The global coordinates for the given location
      */
     auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
     {
         const auto region = regionAtPos(globalPos);
         if (region == Caprock)
            return makeFluidMatrixInteraction(pcKrSwCurve_[Shale]);
         else
            return makeFluidMatrixInteraction(pcKrSwCurve_[OtherLayers]);
     }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        return FluidSystem::phase0Idx;
    }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

    //! Returns the temperature in the domain at the given position
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 273.15 + 22.5 - 50.0/2000*(globalPos[1] + 500); }

private:
    const std::array<std::string, regionNum> regionName = {"Fault", "Caprock", "Reservoir"};
    const std::array<std::string, regionNum> curveName = {"Shale","OtherLayers"};

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

    std::shared_ptr<const CouplingManager> couplingManagerPtr_;
    std::array<Scalar,regionNum> initPermeability_;
    std::array<Scalar,regionNum> initPorosity_;
    std::vector<PcKrSwCurve> pcKrSwCurve_;
    Scalar eps_ = 1e-5;

    std::vector<int> onlyForOutput_;
};

} // end namespace Dumux

#endif
