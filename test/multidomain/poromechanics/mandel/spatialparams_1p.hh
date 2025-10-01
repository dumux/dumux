// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_TEST_MD_MANDEL_ONEP_SPATIALPARAMS_HH
#define DUMUX_TEST_MD_MANDEL_ONEP_SPATIALPARAMS_HH

#include <dumux/discretization/elementsolution.hh>

#include <dumux/porousmediumflow/fvspatialparams1p.hh>
#include <dumux/material/gstatrandomfield.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>

#include "analyticalsolution/analyticalsolution.hh"

namespace Dumux {

template<class GridGeometry, class Scalar, class CouplingManager>
class OnePSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<
    GridGeometry, Scalar, OnePSpatialParams<GridGeometry, Scalar, CouplingManager>
>
{
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using ThisType = OnePSpatialParams<GridGeometry, Scalar, CouplingManager>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;

    using AnalyticalSolution = MandelAnalyticalSolution<Scalar>;
public:
    // export permeability type
    using PermeabilityType = Scalar;

    OnePSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::shared_ptr<CouplingManager> couplingManagerPtr)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManagerPtr)
    , permeability_(getParam<Scalar>("Problem.Permeability"))
    , initPorosity_(getParam<Scalar>("Problem.InitialPorosity"))
    {}

    //! Returns the permeability at a given position.
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPoss) const
    { return permeability_; }

    //! Returns the porosity for a sub-control volume.
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        static constexpr auto poroMechId = CouplingManager::poroMechId;

        const auto& poroMechGridGeom = couplingManager_->problem(poroMechId).gridGeometry();

        // if (this->gridGeometry().elementMapper().index(element) == 0)
        //     {
        //         const auto& curSol = couplingManager_->curSol(poroMechId);
        //         std::cout << "Element: " << this->gridGeometry().elementMapper().index(element)
        //                   << ", scv: " << scv.dofIndex()
        //                   << ", curSol addr: " << &curSol << std::endl;

        //         const auto& prevSol = couplingManager_->problem(poroMechId).spatialParams().prevSolution();
        //         std::cout << " prevSol addr: " << &prevSol << std::endl;
        //     }

        const auto poroMechElemSol = elementSolution(element, couplingManager_->curSol(poroMechId), poroMechGridGeom);


        // evaluate the deformation-dependent porosity at the scv center
        const auto porosity = PorosityDeformation<Scalar>::evaluatePorosity(
            poroMechGridGeom,
            element,
            scv.center(),
            poroMechElemSol,
            analyticalSolution().initialDivU(scv.center()),
            initPorosity_,
            0.0/*min porosity*/
        );

        // if (this->gridGeometry().elementMapper().index(element) == 0)
        //     std::cout << "Element: " << this->gridGeometry().elementMapper().index(element) << ", scv: " << scv.dofIndex() << ", porosity: " << porosity << std::endl;
        return porosity;
    }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    const auto& analyticalSolution() const
    { return *analyticalSolution_; }

    void setAnalyticalSolution(const AnalyticalSolution& sol)
    { analyticalSolution_ = &sol; }

private:
    std::shared_ptr<const CouplingManager> couplingManager_;
    const AnalyticalSolution *analyticalSolution_;
    Scalar permeability_;
    Scalar initPorosity_;
};

} // end namespace Dumux

#endif
