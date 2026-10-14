// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief The spatial parameters of the flow sub-problem of Mandel's problem
 */
#ifndef DUMUX_TEST_MANDEL_ONEP_SPATIALPARAMS_HH
#define DUMUX_TEST_MANDEL_ONEP_SPATIALPARAMS_HH

#include <memory>

#include <dumux/common/parameters.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>
#include <dumux/material/fluidmatrixinteractions/porositylinearporoelasticity.hh>

#include "analyticalsolution/analyticalsolution.hh"
#include "volumetricstrain.hh"

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief The spatial parameters of the flow sub-problem of Mandel's problem
 *
 * The porosity changes with the volumetric strain and the pressure relative to the
 * undrained state at \f$t = 0\f$. As in linear poroelasticity, the fluid density is constant
 * and the compression \f$\phi_0 \Delta p / K_f\f$ of the pore fluid enters the storage term
 * through the porosity, such that the storage term is the fluid content of Biot's theory.
 */
template<class GridGeometry, class Scalar, class CouplingManager>
class OnePSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<
    GridGeometry, Scalar, OnePSpatialParams<GridGeometry, Scalar, CouplingManager>
>
{
    using ThisType = OnePSpatialParams<GridGeometry, Scalar, CouplingManager>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using AnalyticalSolution = MandelAnalyticalSolution<Scalar>;

public:
    using PermeabilityType = Scalar;

    OnePSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::shared_ptr<const CouplingManager> couplingManager,
                      std::shared_ptr<const AnalyticalSolution> analyticalSolution)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , analyticalSolution_(analyticalSolution)
    , porosityLaw_(analyticalSolution->params().initialPorosity(),
                   analyticalSolution->params().biotCoefficient(),
                   analyticalSolution->params().bulkModulus())
    , fluidStorageCoefficient_(analyticalSolution->params().initialPorosity()/analyticalSolution->params().fluidBulkModulus())
    , permeability_(getParam<Scalar>("Problem.Permeability"))
    {}

    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        static constexpr auto poroMechId = CouplingManager::poroMechId;
        const auto& poroMechGridGeometry = couplingManager_->problem(poroMechId).gridGeometry();
        const auto poroMechElemSol = elementSolution(element, couplingManager_->curSol(poroMechId), poroMechGridGeometry);

        const auto& pos = scv.center();
        const Scalar strainChange = volumetricStrain(element, poroMechGridGeometry, poroMechElemSol, pos)
                                    - analyticalSolution_->initialVolumetricStrain(pos);
        const Scalar pressureChange = elemSol[0][0] - analyticalSolution_->initialPressure(pos);
        return porosityLaw_.evaluatePorosity(strainChange, pressureChange) + fluidStorageCoefficient_*pressureChange;
    }

    const AnalyticalSolution& analyticalSolution() const
    { return *analyticalSolution_; }

private:
    std::shared_ptr<const CouplingManager> couplingManager_;
    std::shared_ptr<const AnalyticalSolution> analyticalSolution_;
    PorosityLinearPoroElasticity<Scalar> porosityLaw_;
    Scalar fluidStorageCoefficient_;
    Scalar permeability_;
};

} // end namespace Dumux

#endif
