// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief The spatial parameters of the poroelastic sub-problem of Mandel's problem
 */
#ifndef DUMUX_TEST_MANDEL_POROELASTIC_SPATIALPARAMS_HH
#define DUMUX_TEST_MANDEL_POROELASTIC_SPATIALPARAMS_HH

#include <memory>

#include <dumux/solidmechanics/elastic/lameparams.hh>
#include <dumux/poromechanics/poroelastic/fvspatialparams.hh>
#include <dumux/material/fluidmatrixinteractions/porositylinearporoelasticity.hh>

#include "analyticalsolution/analyticalsolution.hh"
#include "volumetricstrain.hh"

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief The spatial parameters of the poroelastic sub-problem of Mandel's problem
 */
template<class Scalar, class GridGeometry, class CouplingManager>
class PoroElasticSpatialParams
: public FVPoroElasticSpatialParams<
    GridGeometry, Scalar, PoroElasticSpatialParams<Scalar, GridGeometry, CouplingManager>
>
{
    using ThisType = PoroElasticSpatialParams<Scalar, GridGeometry, CouplingManager>;
    using ParentType = FVPoroElasticSpatialParams<GridGeometry, Scalar, ThisType>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using AnalyticalSolution = MandelAnalyticalSolution<Scalar>;

public:
    using LameParams = Dumux::LameParams<Scalar>;

    PoroElasticSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                             std::shared_ptr<const CouplingManager> couplingManager,
                             std::shared_ptr<const AnalyticalSolution> analyticalSolution)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , analyticalSolution_(analyticalSolution)
    , porosityLaw_(analyticalSolution->params().initialPorosity(),
                   analyticalSolution->params().biotCoefficient(),
                   analyticalSolution->params().bulkModulus())
    {
        const auto& params = analyticalSolution->params();
        lameParams_.setMu(params.shearModulus());
        lameParams_.setLambda(params.bulkModulus() - 2.0/3.0*params.shearModulus());
    }

    const LameParams& lameParamsAtPos(const GlobalPosition& globalPos) const
    { return lameParams_; }

    Scalar biotCoefficientAtPos(const GlobalPosition& globalPos) const
    { return analyticalSolution_->params().biotCoefficient(); }

    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        static constexpr auto pmFlowId = CouplingManager::pmFlowId;
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        const Scalar pressure = couplingManager_->curSol(pmFlowId)[eIdx][0];

        const auto& pos = scv.center();
        const Scalar strainChange = volumetricStrain(element, this->gridGeometry(), elemSol, pos)
                                    - analyticalSolution_->initialVolumetricStrain(pos);
        const Scalar pressureChange = pressure - analyticalSolution_->initialPressure(pos);
        return porosityLaw_.evaluatePorosity(strainChange, pressureChange);
    }

    Scalar effectiveFluidDensity(const Element& element,
                                 const SubControlVolume& scv) const
    { return couplingManager_->getPMFlowVolVars(element).density(); }

    template<class ElementVolumeVariables, class FluxVarsCache>
    Scalar effectivePorePressure(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& elemVolVars,
                                 const FluxVarsCache& fluxVarsCache) const
    { return couplingManager_->getPMFlowVolVars(element).pressure(); }

    const AnalyticalSolution& analyticalSolution() const
    { return *analyticalSolution_; }

private:
    std::shared_ptr<const CouplingManager> couplingManager_;
    std::shared_ptr<const AnalyticalSolution> analyticalSolution_;
    PorosityLinearPoroElasticity<Scalar> porosityLaw_;
    LameParams lameParams_;
};

} // end namespace Dumux

#endif
