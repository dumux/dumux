// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Spatial parameters of the mechanical subdomain of the small-deformation Cryer benchmark.
 */
#ifndef DUMUX_TEST_POROELASTIC_CRYER_POROELASTIC_SPATIALPARAMS_HH
#define DUMUX_TEST_POROELASTIC_CRYER_POROELASTIC_SPATIALPARAMS_HH

#include <memory>

#include <dumux/common/parameters.hh>
#include <dumux/poromechanics/poroelastic/fvspatialparams.hh>
#include <dumux/solidmechanics/elastic/lameparams.hh>

#include "volumetricstrain.hh"

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief Spatial parameters of the mechanical subdomain: Lamé parameters from the drained bulk
 *        modulus \f$ K \f$ and the shear modulus \f$ G \f$, and the Biot coefficient.
 */
template<class GridDiscretization, class Scalar, class CouplingManager>
class CryerPoroElasticSpatialParams
: public FVPoroElasticSpatialParams<
    GridDiscretization, Scalar, CryerPoroElasticSpatialParams<GridDiscretization, Scalar, CouplingManager>
>
{
    using ThisType = CryerPoroElasticSpatialParams<GridDiscretization, Scalar, CouplingManager>;
    using ParentType = FVPoroElasticSpatialParams<GridDiscretization, Scalar, ThisType>;
    using SubControlVolume = typename GridDiscretization::SubControlVolume;
    using FVElementGeometry = typename GridDiscretization::LocalView;
    using Element = typename GridDiscretization::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using LameParams = Dumux::LameParams<Scalar>;

    CryerPoroElasticSpatialParams(std::shared_ptr<const GridDiscretization> gridDiscretization,
                                  std::shared_ptr<const CouplingManager> couplingManager)
    : ParentType(gridDiscretization)
    , couplingManager_(couplingManager)
    , porosity0_(getParam<Scalar>("SpatialParams.InitialPorosity"))
    , biotCoefficient_(getParam<Scalar>("SpatialParams.BiotCoefficient", 1.0))
    , storageCoefficient_(getParam<Scalar>("SpatialParams.StorageCoefficient", 0.0))
    {
        const auto shearModulus = getParam<Scalar>("SpatialParams.ShearModulus");
        const auto bulkModulus = getParam<Scalar>("SpatialParams.BulkModulus");
        lameParams_.setMu(shearModulus);
        lameParams_.setLambda(bulkModulus - 2.0/3.0*shearModulus);
    }

    const LameParams& lameParamsAtPos(const GlobalPosition&) const
    { return lameParams_; }

    Scalar biotCoefficientAtPos(const GlobalPosition&) const
    { return biotCoefficient_; }

    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        static constexpr auto pmFlowId = CouplingManager::pmFlowId;
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        const Scalar pressure = couplingManager_->curSol(pmFlowId)[eIdx][0];
        const Scalar divU = volumetricStrain(element, this->gridGeometry(), elemSol, scv.center());
        return porosity0_ + biotCoefficient_*divU + storageCoefficient_*pressure;
    }

    Scalar effectiveFluidDensity(const Element& element, const SubControlVolume&) const
    { return couplingManager_->getPMFlowVolVars(element).density(); }

    template<class ElementVolumeVariables, class FluxVarsCache>
    Scalar effectivePorePressure(const Element& element,
                                 const FVElementGeometry&,
                                 const ElementVolumeVariables&,
                                 const FluxVarsCache&) const
    { return couplingManager_->getPMFlowVolVars(element).pressure(); }

private:
    std::shared_ptr<const CouplingManager> couplingManager_;
    Scalar porosity0_;
    Scalar biotCoefficient_;
    Scalar storageCoefficient_;
    LameParams lameParams_;
};

} // end namespace Dumux

#endif
