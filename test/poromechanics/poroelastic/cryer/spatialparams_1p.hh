// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Spatial parameters of the flow subdomain of the small-deformation Cryer benchmark.
 */
#ifndef DUMUX_TEST_POROELASTIC_CRYER_ONEP_SPATIALPARAMS_HH
#define DUMUX_TEST_POROELASTIC_CRYER_ONEP_SPATIALPARAMS_HH

#include <memory>

#include <dumux/common/parameters.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

#include "volumetricstrain.hh"

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief Spatial parameters of the flow subdomain.
 *
 * The porosity \f$ \phi = \phi_0 + \alpha_B \nabla \cdot \mathbf{u} + S_p p \f$ gives the storage
 * \f$ \alpha_B \partial_t \nabla \cdot \mathbf{u} + S_p \partial_t p \f$ of linear poroelasticity.
 */
template<class GridDiscretization, class Scalar, class CouplingManager>
class CryerOnePSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<
    GridDiscretization, Scalar, CryerOnePSpatialParams<GridDiscretization, Scalar, CouplingManager>
>
{
    using ThisType = CryerOnePSpatialParams<GridDiscretization, Scalar, CouplingManager>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridDiscretization, Scalar, ThisType>;
    using SubControlVolume = typename GridDiscretization::SubControlVolume;
    using Element = typename GridDiscretization::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;

    CryerOnePSpatialParams(std::shared_ptr<const GridDiscretization> gridDiscretization,
                           std::shared_ptr<const CouplingManager> couplingManager)
    : ParentType(gridDiscretization)
    , couplingManager_(couplingManager)
    , porosity0_(getParam<Scalar>("SpatialParams.InitialPorosity"))
    , biotCoefficient_(getParam<Scalar>("SpatialParams.BiotCoefficient", 1.0))
    , storageCoefficient_(getParam<Scalar>("SpatialParams.StorageCoefficient", 0.0))
    , permeability_(getParam<Scalar>("SpatialParams.InitialPermeability"))
    {}

    PermeabilityType permeabilityAtPos(const GlobalPosition&) const
    { return permeability_; }

    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        static constexpr auto poroMechId = CouplingManager::poroMechId;
        const auto& poroMechGridDiscretization = couplingManager_->problem(poroMechId).gridGeometry();
        const auto& poroMechSol = couplingManager_->curSol(poroMechId);
        const auto poroMechElemSol = elementSolution(element, poroMechSol, poroMechGridDiscretization);
        const Scalar divU = volumetricStrain(element, poroMechGridDiscretization, poroMechElemSol, scv.center());
        return porosity0_ + biotCoefficient_*divU + storageCoefficient_*elemSol[0][0];
    }

private:
    std::shared_ptr<const CouplingManager> couplingManager_;
    Scalar porosity0_;
    Scalar biotCoefficient_;
    Scalar storageCoefficient_;
    Scalar permeability_;
};

} // end namespace Dumux

#endif
