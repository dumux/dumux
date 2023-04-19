// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief The spatial parameters class for the test problem using the
 *        1p box model.
 */

#ifndef DUMUX_1P_TEST_SPATIALPARAMS_HH
#define DUMUX_1P_TEST_SPATIALPARAMS_HH

#include <dumux/discretization/elementsolution.hh>

#include <dumux/porousmediumflow/fvspatialparams1p.hh>
#include <dumux/material/gstatrandomfield.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief The spatial parameters class for the test problem using the
 *        1p box model.
 */
template<class GridGeometry, class Scalar, class CouplingManager>
class OnePSpatialParams : public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                                                 OnePSpatialParams<GridGeometry, Scalar, CouplingManager>>
{
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using ThisType = OnePSpatialParams<GridGeometry, Scalar, CouplingManager>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    OnePSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::shared_ptr<CouplingManager> couplingManagerPtr)
    : ParentType(gridGeometry)
    , couplingManagerPtr_(couplingManagerPtr)
    , permeability_(getParam<Scalar>("SpatialParams.Permeability"))
    , initPorosity_(getParam<Scalar>("SpatialParams.InitialPorosity"))
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

        const auto& poroMechGridGeom = couplingManagerPtr_->problem(poroMechId).gridGeometry();
        const auto poroMechElemSol = elementSolution(element, couplingManagerPtr_->curSol(poroMechId), poroMechGridGeom);

        // evaluate the deformation-dependent porosity at the scv center
        return PorosityDeformation<Scalar>::evaluatePorosity(poroMechGridGeom, element, scv.center(), poroMechElemSol, initPorosity_);
    }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<const CouplingManager> couplingManagerPtr_;
    Scalar permeability_;
    Scalar initPorosity_;
};

} // end namespace Dumux

#endif
