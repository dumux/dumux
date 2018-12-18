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
 * \brief The spatial parameters class for the test problem using the
 *        1p box model.
 */

#ifndef DUMUX_1P_TEST_SPATIALPARAMS_HH
#define DUMUX_1P_TEST_SPATIALPARAMS_HH

#include <dumux/discretization/elementsolution.hh>

#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/material/spatialparams/gstatrandomfield.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief The spatial parameters class for the test problem using the
 *        1p box model.
 */
template<class FVGridGeometry, class Scalar, class CouplingManager>
class OnePSpatialParams : public FVSpatialParamsOneP<FVGridGeometry, Scalar,
                                                     OnePSpatialParams<FVGridGeometry, Scalar, CouplingManager>>
{
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using ThisType = OnePSpatialParams<FVGridGeometry, Scalar, CouplingManager>;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar, ThisType>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    OnePSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                      std::shared_ptr<CouplingManager> couplingManagerPtr)
    : ParentType(fvGridGeometry)
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

        const auto& poroMechGridGeom = couplingManagerPtr_->template problem<poroMechId>().fvGridGeometry();
        const auto poroMechElemSol = elementSolution(element, couplingManagerPtr_->curSol()[poroMechId], poroMechGridGeom);

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
