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
 * \brief The spatial parameters class for the bulk flow domain in the
 *        elastic single-phase facet coupling test.
 */
#ifndef DUMUX_ANALYTIC_CRACK_BULK_FLOW_SPATIALPARAMS_HH
#define DUMUX_ANALYTIC_CRACK_BULK_FLOW_SPATIALPARAMS_HH

#include <dumux/discretization/elementsolution.hh>

#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>

namespace Dumux {

/*!
 * \brief The spatial parameters for the bulk flow domain in the single-phase facet coupling test.
 */
template<class FVGridGeometry, class Scalar, class CouplingManager>
class OnePBulkSpatialParams : public FVSpatialParamsOneP<FVGridGeometry, Scalar,
                                                         OnePBulkSpatialParams<FVGridGeometry, Scalar, CouplingManager>>
{
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;

    using ThisType = OnePBulkSpatialParams<FVGridGeometry, Scalar, CouplingManager>;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar, ThisType>;

public:
    //! Export the type used for permeability
    using PermeabilityType = Scalar;

    OnePBulkSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                          std::shared_ptr<const CouplingManager> couplingManagerPtr,
                          const std::string& paramGroup = "")
    : ParentType(fvGridGeometry)
    , couplingManagerPtr_(couplingManagerPtr)
    , initialPorosity_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.InitialPorosity"))
    , initialPermeability_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.InitialPermeability"))
    , useConstantPorosity_(getParamFromGroup<bool>(paramGroup, "SpatialParams.UseConstantPorosity"))
    {}

    //! Function for defining the (intrinsic) permeability \f$[m^2]\f$.
    template<class ElementSolution>
    Scalar permeability(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    { return initialPermeability_; }

    //! Returns the porosity for a sub-control volume.
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        if (useConstantPorosity_)
            return initialPorosity_;

        static constexpr auto mechId = CouplingManager::poroMechId;
        static constexpr auto matrixFlowId = CouplingManager::matrixFlowId;

        const auto eIdx = this->fvGridGeometry().elementMapper().index(element);
        const auto poroMechElementIdx = couplingManagerPtr_->bulkIndexMap().map(matrixFlowId, eIdx);
        const auto& poroMechGridGeom = couplingManagerPtr_->problem(mechId).gridGeometry();
        const auto poroMechElement = poroMechGridGeom.element(poroMechElementIdx);

        auto poroMechElemSol = elementSolution(poroMechElement, couplingManagerPtr_->curSol()[mechId], poroMechGridGeom);
        // evaluate the deformation-dependent porosity at the scv center
        return PorosityDeformation<Scalar>::evaluatePorosity(poroMechGridGeom, poroMechElement, scv.center(), poroMechElemSol, initialPorosity_);
    }

private:
    std::shared_ptr<const CouplingManager> couplingManagerPtr_;
    Scalar initialPorosity_;
    Scalar initialPermeability_;
    bool useConstantPorosity_;
};

} // end namespace Dumux

#endif // DUMUX_TEST_FACETCOUPLING_ELONEP_ONEP_BULK_FLOW_SPATIALPARAMS_HH
