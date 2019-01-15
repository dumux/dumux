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
 * \ingroup FacetTests
 * \brief The spatial parameters class for the facet flow domain in the
 *        elastic single-phase facet coupling test.
 */
#ifndef DUMUX_TEST_FACETCOUPLING_ELONEP_ONEP_FACET_FLOW_SPATIALPARAMS_HH
#define DUMUX_TEST_FACETCOUPLING_ELONEP_ONEP_FACET_FLOW_SPATIALPARAMS_HH

#include <dumux/discretization/elementsolution.hh>
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief The spatial parameters for the single-phase facet coupling test.
 */
template<class FVGridGeometry, class Scalar, class CouplingManager>
class OnePFacetSpatialParams : public FVSpatialParamsOneP<FVGridGeometry, Scalar,
                                                          OnePFacetSpatialParams<FVGridGeometry, Scalar, CouplingManager>>
{
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;

    using ThisType = OnePFacetSpatialParams<FVGridGeometry, Scalar, CouplingManager>;
    using ParentType = FVSpatialParamsOneP<FVGridGeometry, Scalar, ThisType>;

public:
    //! Export the type used for permeability
    using PermeabilityType = Scalar;

    OnePFacetSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                           std::shared_ptr<CouplingManager> couplingManagerPtr,
                           const std::string& paramGroup = "")
    : ParentType(fvGridGeometry)
    , couplingManagerPtr_(couplingManagerPtr)
    , initialAperture_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.InitialAperture"))
    {}

    //! Function for defining the (intrinsic) permeability \f$[m^2]\f$.
    template<class ElementSolution>
    Scalar permeability(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolution& elemSol) const
    {
        const auto a = couplingManagerPtr_->computeAperture(element, scv, initialAperture_);
        return a*a/12.0;
    }

    //! Returns the porosity for a sub-control volume.
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

private:
    std::shared_ptr<const CouplingManager> couplingManagerPtr_;
    Scalar initialAperture_;
};

} // end namespace Dumux

#endif // DUMUX_TEST_FACETCOUPLING_ELONEP_ONEP_FACET_FLOW_SPATIALPARAMS_HH
