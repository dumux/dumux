// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief The problem for the (d-1)-dimensional facet domain in the single-phase
 *        facet coupling test involving three domains.
 */
#ifndef DUMUX_TEST_FACETCOUPLING_THREEDOMAIN_ONEP_FACETPROBLEM_HH
#define DUMUX_TEST_FACETCOUPLING_THREEDOMAIN_ONEP_FACETPROBLEM_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include "spatialparams.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePFacetProblem;

namespace Properties {
// create the type tag nodes
NEW_TYPE_TAG(OnePFacet, INHERITS_FROM(OneP));
NEW_TYPE_TAG(OnePFacetTpfa, INHERITS_FROM(OnePFacet, CCTpfaFacetCouplingModel));

// Set the grid type
SET_TYPE_PROP(OnePFacet, Grid, Dune::FoamGrid<2, 3>);
// Set the problem type
SET_TYPE_PROP(OnePFacet, Problem, OnePFacetProblem<TypeTag>);
// set the spatial params
SET_TYPE_PROP(OnePFacet, SpatialParams, OnePSpatialParams<TypeTag>);

// the fluid system
SET_PROP(OnePFacet, FluidSystem)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
};

} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief The (d-1)-dimensional test problem for the incompressible
 *        one-phase model with coupling across the bulk grid facets
 */
template<class TypeTag>
class OnePFacetProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

public:
    //! The constructor
    OnePFacetProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                     std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                     const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    , aperture_(getParam<Scalar>("Extrusion.Aperture"))
    {}

    //! Specifies the kind of boundary condition at a boundary position
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Specifies which kind of interior boundary condition should be
     *        used for which equation on a given sub-control volume face
     *        that couples to a facet element.
     *
     * \param element The finite element the scvf is embedded in
     * \param scvf The sub-control volume face
     */
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     */
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        // forward to solution independent, fully-implicit specific interface
        auto source = couplingManagerPtr_->evalSourcesFromBulk(element, fvGeometry, elemVolVars, scv);
        source /= scv.volume()*elemVolVars[scv].extrusionFactor();
        return source;
    }

    //! Set the aperture as extrusion factor.
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    { return aperture_; }

    //! Evaluate the initial conditions
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(1.0); }

    //! Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
    Scalar temperature() const
    { return 283.15; /*10°*/ }

    //! Return const reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

    //! sets the pointer to the coupling manager.
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManagerPtr_ = cm; }

private:
    Scalar aperture_;
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
};

} // end namespace Dumux

#endif
