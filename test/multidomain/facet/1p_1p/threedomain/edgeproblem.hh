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
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief The problem for the (d-2)-dimensional edge domain in the single-phase
 *        facet coupling test involving three domains.
 */
#ifndef DUMUX_TEST_FACETCOUPLING_THREEDOMAIN_ONEP_EDGEPROBLEM_HH
#define DUMUX_TEST_FACETCOUPLING_THREEDOMAIN_ONEP_EDGEPROBLEM_HH

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include "spatialparams.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePEdgeProblem;

namespace Properties {
// create the type tag nodes
NEW_TYPE_TAG(OnePEdge, INHERITS_FROM(OneP));
NEW_TYPE_TAG(OnePEdgeTpfa, INHERITS_FROM(CCTpfaModel, OnePEdge));

// Set the grid type
SET_TYPE_PROP(OnePEdge, Grid, Dune::FoamGrid<1, 3>);
// Set the problem type
SET_TYPE_PROP(OnePEdge, Problem, OnePEdgeProblem<TypeTag>);
// set the spatial params
SET_TYPE_PROP(OnePEdge, SpatialParams, OnePSpatialParams<TypeTag>);

// the fluid system
SET_PROP(OnePEdge, FluidSystem)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
};

// Disable caching (for testing purposes)
SET_BOOL_PROP(OnePEdge, EnableGridVolumeVariablesCache, false);
SET_BOOL_PROP(OnePEdge, EnableGridFluxVariablesCache, false);
SET_BOOL_PROP(OnePEdge, EnableFVGridGeometryCache, false);

} // end namespace Properties
/*!
 * \ingroup OnePTests
 * \brief The (d-2)-dimensional test problem for the incompressible
 *        one-phase model with coupling across the bulk grid facets
 */
template<class TypeTag>
class OnePEdgeProblem : public PorousMediumFlowProblem<TypeTag>
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
    OnePEdgeProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                    std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                    const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    {
        const auto a = getParam<Scalar>("Extrusion.Aperture");
        exFactor_ = a*a;
    }

    //!Specifies the type of boundary condition on a boundary position
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (globalPos[1] > 1.0 - 1e-6)
            values.setAllDirichlet();
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

    //! Evaluate the Dirichlet boundary conditions for a boundary position
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return PrimaryVariables({1.0}); }

    //! Set the aperture squared as extrusion factor.
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    { return exFactor_; }

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
    Scalar exFactor_;
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
};

} // end namespace Dumux

#endif
