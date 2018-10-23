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
 * \brief The problem for the lower-dimensional domain in the single-phase two-components facet coupling test
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEPNC_LOWDIMPROBLEM_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEPNC_LOWDIMPROBLEM_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include "spatialparams.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePNCLowDimProblem;

namespace Properties {
// create the type tag nodes
NEW_TYPE_TAG(OnePNCLowDim, INHERITS_FROM(OnePNC));
NEW_TYPE_TAG(OnePNCNILowDim, INHERITS_FROM(OnePNCLowDim, OnePNCNI));
NEW_TYPE_TAG(OnePNCLowDimTpfa, INHERITS_FROM(CCTpfaModel, OnePNCLowDim));
NEW_TYPE_TAG(OnePNCLowDimBox, INHERITS_FROM(BoxModel, OnePNCLowDim));
NEW_TYPE_TAG(OnePNCNILowDimTpfa, INHERITS_FROM(CCTpfaModel, OnePNCNILowDim));
NEW_TYPE_TAG(OnePNCNILowDimBox, INHERITS_FROM(BoxModel, OnePNCNILowDim));

// Set the grid type
SET_TYPE_PROP(OnePNCLowDim, Grid, Dune::FoamGrid<1, 2>);
// Set the problem type
SET_TYPE_PROP(OnePNCLowDim, Problem, OnePNCLowDimProblem<TypeTag>);
// set the spatial params
SET_TYPE_PROP(OnePNCLowDim, SpatialParams, OnePSpatialParams< typename GET_PROP_TYPE(TypeTag, FVGridGeometry),
                                                              typename GET_PROP_TYPE(TypeTag, Scalar) >);

// the fluid system
SET_TYPE_PROP(OnePNCLowDim,
              FluidSystem,
              FluidSystems::OnePAdapter< FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), FluidSystems::H2ON2DefaultPolicy<true>> >);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(OnePNCLowDim, UseMoles, true);
} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief The lower-dimensional test problem for the incompressible
 *        one-phase model with coupling across the bulk grid facets
 */
template<class TypeTag>
class OnePNCLowDimProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    // copy some indices for convenience
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    enum
    {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,

        // component indices
        H2OIdx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx),
        N2Idx = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::N2Idx),
    };

    static constexpr bool isNonIsothermal = GET_PROP_TYPE(TypeTag, ModelTraits)::enableEnergyBalance();

public:
    OnePNCLowDimProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                        std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                        const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    , aperture_(getParam<Scalar>("Problem.FractureAperture"))
    {}

    //! Specifies the type of boundary condition at a given position
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    //! Evaluate the source term at a given position
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        // evaluate sources from bulk domain
        auto source = couplingManagerPtr_->evalSourcesFromBulk(element, fvGeometry, elemVolVars, scv);
        source /= scv.volume()*elemVolVars[scv].extrusionFactor();
        return source;
    }

    //! evaluates the Dirichlet boundary condition for a given position
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(initialAtPos(globalPos)); }

    //! Set the aperture as extrusion factor.
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    { return aperture_; }

    //! evaluate the initial conditions (isothermal case)
    template<bool ni = isNonIsothermal, std::enable_if_t<!ni, int> = 0>
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        values[pressureIdx] = 1.0e5;
        values[N2Idx] = 0.0;
        return values;
    }

    //! evaluate the initial conditions (non-isothermal case)
    template<bool ni = isNonIsothermal, std::enable_if_t<ni, int> = 0>
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        values[pressureIdx] = 1.0e5;
        values[N2Idx] = 0.0;
        values[Indices::temperatureIdx] = temperature();
        return values;
    }

    //! returns the temperature in \f$\mathrm{[K]}\f$ in the domain
    Scalar temperature() const
    { return 283.15; /*10°*/ }

    //! sets the pointer to the coupling manager.
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManagerPtr_ = cm; }

    //! returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    Scalar aperture_;
};

} // end namespace Dumux

#endif
