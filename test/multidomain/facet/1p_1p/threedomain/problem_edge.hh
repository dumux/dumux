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
 * \brief The problem for the (d-2)-dimensional edge domain in the single-phase
 *        facet coupling test involving three domains.
 */

#ifndef DUMUX_TEST_FACETCOUPLING_THREEDOMAIN_ONEP_EDGEPROBLEM_HH
#define DUMUX_TEST_FACETCOUPLING_THREEDOMAIN_ONEP_EDGEPROBLEM_HH

#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include "spatialparams.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePEdgeProblem;

namespace Properties {
// create the type tag nodes
// Create new type tags
namespace TTag {
struct OnePEdge { using InheritsFrom = std::tuple<OneP>; };
struct OnePEdgeTpfa { using InheritsFrom = std::tuple<OnePEdge, CCTpfaModel>; };
struct OnePEdgeMpfa { using InheritsFrom = std::tuple<OnePEdge, CCTpfaModel>; };
struct OnePEdgeBox { using InheritsFrom = std::tuple<OnePEdge, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePEdge> { using type = Dune::FoamGrid<1, 3>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePEdge> { using type = OnePEdgeProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePEdge>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePSpatialParams<GridGeometry, Scalar>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePEdge>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
};

} // end namespace Properties

/*!
 * \ingroup FacetTests
 * \brief The (d-2)-dimensional test problem for the incompressible
 *        one-phase model with coupling across the bulk grid facets.
 */
template<class TypeTag>
class OnePEdgeProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    OnePEdgeProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                    std::shared_ptr<CouplingManager> couplingManagerPtr,
                    const std::string& paramGroup = "Edge")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManagerPtr)
    {
        const auto a = getParam<Scalar>("Extrusion.Aperture");
        exFactor_ = a*a;
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

    //! Specifies the type of boundary condition on a boundary position.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub-control volume.
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

    //! Sets the aperture squared as extrusion factor.
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    { return exFactor_; }

    //! Evaluates the initial conditions.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(1.0); }

    //! Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
    Scalar temperature() const
    { return 283.15; /*10°*/ }

    //! Returns const reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    std::string problemName_;
    Scalar exFactor_;
};

} // end namespace Dumux

#endif
