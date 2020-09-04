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

#ifndef DUMUX_TEST_FREEFLOW_PIPE_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_PIPE_PROBLEM_HH

#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

namespace Dumux {

template<class TypeTag>
class FreeFlowPipeProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct PipeFlow {};
struct PipeFlowMomentum { using InheritsFrom = std::tuple<PipeFlow, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct PipeFlowMass { using InheritsFrom = std::tuple<PipeFlow, NavierStokesMassOneP, CCTpfaModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PipeFlow>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::Constant<1, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PipeFlow>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PipeFlow>
{ using type = FreeFlowPipeProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::PipeFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::PipeFlow> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PipeFlow> { static constexpr bool value = true; };

// rotation-symmetric grid geometry forming a cylinder channel
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PipeFlowMomentum>
{
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;

    struct GGTraits : public FaceCenteredStaggeredDefaultGridGeometryTraits<GridView>
    { using Extrusion = RotationalExtrusion<0>; };

    using type = FaceCenteredStaggeredFVGridGeometry<GridView, enableCache, GGTraits>;
};

// rotation-symmetric grid geometry forming a cylinder channel
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PipeFlowMass>
{
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;

    struct GGTraits : public CCTpfaDefaultGridGeometryTraits<GridView>
    { using Extrusion = RotationalExtrusion<0>; };

    using type = CCTpfaFVGridGeometry<GridView, enableCache, GGTraits>;
};

} // end namespace Dumux::Properties

/*!
 * \brief Freeflow problem for pipe flow
 * Simulation of a radially-symmetric pipe flow with circular cross-section
 */
template <class TypeTag>
class FreeFlowPipeProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PrimaryVariables = typename ParentType::PrimaryVariables;
    using NumEqVector = typename ParentType::NumEqVector;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using Indices = typename ModelTraits::Indices;

    FreeFlowPipeProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)

    {
        name_ = getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        meanInletVelocity_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.MeanInletVelocity");
        mu_ = getParam<Scalar>("Component.LiquidKinematicViscosity")*getParam<Scalar>("Component.LiquidDensity");

        pipeRadius_ = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
        pipeLength_ = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];
        eps_ = 1e-7*pipeRadius_;

        std::cout << "-- Reynolds number: " << 2*pipeRadius_*meanInletVelocity_/getParam<Scalar>("Component.LiquidKinematicViscosity") << std::endl;
    }

    const std::string& name() const
    {
        return name_;
    }

    Scalar temperature() const
    { return 293.15; }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            if (onLowerBoundary_(globalPos) || onOuterBoundary_(globalPos))
                values.setAllDirichlet();
            else if (onInnerBoundary_(globalPos))
            {
                values.setDirichlet(Indices::velocityXIdx);
                values.setNeumann(Indices::momentumYBalanceIdx);
            }
            else
                values.setAllNeumann();
        }
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto& globalPos = scvf.ipGlobal();

        if constexpr (ParentType::isMomentumProblem())
        {
            if (onInnerBoundary_(globalPos))
                return values; // zero shear stress at symmetry axis
            if (onUpperBoundary_(globalPos))
                values = NavierStokesMomentumBoundaryFluxHelper::fixedPressureMomentumFlux(*this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, analyticalPressure(globalPos), true /*zeroNormalVelocityGradient*/);
        }
        else
        {
            using FluxHelper = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>;
            if (onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos))
                values = FluxHelper::scalarOutflowFlux(*this, element, fvGeometry, scvf, elemVolVars);
        }

        return values;
    }

    /*!
      * \brief Evaluates the initial value for a control volume.
      *
      * \param globalPos The global position
      */
     PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
     {
         return analyticalSolution(globalPos);
     }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return analyticalSolution(globalPos); }

    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);

        // paraboloid velocity profile
        if constexpr (ParentType::isMomentumProblem())
        {
            const auto r = globalPos[0] - this->gridGeometry().bBoxMin()[0];
            values[Indices::velocityXIdx] = 0.0;
            values[Indices::velocityYIdx] = 2.0*meanInletVelocity_*(1.0 - r*r/(pipeRadius_*pipeRadius_));
        }
        else
            values[Indices::pressureIdx] = analyticalPressure(globalPos);

        return values;
    }

    Scalar analyticalPressure(const GlobalPosition& globalPos) const
    {
        const auto y = globalPos[1];
        return (pipeLength_-y)*meanInletVelocity_*8.0*mu_/(pipeRadius_*pipeRadius_);
    }

private:
    bool onInnerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onOuterBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    std::string name_;
    Scalar meanInletVelocity_;
    Scalar mu_;
    Scalar pipeRadius_, pipeLength_;
    Scalar eps_;
};

} // end namespace Dumux

#endif
