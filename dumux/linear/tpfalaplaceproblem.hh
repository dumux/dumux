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
 * \ingroup Linear
 * \brief Laplace problem for preconditioners
 */
#ifndef DUMUX_LINEAR_TPFA_LAPLACE_PROBLEM
#define DUMUX_LINEAR_TPFA_LAPLACE_PROBLEM

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief The spatial parameters class for the test problem using the
 *        incompressible 1p model
 */
template<class GridGeometry, class Scalar>
class OnePTestSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                             OnePTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using ThisType = OnePTestSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using PermeabilityType = Scalar;

    OnePTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        alpha_ = getParam<Scalar>("Problem.Alpha");
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param element The element
     * \param scv The sub-control volume
     * \param elemSol The element solution vector
     * \return The intrinsic permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        return alpha_;
    }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    /*!
     * \brief Define the temperature in the domain in Kelvin.
     * \param globalPos The global position
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 283.15; }

private:
    Scalar alpha_;

};

} // end namespace Dumux

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {
/*!
 * \ingroup OnePTests
 * \brief  Test problem for the incompressible one-phase model
 */
template<class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

public:
    OnePTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        beta_ = getParam<Scalar>("Problem.Beta");
        viscosity_ = getParam<Scalar>("Component.LiquidDynamicViscosity");
        useDirichlet_ = getParam<bool>("Problem.UseDirichletInLaplacian", false);

        // gradP is 1/m
        deltaP_ = (this->gridGeometry().bBoxMax()[0]-this->gridGeometry().bBoxMin()[0]);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if (useDirichlet_ && (isOutlet_(globalPos) || isInlet_(globalPos)))
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return { isInlet_(globalPos) ? deltaP_ : 0.0 };
    }

    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        // this creates the mass matrix contribution
        return { -beta_ / viscosity_ * elemVolVars[scv].pressure() };
    }

private:
    Scalar beta_, viscosity_;
    bool useDirichlet_;

    bool isInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    static constexpr Scalar eps_ = 1e-8;
    Scalar deltaP_;
};

} // end namespace Dumux

#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dune/grid/uggrid.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
template<class T> struct OnePIncompressible { using InheritsFrom = std::tuple<OneP>; };
template<class T> struct OnePIncompressibleTpfa { using InheritsFrom = std::tuple<OnePIncompressible<T>, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag, class T>
struct Grid<TypeTag, TTag::OnePIncompressible<T>>
{
    using type = typename T::Grid;
};

// Set the problem type
template<class TypeTag, class T>
struct Problem<TypeTag, TTag::OnePIncompressible<T>>
{ using type = OnePTestProblem<TypeTag>; };

// set the spatial params
template<class TypeTag, class T>
struct SpatialParams<TypeTag, TTag::OnePIncompressible<T>>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePTestSpatialParams<GridGeometry, Scalar>;
};

// use the incompressible local residual (provides analytic jacobian)
template<class TypeTag, class T>
struct LocalResidual<TypeTag, TTag::OnePIncompressible<T>>
{ using type = OnePIncompressibleLocalResidual<TypeTag>; };

// the fluid system
template<class TypeTag, class T>
struct FluidSystem<TypeTag, TTag::OnePIncompressible<T>>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<0, Scalar> >;
};

// Enable caching
template<class TypeTag, class T>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePIncompressible<T>> { static constexpr bool value = true; };
template<class TypeTag, class T>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePIncompressible<T>> { static constexpr bool value = true; };
template<class TypeTag, class T>
struct EnableGridGeometryCache<TypeTag, TTag::OnePIncompressible<T>> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#include <dumux/assembly/fvassembler.hh>

namespace Dumux {

template<class GridView>
auto makeTpfaLaplaceMatrix(const GridView& leafGridView)
{
    using TypeTag = Properties::TTag::OnePIncompressibleTpfa<GridView>;

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the grid variables
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    SolutionVector x(gridGeometry->numDofs());
    gridVariables->init(x);

    // create assembler & linear solver
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using A = typename Assembler::JacobianMatrix;
    using V = typename Assembler::ResidualType;
    auto jacobian = std::make_shared<A>();
    auto residual = std::make_shared<V>();
    assembler->setLinearSystem(jacobian, residual);
    assembler->assembleJacobianAndResidual(x);
    return jacobian;
}

template<class LinearOperator, class GridView>
std::shared_ptr<LinearOperator> makeTpfaLaplaceOperator(const GridView& leafGridView)
{
    return std::make_shared<LinearOperator>(makeTpfaLaplaceMatrix(leafGridView));
}

template<class LinearOperator, class GridView, class Comm>
std::shared_ptr<LinearOperator> makeTpfaLaplaceOperator(const GridView& leafGridView, const Comm& comm)
{
    return std::make_shared<LinearOperator>(makeTpfaLaplaceMatrix(leafGridView), comm);
}

} // end namespace Dumux

#endif
