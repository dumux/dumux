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
 * \brief Definition of the spatial parameters for the two-phase flow
 *        sub-problem in the coupled poro-mechanical elp problem.
 */

#ifndef DUMUX_2P_SUB_PROBLEM_HH
#define DUMUX_2P_SUB_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/brineco2.hh>

#include "spatialparams_2p.hh"
#include "co2tables_el2p.hh"

namespace Dumux {

// forward declaration of the problem class
template <class TypeTag>
class TwoPSubProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct TwoPSub { using InheritsFrom = std::tuple<TwoP, CCTpfaModel>; };
} // end namespace TTag

// Set the fluid system for TwoPSubProblem
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TwoPSub>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::BrineCO2<Scalar, El2P::CO2Tables>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPSub> { using type = Dune::YaspGrid<3>; };
// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPSub> { using type = TwoPSubProblem<TypeTag> ; };
// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TwoPSub>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using type = TwoPSpatialParams<GridGeometry, Scalar, CouplingManager>;
};
} // end namespace Properties

/*!
 * \ingroup PoromechanicsTests
 * \brief The two-phase sub problem in the el2p coupled problem.
 */
template <class TypeTag>
class TwoPSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    // copy pressure index for convenience
    enum {
          pressureIdx = GetPropType<TypeTag, Properties::ModelTraits>::Indices::pressureIdx,
          saturationNIdx = GetPropType<TypeTag, Properties::ModelTraits>::Indices::saturationIdx,
          waterPhaseIdx = FluidSystem::phase0Idx,
          gasPhaseIdx = FluidSystem::phase1Idx,
          dimWorld = GridView::dimensionworld
    };

    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

public:
    TwoPSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                   std::shared_ptr<GetPropType<TypeTag, Properties::SpatialParams>> spatialParams,
                   const std::string& paramGroup = "TwoP")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    {
        FluidSystem::init();
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

    //! Returns the temperature within the domain in [K].
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    //! Evaluates the boundary conditions for a Dirichlet boundary segment.
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    { return initialAtPos(globalPos); }

    //! Evaluates the initial value for a control volume.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
      PrimaryVariables values;

      values[pressureIdx] = 1.5e7;
      values[saturationNIdx] = 0.0;
      return values;
    }

    //! Evaluates source terms.
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);

        static const Scalar sourceG = getParam<Scalar>("Problem.InjectionRateGas");
        static const Scalar sourceW = getParam<Scalar>("Problem.InjectionRateWater");
        if(globalPos[0] > 250 + eps_ && globalPos[0] < 750 - eps_
           && globalPos[1] > 250 + eps_ && globalPos[1] < 750 - eps_
           && globalPos[dimWorld-1] > 250 + eps_ && globalPos[dimWorld-1] < 750 - eps_)
        {
            values[gasPhaseIdx] = sourceG;
            values[waterPhaseIdx] = sourceW;
        }
        return values;
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[dimWorld-1] < eps_)
            values.setAllNeumann();
        else
            values.setAllDirichlet();
        return values;
    }

private:
    static constexpr Scalar eps_ = 1.0e-6;
    std::string problemName_;
};

} // end namespace Dumux

#endif
