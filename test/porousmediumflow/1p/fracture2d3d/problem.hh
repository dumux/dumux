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
 * \ingroup OnePTests
 * \brief A discrete fracture network embedded in an impermeable matrix.
 *
 * The fracture is a 2D network embedded in 3D.
 */

#ifndef DUMUX_ONEP_FRACTURE_TEST_PROBLEM_HH
#define DUMUX_ONEP_FRACTURE_TEST_PROBLEM_HH

#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>

#include "spatialparams.hh"

namespace Dumux {

template <class TypeTag>
class FractureProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct Fracture { using InheritsFrom = std::tuple<OneP>; };
struct FractureBox { using InheritsFrom = std::tuple<Fracture, BoxModel>; };
struct FractureCCTpfa { using InheritsFrom = std::tuple<Fracture, CCTpfaModel>; };
struct FractureCCMpfa { using InheritsFrom = std::tuple<Fracture, CCMpfaModel>; };
} // end namespace TTag

//! Enable caching (more memory, but faster runtime)
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::Fracture> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::Fracture> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::Fracture> { static constexpr bool value = true; };

//! The grid type
#if HAVE_DUNE_FOAMGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::Fracture> { using type = Dune::FoamGrid<2, 3>; };
#endif

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Fracture> { using type = Dumux::FractureProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Fracture>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};
// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::Fracture>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FractureSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief A discrete fracture network embedded in an impermeable matrix.
 *
 * The fracture is a 2D network embedded in 3D.
 * This problem uses the \ref OnePModel.
 */
template <class TypeTag>
class FractureProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum { dimWorld = GridView::dimensionworld };

    enum {
        // index of the primary variable
        pressureIdx = Indices::pressureIdx
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

public:
    FractureProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    {
        return name_;
    }

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * This problem assumes a uniform temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setAllNeumann();
        const auto& gg = this->gridGeometry();
        if (globalPos[0] > gg.bBoxMax()[0] - eps_ || globalPos[0] < gg.bBoxMin()[0] + eps_)
            values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initialAtPos(globalPos);
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{


    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        const auto& gg = this->gridGeometry();
        values[pressureIdx] = 1.0e5*(globalPos[0] - gg.bBoxMin()[0])/(gg.bBoxMax()[0] - gg.bBoxMin()[0]) + 1.0e5;
        return values;
    }
    // \}

private:
    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;
};

} // end namespace Dumux

#endif
