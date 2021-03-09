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
 * \ingroup RichardsTests
 * \brief A water infiltration problem with a low-permeability lens
 *        embedded into a high-permeability domain which uses the
 *        Richards box model.
 */

#ifndef DUMUX_RICHARDS_LENSPROBLEM_HH
#define DUMUX_RICHARDS_LENSPROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "spatialparams.hh"

namespace Dumux {

template <class TypeTag>
class RichardsLensProblem;

// Specify the properties for the lens problem
namespace Properties {
// Create new type tags
namespace TTag {
struct RichardsLens { using InheritsFrom = std::tuple<Richards>; };
struct RichardsLensBox { using InheritsFrom = std::tuple<RichardsLens, BoxModel>; };
struct RichardsLensCC { using InheritsFrom = std::tuple<RichardsLens, CCTpfaModel>; };
struct RichardsLensCCMpfa { using InheritsFrom = std::tuple<RichardsLens, CCMpfaModel>; };
} // end namespace TTag

#ifndef GRIDTYPE
// Use 2d YaspGrid
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsLens> { using type = Dune::YaspGrid<2>; };
#else
// Use GRIDTYPE from CMakeLists.txt
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsLens> { using type = GRIDTYPE; };
#endif

// Set the physical problem to be solved
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsLens> { using type = RichardsLensProblem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsLens>
{
    using type = RichardsLensSpatialParams<GetPropType<TypeTag, Properties::GridGeometry>,
                                           GetPropType<TypeTag, Properties::Scalar>>;
};
} // end namespace Dumux

/*!
 * \ingroup RichardsTests
 *
 * \brief A water infiltration problem with a low-permeability lens
 *        embedded into a high-permeability domain which uses the
 *        Richards model.
 *
 * The domain is box shaped. Left and right boundaries are Dirichlet
 * boundaries with fixed water pressure (fixed Saturation \f$S_w = 0\f$),
 * bottom boundary is closed (Neumann 0 boundary), the top boundary
 * (Neumann 0 boundary) is also closed except for infiltration
 * section, where water is infiltrating into an initially unsaturated
 * porous medium. This problem is very similar to the LensProblem
 * which uses the TwoPBoxModel, with the main difference being that
 * the domain is initially fully saturated by gas instead of water and
 * water instead of a %DNAPL infiltrates from the top.
 *
 * This problem uses the \ref RichardsModel
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_boxrichards -parameterFile test_boxrichards.input -TimeManager.TEnd 10000000</tt>
 * <tt>./test_ccrichards -parameterFile test_ccrichards.input -TimeManager.TEnd 10000000</tt>
 *
 * where the initial time step is 100 seconds, and the end of the
 * simulation time is 10,000,000 seconds (115.7 days)
 */
template <class TypeTag>
class RichardsLensProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    enum {
        // copy some indices for convenience
        pressureIdx = Indices::pressureIdx,
        conti0EqIdx = Indices::conti0EqIdx,
        bothPhases = Indices::bothPhases,

        // world dimension
        dimWorld = GridView::dimensionworld
    };
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    RichardsLensProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Returns the temperature [K] within a finite volume
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; }; // -> 10°C

    /*!
     * \brief Returns the reference pressure [Pa] of the nonwetting
     *        fluid phase within a finite volume
     *
     * This problem assumes a constant reference pressure of 1 bar.
     */
    Scalar nonwettingReferencePressure() const
    { return 1.0e5; };

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the boundary type is set
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The position for which the Dirichlet value is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     *
     * \param globalPos The position for which the Neumann value is set
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);
        if (onInlet_(globalPos))
            values[conti0EqIdx] = -0.04; // kg/(m*s)
        return values;
    }

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     *
     * \param globalPos The position for which the boundary type is set
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return initial_(globalPos); };

    // \}

private:
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        const Scalar sw = 0.0;
        const Scalar pc = this->spatialParams().fluidMatrixInteractionAtPos(globalPos).pc(sw);
        values[pressureIdx] = nonwettingReferencePressure() - pc;
        values.setState(bothPhases);
        return values;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0];
        Scalar lambda = (this->gridGeometry().bBoxMax()[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda + eps_ && lambda < 2.0/3.0 + eps_;
    }

    static constexpr Scalar eps_ = 1.5e-7;

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;
    std::string name_;
};

} // end namespace Dumux

#endif
