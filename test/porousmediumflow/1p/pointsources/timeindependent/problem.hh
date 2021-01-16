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
 * \brief A test problem for the one-phase model:
 * Water is injected in one single point in the middle of the domain.
 */
#ifndef DUMUX_1P_SINGULARITY_PROBLEM_HH
#define DUMUX_1P_SINGULARITY_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/spatialparams/fv1pconstant.hh>

#ifndef GRIDTYPE // default to yasp grid if not provided by CMake
#define GRIDTYPE Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >
#endif

namespace Dumux {
template <class TypeTag>
class OnePSingularityProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct OnePSingularity { using InheritsFrom = std::tuple<OneP>; };
struct OnePSingularityBox { using InheritsFrom = std::tuple<OnePSingularity, BoxModel>; };
struct OnePSingularityCCTpfa { using InheritsFrom = std::tuple<OnePSingularity, CCTpfaModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePSingularity>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePSingularity> { using type = GRIDTYPE; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePSingularity> { using type = OnePSingularityProblem<TypeTag> ; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePSingularity>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FVSpatialParamsOnePConstant<GridGeometry, Scalar>;
};
} // end namespace Dumux

/*!
 * \ingroup OnePTests
 * \brief  Test problem for the one-phase model:
 * Water is injected in a single point in the middle of the domain.
 *
 * The domain is box shaped. All sides have Dirichlet boundary conditions.
 *
 * In the middle of the domain, water is injected simulating a small diameter well.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box1p_pointsources</tt> or
 * <tt>./test_cc1p_pointsources</tt>
 */
template <class TypeTag>
class OnePSingularityProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    // copy some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum { dimWorld = GridView::dimensionworld };
    enum {
        // index of the primary variable
        pressureIdx = Indices::pressureIdx
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    OnePSingularityProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    {
        return name_;
    }

    /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllDirichlet();
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
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
     * \brief Applies a vector of point sources which
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in untis
     * \f$ [ \textnormal{unit of conserved quantity} / s ] \f$.
     * Positive values mean that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    {
        // inject 10 kg/s water at origin;
        pointSources.push_back(PointSource(GlobalPosition(0.0), {10}));
    }

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars[pressureIdx] = 1.0e5;
        return priVars;
    }

    // \}

private:
    std::string name_;
};

} // end namespace Dumux

#endif
