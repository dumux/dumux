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
 * \ingroup BoundaryTests
 * \brief The convergence test
 */

#ifndef DUMUX_CONVERGENCE_TEST_ONEP_PROBLEM_HH
#define DUMUX_CONVERGENCE_TEST_ONEP_PROBLEM_HH

#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/ccwmpfa.hh>
#include <dumux/discretization/cellcentered/wmpfa/methods.hh>
#include <dumux/common/boundarytypes.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "spatialparams.hh"

namespace Dumux {
template <class TypeTag>
class ConvergenceProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct OnePConvergence { using InheritsFrom = std::tuple<OneP, CCWMpfaModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePConvergence> { using type = Dumux::ConvergenceProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePConvergence>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::Constant<1, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePConvergence> { using type = GRIDTYPE; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePConvergence>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ConvergenceTestSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePConvergence> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePConvergence> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePConvergence> { static constexpr bool value = true; };

template<class TypeTag>
struct DiscretizationSubmethod<TypeTag, TTag::OnePConvergence> { static constexpr WMpfaMethod value = WMpfaMethod::nltpfa; };

} // end namespace Properties

/*!
 * \ingroup BoundaryTests
 * \brief The convergence test
 */
template <class TypeTag>
class ConvergenceProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto velocityXIdx = 0;
    static constexpr auto velocityYIdx = 1;
    static constexpr auto pressureIdx = 2;

public:
    //! export the Indices
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    ConvergenceProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain in [K].
     *
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
      * \brief Specifies which kind of boundary condition should be
      *        used for which equation on a given boundary control volume.
      *
      * \param element The element
      * \param scvf The boundary sub control volume face
      */
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;

        values.setAllDirichlet();

        return values;
    }

        /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param element The element for which the Dirichlet boundary condition is set
     * \param scvf The boundary subcontrolvolumeface
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        const auto p = analyticalSolution(scvf.center())[pressureIdx];
        return PrimaryVariables(p);
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub control volume.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        using std::exp; using std::sin; using std::cos;
        static constexpr Scalar omega = M_PI;
        static constexpr Scalar c = 0.9;
        const Scalar cosOmegaX = cos(omega*x);
        static const Scalar expTwo = exp(2);
        const Scalar expYPlusOne = exp(y+1);

        const Scalar result = (-(c*cosOmegaX + 1)*exp(y - 1)
                                             + 1.5*c*expYPlusOne*cosOmegaX
                                             + omega*omega*(expYPlusOne - expTwo + 2))
                                             *sin(omega*x);
        return NumEqVector(result);
    }

    // \}

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param element The element
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initial(const Element &element) const
    {
        return PrimaryVariables(0.0);
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     */
    auto analyticalSolution(const GlobalPosition& globalPos) const
    {
        Dune::FieldVector<Scalar, 3> sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        static constexpr Scalar omega = M_PI;
        static constexpr Scalar c = 0.9;
        using std::exp; using std::sin; using std::cos;
        const Scalar sinOmegaX = sin(omega*x);
        const Scalar cosOmegaX = cos(omega*x);
        static const Scalar expTwo = exp(2);
        const Scalar expYPlusOne = exp(y+1);

        sol[pressureIdx] = (expYPlusOne + 2 - expTwo)*sinOmegaX + 10.0;
        sol[velocityXIdx] = c/(2*omega)*expYPlusOne*sinOmegaX*sinOmegaX
                            -omega*(expYPlusOne + 2 - expTwo)*cosOmegaX;
        sol[velocityYIdx] = (0.5*c*(expYPlusOne + 2 - expTwo)*cosOmegaX
                            -(c*cosOmegaX + 1)*exp(y-1))*sinOmegaX;

        return sol;
    }

    // \}

private:
    static constexpr Scalar eps_ = 1e-7;
};
} // end namespace Dumux

#endif //DUMUX_DARCY_SUBPROBLEM_HH
