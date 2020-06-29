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
 * \brief A test for internal Dirichlet constraints
 */

#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_INTERNAL_DIRICHLET_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_INTERNAL_DIRICHLET_HH

#include <bitset>
#include <dumux/common/boundarytypes.hh>
#include <test/porousmediumflow/1p/implicit/incompressible/problem.hh>

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePTestProblemInternalDirichlet;

namespace Properties {
// Create new type tags
namespace TTag {
struct OnePInternalDirichlet {};
struct OnePInternalDirichletTpfa { using InheritsFrom = std::tuple<OnePInternalDirichlet, OnePIncompressibleTpfa>; };
struct OnePInternalDirichletBox { using InheritsFrom = std::tuple<OnePInternalDirichlet, OnePIncompressibleBox>; };
} // end namespace TTag

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePInternalDirichlet>
{ using type = OnePTestProblemInternalDirichlet<TypeTag>; };

} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief A test for internal Dirichlet constraints
 */
template<class TypeTag>
class OnePTestProblemInternalDirichlet : public OnePTestProblem<TypeTag>
{
    using ParentType = OnePTestProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NeumannValues = GetPropType<TypeTag, Properties::NumEqVector>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using BoundaryTypes = Dumux::BoundaryTypes<ModelTraits::numEq()>;
    using Indices = typename ModelTraits::Indices;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto numEq = ModelTraits::numEq();

public:
    OnePTestProblemInternalDirichlet(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param globalPos The position of the boundary face's integration point in global coordinates
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would be the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    NeumannValues neumannAtPos(const GlobalPosition& globalPos) const
    {
        const auto& gg = this->gridGeometry();
        if (globalPos[0] < gg.bBoxMin()[0] + eps_)
            return NeumannValues(1e3);
        else if (globalPos[1] < gg.bBoxMin()[1] + eps_)
            return NeumannValues(-1e3);
        else
            return NeumannValues(0.0);
    }

    //! Enable internal Dirichlet constraints
    static constexpr bool enableInternalDirichletConstraints()
    { return true; }

    /*!
     * \brief Tag a degree of freedom to carry internal Dirichlet constraints.
     *        If true is returned for a dof, the equation for this dof is replaced
     *        by the constraint that its primary variable values must match the
     *        user-defined values obtained from the function internalDirichlet(),
     *        which must be defined in the problem.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     */
    std::bitset<numEq> hasInternalDirichletConstraint(const Element& element,
                                                      const SubControlVolume& scv) const
    {
        // the pure Neumann problem is only defined up to a constant
        // we create a well-posed problem by fixing the pressure at one dof in the middle of the domain
        std::bitset<numEq> values;
        if (scv.dofIndex() == static_cast<std::size_t>(this->gridGeometry().numDofs()/2))
            values.set(Indices::pressureIdx);
        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(1e5); }

private:
    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux

#endif
