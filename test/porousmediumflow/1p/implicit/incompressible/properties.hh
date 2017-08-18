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
 * \brief The properties for the incompressible test
 */
#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROPERTIES_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROPERTIES_HH

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/porousmediumflow/1p/implicit/propertydefaults.hh>

namespace Dumux
{

template<class TypeTag>
class MockProblem
{
    using ElementMapper = typename GET_PROP_TYPE(TypeTag, DofMapper);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
public:
    MockProblem(const GridView& gridView) : mapper_(gridView) {}

    const ElementMapper& elementMapper() const
    { return mapper_; }

    template<class Element, class Intersection>
    bool isInteriorBoundary(const Element& e, const Intersection& i) const
    { return false; }

    std::vector<unsigned int> getAdditionalDofDependencies(unsigned int index) const
    { return std::vector<unsigned int>(); }
private:
    ElementMapper mapper_;
};

namespace Properties
{

NEW_PROP_TAG(EnableFVGridGeometryCache);
NEW_PROP_TAG(FVGridGeometry);

NEW_TYPE_TAG(IncompressibleTestProblem, INHERITS_FROM(OneP, CCTpfaModel));

// Set the grid type
SET_TYPE_PROP(IncompressibleTestProblem, Grid, Dune::YaspGrid<2>);

// Set the finite volume grid geometry
SET_TYPE_PROP(IncompressibleTestProblem, FVGridGeometry, CCTpfaFVGridGeometry<TypeTag, false>);

// Set the problem type
SET_TYPE_PROP(IncompressibleTestProblem, Problem, Dumux::MockProblem<TypeTag>);

} // end namespace Properties

} // end namespace Dumux

#endif
