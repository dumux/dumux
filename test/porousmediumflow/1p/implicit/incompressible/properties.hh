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

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/implicit/gridvariables.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/porousmediumflow/1p/implicit/propertydefaults.hh>

namespace Dumux
{

template<class TypeTag>
class MockSpatialParams
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimension>;

public:
    using PermeabilityType = Scalar;
    PermeabilityType permeability(const Element &element,
                        const SubControlVolume &scv,
                        const ElementSolutionVector &elemSol) const
    { return 1e-12; }

    PermeabilityType permeabilityAtPos(const GlobalPosition &globalPos) const
    { return 1e-12; }

    Scalar porosity(const Element &element,
                        const SubControlVolume &scv,
                        const ElementSolutionVector &elemSol) const
    { return 0.2; }
};


template<class TypeTag>
class MockProblem
{
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimension>;

public:
    MockProblem(const GridView& gridView) {}

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
        return PrimaryVariables(0.0);
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at a given global position.
     *
     * This is not specific to the discretization. By default it just
     * calls temperature().
     *
     * \param globalPos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    { return 273.15 + 10; }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This is discretization independent interface. By default it
     * just calls gravity().
     */
    GlobalPosition gravityAtPos(const GlobalPosition &pos) const
    { return GlobalPosition({0.0, 0.0, -9.81}); }

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolutionVector &elemSol) const
    {
        return 1.0;
    }

    const MockSpatialParams<TypeTag>& spatialParams() const
    { return spatialParams_; }

private:
    MockSpatialParams<TypeTag> spatialParams_;


};

namespace Properties
{

NEW_PROP_TAG(EnableFVGridGeometryCache);
NEW_PROP_TAG(FVGridGeometry);

NEW_TYPE_TAG(IncompressibleTestProblem, INHERITS_FROM(CCTpfaModel, OneP));

// Set the grid type
SET_TYPE_PROP(IncompressibleTestProblem, Grid, Dune::YaspGrid<2>);

// Set the finite volume grid geometry
SET_TYPE_PROP(IncompressibleTestProblem, FVGridGeometry, CCTpfaFVGridGeometry<TypeTag, true>);

// Set the problem type
SET_TYPE_PROP(IncompressibleTestProblem, Problem, MockProblem<TypeTag>);
SET_TYPE_PROP(IncompressibleTestProblem, SpatialParams, MockSpatialParams<TypeTag>);

// the grid variables
SET_TYPE_PROP(IncompressibleTestProblem, GridVariables, GridVariables<TypeTag>);

// the fluid system
SET_PROP(IncompressibleTestProblem, Fluid)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar> >;
};

// Enable caching
SET_BOOL_PROP(IncompressibleTestProblem, EnableGlobalVolumeVariablesCache, true);
SET_BOOL_PROP(IncompressibleTestProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(IncompressibleTestProblem, EnableGlobalFVGeometryCache, true);
SET_BOOL_PROP(IncompressibleTestProblem, EnableFVGridGeometryCache, true);

} // end namespace Properties

} // end namespace Dumux

#endif
