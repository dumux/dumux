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
 * \ingroup SequentialTwoPTests
 * \brief test problem for the explicit transport model
 */
#ifndef DUMUX_TEST_TRANSPORT_PROBLEM_HH
#define DUMUX_TEST_TRANSPORT_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/porousmediumflow/2p/sequential/transport/cellcentered/properties.hh>
#include <dumux/porousmediumflow/2p/sequential/transport/problem.hh>

#include "test_transportspatialparams.hh"

namespace Dumux
{
/*!
 * \ingroup SequentialTwoPTests
 * \brief test problem for the explicit transport model
 */
template<class TypeTag>
class TestTransportProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(TransportTestTypeTag, INHERITS_FROM(FVTransportTwoP, TestTransportSpatialParams));

// Set the grid type
SET_TYPE_PROP(TransportTestTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(TransportTestTypeTag, Problem, TestTransportProblem<TypeTag>);

// Set the fluid system
SET_PROP(TransportTestTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

SET_INT_PROP(TransportTestTypeTag, VelocityFormulation, SequentialTwoPCommonIndices::velocityTotal);
}

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for the explicit transport model
 *
 * A unit "fluid" is injected from the left side into a rectangular 2D
 * domain also this testing fluid. Upper and lower boundary are closed (Neumann = 0),
 * and there is free outflow on the right side.
 *
 * This test solely applies the 2p transport on a given velocity field, without a
 * pressure field being solved.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_transport -parameterFile ./test_transport.input</tt>,
 * where the arguments define the parameter file.
 */
template<class TypeTag>
class TestTransportProblem: public TransportProblem2P<TypeTag>
{
    using ParentType = TransportProblem2P<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Grid = typename GridView::Grid;

    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;

    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using SolutionTypes = typename GET_PROP(TypeTag, SolutionTypes);
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;

    using CellData = typename GET_PROP_TYPE(TypeTag, CellData);

    enum
    {
        dimWorld = GridView::dimensionworld
    };

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;
public:
    TestTransportProblem(TimeManager& timeManager, Grid& grid) :
        ParentType(timeManager, grid)
    {}


    //!set initial velocity field -> v_total is assumed to be constant in this test!
    void init()
    {
        ParentType::init();

        VelocityVector vel(0);
        vel[0] = 1e-5;

        // compute update vector
        for (const auto& element : elements(this->gridView()))
        {
            // cell index
            int eIdxGlobal = this->elementMapper().index(element);

            CellData& cellData = this->variables().cellData(eIdxGlobal);

            // run through all intersections with neighbors and boundary
            for (const auto& intersection : intersections(this->gridView(), element))
            {
                // local number of facet
                int indexInInside = intersection.indexInInside();

                cellData.fluxData().setVelocity(wPhaseIdx, indexInInside, vel);

                const GlobalPosition& unitOuterNormal = intersection.centerUnitOuterNormal();

                Scalar pot = vel * unitOuterNormal;

                cellData.fluxData().setUpwindPotential(wPhaseIdx, indexInInside, pot);
                cellData.fluxData().setUpwindPotential(nPhaseIdx, indexInInside, pot);
            }
        }
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
    std::string name() const
    {
        return "test_transport";
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return 273.15 + 10; // -> 10°C
    }

    // \}


    //! Returns the reference pressure for evaluation of constitutive relations
    Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
    {
        return 1e5; // -> 10°C
    }

    void sourceAtPos(PrimaryVariables &values,const GlobalPosition& globalPos) const
    {
        values = 0;
    }

    /*!
    * \brief Returns the type of boundary condition.
    *
    *
    * BC for saturation equation can be dirichlet (saturation), neumann (flux), or outflow.
    */
    void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
    {
            if (globalPos[0] < eps_)
            {
                bcTypes.setAllDirichlet();
            }
            else if (globalPos[0] > this->bBoxMax()[0] - eps_)
            {
                bcTypes.setAllOutflow();
            }
            // all other boundaries
            else
            {
                bcTypes.setAllNeumann();
            }
    }

    //! set dirichlet condition  (saturation [-])
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = 0;
        if (globalPos[0] < eps_)
        {
            values = 1.0;
        }
    }

    //! set neumann condition for phases (flux, [kg/(m^2 s)])
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = 0;
    }

    //! return initial solution
    void initialAtPos(PrimaryVariables &values,
            const GlobalPosition &globalPos) const
    {
        values = 0;
    }

private:
    static constexpr Scalar eps_ = 1e-6;
};
} //end namespace

#endif
