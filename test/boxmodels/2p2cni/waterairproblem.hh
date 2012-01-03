// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Klaus Mosthaf                                     *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium.
 */
#ifndef DUMUX_WATER_AIR_PROBLEM_HH
#define DUMUX_WATER_AIR_PROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>

#include <dumux/boxmodels/2p2cni/2p2cnimodel.hh>

#include "waterairspatialparameters.hh"

#define ISOTHERMAL 0

namespace Dumux
{
template <class TypeTag>
class WaterAirProblem;

namespace Properties
{
NEW_TYPE_TAG(WaterAirProblem, INHERITS_FROM(BoxTwoPTwoCNI, WaterAirSpatialParameters));

// Set the grid type
SET_PROP(WaterAirProblem, Grid)
{
    typedef Dune::YaspGrid<2> type;
};

// Set the problem property
SET_PROP(WaterAirProblem, Problem)
{
    typedef Dumux::WaterAirProblem<TypeTag> type;
};

// Set the wetting phase
SET_TYPE_PROP(WaterAirProblem, FluidSystem, Dumux::FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)), false>);

// Enable gravity
SET_BOOL_PROP(WaterAirProblem, EnableGravity, true);

// Use forward differences instead of central differences
SET_INT_PROP(WaterAirProblem, NumericDifferenceMethod, +1);

// Write newton convergence
SET_BOOL_PROP(WaterAirProblem, NewtonWriteConvergence, false);
}


/*!
 * \ingroup TwoPTwoCNIModel
 * \ingroup BoxTestProblems
 * \brief Non-isothermal gas injection problem where a gas (e.g. air)
 *        is injected into a fully water saturated medium. During
 *        buoyancy driven upward migration the gas passes a high
 *        temperature area.
 *
 * The domain is sized 40 m times 40 m in a depth of 1000 m. The rectangular area
 * with the increased temperature (380 K) starts at (20 m, 5 m) and ends at
 * (30 m, 35 m).
 *
 * For the mass conservation equation neumann boundary conditions are used on
 * the top and on the bottom of the domain, while dirichlet conditions
 * apply on the left and the right boundary.
 * For the energy conservation equation dirichlet boundary conditions are applied
 * on all boundaries.
 *
 * Gas is injected at the bottom boundary from 15 m to 25 m at a rate of
 * 0.001 kg/(s m), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * At the dirichlet boundaries a hydrostatic pressure, a gas saturation of zero and
 * a geothermal temperature gradient of 0.03 K/m are applied.
 *
 * This problem uses the \ref TwoPTwoCNIModel.
 *
 * This problem should typically be simulated for 300000 s.
 * A good choice for the initial time step size is 1000 s.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_2p2cni ./grids/test_2p2cni.dgf 300000 1000</tt>
 *  */
template <class TypeTag >
class WaterAirProblem : public TwoPTwoCNIProblem<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GridView::Grid Grid;

    typedef WaterAirProblem<TypeTag> ThisType;
    typedef TwoPTwoCNIProblem<TypeTag> ParentType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;
    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,
#if !ISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif

        // Phase State
        lPhaseOnly = Indices::lPhaseOnly,
        gPhaseOnly = Indices::gPhaseOnly,
        bothPhases = Indices::bothPhases,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    WaterAirProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        maxDepth_ = 1000.0; // [m]
        eps_ = 1e-6;

        FluidSystem::init();
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
    const char *name() const
    { return "waterair"; }

#if ISOTHERMAL
    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param element The element
     * \param fvElemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local vertex index (SCV index)
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    {
        return 273.15 + 10; // -> 10°C
    };
#endif

    void sourceAtPos(PrimaryVariables &values,
                const GlobalPosition &globalPos) const
    {
        values = 0;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param vertex The vertex for which the boundary type is set
     */
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        if(globalPos[0] > 40 - eps_ || globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

#if !ISOTHERMAL
        values.setDirichlet(temperatureIdx, energyEqIdx);
#endif
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex for which the boundary type is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        initial_(values, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations
     * \param element The finite element
     * \param fvElemGeom The finite-volume geometry in the box scheme
     * \param is The intersection between element and boundary
     * \param scvIdx The local vertex index
     * \param boundaryFaceIdx The index of the boundary face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 const Intersection &is,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);
        values = 0;

        // negative values for injection
        if (globalPos[0] > 15 && globalPos[0] < 25 &&
            globalPos[1] < eps_)
        {
            values[Indices::contiGEqIdx] = -1e-3;
        }
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param element The finite element
     * \param fvElemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local vertex index
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 int scvIdx) const
    {
           const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        initial_(values, globalPos);

#if !ISOTHERMAL
            if (globalPos[0] > 20 && globalPos[0] < 30 && globalPos[1] < 30)
               values[temperatureIdx] = 380;
#endif
    }

    /*!
     * \brief Return the initial phase state inside a control volume.
     *
     * \param vert The vertex
     * \param globalIdx The index of the global vertex
     * \param globalPos The global position
     */
    int initialPhasePresence(const Vertex &vert,
                             int &globalIdx,
                             const GlobalPosition &globalPos) const
    {
        return lPhaseOnly;
    }

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 + (maxDepth_ - globalPos[1])*densityW*9.81;
        values[switchIdx] = 0.0;
#if !ISOTHERMAL
        values[temperatureIdx] = 283.0 + (maxDepth_ - globalPos[1])*0.03;
#endif
    }

    Scalar maxDepth_;
    Scalar eps_;
};
} //end namespace

#endif
