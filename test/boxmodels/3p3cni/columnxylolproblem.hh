// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Holger Class                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later vesion.                                     *
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
 *        is injected into a fully water saturated medium with a residually
 *        trapped NAPL contamination.
 */
#ifndef DUMUX_COLUMNXYLOLPROBLEM_HH
#define DUMUX_COLUMNXYLOLPROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#warning "TODO: cleanup the h2o-air-xylene fluid system. we use H2OAirMesitylene in the mean time!"
//#include <dumux/material/fluidsystems/h2o_air_xylene_system.hh>
#include <dumux/material/fluidsystems/h2oairmesitylenefluidsystem.hh>

#include <dumux/boxmodels/3p3cni/3p3cnimodel.hh>

#include "columnxylolspatialparameters.hh"

#define ISOTHERMAL 0

namespace Dumux
{
template <class TypeTag>
class ColumnProblem;

namespace Properties
{
NEW_TYPE_TAG(ColumnProblem, INHERITS_FROM(BoxThreePThreeCNI, ColumnSpatialParameters));

// Set the grid type
SET_PROP(ColumnProblem, Grid)
{
    typedef Dune::YaspGrid<2> type;
};

// Set the problem property
SET_PROP(ColumnProblem, Problem)
{
    typedef Dumux::ColumnProblem<TypeTag> type;
};

// Set the fluid system
SET_TYPE_PROP(ColumnProblem, 
              FluidSystem,
              Dumux::FluidSystems::H2OAirMesitylene<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Enable gravity
SET_BOOL_PROP(ColumnProblem, EnableGravity, true);

// Use forward differences instead of central differences
SET_INT_PROP(ColumnProblem, NumericDifferenceMethod, 0);

// Write newton convergence
SET_BOOL_PROP(ColumnProblem, NewtonWriteConvergence, true);

//! Set the formulation
SET_INT_PROP(ColumnProblem, Formulation, ThreePThreeCFormulation::pgSwSn);

// Set the maximum time step
SET_SCALAR_PROP(ColumnProblem, MaxTimeStepSize, 5.);
}


/*!
 * \ingroup ThreePThreeCNIBoxModel
 *
 *  */
template <class TypeTag >
class ColumnProblem : public ThreePThreeCNIProblem<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GridView::Grid Grid;

    typedef ColumnProblem<TypeTag> ThisType;
    typedef ThreePThreeCNIProblem<TypeTag> ParentType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, ThreePThreeCIndices) Indices;
    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),

        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,

        // Phase State
        threePhases = Indices::threePhases,
        wPhaseOnly  = Indices::wPhaseOnly,
        gnPhaseOnly = Indices::gnPhaseOnly,
        wnPhaseOnly = Indices::wnPhaseOnly,
        gPhaseOnly  = Indices::gPhaseOnly,
        wgPhaseOnly = Indices::wgPhaseOnly,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    ColumnProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
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
    { return "columnxylol"; }


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

        if(globalPos[1] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
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
        if (globalPos[1] > 1.2 - eps_)
        {
            values[Indices::contiWEqIdx] = -0.395710;
            values[Indices::contiAEqIdx] = -0.000001;
            values[Indices::contiCEqIdx] = -0.00;
            values[Indices::energyEqIdx] = -17452.97;
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
        return threePhases;
    }

    /*
      bool shouldWriteOutput() const
      {
      return
      this->timeManager().timeStepIndex() == 0 ||
      this->timeManager().timeStepIndex() == 10 ||
      this->timeManager().timeStepIndex() == 20 ||
      this->timeManager().timeStepIndex() == 50 ||
      this->timeManager().timeStepIndex() == 100 ||
      this->timeManager().timeStepIndex() == 200 ||
      this->timeManager().timeStepIndex() == 300 ||
      this->timeManager().timeStepIndex() == 400 ||
      this->timeManager().timeStepIndex() == 500 ||
      this->timeManager().timeStepIndex() == 600 ||
      this->timeManager().timeStepIndex() == 700 ||
      this->timeManager().timeStepIndex() == 800 ||
      this->timeManager().timeStepIndex() == 900 ||
      this->timeManager().timeStepIndex() == 1000 ||
      this->timeManager().timeStepIndex() == 1100 ||
      this->timeManager().timeStepIndex() == 1200 ||
      this->timeManager().timeStepIndex() == 1300 ||
      this->timeManager().timeStepIndex() == 1400 ||
      this->timeManager().timeStepIndex() == 1500 ||
      this->timeManager().timeStepIndex() == 1600 ||
      this->timeManager().timeStepIndex() == 1700 ||
      this->timeManager().timeStepIndex() == 1800 ||
      this->timeManager().timeStepIndex() == 1900 ||
      this->timeManager().timeStepIndex() == 2000 ||
      this->timeManager().timeStepIndex() == 2100 ||
      this->timeManager().timeStepIndex() == 2200 ||
      this->timeManager().timeStepIndex() == 2300 ||
      this->timeManager().timeStepIndex() == 2400 ||
      this->timeManager().timeStepIndex() == 2500 ||
      this->timeManager().timeStepIndex() == 2600 ||
      this->timeManager().timeStepIndex() == 2700 ||
      this->timeManager().timeStepIndex() == 2800 ||
      this->timeManager().timeStepIndex() == 2900 ||
      this->timeManager().willBeFinished();
      }
    */


private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar y = globalPos[1];

        values[temperatureIdx] = 296.15;
        values[pressureIdx] = 1.e5;

        if(y > 1.2-eps_){
            values[switch2Idx] = 0.112; // almost no contaminant component
            values[switch1Idx] = 0.005;
        } else if(y < 1.2-0.3){ // extended domain
            values[switch2Idx] = 1.e-4; // almost no contaminant component
            values[switch1Idx] = 0.005;
        } else {
            values[switch1Idx] = 0.005;
            if((y<=1.2-0.001)&&(y>=1.2-0.0148)) values[switch2Idx] = 0+((1.2-y)/0.0148)*0.112;
            if((y<1.2-0.0148)&&(y>=1.2-0.0296)) values[switch2Idx] = 0.112+(((1.2-y)-0.0148)/0.0148)*(0.120-0.112);
            if((y<1.2-0.0296)&&(y>=1.2-0.0444)) values[switch2Idx] = 0.120+(((1.2-y)-0.0296)/0.0148)*(0.125-0.120);
            if((y<1.2-0.0444)&&(y>=1.2-0.0592)) values[switch2Idx] = 0.125+(((1.2-y)-0.0444)/0.0148)*(0.137-0.125);
            if((y<1.2-0.0592)&&(y>=1.2-0.0740)) values[switch2Idx] = 0.137+(((1.2-y)-0.0592)/0.0148)*(0.150-0.137);
            if((y<1.2-0.0740)&&(y>=1.2-0.0888)) values[switch2Idx] = 0.150+(((1.2-y)-0.0740)/0.0148)*(0.165-0.150);
            if((y<1.2-0.0888)&&(y>=1.2-0.1036)) values[switch2Idx] = 0.165+(((1.2-y)-0.0888)/0.0148)*(0.182-0.165);
            if((y<1.2-0.1036)&&(y>=1.2-0.1184)) values[switch2Idx] = 0.182+(((1.2-y)-0.1036)/0.0148)*(0.202-0.182);
            if((y<1.2-0.1184)&&(y>=1.2-0.1332)) values[switch2Idx] = 0.202+(((1.2-y)-0.1184)/0.0148)*(0.226-0.202);
            if((y<1.2-0.1332)&&(y>=1.2-0.1480)) values[switch2Idx] = 0.226+(((1.2-y)-0.1332)/0.0148)*(0.257-0.226);
            if((y<1.2-0.1480)&&(y>=1.2-0.1628)) values[switch2Idx] = 0.257+(((1.2-y)-0.1480)/0.0148)*(0.297-0.257);
            if((y<1.2-0.1628)&&(y>=1.2-0.1776)) values[switch2Idx] = 0.297+(((1.2-y)-0.1628)/0.0148)*(0.352-0.297);
            if((y<1.2-0.1776)&&(y>=1.2-0.1924)) values[switch2Idx] = 0.352+(((1.2-y)-0.1776)/0.0148)*(0.426-0.352);
            if((y<1.2-0.1924)&&(y>=1.2-0.2072)) values[switch2Idx] = 0.426+(((1.2-y)-0.1924)/0.0148)*(0.522-0.426);
            if((y<1.2-0.2072)&&(y>=1.2-0.2220)) values[switch2Idx] = 0.522+(((1.2-y)-0.2072)/0.0148)*(0.640-0.522);
            if((y<1.2-0.2220)&&(y>=1.2-0.2368)) values[switch2Idx] = 0.640+(((1.2-y)-0.2220)/0.0148)*(0.767-0.640);
            if((y<1.2-0.2368)&&(y>=1.2-0.2516)) values[switch2Idx] = 0.767+(((1.2-y)-0.2368)/0.0148)*(0.878-0.767);
            if((y<1.2-0.2516)&&(y>=1.2-0.2664)) values[switch2Idx] = 0.878+(((1.2-y)-0.2516)/0.0148)*(0.953-0.878);
            if((y<1.2-0.2664)&&(y>=1.2-0.2812)) values[switch2Idx] = 0.953+(((1.2-y)-0.2664)/0.0148)*(0.988-0.953);
            if((y<1.2-0.2812)&&(y>=1.2-0.3000)) values[switch2Idx] = 0.988;
        }

    }

    static constexpr Scalar eps_ = 1e-6;
};
} //end namespace

#endif
