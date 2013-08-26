// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *          See the file COPYING for full copying permissions.               *
 *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
 * \brief Isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 */
#ifndef DUMUX_NAPLINFILTRATIONPROBLEM3P_HH
#define DUMUX_NAPLINFILTRATIONPROBLEM3P_HH

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>

#include <dumux/implicit/3p/3pmodel.hh>
#include <dumux/implicit/common/implicitporousmediaproblem.hh>

#include <dumux/material/fluidsystems/h2oairmesitylenefluidsystem.hh>

#include "naplinfiltrationspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class InfiltrationProblem;

namespace Properties
{
NEW_TYPE_TAG(InfiltrationProblem, INHERITS_FROM(BoxThreeP, InfiltrationSpatialParams));

// Set the grid type
SET_TYPE_PROP(InfiltrationProblem, Grid, Dune::SGrid<2,2>);

// Set the problem property
SET_TYPE_PROP(InfiltrationProblem, Problem, Dumux::InfiltrationProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(InfiltrationProblem,
              FluidSystem,
              Dumux::FluidSystems::H2OAirMesitylene<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Enable gravity?
SET_BOOL_PROP(InfiltrationProblem, ProblemEnableGravity, true);

// Write newton convergence?
SET_BOOL_PROP(InfiltrationProblem, NewtonWriteConvergence, false);

// Maximum tolerated relative error in the Newton method
SET_SCALAR_PROP(InfiltrationProblem, NewtonRelTolerance, 1e-4);

// -1 backward differences, 0: central differences, +1: forward differences
SET_INT_PROP(InfiltrationProblem, ImplicitNumericDifferenceMethod, 1);
}

/*!
 * \ingroup ThreePBoxModel
 * \ingroup BoxTestProblems
 * \brief Isothermal NAPL infiltration problem: LNAPL contaminates
 *        the unsaturated and the saturated groundwater zone.
 *
 * The 2D domain of this test problem is 500 m long and 10 m deep, where
 * the lower part represents a slightly inclined groundwater table, and the
 * upper part is the vadose zone. 
 * A LNAPL (Non-Aqueous Phase Liquid which is lighter than water) infiltrates
 * (modelled with a Neumann boundary condition) into the vadose zone. Upon
 * reaching the water table, it spreads (since lighter than water) and migrates
 * on top of the water table in the direction of the slope.
 * On its way through the vadose zone, it leaves a trace of residually trapped
 * immobile NAPL, which can in the following dissolve and evaporate slowly,
 * and eventually be transported by advection and diffusion.
 *
 * Left and right boundaries are constant head boundaries (Dirichlet),
 * Top and bottom are Neumann boundaries, all no-flow except for the small
 * infiltration zone in the upper left part.
 *
 * This problem uses the \ref ThreePModel.
 *
 * This problem should typically be simulated for 30 days.
 * A good choice for the initial time step size is 60 s.
 * To adjust the simulation time it is necessary to edit the file naplinfiltrationexercise.input
 *
 * To run the simulation execute the following line in shell:
 * <tt>./naplinfiltrationexercise -parameterFile naplinfiltrationexercise.input</tt>
 *  */
template <class TypeTag >
class InfiltrationProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        pressureIdx = Indices::pressureIdx,
        swIdx = Indices::swIdx,
        snIdx = Indices::snIdx,

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

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    InfiltrationProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
        , eps_(1e-6)
    {
        temperature_ = 273.15 + 10.0; // -> 10 degrees Celsius
        FluidSystem::init( 272.15, 274.15, 3, 8e4, 3e5, 200);
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
    { return "infiltration"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param element The element
     * \param fvElemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local vertex index (SCV index)
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar boxTemperature(const Element &element,
                          const FVElementGeometry &fvElemGeom,
                          int scvIdx) const
    {
        return temperature_;
    };

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

        if(globalPos[0] > 500. - eps_)
            values.setAllDirichlet();
        else if(globalPos[0] < eps_)
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

        Scalar y = globalPos[1];
        Scalar x = globalPos[0];
        Scalar Sw, Swr=0.12, Sgr=0.03;

        if(y >(-1.E-3*x+5) )
        {
            Scalar pc = 9.81 * 1000.0 * (y - (-5E-4*x+5));
            if (pc < 0.0) pc = 0.0;

            Sw = invertPCGW_(pc,
                             this->spatialParams().materialLawParams());
            if (Sw < Swr) Sw = Swr;
            if (Sw > 1.-Sgr) Sw = 1.-Sgr;

            values[pressureIdx] = 1e5 ;
            values[swIdx] = Sw;
            values[snIdx] = 0.;
        }else {
            values[pressureIdx] = 1e5 + 9.81 * 1000.0 * ((-5E-4*x+5) - y);
            values[swIdx] = 1.-Sgr;
            values[snIdx] = 0.;
        }

        //initial_(values, globalPos, element);
        //const MaterialLawParams& materialParams = this->spatialParams().materialLawParams();;
        //MaterialLaw::pCGW(materialParams, 1.0);
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
        if (this->timeManager().time()<2592000.)
        {
            if ((globalPos[0] <= 175.+eps_) && (globalPos[0] >= 150.+eps_) && (globalPos[1] >= 10.-eps_))
            {
                values[Indices::contiWEqIdx] = -0.0;
                values[Indices::contiNEqIdx] = -0.001, // /*Molfluss, umr. über M(Mesit.)=0,120 kg/mol --> 1.2e-4  kg/(sm)
                values[Indices::contiGEqIdx] = -0.0;
            }
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

        initial_(values, globalPos, element);
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a uniform temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; };

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos, const Element &element) const
    {
        Scalar y = globalPos[1];
        Scalar x = globalPos[0];
        Scalar Sw, Swr=0.12, Sgr=0.03;

        if(y >(-1.E-3*x+5) )
        {
            Scalar pc = 9.81 * 1000.0 * (y - (-5E-4*x+5));
            if (pc < 0.0) pc = 0.0;

            Sw = invertPCGW_(pc,
                             this->spatialParams().materialLawParams());
            if (Sw < Swr) Sw = Swr;
            if (Sw > 1.-Sgr) Sw = 1.-Sgr;

            values[pressureIdx] = 1e5 ;
            values[swIdx] = Sw;
            values[snIdx] = 0.;
        }else {
            values[pressureIdx] = 1e5 + 9.81 * 1000.0 * ((-5E-4*x+5) - y);
            values[swIdx] = 1.-Sgr;
            values[snIdx] = 0.;
        }
    }

    static Scalar invertPCGW_(Scalar pcIn, const MaterialLawParams &pcParams)
    {
        Scalar lower,upper;
        int k;
        int maxIt = 50;
        Scalar bisLimit = 1.;
        Scalar Sw, pcGW;
        lower=0.0; upper=1.0;
        for (k=1; k<=25; k++)
        {
            Sw = 0.5*(upper+lower);
            pcGW = MaterialLaw::pCGW(pcParams, Sw);
            Scalar delta = pcGW-pcIn;
            if (delta<0.) delta*=-1.;
            if (delta<bisLimit)
            {
                return(Sw);
            }
            if (k==maxIt) {
                return(Sw);
            }
            if (pcGW>pcIn) lower=Sw;
            else upper=Sw;
        }
        return(Sw);
    }

    Scalar temperature_;
    Scalar eps_;
};
} //end namespace

#endif
