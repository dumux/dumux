/*****************************************************************************
 *   Copyright (C) 2011 by Holger Class                                    *
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
#ifndef DUMUX_INFILTRATIONPROBLEM_HH
#define DUMUX_INFILTRATIONPROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidsystems/holle_h20_air_mesitylene_system.hh>

#include <dumux/boxmodels/holle3p3c/holle3p3cmodel.hh>

#include "infiltrationspatialparameters.hh"

#define ISOTHERMAL 0

namespace Dumux
{
template <class TypeTag>
class InfiltrationProblem;

namespace Properties
{
NEW_TYPE_TAG(InfiltrationProblem, INHERITS_FROM(BoxHolleThreePThreeC));

// Set the grid type
SET_PROP(InfiltrationProblem, Grid)
{
    typedef Dune::YaspGrid<2> type;
};

#if HAVE_DUNE_PDELAB
SET_PROP(InfiltrationProblem, LocalFEMSpace)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum{dim = GridView::dimension};

public:
    typedef Dune::PDELab::Q1LocalFiniteElementMap<Scalar,Scalar,dim> type; // for cubes
//    typedef Dune::PDELab::P1LocalFiniteElementMap<Scalar,Scalar,dim> type; // for simplices
};
#endif // HAVE_DUNE_PDELAB

// Set the problem property
SET_PROP(InfiltrationProblem, Problem)
{
    typedef Dumux::InfiltrationProblem<TypeTag> type;
};

// Set the wetting phase
SET_TYPE_PROP(InfiltrationProblem, FluidSystem, Dumux::Holle_H2O_Air_Mesitylene_System<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(InfiltrationProblem,
              SpatialParameters,
              Dumux::InfiltrationSpatialParameters<TypeTag>);

// Enable gravity
SET_BOOL_PROP(InfiltrationProblem, EnableGravity, true);

// Use forward differences instead of central differences
SET_INT_PROP(InfiltrationProblem, NumericDifferenceMethod, 0);

// Write newton convergence
SET_BOOL_PROP(InfiltrationProblem, NewtonWriteConvergence, true);

//! Set the formulation 
SET_INT_PROP(InfiltrationProblem, Formulation, HolleThreePThreeCFormulation::pgSwSn);
}


/*!
 * \ingroup HolleTHreePThreeCBoxModel
 *
 *  */
template <class TypeTag >
class InfiltrationProblem : public HolleThreePThreeCProblem<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GridView::Grid Grid;

    typedef InfiltrationProblem<TypeTag> ThisType;
    typedef HolleThreePThreeCProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLawParams)) MaterialLawParams;


    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(HolleThreePThreeCIndices)) Indices;
    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        pressureIdx = Indices::pressureIdx,
        switch1Idx = Indices::switch1Idx,
        switch2Idx = Indices::switch2Idx,
// #if !ISOTHERMAL
//        temperatureIdx = Indices::temperatureIdx,
//        energyEqIdx = Indices::energyEqIdx,
// #endif

        // Phase State
        threePhases = Indices::threePhases,
        wPhaseOnly  = Indices::wPhaseOnly,
        gnPhaseOnly = Indices::gnPhaseOnly,
        wnPhaseOnly = Indices::wnPhaseOnly,
        gPhaseOnly  = Indices::gPhaseOnly,
        wgPhaseOnly = Indices::wgPhaseOnly,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
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
    InfiltrationProblem(TimeManager &timeManager, const GridView &gridView)
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
    Scalar temperature(const Element &element,
                       const FVElementGeometry &fvElemGeom,
                       int scvIdx) const
    {
        return (273.15 + 10.0); // -> Temperatur 10°C
    };

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

        if(globalPos[0] > 125. - eps_)
            values.setAllDirichlet();
        else if(globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

//#if !ISOTHERMAL
//        values.setDirichlet(temperatureIdx, energyEqIdx);
//#endif
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
        Scalar Sw, Swr=0.12, Sgr=0.03;

        if(y >5. )
        {
          Scalar pc = 9.81 * 1000.0 * (y - 5);         /*capillary Eq.*/
                  if (pc < 0.0) pc = 0.0;

           Sw = invertPCGW_(pc,
                this->spatialParameters().materialLawParams());
           if (Sw < Swr) Sw = Swr;
           if (Sw > 1.-Sgr) Sw = 1.-Sgr;

          values[pressureIdx] = 1e5 ;
          values[switch1Idx] = Sw;
          values[switch2Idx] = 0.001;
        }else {
          values[pressureIdx] = 1e5 + 9.81 * 1000.0 * (5 - y);
          values[switch1Idx] = 1.-Sgr;
          values[switch2Idx] = 0.001;
        }

      //initial_(values, globalPos, element);
	//const MaterialLawParams& materialParams = this->spatialParameters().materialLawParams();;
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
    using ParentType::neumann;
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
        if ((globalPos[0] <= 80.+eps_) && (globalPos[0] >= 75.+eps_) && (globalPos[1] >= 10.-eps_))
        {
            values[Indices::contiWEqIdx] = -0.0;
            values[Indices::contiCEqIdx] = -1.2e-4;
            values[Indices::contiAEqIdx] = -0.0;
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

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos, const Element &element) const
    {
        Scalar y = globalPos[1];
        Scalar Sw, Swr=0.12, Sgr=0.03;
         
        if(y >5. )
        {
          Scalar pc = 9.81 * 1000.0 * (y - 5);         /*capillary Eq.*/
                  if (pc < 0.0) pc = 0.0;

           Sw = invertPCGW_(pc,
                this->spatialParameters().materialLawParams());
           if (Sw < Swr) Sw = Swr;
           if (Sw > 1.-Sgr) Sw = 1.-Sgr;

          values[pressureIdx] = 1e5 ;
          values[switch1Idx] = Sw;
          values[switch2Idx] = 0.001;
        }else {
          values[pressureIdx] = 1e5 + 9.81 * 1000.0 * (5 - y);
          values[switch1Idx] = 1.-Sgr;
          values[switch2Idx] = 0.001;
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

    static const Scalar eps_ = 1e-6;
};
} //end namespace

#endif
