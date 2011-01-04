// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Melanie Darcis                                    *
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
 * \brief Nonisothermal gas injection problem where a gas (e.g. air) is injected into a fully
 *        water saturated medium. During buoyancy driven upward migration the gas
 *        passes a high temperature area.
 */

#ifndef DUMUX_INJECTIONPROBLEM2PNI_HH
#define DUMUX_INJECTIONPROBLEM2PNI_HH

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/boxmodels/2pni/2pnimodel.hh>

#include <dumux/material/fluidsystems/h2o_n2_system.hh>

// use the same spatial parameters as the injection problem of the
// 2p2c test program
#include "../2p2c/injectionspatialparameters.hh"

#define ISOTHERMAL 0

namespace Dumux {

template <class TypeTag>
class InjectionProblem2PNI;

namespace Properties
{
#if !ISOTHERMAL
NEW_TYPE_TAG(InjectionProblem2PNI, INHERITS_FROM(BoxTwoPNI));
#else
NEW_TYPE_TAG(InjectionProblem2PNI, INHERITS_FROM(BoxTwoP));
#endif

// Set the grid type
SET_PROP(InjectionProblem2PNI, Grid)
{
#if HAVE_UG
    typedef Dune::UGGrid<2> type;
#else
    typedef Dune::SGrid<2, 2> type;
    //typedef Dune::YaspGrid<2> type;
#endif
};

#if HAVE_DUNE_PDELAB
SET_PROP(InjectionProblem2PNI, LocalFEMSpace)
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
SET_PROP(InjectionProblem2PNI, Problem)
{
    typedef Dumux::InjectionProblem2PNI<TypeTag> type;
};

// Set the spatial parameters. we use the same spatial parameters as the
// 2p2c injection problem
SET_PROP(InjectionProblem2PNI, SpatialParameters)
{ typedef InjectionSpatialParameters<TypeTag> type; };

#if 1
// Use the same fluid system as the 2p2c injection problem
SET_PROP(InjectionProblem2PNI, FluidSystem)
{
    typedef H2O_N2_System<TypeTag> type;
};
#else
// Set the wetting phase
SET_PROP(InjectionProblem2PNI, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(InjectionProblem2PNI, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::GasPhase<Scalar, Dumux::N2<Scalar> > type;
};
#endif

// Enable gravity
SET_BOOL_PROP(InjectionProblem2PNI, EnableGravity, true);

// write convergence behaviour to disk?
SET_BOOL_PROP(InjectionProblem2PNI, NewtonWriteConvergence, true);
}

/*!
 * \ingroup TwoPNIBoxModel
 * \brief Nonisothermal gas injection problem where a gas (e.g. air) is injected into a fully
 *        water saturated medium. During buoyancy driven upward migration the gas
 *        passes a high temperature area.
 *
 * The domain is sized 40 m times 40 m. The rectangular area with the increased temperature (380 K)
 * starts at (20 m, 5 m) and ends at (30 m, 35 m)
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
 * This problem uses the \ref TwoPNIBoxModel.
 *
 * This problem should typically be simulated for 300000 seconds.
 * A good choice for the initial time step size is 1000 seconds.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_2pni ./grids/test_2pni.dgf 300000 1000</tt>
 */
template<class TypeTag>
class InjectionProblem2PNI
#if !ISOTHERMAL
  : public TwoPNIProblem<TypeTag>
#else
  : public TwoPNIProblem<TypeTag>
#endif
{
    typedef InjectionProblem2PNI<TypeTag> ThisType;
    typedef TwoPNIProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

#if ISOTHERMAL
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;
#else
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPNIIndices)) Indices;
#endif
    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,

#if !ISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif

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
    InjectionProblem2PNI(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        // initialize the tables of the fluid system
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
    { return "injection2pni"; }

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
    Scalar temperature(const Element &element,
                       const FVElementGeometry &fvElemGeom,
                       int scvIdx) const
    {
        return 273.15 + 30; // [K]
    };
#endif

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

        if (globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

#if !ISOTHERMAL
        // set a dirichlet value for the temperature, use the energy
        // equation to set the value
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

        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 + (depthBOR_ - globalPos[1])*densityW*9.81;
        values[saturationIdx] = 0.0;
#if !ISOTHERMAL
        values[temperatureIdx] = 283.0 + (depthBOR_ - globalPos[1])*0.03;
#endif
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
        if (globalPos[1] < 15 && globalPos[1] > 5) {
            // inject air. negative values mean injection
            values[contiNEqIdx] = -1e-3; // kg/(s*m^2)
        }
    }

    // \}


    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param values The source and sink values for the conservation equations
     * \param element The finite element
     * \param fvElemGeom The finite-volume geometry in the box scheme
     * \param scvIdx The local vertex index
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvElemGeom,
                int scvIdx) const
    {
           values = Scalar(0.0);
    }

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

        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 + (depthBOR_ - globalPos[1])*densityW*9.81;
        values[saturationIdx] = 0.0;

#if !ISOTHERMAL
        values[temperatureIdx] = 283.0 + (depthBOR_ - globalPos[1])*0.03;
        if (globalPos[0] > 20 && globalPos[0] < 30 && globalPos[1] > 5 && globalPos[1] < 35)
            values[temperatureIdx] = 380;
#endif // !ISOTHERMAL
    }
    // \}

private:
    static const Scalar depthBOR_ = 2700.0; // [m]
    static const Scalar eps_ = 1e-6;
};
} //end namespace

#endif
