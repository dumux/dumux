/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
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
 * \brief A water infiltration problem with a low-permeability lens
 *        embedded into a high- permeability domain which uses the
 *        Richards box model.
 */
#ifndef DUMUX_RICHARDS_LENSPROBLEM_HH
#define DUMUX_RICHARDS_LENSPROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/boxmodels/richards/richardsmodel.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include "richardslensspatialparameters.hh"

namespace Dumux
{

template <class TypeTag>
class RichardsLensProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
NEW_TYPE_TAG(RichardsLensProblem, INHERITS_FROM(BoxRichards));

// Set the grid type. Use UG if available, else SGrid
#if HAVE_UG
SET_TYPE_PROP(RichardsLensProblem, Grid, Dune::UGGrid<2>);
#else
SET_PROP(RichardsLensProblem, Grid) { typedef Dune::SGrid<2, 2> type; };
//SET_TYPE_PROP(RichardsLensProblem, Grid, Dune::YaspGrid<2>);
#endif

#if HAVE_DUNE_PDELAB
// set the local finite element space to be used to calculate the
// gradients in the flux calculation
SET_PROP(RichardsLensProblem, LocalFEMSpace)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum{dim = GridView::dimension};

public:
    typedef Dune::PDELab::Q1LocalFiniteElementMap<Scalar,Scalar,dim> type; // for cubes
//    typedef Dune::PDELab::P1LocalFiniteElementMap<Scalar,Scalar,dim> type; // for simplices
};
#endif // HAVE_DUNE_PDELAB

// Set the phsical problem to be solved
SET_PROP(RichardsLensProblem, Problem)
{ typedef Dumux::RichardsLensProblem<TypeTag> type; };

// Set the wetting phase
SET_PROP(RichardsLensProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the spatial parameters
SET_PROP(RichardsLensProblem, SpatialParameters)
{ typedef Dumux::RichardsLensSpatialParameters<TypeTag> type; };

// Enable gravity?
SET_BOOL_PROP(RichardsLensProblem, EnableGravity, true);

// Enable partial reassembly of the Jacobian matrix?
SET_BOOL_PROP(RichardsLensProblem, EnablePartialReassemble, true);

// Use forward diffferences to approximate the Jacobian matrix
SET_BOOL_PROP(RichardsLensProblem, NumericDifferenceMethod, +1);

// Write the intermediate results of the newton method?
SET_BOOL_PROP(RichardsLensProblem, NewtonWriteConvergence, false);
}

/*!
 * \ingroup RichardsModel
 *
 * \brief A water infiltration problem with a low-permeability lens
 *        embedded into a high- permeability domain which uses the
 *        Richards box model.
 *
 * The domain is box shaped. Left and right boundaries are Dirichlet
 * boundaries with fixed water pressure (fixed Saturation \f$S_w = 0\f$),
 * bottom boundary is closed (Neumann 0 boundary), the top boundary
 * (Neumann 0 boundary) is also closed except for infiltration
 * section, where water is infiltrating into an initially unsaturated
 * porous medium. This problem is very similar the the LensProblem
 * which uses the TwoPBoxModel, with the main difference being that
 * the domain is initally fully saturated by gas instead of water and
 * water instead of a %DNAPL infiltrates from the top.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_richards grids/richardslens.dgf 10e6 100</tt>
 *
 * where the initial time step is 100 seconds, and the end of the
 * simulation time is 10,000,000 seconds (115.7 days)
 */
template <class TypeTag = TTAG(RichardsLensProblem) >
class RichardsLensProblem : public RichardsBoxProblem<TypeTag>
{
    typedef RichardsLensProblem<TypeTag> ThisType;
    typedef RichardsBoxProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(RichardsIndices)) Indices;
    enum {
        // copy some indices for convenience
        pwIdx = Indices::pwIdx,
        contiEqIdx = Indices::contiEqIdx,

        // Grid and world dimension
        dim = GridView::dimensionworld,
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dim> GlobalPosition;

public:
    /*!
     * \brief Constructor
     *
     * \param timeManager The Dumux TimeManager for simulation management.
     * \param gridView The grid view on the spatial domain of the problem
     * \param lensLowerLeft The lower left coordinate of the
     *                      low-permeability lens
     * \param lensUpperRight The upper right coordinate of the
     *                       low-permeability lens
     */
    RichardsLensProblem(TimeManager &timeManager,
                        const GridView &gridView,
                        const GlobalPosition &lensLowerLeft,
                        const GlobalPosition &lensUpperRight)
        : ParentType(timeManager, gridView)
    {
        lensLowerLeft_=lensLowerLeft;
        lensUpperRight_=lensUpperRight;
        pnRef_ = 1e5;
        this->spatialParameters().setLensCoords(lensLowerLeft_, lensUpperRight_);
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    { return "richardslens"; }

    /*!
     * \brief Returns the temperature [K] within a finite volume
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     *
     * \param element The DUNE Codim<0> entity which intersects with
     *                the finite volume in question
     * \param fvElemGeom The finite volume geometry of the element
     * \param scvIdx The sub control volume index inside the finite
     *               volume geometry
     */
    Scalar temperature(const Element &element,
                       const FVElementGeometry &fvElemGeom,
                       int scvIdx) const
    { return 273.15 + 10; }; // -> 10°C

    /*!
     * \brief Returns the reference pressure [Pa] of the non-wetting
     *        fluid phase within a finite volume
     *
     * This problem assumes a constant reference pressure of 1 bar.
     *
     * \param element The DUNE Codim<0> entity which intersects with
     *                the finite volume in question
     * \param fvElemGeom The finite volume geometry of the element
     * \param scvIdx The sub control volume index inside the finite
     *               volume geometry
     */
    Scalar referencePressure(const Element &element,
                             const FVElementGeometry &fvElemGeom,
                             int scvIdx) const
    { return pnRef_; };

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

        if (onLeftBoundary_(globalPos) ||
            onRightBoundary_(globalPos))
        {
            values.setAllDirichlet();
        }
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
        // use initial values as boundary conditions
        initial_(values, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     *
     * \param values The neumann values for the conservation equations
     * \param element The DUNE Codim<0> entity which intersects with
     *                the finite volume in question
     * \param fvElemGeom The finite volume geometry of the element
     * \param is The DUNE boundary intersection of the boundary segment
     * \param scvIdx The sub control volume index of the finite
     *               volume geometry
     * \param boundaryFaceIdx The index of the boundary face of the
     *                        finite volume geometry
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 const Intersection &is,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {

        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        values = 0.0;
        if (onInlet_(globalPos)) {
            // inflow of water
            values[contiEqIdx] = -0.04; // kg / (m * s)
        }
    }

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Specify the source term [kg/m^3] for the wetting phase
     *        within a given sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate
     * wetting phase mass is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     *
     * \param values Storage for all values for the source terms
     * \param element The DUNE Codim<0> entity which intersects with
     *                the finite volume in question
     * \param fvElemGeom The finite volume geometry of the element
     * \param scvIdx The sub control volume index of the finite
     *               volume geometry
     */
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvElemGeom,
                int scvIdx) const
    {
        values = Scalar(0.0);
    }

    /*!
     * \brief Evaluate the initial values for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     *
     * \param values Storage for all primary variables of the initial condition
     * \param element The DUNE Codim<0> entity which intersects with
     *                the finite volume in question
     * \param fvElemGeom The finite volume geometry of the element
     * \param scvIdx The sub control volume index of the finite
     *               volume geometry
     */
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 int scvIdx) const
    {
        const GlobalPosition pos = element.geometry().corner(scvIdx);

        initial_(values, pos);
    };

    // \}

private:
    void initial_(PrimaryVariables &values, const GlobalPosition &pos) const
    {
        Scalar Sw = 0.0;
        Scalar pc =
            MaterialLaw::pC(this->spatialParameters().materialLawParams(pos),
                            Sw);
        values[pwIdx] = pnRef_ - pc;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->bboxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->bboxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->bboxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->bboxMax()[1] - eps_;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->bboxMax()[0] - this->bboxMin()[0];
        Scalar lambda = (this->bboxMax()[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda && lambda < 2.0/3.0;
    }

    static const Scalar eps_ = 3e-6;
    Scalar pnRef_;
    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;
};
} //end namespace

#endif
