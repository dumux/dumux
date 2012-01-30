/*****************************************************************************
 *   Copyright (C) 2009-2011 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
/**
 * @file
 * @brief  Definition of a simple Stokes problem
 * @author Klaus Mosthaf, Andreas Lauser, Bernd Flemisch
 */
#ifndef DUMUX_STOKES2CNITESTPROBLEM_HH
#define DUMUX_STOKES2CNITESTPROBLEM_HH

#if HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>

#include <dumux/freeflow/stokes2cni/stokes2cnimodel.hh>

namespace Dumux
{

template <class TypeTag>
class Stokes2cniTestProblem;

//////////
// Specify the properties for the stokes problem
//////////
namespace Properties
{
NEW_TYPE_TAG(Stokes2cniTestProblem, INHERITS_FROM(BoxStokes2cni));

// Set the grid type
SET_TYPE_PROP(Stokes2cniTestProblem, Grid, Dune::SGrid<2, 2>);

// Set the problem property
SET_TYPE_PROP(Stokes2cniTestProblem, Problem, Stokes2cniTestProblem<TypeTag>);

//! Select the fluid system
SET_PROP(BoxStokes2cni, FluidSystem)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dumux::FluidSystems::H2ON2<Scalar> type;
};

//! Scalar is set to type long double for higher accuracy
//SET_TYPE_PROP(BoxStokes, Scalar, long double); // for a higher accuracy

//! a stabilization factor. Set to zero for no stabilization
SET_SCALAR_PROP(BoxStokes2cni, StabilizationAlpha, -1.0); // -1.0

//! stabilization at the boundaries
SET_SCALAR_PROP(BoxStokes2cni, StabilizationBeta, 0.0);

// Enable gravity
SET_BOOL_PROP(Stokes2cniTestProblem, EnableGravity, true);
}

/*!
 * \ingroup BoxStokes2cniModel
 * \ingroup BoxTestProblems
 * \brief Stokes2cni problem with air (N2) flowing
 *        from the left to the right.
 *
 * The domain is sized 1m times 1m. The boundary conditions for the momentum balances
 * are all set to Dirichlet. The mass balance has outflow boundary conditions, which are
 * replaced in the localresidual by the sum
 * of the two momentum balances. In the middle of the right boundary,
 * one vertex obtains Dirichlet bcs to fix the pressure at one point.
 *
 * This problem uses the \ref BoxStokes2cniModel.
 *
 * This problem is non-stationary and can be simulated until \f$t_{\text{end}} =
 * 100\;s\f$ is reached. A good choice for the initial time step size
 * is \f$t_{\text{inital}} = 1\;s\f$.
 * To run the simulation execute the following line in shell:
 * <tt>./test_stokes2cni grids/test_stokes2cni.dgf 100 1</tt>
 */
template <class TypeTag>
class Stokes2cniTestProblem : public StokesProblem<TypeTag>
{
    typedef StokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Stokes2cniIndices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { // Number of equations and grid dimension
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension
    };
    enum { // copy some indices for convenience
        massBalanceIdx = Indices::massBalanceIdx,
        momentumXIdx = Indices::momentumXIdx, //!< Index of the x-component of the momentum balance
        momentumYIdx = Indices::momentumYIdx, //!< Index of the y-component of the momentum balance
        momentumZIdx = Indices::momentumZIdx, //!< Index of the z-component of the momentum balance
        transportIdx = Indices::transportIdx, //!< Index of the transport equation (massfraction)
        energyIdx =    Indices::energyIdx     //!< Index of the energy equation (temperature)
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<CoordScalar, dim> GlobalPosition;

public:
    Stokes2cniTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 1e-6;

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
    { return "stokes2cni"; }

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
     * \param vertex The vertex on the boundary for which the
     *               conditions needs to be specified
     */
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        values.setAllDirichlet();

//        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
//        {
//            values.setOutflow(transportIdx);
//            values.setOutflow(energyIdx);
//        }

        if (onRightBoundary_(globalPos))
            values.setAllOutflow();

        // set all corner points to dirichlet
        if ((onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)) &&
                (onUpperBoundary_(globalPos) || onLowerBoundary_(globalPos)))
            values.setAllDirichlet();

        // the mass balance has to be of type outflow
        values.setOutflow(massBalanceIdx);

        // set pressure at one point
        const Scalar middle = (this->bboxMax()[1] - this->bboxMin()[1])/2;
//        const Scalar middle = this->bboxMin()[1] + (this->bboxMax()[1] - this->bboxMin()[1])/2;

        if (onLeftBoundary_(globalPos) &&
                globalPos[1] > middle - eps_ && globalPos[1] < middle + eps_)
            values.setDirichlet(massBalanceIdx);
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex representing the "half volume on the boundary"
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        initial_(values, globalPos);

        if(onLeftBoundary_(globalPos)
                && globalPos[1]<0.75 && globalPos[1]>0.25)
        {
            values[transportIdx] = 0.9e-4;
            values[energyIdx] = 284.15;
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
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
//        const GlobalPosition &globalPos
//            = is.geometry().center();

        values = 0.0;
    }

    /*!
     * \brief Evaluate the Beavers-Joseph coefficient
     *        at the center of a given intersection
     *
     * \return Beavers-Joseph coefficient
     */
    Scalar beaversJosephCoeff(const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 const Intersection &is,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
//        const GlobalPosition &globalPos = is.geometry().center();

//        if (onLowerBoundary_(globalPos))
//                && globalPos[0] > this->bboxMin()[0]+eps_ &&  globalPos[0] < this->bboxMax()[0]-eps_)
//            return 1.0;
//        else
            return 0.0;
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
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &,
                int subControlVolumeIdx) const
    {
        // ATTENTION: The source term of the mass balance has to be chosen as
        // div (q_momentum) in the problem file
        values = Scalar(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 int scvIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        initial_(values, globalPos);
    }

    /*!
    * \brief Evaluate the intrinsic permeability
    *        at the corner of a given element
    *
    * \return permeability in x-direction
    */
   Scalar permeability(const Element &element,
                const FVElementGeometry &fvElemGeom,
                int scvIdx) const
   {
       return 1e-8;
   }
   // \}

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
//        const Scalar v0 = 1.0;
        const Scalar v0 = 2.0;

        values[momentumXIdx] = 0.0;//v0*(globalPos[1]-1.0) + 1e-4;
        values[momentumYIdx] = 0.0;
        if (globalPos[1] < this->bboxMax()[1] && globalPos[1] > this->bboxMin()[1])
            values[momentumXIdx] = v0*(this->bboxMax()[1]-globalPos[1])*(globalPos[1]-this->bboxMin()[1]);
//        if (//onUpperBoundary_(globalPos) &&
//               globalPos[0]<0.75 && globalPos[0]>0.25)
//            values[momentumYIdx] = v1*(0.75-globalPos[0])*(globalPos[0]-0.25);

        values[massBalanceIdx] = 1e5;
        values[transportIdx] = 1e-4;
        values[energyIdx] = 283.15;
    }
    // \}

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->bboxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->bboxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->bboxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->bboxMax()[1] - eps_; }

    bool onBoundary_(const GlobalPosition &globalPos) const
    {
        return (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)
                || onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos));
    }

    Scalar eps_;
};
} //end namespace

#endif
