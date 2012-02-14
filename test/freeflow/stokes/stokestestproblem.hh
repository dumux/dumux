// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2010 by Klaus Mosthaf                                *
 *   Copyright (C) 2010 by Katherina Baber                                   *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
/**
 * @file
 * @brief  Definition of a simple Stokes problem
 * @author Klaus Mosthaf, Andreas Lauser, Bernd Flemisch
 */
#ifndef DUMUX_STOKESTESTPROBLEM_HH
#define DUMUX_STOKESTESTPROBLEM_HH

#if HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>
#include <dumux/material/fluidsystems/gasphase.hh>

#include <dumux/freeflow/stokes/stokesmodel.hh>

namespace Dumux
{

template <class TypeTag>
class StokesTestProblem;

//////////
// Specify the properties for the stokes problem
//////////
namespace Properties
{
NEW_TYPE_TAG(StokesTestProblem, INHERITS_FROM(BoxStokes));

// Set the grid type
SET_TYPE_PROP(StokesTestProblem, Grid, Dune::SGrid<2, 2>);

// Set the problem property
SET_TYPE_PROP(StokesTestProblem, Problem, Dumux::StokesTestProblem<TypeTag>);

SET_PROP(StokesTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::GasPhase<Scalar, Dumux::N2<Scalar> > type;
};
//! Scalar is set to type long double for higher accuracy
SET_TYPE_PROP(BoxStokes, Scalar, double);
//SET_TYPE_PROP(BoxStokes, Scalar, long double);

//! A stabilization factor. Set negative for stabilization and to zero for no stabilization
SET_SCALAR_PROP(BoxStokes, StabilizationAlpha, -1.0);

//! Stabilization factor for the boundaries
SET_SCALAR_PROP(BoxStokes, StabilizationBeta, 0.0);

// Enable gravity
SET_BOOL_PROP(StokesTestProblem, EnableGravity, false);
}

/*!
 * \ingroup BoxStokesModel
 * \ingroup BoxTestProblems
 * \brief Stokes flow problem with nitrogen (N2) flowing
 *        from the left to the right.
 *
 * The domain is sized 1m times 1m. The boundary conditions for the momentum balances
 * are set to Dirichlet with outflow on the right boundary. The mass balance has
 * outflow bcs, which are replaced in the localresidual by the sum
 * of the two momentum balances. In the middle of the right boundary,
 * one vertex receives Dirichlet bcs to set the pressure level.
 *
 * This problem uses the \ref BoxStokesModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_stokes -parameterFile ./test_stokes.input</tt>
 */
template <class TypeTag>
class StokesTestProblem : public StokesProblem<TypeTag>
{
    typedef StokesProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, StokesIndices) Indices;
    enum {

        // Number of equations and grid dimension
        dim = GridView::dimension,

        // copy some indices for convenience
        massBalanceIdx = Indices::massBalanceIdx,
        momentumXIdx = Indices::momentumXIdx, //!< Index of the x-component of the momentum balance
        momentumYIdx = Indices::momentumYIdx //!< Index of the y-component of the momentum balance
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<CoordScalar, dim> GlobalPosition;

public:
    StokesTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 1e-6;
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
    { return "stokes"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a constant temperature of 10 degrees Celsius.
     */
    Scalar boxTemperature(const Element &element,
                       const FVElementGeometry &fvElemGeom,
                       int scvIdx) const
    {
        return 273.15 + 10; // -> 10C
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
     * \param vertex The vertex on the boundary for which the
     *               conditions needs to be specified
     */
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        values.setAllDirichlet();

        // the mass balance has to be of type outflow
        values.setOutflow(massBalanceIdx);

        if(onRightBoundary_(globalPos) &&
                globalPos[1] < this->bboxMax()[1]-eps_ && globalPos[1] > this->bboxMin()[1]+eps_)
            values.setAllOutflow();

        // set pressure at one point
        const Scalar middle = (this->bboxMax()[1] - this->bboxMin()[1])/2;
        if (onRightBoundary_(globalPos) &&
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

        const Scalar v0 = 1.0;
        // parabolic velocity profile
        values[momentumXIdx] =  v0*(globalPos[1] - this->bboxMin()[1])*(this->bboxMax()[1] - globalPos[1])
                               / (0.25*(this->bboxMax()[1] - this->bboxMin()[1])*(this->bboxMax()[1] - this->bboxMin()[1]));
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     *
     * A NEUMANN condition for the Stokes equation corresponds to:
     * \f[ -\mu \nabla {\bf v} \cdot {\bf n} + p \cdot {\bf n} = q_N \f]
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 const Intersection &is,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        values = 0.0;
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
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        initial_(values, globalPos);
    }
    // \}

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        values = 0.0;
        values[massBalanceIdx] = 1e5;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->bboxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->bboxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->bboxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->bboxMax()[1] - eps_; }

    Scalar eps_;
};

} //end namespace

#endif
