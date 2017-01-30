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
 *
 * \brief A discrete fracture network embedded in an impermeable matrix.
 *        The fracture is a 2D network embedded in 3D.
 */
#ifndef DUMUX_TWOP_FRACTURE_TEST_PROBLEM_HH
#define DUMUX_TWOP_FRACTURE_TEST_PROBLEM_HH

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>
#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/implicit/cellcentered/propertydefaults.hh>

#include "fracturespatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class FractureProblem;

namespace Properties
{
NEW_TYPE_TAG(FractureProblem, INHERITS_FROM(TwoP, FractureSpatialParams));
NEW_TYPE_TAG(FractureBoxProblem, INHERITS_FROM(BoxModel, FractureProblem));
NEW_TYPE_TAG(FractureCCProblem, INHERITS_FROM(CCTpfaModel, FractureProblem));
NEW_TYPE_TAG(FractureCCMpfaProblem, INHERITS_FROM(CCMpfaModel, FractureProblem));

SET_BOOL_PROP(FractureProblem, EnableGlobalFVGeometryCache, true);
SET_BOOL_PROP(FractureProblem, EnableGlobalVolumeVariablesCache, true);
SET_BOOL_PROP(FractureProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(FractureProblem, SolutionDependentAdvection, false);

#if HAVE_DUNE_FOAMGRID
SET_TYPE_PROP(FractureProblem, Grid, Dune::FoamGrid<2, 3>);
#else
SET_TYPE_PROP(FractureProblem, Grid, Dune::YaspGrid<3>);
#endif

// Set the problem property
SET_TYPE_PROP(FractureProblem, Problem, Dumux::FractureProblem<TypeTag>);

// Set the wetting phase
SET_PROP(FractureProblem, WettingPhase)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar>>;
};

// Set the non-wetting phase
SET_PROP(FractureProblem, NonwettingPhase)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::LiquidPhase<Scalar, DNAPL<Scalar>>;
};

// Linear solver settings
SET_TYPE_PROP(FractureProblem, LinearSolver, Dumux::ILU0BiCGSTABBackend<TypeTag>);

SET_BOOL_PROP(FractureProblem, ProblemEnableGravity, false);
}

/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeablility which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability. Note that
 * this problem is discretized using only two dimensions, so from the
 * point of view of the two-phase model, the depth of the domain
 * implicitly is 1 m everywhere.
 *
 * On the top and the bottom of the domain neumann boundary conditions
 * are used, while dirichlet conditions apply on the left and right
 * boundaries.
 *
 * DNAPL is injected at the top boundary from 3m to 4m at a rate of
 * 0.04 kg/(s m^2), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * The dirichlet boundaries on the left boundary is the hydrostatic
 * pressure scaled by a factor of 1.125, while on the right side it is
 * just the hydrostatic pressure. The DNAPL saturation on both sides
 * is zero.
 *
 * This problem uses the \ref TwoPModel.
 *
 * This problem should typically be simulated until \f$t_{\text{end}}
 * \approx 20\,000\;s\f$ is reached. A good choice for the initial time step
 * size is \f$t_{\text{inital}} = 250\;s\f$.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2p -parameterFile test_box2p.input</tt> or
 * <tt>./test_cc2p -parameterFile test_cc2p.input</tt>
 */
template <class TypeTag>
class FractureProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using WettingPhase = typename GET_PROP_TYPE(TypeTag, WettingPhase);
    using NonwettingPhase = typename GET_PROP_TYPE(TypeTag, NonwettingPhase);

    enum {

        // primary variable indices
        pwIdx = Indices::pwIdx,
        snIdx = Indices::snIdx,

        // equation indices
        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,


        // world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };


    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using ReferenceElements = typename Dune::ReferenceElements<Scalar, dim>;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    FractureProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    {
        return name_;
    }

    /*!
     * \brief User defined output after the time integration
     *
     * Will be called diretly after the time integration.
     */
    void postTimeStep()
    {
        // Calculate storage terms
        PrimaryVariables storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout<<"Storage: " << storage << std::endl;
        }
    }

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * This problem assumes a uniform temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; }

    /*!
     * \brief Returns the source term
     *
     * \param values Stores the source values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param globalPos The global position
     */
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        return PrimaryVariables(0.0);
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
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setAllDirichlet();
        if (onInlet_(globalPos))
            values.setAllNeumann();
        if (globalPos[2] > 1.0 - eps_ || globalPos[2] < eps_)
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values[pwIdx] = 1e5;
        values[snIdx] = 0.0;
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        if (onInlet_(globalPos)) {
            values[contiNEqIdx] = -0.04; // kg / (m * s)
        }
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{


    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return dirichletAtPos(globalPos);
    }
    // \}

private:
    bool onInlet_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < eps_ && globalPos[1] > -0.5 - eps_;
    }

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;
};

} //end namespace Dumux

#endif
