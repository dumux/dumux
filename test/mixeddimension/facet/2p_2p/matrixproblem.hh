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
 * \brief A test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 */
#ifndef DUMUX_2P_MATRIX_PROBLEM_HH
#define DUMUX_2P_MATRIX_PROBLEM_HH

#include <dumux/mixeddimension/facet/mpfa/properties.hh>
#include <dumux/mixeddimension/subproblemproperties.hh>

#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/dnapl.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include "matrixspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class TwoPMpfaMatrixProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoPMatrixProblem, INHERITS_FROM(TwoP));
NEW_TYPE_TAG(TwoPNIMatrixProblem, INHERITS_FROM(TwoPNI));
NEW_TYPE_TAG(TwoPCCMpfaMatrixProblem, INHERITS_FROM(FacetCouplingBulkMpfaModel, TwoPMatrixProblem, MatrixSpatialParams));
NEW_TYPE_TAG(TwoPNICCMpfaMatrixProblem, INHERITS_FROM(FacetCouplingBulkMpfaModel, TwoPNIMatrixProblem, MatrixSpatialParams));

// Set fluid configuration
SET_PROP(TwoPMatrixProblem, FluidSystem)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::TwoPImmiscible<Scalar,
                                              FluidSystems::LiquidPhase<Scalar, H2O<Scalar>>,
                                              FluidSystems::LiquidPhase<Scalar, DNAPL<Scalar>>>;
};
SET_PROP(TwoPNIMatrixProblem, FluidSystem)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
public:
    using type = FluidSystems::TwoPImmiscible<Scalar,
                                              FluidSystems::LiquidPhase<Scalar, H2O<Scalar>>,
                                              FluidSystems::LiquidPhase<Scalar, DNAPL<Scalar>>>;
};

// Set the problem property
SET_TYPE_PROP(TwoPMatrixProblem, Problem, TwoPMpfaMatrixProblem<TypeTag>);
SET_TYPE_PROP(TwoPNIMatrixProblem, Problem, TwoPMpfaMatrixProblem<TypeTag>);

// Linear solver settings
SET_TYPE_PROP(TwoPMatrixProblem, LinearSolver, SuperLUBackend<TypeTag>);
SET_TYPE_PROP(TwoPNIMatrixProblem, LinearSolver, SuperLUBackend<TypeTag>);

// Enable gravity
SET_BOOL_PROP(TwoPMatrixProblem, ProblemEnableGravity, true);
SET_BOOL_PROP(TwoPNIMatrixProblem, ProblemEnableGravity, true);

// Solution-independent tensors
SET_BOOL_PROP(TwoPCCMpfaMatrixProblem, SolutionDependentAdvection, false);
SET_BOOL_PROP(TwoPNICCMpfaMatrixProblem, SolutionDependentAdvection, false);
}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the one-phase model
 */

template <class TypeTag>
class TwoPMpfaMatrixProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using GlobalProblemTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalProblemTypeTag, CouplingManager);

    // copy some indices for convenience
    enum
    {
        // primary variable indices
        pwIdx = Indices::pwIdx,
        snIdx = Indices::snIdx,

        // equation indices
        contiWEqIdx = Indices::contiWEqIdx,
        contiNEqIdx = Indices::contiNEqIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
    };

    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;


public:
    TwoPMpfaMatrixProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "_matrix";
        eps_ = 1e-6;
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; }

    /*!
     * \brief Return the sources within the domain.
     */
    PrimaryVariables sourceAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     */
    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        const auto globalPos = scvf.ipGlobal();

        values.setAllNeumann();
        if (globalPos[0] < eps_ || globalPos[0] > this->bBoxMax()[0] - eps_)
            values.setAllDirichlet();

        if (couplingManager().isInteriorBoundary(element, scvf))
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Specifies if a given intersection is on an interior boundary
     */
    bool isInteriorBoundary(const Element& element, const Intersection& is) const
    { return couplingManager().isInteriorBoundary(element, is); }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume (isothermal case).
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        auto values = initialAtPos(globalPos);

        if (globalPos[0] < eps_ && globalPos[2] > 120 && globalPos[1] < 355 && globalPos[1] > 305)
            values[snIdx] = 0.5;

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     */
    PrimaryVariables neumann(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf) const
    { return PrimaryVariables(0.0); }

    /*!
     * \brief Evaluate the initial value for a control volume (isothermal case)
     */
    template<class T = TypeTag>
    typename std::enable_if<std::is_same<T, TTAG(TwoPCCMpfaMatrixProblem)>::value, PrimaryVariables>::type
    initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 8e5 - 7e5*globalPos[0]/this->bBoxMax()[0] + (this->bBoxMax()[2] - globalPos[2])*9.81*1000;
        return values;
    }

    /*!
     * \brief Evaluate the initial value for a control volume (non-isothermal case)
     */
    template<class T = TypeTag>
    typename std::enable_if<std::is_same<T, TTAG(TwoPNICCMpfaMatrixProblem)>::value, PrimaryVariables>::type
    initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 8e5 - 7e5*globalPos[0]/this->bBoxMax()[0] + (this->bBoxMax()[2] - globalPos[2])*9.81*1000;
        values[Indices::temperatureIdx] = temperature();
        return values;
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::string name_;
    Scalar eps_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} //end namespace

#endif
