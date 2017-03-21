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
 * \brief A test problem for the one-dimensional two-phase fracture model
 */
#ifndef DUMUX_2P_FRACTURE_PROBLEM_HH
#define DUMUX_2P_FRACTURE_PROBLEM_HH

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/mixeddimension/subproblemproperties.hh>

#include <dumux/porousmediumflow/1p2c/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>
#include "fracturespatialparams.hh"

namespace Dumux
{

//! Forward declaration of the problem class
template <class TypeTag>
class TwoPFractureProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoPFractureProblem, INHERITS_FROM(TwoP));
NEW_TYPE_TAG(TwoPNIFractureProblem, INHERITS_FROM(TwoPNI));
NEW_TYPE_TAG(TwoPCCMpfaFractureProblem, INHERITS_FROM(CCMpfaModel, TwoPFractureProblem, FractureSpatialParams));
NEW_TYPE_TAG(TwoPNICCMpfaFractureProblem, INHERITS_FROM(CCTpfaModel, TwoPNIFractureProblem, FractureSpatialParams));

// Set the wetting phase
SET_PROP(TwoPFractureProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(TwoPFractureProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, DNAPL<Scalar> > type;
};

// Set the problem property
SET_TYPE_PROP(TwoPFractureProblem, Problem, TwoPFractureProblem<TypeTag>);
SET_TYPE_PROP(TwoPNIFractureProblem, Problem, TwoPFractureProblem<TypeTag>);

// Linear solver settings
SET_TYPE_PROP(TwoPFractureProblem, LinearSolver, SuperLUBackend<TypeTag>);
SET_TYPE_PROP(TwoPNIFractureProblem, LinearSolver, SuperLUBackend<TypeTag>);

// Enable gravity
SET_BOOL_PROP(TwoPFractureProblem, ProblemEnableGravity, true);
SET_BOOL_PROP(TwoPNIFractureProblem, ProblemEnableGravity, true);

// Solution-independent tensors
SET_BOOL_PROP(TwoPCCMpfaFractureProblem, SolutionDependentAdvection, false);
SET_BOOL_PROP(TwoPNICCMpfaFractureProblem, SolutionDependentAdvection, false);

SET_BOOL_PROP(TwoPCCMpfaFractureProblem, EnableGlobalFVGeometryCache, true);
SET_BOOL_PROP(TwoPNICCMpfaFractureProblem, EnableGlobalFVGeometryCache, true);
SET_BOOL_PROP(TwoPCCMpfaFractureProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(TwoPNICCMpfaFractureProblem, EnableGlobalFluxVariablesCache, true);
}

template <class TypeTag>
class TwoPFractureProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
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
    TwoPFractureProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "_fracture";
        eps_ = 1e-6;
    }

    /*!
     * \brief The problem name.
     *        This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 50C

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     */
    Scalar extrusionFactorAtPos(const GlobalPosition &globalPos) const
    {
        //! return the user-specified fracture aperture
        static const Scalar a = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, FractureAperture);
        return a;
    }


    /*!
     * \brief Return the sources within the domain.
     */
    PrimaryVariables source(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    {
        // we have only sources coming from the bulk domain
        return couplingManager().evalSourcesFromBulk(element, fvGeometry, elemVolVars, scv);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be used.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume (isothermal case).
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }

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
    typename std::enable_if<std::is_same<T, TTAG(TwoPCCMpfaFractureProblem)>::value, PrimaryVariables>::type
    initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 2e5 + (this->bBoxMax()[2] - globalPos[2])*9.81*1000;
        return values;
    }

    /*!
     * \brief Evaluate the initial value for a control volume (non-isothermal case)
     */
    template<class T = TypeTag>
    typename std::enable_if<std::is_same<T, TTAG(TwoPNICCMpfaFractureProblem)>::value, PrimaryVariables>::type
    initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 2e5 + (170 - globalPos[2])*9.81*1000;
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
