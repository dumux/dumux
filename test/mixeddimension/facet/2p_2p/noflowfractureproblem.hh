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
 * \brief A test problem for the one-dimensional single-phase fracture model
 */
#ifndef DUMUX_1P2C_FRACTURE_PROBLEM_HH
#define DUMUX_1P2C_FRACTURE_PROBLEM_HH

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/mixeddimension/subproblemproperties.hh>

#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>
#include "fracturespatialparams.hh"

#define COMPLEXPHYSICS 0

namespace Dumux
{

//! Forward declaration of the problem class
template <class TypeTag>
class NoFlowDomainFractureProblem;

namespace Properties
{
NEW_TYPE_TAG(NoFlowDomainFractureProblem, INHERITS_FROM(CCTpfaModel, TwoP, FractureSpatialParams));

// Set fluid configuration
#if COMPLEXPHYSICS
SET_TYPE_PROP(NoFlowDomainFractureProblem, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), true>);
#else
SET_TYPE_PROP(NoFlowDomainFractureProblem, FluidSystem, FluidSystems::TwoPImmiscible<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                                                     FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar), SimpleH2O<typename GET_PROP_TYPE(TypeTag, Scalar)>>,
                                                                                     FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar), DNAPL<typename GET_PROP_TYPE(TypeTag, Scalar)>>>);
#endif

// Set the problem property
SET_TYPE_PROP(NoFlowDomainFractureProblem, Problem, NoFlowDomainFractureProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(NoFlowDomainFractureProblem, SpatialParams, TwoPFractureSpatialParams<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(NoFlowDomainFractureProblem, UseMoles, true);

// Linear solver settings
SET_TYPE_PROP(NoFlowDomainFractureProblem, LinearSolver, SuperLUBackend<TypeTag>);

// Enable gravity
SET_BOOL_PROP(NoFlowDomainFractureProblem, ProblemEnableGravity, false);

// Solution-independent tensors
SET_BOOL_PROP(NoFlowDomainFractureProblem, SolutionDependentAdvection, false);
}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the one-phase model:
 */
template <class TypeTag>
class NoFlowDomainFractureProblem : public ImplicitPorousMediaProblem<TypeTag>
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

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    NoFlowDomainFractureProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "_fracture";
        eps_ = 1e-6;

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions"<<std::endl;
        else
            std::cout<<"problem uses mass fractions"<<std::endl;
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
        auto sources = couplingManager().evalSourcesFromBulk(element, fvGeometry, elemVolVars, scv);
        sources /= scv.volume()*elemVolVars[scv].extrusionFactor();
        return sources;
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
     * \brief Evaluate the boundary conditions for a dirichlet control volume.
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
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pwIdx] = 1e5;
        return values;
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    //! Get the coupling manager
    CouplingManager& couplingManager()
    { return *couplingManager_; }

private:
    std::string name_;
    Scalar eps_;
    std::shared_ptr<CouplingManager> couplingManager_;
};
} //end namespace

#endif
