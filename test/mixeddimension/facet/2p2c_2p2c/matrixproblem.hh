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

#ifndef DUMUX_2P2C_MATRIX_PROBLEM_HH
#define DUMUX_2P2C_MATRIX_PROBLEM_HH

#include <dumux/mixeddimension/facet/mpfa/properties.hh>
#include <dumux/mixeddimension/subproblemproperties.hh>

#include <dumux/porousmediumflow/2p2c/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include "matrixspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class TwoPTwoCMpfaMatrixProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoPTwoCIMatrixProblem, INHERITS_FROM(TwoPTwoC));
NEW_TYPE_TAG(TwoPTwoCNIMatrixProblem, INHERITS_FROM(TwoPTwoCNI));
NEW_TYPE_TAG(TwoPTwoCICCMpfaMatrixProblem, INHERITS_FROM(FacetCouplingBulkMpfaModel, TwoPTwoCIMatrixProblem, TwoPTwoCMatrixSpatialParams));
NEW_TYPE_TAG(TwoPTwoCNICCMpfaMatrixProblem, INHERITS_FROM(FacetCouplingBulkMpfaModel, TwoPTwoCNIMatrixProblem, TwoPTwoCMatrixSpatialParams));

// Set fluid configuration
SET_TYPE_PROP(TwoPTwoCIMatrixProblem, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), true>);
SET_TYPE_PROP(TwoPTwoCNIMatrixProblem, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), true>);

// Set the problem property
SET_TYPE_PROP(TwoPTwoCIMatrixProblem, Problem, TwoPTwoCMpfaMatrixProblem<TypeTag>);
SET_TYPE_PROP(TwoPTwoCNIMatrixProblem, Problem, TwoPTwoCMpfaMatrixProblem<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(TwoPTwoCIMatrixProblem, UseMoles, true);
SET_BOOL_PROP(TwoPTwoCNIMatrixProblem, UseMoles, true);

// Linear solver settings
SET_TYPE_PROP(TwoPTwoCIMatrixProblem, LinearSolver, SuperLUBackend<TypeTag>);
SET_TYPE_PROP(TwoPTwoCNIMatrixProblem, LinearSolver, SuperLUBackend<TypeTag>);

// Enable gravity
SET_BOOL_PROP(TwoPTwoCIMatrixProblem, ProblemEnableGravity, false);
SET_BOOL_PROP(TwoPTwoCNIMatrixProblem, ProblemEnableGravity, false);

// Solution-independent tensors
SET_BOOL_PROP(TwoPTwoCICCMpfaMatrixProblem, SolutionDependentAdvection, false);
SET_BOOL_PROP(TwoPTwoCNICCMpfaMatrixProblem, SolutionDependentAdvection, false);

// change mpfa method
// SET_PROP(TwoPTwoCICCMpfaMatrixProblem, MpfaMethod) { static const MpfaMethods value = MpfaMethods::lMethod; };
// SET_PROP(TwoPTwoCNICCMpfaMatrixProblem, MpfaMethod) { static const MpfaMethods value = MpfaMethods::lMethod; };
}

template <class TypeTag>
class TwoPTwoCMpfaMatrixProblem : public ImplicitPorousMediaProblem<TypeTag>
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
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        contiH2OEqIdx = Indices::contiWEqIdx,
        contiN2EqIdx = Indices::contiNEqIdx
    };

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;


public:
    TwoPTwoCMpfaMatrixProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name) + "_matrix";
        eps_ = 1e-6;

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions"<<std::endl;
        else
            std::cout<<"problem uses mass fractions"<<std::endl;

        // initialize the fracture segments for the given grid
        openFractureSegments_.resize(4);
        openFractureSegments_[0] = std::make_pair(GlobalPosition({0.0, 630.9}),
                                                  GlobalPosition({316.7, 592.4}) - GlobalPosition({0.0, 630.9}));
        openFractureSegments_[1] = std::make_pair(GlobalPosition({316.7, 592.4}),
                                                  GlobalPosition({687.6, 535.6}) - GlobalPosition({316.7, 592.4}));
        openFractureSegments_[2] = std::make_pair(GlobalPosition({687.6, 535.6}),
                                                  GlobalPosition({878.3, 506.4}) - GlobalPosition({687.6, 535.6}));
        openFractureSegments_[3] = std::make_pair(GlobalPosition({878.3, 506.4}),
                                                  GlobalPosition({878.3, 506.4}) - GlobalPosition({314.7789, 95.6516}));

        barrierSegments_.resize(5);
        barrierSegments_[0] = std::make_pair(GlobalPosition({352.6, 900}),
                                             GlobalPosition({316.7, 592.4}) - GlobalPosition({352.6, 900}));
        barrierSegments_[1] = std::make_pair(GlobalPosition({316.7, 592.4}),
                                             GlobalPosition({304.3, 481.8}) - GlobalPosition({316.7, 592.4}));
        barrierSegments_[2] = std::make_pair(GlobalPosition({748, 716.8}),
                                             GlobalPosition({687.6, 535.6}) - GlobalPosition({748, 716.8}));
        barrierSegments_[3] = std::make_pair(GlobalPosition({687.6, 535.6}),
                                             GlobalPosition({596.2, 261.7}) - GlobalPosition({687.6, 535.6}));
        barrierSegments_[4] = std::make_pair(GlobalPosition({596.2, 261.7}),
                                             GlobalPosition({538.8, 80}) - GlobalPosition({596.2, 261.7}));
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
    { return 273.15 + 20; }

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
        if (globalPos[0] > this->bBoxMax()[0] - eps_)
            values.setAllDirichlet();

        if (couplingManager().isInteriorBoundary(element, scvf))
        {
          if (isOpenFracture(scvf.ipGlobal()))
             values.setAllDirichlet();
          else
             values.setAllNeumann();
        }

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
    template<class T = TypeTag>
    typename std::enable_if<std::is_same<T, TTAG(TwoPTwoCICCMpfaMatrixProblem)>::value, PrimaryVariables>::type
    dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume (non-isothermal case).
     */
    template<class T = TypeTag>
    typename std::enable_if<std::is_same<T, TTAG(TwoPTwoCNICCMpfaMatrixProblem)>::value, PrimaryVariables>::type
    dirichletAtPos(const GlobalPosition& globalPos) const
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
    typename std::enable_if<std::is_same<T, TTAG(TwoPTwoCICCMpfaMatrixProblem)>::value, PrimaryVariables>::type
    initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = 1e5;
        return values;
    }

    /*!
     * \brief Evaluate the initial value for a control volume (non-isothermal case)
     */
    template<class T = TypeTag>
    typename std::enable_if<std::is_same<T, TTAG(TwoPTwoCNICCMpfaMatrixProblem)>::value, PrimaryVariables>::type
    initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = 1e5;
        values[Indices::temperatureIdx] = temperature();
        return values;
    }

    /*!
     * \brief Returns the initial phase state for a control volume.
     *
     * \param scv The sub-control volume
     */
    int initialPhasePresence(const SubControlVolume& scv) const
    { return Indices::wPhaseOnly; }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    bool isOpenFracture(const GlobalPosition& globalPos) const
    {
        // for (const auto& segment : openFractureSegments_)
        //     if (isOnSegment_(segment, globalPos))
        //         return true;
        return false;
    }

    bool isBarrier(const GlobalPosition& globalPos) const
    {
        // for (const auto& segment : barrierSegments_)
        //     if (isOnSegment_(segment, globalPos))
        //         return true;
        return false;
    }

private:

    bool isOnSegment_(const std::pair<GlobalPosition, GlobalPosition>& segment, const GlobalPosition& globalPos) const
    {
        const auto v = globalPos-segment.first;
        const auto vNorm = v.two_norm();
        auto vUnit = v;
        vUnit /= vNorm;

        const auto dNorm = segment.second.two_norm();
        auto dUnit = segment.second;
        dUnit /= dNorm;

        using std::abs;
        if (abs(vUnit*dUnit - 1.0) < eps_ && vNorm <= dNorm)
            return true;
        return false;
    }

    std::string name_;
    Scalar eps_;
    std::shared_ptr<CouplingManager> couplingManager_;
    std::vector<std::pair<GlobalPosition, GlobalPosition>> openFractureSegments_;
    std::vector<std::pair<GlobalPosition, GlobalPosition>> barrierSegments_;
};
} //end namespace

#endif
