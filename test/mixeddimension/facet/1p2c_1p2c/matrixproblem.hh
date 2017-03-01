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
#ifndef DUMUX_1P2C_MATRIX_PROBLEM_HH
#define DUMUX_1P2C_MATRIX_PROBLEM_HH

#include <dumux/mixeddimension/facet/mpfa/properties.hh>
#include <dumux/mixeddimension/subproblemproperties.hh>

#include <dumux/porousmediumflow/1p2c/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include "matrixspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class OnePTwoCMpfaMatrixProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePTwoCIMatrixProblem, INHERITS_FROM(OnePTwoC));
NEW_TYPE_TAG(OnePTwoCNIMatrixProblem, INHERITS_FROM(OnePTwoCNI));
NEW_TYPE_TAG(OnePTwoCICCMpfaMatrixProblem, INHERITS_FROM(FacetCouplingBulkMpfaModel, OnePTwoCIMatrixProblem));
NEW_TYPE_TAG(OnePTwoCNICCMpfaMatrixProblem, INHERITS_FROM(FacetCouplingBulkMpfaModel, OnePTwoCNIMatrixProblem));

// Set fluid configuration
SET_TYPE_PROP(OnePTwoCIMatrixProblem, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), true>);
SET_TYPE_PROP(OnePTwoCNIMatrixProblem, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), true>);

// Set the problem property
SET_TYPE_PROP(OnePTwoCIMatrixProblem, Problem, OnePTwoCMpfaMatrixProblem<TypeTag>);
SET_TYPE_PROP(OnePTwoCNIMatrixProblem, Problem, OnePTwoCMpfaMatrixProblem<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(OnePTwoCIMatrixProblem, UseMoles, true);
SET_BOOL_PROP(OnePTwoCNIMatrixProblem, UseMoles, true);

// Set the spatial parameters
SET_TYPE_PROP(OnePTwoCIMatrixProblem, SpatialParams, OnePMatrixSpatialParams<TypeTag>);
SET_TYPE_PROP(OnePTwoCNIMatrixProblem, SpatialParams, OnePMatrixSpatialParams<TypeTag>);

// Linear solver settings
SET_TYPE_PROP(OnePTwoCIMatrixProblem, LinearSolver, SuperLUBackend<TypeTag>);
SET_TYPE_PROP(OnePTwoCNIMatrixProblem, LinearSolver, SuperLUBackend<TypeTag>);

// Enable gravity
SET_BOOL_PROP(OnePTwoCIMatrixProblem, ProblemEnableGravity, false);
SET_BOOL_PROP(OnePTwoCNIMatrixProblem, ProblemEnableGravity, false);

// Solution-independent tensors
SET_BOOL_PROP(OnePTwoCICCMpfaMatrixProblem, SolutionDependentAdvection, false);
SET_BOOL_PROP(OnePTwoCNICCMpfaMatrixProblem, SolutionDependentAdvection, false);

// change mpfa method
// SET_PROP(OnePTwoCICCMpfaMatrixProblem, MpfaMethod) { static const MpfaMethods value = MpfaMethods::lMethod; };
// SET_PROP(OnePTwoCNICCMpfaMatrixProblem, MpfaMethod) { static const MpfaMethods value = MpfaMethods::lMethod; };
}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the one-phase model
 */

template <class TypeTag>
class OnePTwoCMpfaMatrixProblem : public ImplicitPorousMediaProblem<TypeTag>
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
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,

        // indices of the equations
        conti0EqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx
    };

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;


public:
    OnePTwoCMpfaMatrixProblem(TimeManager &timeManager, const GridView &gridView)
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
        openFractureSegments_.resize(39);
        openFractureSegments_[0] = std::make_pair(GlobalPosition({0, 333.7332}),
                                                  GlobalPosition({32.2042, 309.3757}) - GlobalPosition({0, 333.7332}));
        openFractureSegments_[1] = std::make_pair(GlobalPosition({32.2042, 309.3757}),
                                                  GlobalPosition({53.0462, 293.612}) - GlobalPosition({32.2042, 309.3757}));
        openFractureSegments_[2] = std::make_pair(GlobalPosition({53.0462, 293.612}),
                                                  GlobalPosition({314.7789, 95.6516}) - GlobalPosition({53.0462, 293.612}));
        openFractureSegments_[3] = std::make_pair(GlobalPosition({314.7789, 95.6516}),
                                                  GlobalPosition({387.6907, 40.5051}) - GlobalPosition({314.7789, 95.6516}));
        openFractureSegments_[4] = std::make_pair(GlobalPosition({387.6907, 40.5051}),
                                                  GlobalPosition({441.2444, 0}) - GlobalPosition({387.6907, 40.5051}));
        openFractureSegments_[5] = std::make_pair(GlobalPosition({32.2042, 309.3757}),
                                                  GlobalPosition({51.6592, 371.6782}) - GlobalPosition({32.2042, 309.3757}));
        openFractureSegments_[6] = std::make_pair(GlobalPosition({51.6592, 371.6782}),
                                                  GlobalPosition({68.2924, 424.9442}) - GlobalPosition({51.6592, 371.6782}));
        openFractureSegments_[7] = std::make_pair(GlobalPosition({68.2924, 424.9442}),
                                                  GlobalPosition({82.9736, 471.959}) - GlobalPosition({68.2924, 424.9442}));
        openFractureSegments_[8] = std::make_pair(GlobalPosition({82.9736, 471.959}),
                                                  GlobalPosition({92.0884, 501.1482}) - GlobalPosition({82.9736, 471.959}));
        openFractureSegments_[9] = std::make_pair(GlobalPosition({92.0884, 501.1482}),
                                                  GlobalPosition({96.36, 514.8274}) - GlobalPosition({92.0884, 501.1482}));
        openFractureSegments_[10] = std::make_pair(GlobalPosition({92.0884, 501.1482}),
                                                   GlobalPosition({153.1981, 465.6789}) - GlobalPosition({92.0884, 501.1482}));
        openFractureSegments_[11] = std::make_pair(GlobalPosition({153.1981, 465.6789}),
                                                   GlobalPosition({165.8919, 458.3112}) - GlobalPosition({153.1981, 465.6789}));
        openFractureSegments_[12] = std::make_pair(GlobalPosition({165.8919, 458.3112}),
                                                   GlobalPosition({184.3938, 447.5723}) - GlobalPosition({165.8919, 458.3112}));
        openFractureSegments_[13] = std::make_pair(GlobalPosition({184.3938, 447.5723}),
                                                   GlobalPosition({270.148, 397.7989}) - GlobalPosition({184.3938, 447.5723}));
        openFractureSegments_[14] = std::make_pair(GlobalPosition({270.148, 397.7989}),
                                                   GlobalPosition({356.0905, 347.6924}) - GlobalPosition({270.148, 397.7989}));
        openFractureSegments_[15] = std::make_pair(GlobalPosition({13.4389, 342.6924}),
                                                   GlobalPosition({51.6592, 371.6782}) - GlobalPosition({13.4389, 342.6924}));
        openFractureSegments_[16] = std::make_pair(GlobalPosition({51.6592, 371.6782}),
                                                   GlobalPosition({99.3531, 407.8488}) - GlobalPosition({51.6592, 371.6782}));
        openFractureSegments_[17] = std::make_pair(GlobalPosition({99.3531, 407.8488}),
                                                   GlobalPosition({138.8863, 437.8304}) - GlobalPosition({99.3531, 407.8488}));
        openFractureSegments_[18] = std::make_pair(GlobalPosition({138.8863, 437.8304}),
                                                   GlobalPosition({153.037, 448.5621}) - GlobalPosition({138.8863, 437.8304}));
        openFractureSegments_[19] = std::make_pair(GlobalPosition({153.037, 448.5621}),
                                                   GlobalPosition({165.8919, 458.3112}) - GlobalPosition({153.037, 448.5621}));
        openFractureSegments_[20] = std::make_pair(GlobalPosition({165.8919, 458.3112}),
                                                   GlobalPosition({182.9367, 471.2377}) - GlobalPosition({165.8919, 458.3112}));
        openFractureSegments_[21] = std::make_pair(GlobalPosition({182.9367, 471.2377}),
                                                   GlobalPosition({338.3707, 589.1173}) - GlobalPosition({182.9367, 471.2377}));
        openFractureSegments_[22] = std::make_pair(GlobalPosition({338.3707, 589.1173}),
                                                   GlobalPosition({347.1719, 595.792}) - GlobalPosition({338.3707, 589.1173}));
        openFractureSegments_[23] = std::make_pair(GlobalPosition({337.5811, 600}),
                                                   GlobalPosition({338.3707, 589.1173}) - GlobalPosition({337.5811, 600}));
        openFractureSegments_[24] = std::make_pair(GlobalPosition({338.3707, 589.1173}),
                                                   GlobalPosition({346.1441, 481.9792}) - GlobalPosition({338.3707, 589.1173}));
        openFractureSegments_[25] = std::make_pair(GlobalPosition({346.1441, 481.9792}),
                                                   GlobalPosition({351.2625, 411.4348}) - GlobalPosition({346.1441, 481.9792}));
        openFractureSegments_[26] = std::make_pair(GlobalPosition({351.2625, 411.4348}),
                                                   GlobalPosition({354.1847, 371.1592}) - GlobalPosition({351.2625, 411.4348}));
        openFractureSegments_[27] = std::make_pair(GlobalPosition({354.1847, 371.1592}),
                                                   GlobalPosition({354.8302, 360.4798}) - GlobalPosition({354.1847, 371.1592}));
        openFractureSegments_[28] = std::make_pair(GlobalPosition({347.1719, 374.05}),
                                                   GlobalPosition({354.1847, 371.1592}) - GlobalPosition({347.1719, 374.05}));
        openFractureSegments_[29] = std::make_pair(GlobalPosition({354.1847, 371.1592}),
                                                   GlobalPosition({372.7787, 363.4945}) - GlobalPosition({354.1847, 371.1592}));
        openFractureSegments_[30] = std::make_pair(GlobalPosition({372.7787, 363.4945}),
                                                   GlobalPosition({455.4811, 329.4034}) - GlobalPosition({372.7787, 363.4945}));
        openFractureSegments_[31] = std::make_pair(GlobalPosition({455.4811, 329.4034}),
                                                   GlobalPosition({472.2415, 322.4945}) - GlobalPosition({455.4811, 329.4034}));
        openFractureSegments_[32] = std::make_pair(GlobalPosition({472.2415, 322.4945}),
                                                   GlobalPosition({484.9833, 317.2422}) - GlobalPosition({472.2415, 322.4945}));
        openFractureSegments_[33] = std::make_pair(GlobalPosition({484.9833, 317.2422}),
                                                   GlobalPosition({493.6741, 313.6597}) - GlobalPosition({484.9833, 317.2422}));
        openFractureSegments_[34] = std::make_pair(GlobalPosition({493.6741, 313.6597}),
                                                   GlobalPosition({534.3413, 296.8961}) - GlobalPosition({493.6741, 313.6597}));
        openFractureSegments_[35] = std::make_pair(GlobalPosition({534.3413, 296.8961}),
                                                   GlobalPosition({558.9313, 286.7597}) - GlobalPosition({534.3413, 296.8961}));
        openFractureSegments_[36] = std::make_pair(GlobalPosition({558.9313, 286.7597}),
                                                   GlobalPosition({566.1566, 283.7814}) - GlobalPosition({558.9313, 286.7597}));
        openFractureSegments_[37] = std::make_pair(GlobalPosition({566.1566, 283.7814}),
                                                   GlobalPosition({640.5882, 253.0996}) - GlobalPosition({566.1566, 283.7814}));
        openFractureSegments_[38] = std::make_pair(GlobalPosition({44.7964, 528.5974}),
                                                   GlobalPosition({92.0884, 501.1482}) - GlobalPosition({44.7964, 528.5974}));

        barrierSegments_.resize(5);
        barrierSegments_[0] = std::make_pair(GlobalPosition({269.1347, 147.7814}),
                                             GlobalPosition({269.8085, 314.0322}) - GlobalPosition({269.1347, 147.7814}));
        barrierSegments_[1] = std::make_pair(GlobalPosition({269.8085, 314.0322}),
                                             GlobalPosition({270.148, 397.7989}) - GlobalPosition({269.8085, 314.0322}));
        barrierSegments_[2] = std::make_pair(GlobalPosition({270.148, 397.7989}),
                                             GlobalPosition({270.4584, 474.3837}) - GlobalPosition({270.148, 397.7989}));
        barrierSegments_[3] = std::make_pair(GlobalPosition({270.4584, 474.3837}),
                                             GlobalPosition({270.5362, 493.5848}) - GlobalPosition({270.4584, 474.3837}));
        barrierSegments_[4] = std::make_pair(GlobalPosition({270.5362, 493.5848}),
                                             GlobalPosition({270.6623, 524.7026}) - GlobalPosition({270.5362, 493.5848}));
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
    template<class T = TypeTag>
    typename std::enable_if<std::is_same<T, TTAG(OnePTwoCICCMpfaMatrixProblem)>::value, PrimaryVariables>::type
    dirichletAtPos(const GlobalPosition& globalPos) const
    {
        auto values = initialAtPos(globalPos);
        if (globalPos[0] < eps_ && globalPos[1] > 250 && globalPos[1] < 400)
            values[massOrMoleFracIdx] = 2e-5;
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume (non-isothermal case).
     */
    template<class T = TypeTag>
    typename std::enable_if<std::is_same<T, TTAG(OnePTwoCNICCMpfaMatrixProblem)>::value, PrimaryVariables>::type
    dirichletAtPos(const GlobalPosition& globalPos) const
    {
      auto values = initialAtPos(globalPos);
      if (globalPos[0] < eps_ && globalPos[1] > 250 && globalPos[1] < 400)
      {
          values[massOrMoleFracIdx] = 2e-5;
          values[Indices::temperatureIdx] += 30;
      }
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
    typename std::enable_if<std::is_same<T, TTAG(OnePTwoCICCMpfaMatrixProblem)>::value, PrimaryVariables>::type
    initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx] = 10e5 - 8e5*globalPos[0]/this->bBoxMax()[0];
        return values;
    }

    /*!
     * \brief Evaluate the initial value for a control volume (non-isothermal case)
     */
    template<class T = TypeTag>
    typename std::enable_if<std::is_same<T, TTAG(OnePTwoCNICCMpfaMatrixProblem)>::value, PrimaryVariables>::type
    initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx] = 10e5 - 8e5*globalPos[0]/this->bBoxMax()[0];
        values[Indices::temperatureIdx] = temperature();
        return values;
    }

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    bool isOpenFracture(const GlobalPosition& globalPos) const
    {
        for (const auto& segment : openFractureSegments_)
            if (isOnSegment_(segment, globalPos))
                return true;
        return false;
    }

    bool isBarrier(const GlobalPosition& globalPos) const
    {
        for (const auto& segment : barrierSegments_)
            if (isOnSegment_(segment, globalPos))
                return true;
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
