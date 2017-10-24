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
#ifndef DUMUX_1PDMP_PROBLEM_HH
#define DUMUX_1PDMP_PROBLEM_HH

#include <dumux/porousmediumflow/implicit/problem.hh>
#if PROBLEM==1
#include <dumux/porousmediumflow/1p/mimetic/model.hh>
#elif PROBLEM==2
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/1p/implicit/model.hh>
#elif PROBLEM==4
#include <dumux/porousmediumflow/1p/mimetic/model.hh>
#include "wellmimeticlocalresidual.hh"
#elif PROBLEM==5
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include "wellcclocalresidual.hh"
#endif

#include <dumux/material/components/unit.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include "1pdmpspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class OnePDMPProblem;

namespace Capabilities
{
    template<class TypeTag>
    struct isStationary<OnePDMPProblem<TypeTag>>
    { static const bool value = true; };
}

namespace Properties
{
#if PROBLEM==1
NEW_TYPE_TAG(OnePDMPProblem, INHERITS_FROM(OnePMimetic));
#elif PROBLEM==2
NEW_TYPE_TAG(OnePDMPProblem, INHERITS_FROM(CCMpfaModel, OneP));

SET_PROP(OnePDMPProblem, MpfaMethod)
{
    static const MpfaMethods value = MpfaMethods::oMethod;
};
#elif PROBLEM==4
NEW_TYPE_TAG(OnePDMPProblem, INHERITS_FROM(OnePMimetic));

//! The local residual function
SET_TYPE_PROP(OnePDMPProblem, LocalResidual, WellImmiscibleMimeticLocalResidual<TypeTag>);
#elif PROBLEM==5
NEW_TYPE_TAG(OnePDMPProblem, INHERITS_FROM(CCMpfaModel, OneP));

SET_PROP(OnePDMPProblem, MpfaMethod)
{
    static const MpfaMethods value = MpfaMethods::oMethod;
};

//! The local residual function
SET_TYPE_PROP(OnePDMPProblem, LocalResidual, WellImmiscibleLocalResidual<TypeTag>);
#endif


SET_PROP(OnePDMPProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Set the grid type
SET_TYPE_PROP(OnePDMPProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);

// Set the problem property
SET_TYPE_PROP(OnePDMPProblem, Problem, Dumux::OnePDMPProblem<TypeTag> );

// Set the spatial parameters
SET_TYPE_PROP(OnePDMPProblem, SpatialParams, Dumux::OnePDMPSpatialParams<TypeTag> );


SET_BOOL_PROP(OnePDMPProblem, EnableGlobalFVGeometryCache, true);

SET_BOOL_PROP(OnePDMPProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(OnePDMPProblem, EnableGlobalVolumeVariablesCache, true);

// Enable gravity
SET_BOOL_PROP(OnePDMPProblem, ProblemEnableGravity, false);

SET_TYPE_PROP(OnePDMPProblem, LinearSolver, SuperLUBackend<TypeTag> );
}

template <class TypeTag>
class OnePDMPProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    // Grid and world dimension
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
     // indices of the primary variables
     conti0EqIdx = Indices::conti0EqIdx,
     pressureIdx = Indices::pressureIdx,
     //facePressureIdx = Indices::facePressureIdx
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;

public:
    OnePDMPProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);

        testCase_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                                int,
                                                Problem,
                                                TestCase);

        pi_ = 4.0*atan(1.0);
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
    std::string name() const
    {
        return name_;
    }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C


    PrimaryVariables source(const Element &element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolume &scv) const
    {
        PrimaryVariables values(0.0);

        if(testCase_ == 2)
        {
            const GlobalPosition& globalPos = element.geometry().center();
            Scalar p = elemVolVars[scv].pressure(0);

            Scalar x = globalPos[0];
            Scalar y = globalPos[1];

            if(std::abs(x-7.0/22.0) < 1.0e-8 && std::abs(y-0.5) < 1.0e-8)
                values[conti0EqIdx] = 1e12*(p - 0.0);
            else if(std::abs(x-15.0/22.0) < 1.0e-8 && std::abs(y-0.5) < 1.0e-8)
                values[conti0EqIdx] = 1e12*(p - 1.0);

            values[conti0EqIdx] /= scv.volume();
        }

        return values;

    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        if(testCase_ == 2)
            values.setAllNeumann();
        else
            values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        Scalar eps = 1.0e-8;

        if(x < eps || x > 1.0-eps || y < eps || y > 1.0 -eps )
        {
            PrimaryVariables values(0.0);
            return values;
        }
        else
        {
            PrimaryVariables values(1.0e5);
            return values;
        }
    }


    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    PrimaryVariables neumann(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace &scvf) const
    {
        PrimaryVariables values(0.0);
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables priVars(0);
        return priVars;
    }

    bool shouldWriteOutput() const
    {
        return
            this->timeManager().willBeFinished();
    }

    // \}

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->bBoxMax()[1] - eps_;
    }

    std::string name_;
    Scalar pi_;
    unsigned int testCase_;
    static constexpr Scalar eps_ = 3e-6;
};
} //end namespace

#endif
