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
 * \ingroup TracerTests
 * \brief The properties for the incompressible test
 */
#ifndef DUMUX_INCOMPRESSIBLE_TWOP_TEST_PROBLEM_HH
#define DUMUX_INCOMPRESSIBLE_TWOP_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/material/components/trichloroethene.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include "2pncimmiscible.hh"

#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "2ptestspatialparams.hh"
// #include "2ptestspatialparams_randomfield.hh"
// #include "2pnctestlocalresidual.hh"

#ifndef ENABLEINTERFACESOLVER
#define ENABLEINTERFACESOLVER 0
#endif

namespace Dumux
{
/*!
 * \ingroup TracerTests
 * \brief The properties for the incompressible 2p test
 */
// forward declarations
template<class TypeTag>
class TwoPNCTestProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoPNCTestProblem, INHERITS_FROM(TwoPNC));
NEW_TYPE_TAG(TwoPNCTestProblemTpfa, INHERITS_FROM(CCTpfaModel, TwoPNCTestProblem, SpatialParams));
NEW_TYPE_TAG(TwoPNCTestProblemBox, INHERITS_FROM(BoxModel, TwoPNCTestProblem, SpatialParams));

// Set the grid type
SET_TYPE_PROP(TwoPNCTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem type
SET_TYPE_PROP(TwoPNCTestProblem, Problem, TwoPNCTestProblem<TypeTag>);

// the fluid system
SET_PROP(TwoPNCTestProblem, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Trichloroethene<Scalar> >;
    using type = FluidSystems::TwoPNCImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

// Enable caching
SET_BOOL_PROP(TwoPNCTestProblem, EnableGridVolumeVariablesCache, false);
SET_BOOL_PROP(TwoPNCTestProblem, EnableGridFluxVariablesCache, false);
SET_BOOL_PROP(TwoPNCTestProblem, EnableFVGridGeometryCache, false);

// Maybe enable the box-interface solver
SET_BOOL_PROP(TwoPNCTestProblem, EnableBoxInterfaceSolver, ENABLEINTERFACESOLVER);

} // end namespace Properties

/*!
 * \ingroup TracerTests
 * \brief The incompressible 2p test problem.
 */

template<class TypeTag>
class TwoPNCTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    enum {
        pressureH2OIdx = Indices::pressureIdx,
        saturationDNAPLIdx = Indices::switchIdx,
        contiDNAPLEqIdx = Indices::conti0EqIdx + FluidSystem::comp1Idx,
        waterPhaseIdx = FluidSystem::phase0Idx,
        dnaplPhaseIdx = FluidSystem::phase1Idx,

        //the tracercomponent indices
        contiTracer1EqIdx = Indices::conti0EqIdx + FluidSystem::comp2Idx,
        tracer1Idx = Indices::conti0EqIdx + FluidSystem::comp2Idx,
    };

    // static constexpr int dimWorld = GridView::dimensionworld;

public:
    TwoPNCTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        Dune::FMatrixPrecision<>::set_singular_limit(1e-35);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param values Stores the value of the boundary type
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
            values.setAllDirichlet();
        else
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
        values.setState(Indices::firstPhaseOnly);

        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;
        fluidState.setTemperature(temperature());
        fluidState.setPressure(waterPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(dnaplPhaseIdx, /*pressure=*/1e5);

        Scalar densityW = FluidSystem::density(fluidState, waterPhaseIdx);

        // Scalar height = this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1];
        Scalar depth = this->fvGridGeometry().bBoxMax()[1] - globalPos[1];
        // Scalar alpha = 1 + 1.5/height;
        // Scalar width = this->fvGridGeometry().bBoxMax()[0] - this->fvGridGeometry().bBoxMin()[0];
        Scalar factor = 1;

        // hydrostatic pressure scaled by alpha
        values[pressureH2OIdx] = 1e5 - factor*densityW*this->gravity()[1]*depth;
        values[saturationDNAPLIdx] = 0.0;

        //the tracer component's Dirichlet BC
//  VERSION 1
//      if (onUpperBoundary_(globalPos))

//  VERSION 2
        if (onUpperBoundary_(globalPos) ||  onStripe1_(globalPos) || onStripe2_(globalPos) || onStripe3_(globalPos))

//  VERSION 3
//         if (onUpperTrajectoryPoint_(globalPos) || onLeftTrajectroyPoint_(globalPos) || onRightTrajectoryPoint_(globalPos))
        {
            if (useMoles)
                values[tracer1Idx] = 1e-9;
            else
                values[tracer1Idx] = 1e-9*FluidSystem::molarMass(tracer1Idx)/FluidSystem::molarMass(0);
        }

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
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);
        if (onInlet_(globalPos))
        {
            values[contiDNAPLEqIdx] = -0.04; // kg / (m * s)
            values[Indices::conti0EqIdx] = -0.04;
        }

        // in the test with the oil wet lens, use higher injection rate
        if (this->spatialParams().lensIsOilWet())
            values[contiDNAPLEqIdx] *= 10;

        //no tracer is injected in the tracerproblem, so we do not need to set a Neumann BC for it different from 0!

        return values;
    }

    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;

        values.setState(Indices::firstPhaseOnly);

        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;
        fluidState.setTemperature(temperature());
        fluidState.setPressure(waterPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(dnaplPhaseIdx, /*pressure=*/1e5);

        Scalar densityW = FluidSystem::density(fluidState, waterPhaseIdx);

        Scalar depth = this->fvGridGeometry().bBoxMax()[1] - globalPos[1];

        // hydrostatic pressure
        values[pressureH2OIdx] = 1e5 - densityW*this->gravity()[1]*depth;
        values[saturationDNAPLIdx] = 0;

        //the tracer component's initial values
//  VERSION 1
//      if (onUpperBoundary_(globalPos))

//  VERSION 2
        if (onUpperBoundary_(globalPos) ||  onStripe1_(globalPos) || onStripe2_(globalPos) || onStripe3_(globalPos))

//  VERSION 3
//         if (onUpperTrajectoryPoint_(globalPos) || onLeftTrajectroyPoint_(globalPos) || onRightTrajectoryPoint_(globalPos))
        {
            if (useMoles)
                values[tracer1Idx] = 1e-9;
            else
                values[tracer1Idx] = 1e-9*FluidSystem::molarMass(tracer1Idx)/FluidSystem::molarMass(0);
        }

        return values;
    }


    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    {
        return 283.15; // 10°C
    }

    void updateVtkFields(const SolutionVector& curSol)
    { }

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->fvGridGeometry().bBoxMax()[0] - this->fvGridGeometry().bBoxMin()[0];
        Scalar lambda = (this->fvGridGeometry().bBoxMax()[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda && lambda < 2.0/3.0;
    }

    static constexpr Scalar eps_ = 1e-6;

        Scalar yMax_ = this->fvGridGeometry().bBoxMax()[1];
    Scalar xMax_ = this->fvGridGeometry().bBoxMax()[0];

// TODO allgemeiner Abruf der Anzahl von Zellen in Y-Richtung/ X-Richtung
    Scalar yNumCells_ = 32;
    Scalar xNumCells_ = 48;
    Scalar cellHeight_ = yMax_/yNumCells_;
    Scalar cellWidth_ = xMax_/xNumCells_;
    Scalar width_ = xMax_ - this->fvGridGeometry().bBoxMin()[0];

//  VERSION 1
     bool onUpperBoundary_(const GlobalPosition &globalPos) const
     {
         return globalPos[1] > yMax_ - 0.1 - eps_;
     }


//  VERSION 2
     bool onStripe1_(const GlobalPosition &globalPos) const
     {
//          std::cout << "Stripe1 is: " << (( (yMax_ /4.0 - cellHeight_*0.5) <= globalPos[1] ) &&
//                   ( (yMax_/4.0 + cellHeight_*0.5) > globalPos[1] ));
//          std::cout << " at globalPos:" << globalPos[1] << "\n";
//          std::cout << "yMax is: " << yMax_ << "\n" ;

        return  (
            ( (yMax_ /4.0 - cellHeight_*0.5) <= globalPos[1] ) &&
            ( (yMax_/4.0 + cellHeight_*0.5) > globalPos[1] )
        );
    }



    bool onStripe2_(const GlobalPosition &globalPos) const
    {
        return  (
            ( (2.0 * yMax_ /4.0 - cellHeight_*0.5) <= globalPos[1] ) &&
            ( (2.0 * yMax_/4.0 + cellHeight_*0.5) > globalPos[1] )
        );
    }

    bool onStripe3_(const GlobalPosition &globalPos) const
    {
        return  (
            ( (3.0 * yMax_ /4.0 - cellHeight_*0.5) <= globalPos[1] ) &&
            ( (3.0 * yMax_/4.0 + cellHeight_*0.5) > globalPos[1] )
        );
    }

//  VERSION 3

     bool onUpperTrajectoryPoint_(const GlobalPosition &globalPos) const
     {
//          std::cout << "onUpperTrajectoryPoint is: "
//          << (    globalPos[1] > (yMax_ - 2.0*cellHeight_) &&
//                  globalPos[0] > (width_/2.0 - cellWidth_) &&
//                  globalPos[0] < (width_/2.0 + cellWidth_) );
//          std::cout << " at x- globalPos:" << globalPos[0] << "\n";
//          std::cout << " at y- globalPos:" << globalPos[1] << "\n";

         return (
             globalPos[1] > (yMax_ - 2.0*cellHeight_) &&
             globalPos[0] > (width_/2.0 - cellWidth_) &&
             globalPos[0] < (width_/2.0 + cellWidth_)
         );
     }

     bool onLeftTrajectroyPoint_(const GlobalPosition &globalPos) const
     {
         return (
             ( globalPos[1] >= (3.0 * yMax_ /4.0 - cellHeight_) ) &&
             ( globalPos[1] < (3.0 * yMax_/4.0 + cellHeight_) ) &&
             ( globalPos[0] > (xMax_ /6.0 - cellHeight_) ) &&
             ( globalPos[0] < (xMax_ /6.0 + cellHeight_) )
         );
     }

     bool onRightTrajectoryPoint_(const GlobalPosition &globalPos) const
     {
         return (
             ( globalPos[1] >= (3.0 * yMax_ /4.0 - cellHeight_) ) &&
             ( globalPos[1] < (3.0 * yMax_/4.0 + cellHeight_) ) &&
             ( globalPos[0] > (5.0 * xMax_/6.0 - cellHeight_) ) &&
             ( globalPos[0] < (5.0 * xMax_/6.0 + cellHeight_) )
         );
     }
};

} // end namespace Dumux

#endif
