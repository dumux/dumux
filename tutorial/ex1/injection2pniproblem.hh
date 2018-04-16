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
 * \brief The two-phase porousmediumflow problem for exercise 1
 */

#ifndef DUMUX_EX1_INJECTION_PROBLEM_2PNI_HH
#define DUMUX_EX1_INJECTION_PROBLEM_2PNI_HH

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/h2on2.hh>

#include "injection2pspatialparams.hh"

namespace Dumux {

// forward declare problem
template <class TypeTag>
class InjectionProblem2PNI;

namespace Properties
{
 /*!
* TODO:dumux-course-task:
* inherit from the TwoPNI model instead of TwoP here
*/
NEW_TYPE_TAG(Injection2pNITypeTag, INHERITS_FROM(TwoP, InjectionSpatialParamsTypeTag));
NEW_TYPE_TAG(Injection2pNICCTypeTag, INHERITS_FROM(CCTpfaModel, Injection2pNITypeTag));

// Set the grid type
SET_TYPE_PROP(Injection2pNITypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(Injection2pNITypeTag, Problem, InjectionProblem2PNI<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(Injection2pNITypeTag, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false>);
} // end namespace Properties

/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief Gas injection problem where a gas (here  nitrogen) is injected into a fully
 *        water saturated medium. During buoyancy driven upward migration the gas
 *        passes a high temperature area.
 *
 * The domain is sized 60 m times 40 m.
 *
 * For the mass conservation equation neumann boundary conditions are used on
 * the top, on the bottom and on the right of the domain, while dirichlet conditions
 * apply on the left boundary.
 *
 * Gas is injected at the right boundary from 7 m to 15 m at a rate of
 * 0.001 kg/(s m), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * At the dirichlet boundaries a hydrostatic pressure and a gas saturation of zero a
 *
 * This problem uses the \ref TwoPModel model.
 */
template<class TypeTag>
class InjectionProblem2PNI : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    enum { dimWorld = GridView::dimensionworld };
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimension>;

public:
    InjectionProblem2PNI(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        // initialize the tables of the fluid system
        FluidSystem::init(/*tempMin=*/273.15,
                /*tempMax=*/423.15,
                /*numTemp=*/50,
                /*pMin=*/0.0,
                /*pMax=*/30e6,
                /*numP=*/300);

        // name of the problem and output file
        // getParam<TYPE>("GROUPNAME.PARAMNAME") reads and sets parameter PARAMNAME
        // of type TYPE given in the group GROUPNAME from the input file
        name_ = getParam<std::string>("Problem.Name");
        // depth of the aquifer, units: m
        aquiferDepth_ = getParam<Scalar>("Problem.AquiferDepth");
        // the duration of the injection, units: second
        injectionDuration_ = getParam<Scalar>("Problem.InjectionDuration");
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
    std::string name() const
    { return "injection-2pni"; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param bcTypes The boundary types for the conservation equations
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
         BoundaryTypes bcTypes;
        if (globalPos[0] < eps_)
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();

        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initialAtPos(globalPos);
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
        // initialize values to zero, i.e. no-flow Neumann boundary conditions
        PrimaryVariables values(0.0);

        // if we are inside the injection zone set inflow Neumann boundary conditions
        if (time_ < injectionDuration_
            && globalPos[1] < 15 + eps_ && globalPos[1] > 7 - eps_ && globalPos[0] > 0.9*this->fvGridGeometry().bBoxMax()[0])
        {
            // inject nitrogen. negative values mean injection
            // units kg/(s*m^2)
            values[Indices::contiNEqIdx] = -1e-4;
            values[Indices::contiWEqIdx] = 0.0;

         /*!
          * TODO:dumux-course-task:
          * dumux-course-task:
          * set Neumann noflow conditions for the energy equation everywhere
          * hint: use Indices::energyEqIdx) for that
          */
        }

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
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        // get the water density at atmospheric conditions
        const Scalar densityW = FluidSystem::H2O::liquidDensity(283.15, 1.0e5);

        // assume an intially hydrostatic liquid pressure profile
        // note: we subtract rho_w*g*h because g is defined negative
        const Scalar pw = 1.0e5 - densityW*this->gravity()[dimWorld-1]*(aquiferDepth_ - globalPos[dimWorld-1]);

        values[Indices::pressureIdx] = pw;
        values[Indices::saturationIdx] = 0.0;

        /*!
        *  TODO:dumux-course-task:
        * set a temperature gradient of 0.03 K per m beginning at 283 K here.
        * Hint: you can use aquiferDepth_ and the globalPos similar to the pressure gradient
        * use globalPos[0] and globalPos[1] to implement the high temperature lens with 380 K
        * Hint : use Indices::temperatureIdx to address the initial values for temperature
        */
        return values;
    }

    // \}

    //! set the time for the time dependent boundary conditions (called from main)
    void setTime(Scalar time)
    { time_ = time; }

private:
    static constexpr Scalar eps_ = 1e-6;
    std::string name_; //! Problem name
    Scalar aquiferDepth_; //! Depth of the aquifer in m
    Scalar injectionDuration_; //! Duration of the injection in seconds
    Scalar time_;
};

} //end namespace Dumux

#endif
