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
 * \ingroup TwoPNCTest
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 */
#ifndef DUMUX_TWOPNC_DIFFUSION_PROBLEM_HH
#define DUMUX_TWOPNC_DIFFUSION_PROBLEM_HH

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/2pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/h2on2.hh>

#include "2pncdiffusionspatialparams.hh"
#include <dumux/discretization/maxwellstefanslaw.hh>

#ifndef DIFFUSIONTYPE // default to Fick's law if not set through CMake
#define DIFFUSIONTYPE FicksLaw<TypeTag>
#endif

namespace Dumux
{

template <class TypeTag>
class TwoPNCDiffusionProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoPNCDiffusionTypeTag, INHERITS_FROM(TwoPNC, TwoPNCDiffusionSpatialParams));
NEW_TYPE_TAG(TwoPNCDiffusionCCTypeTag, INHERITS_FROM(CCTpfaModel, TwoPNCDiffusionTypeTag));

// Set the grid type
SET_TYPE_PROP(TwoPNCDiffusionTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(TwoPNCDiffusionTypeTag, Problem, TwoPNCDiffusionProblem<TypeTag>);

// // Set fluid configuration
SET_TYPE_PROP(TwoPNCDiffusionTypeTag,
              FluidSystem,
              FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false /*useComplexRelations*/>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(TwoPNCDiffusionTypeTag, UseMoles, true);

//! Here we set FicksLaw or TwoPNCDiffusionsLaw
SET_TYPE_PROP(TwoPNCDiffusionTypeTag, MolecularDiffusionType, DIFFUSIONTYPE);

//! Set the default formulation to pw-Sn: This can be over written in the problem.
SET_INT_PROP(TwoPNCDiffusionTypeTag, Formulation, TwoPNCFormulation::pwsn);
}


/*!
 * \ingroup TwoPNCTest
 * \brief DOC_ME
 */
template <class TypeTag>
class TwoPNCDiffusionProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    // copy some indices for convenience
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,

        wPhaseOnly = Indices::wPhaseOnly,

        contiH2OEqIdx = Indices::contiWEqIdx,
        contiN2EqIdx = Indices::contiNEqIdx
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    TwoPNCDiffusionProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        // initialize the tables of the fluid system
        FluidSystem::init();

        name_ = getParam<std::string>("Problem.Name");

        //stateing in the console whether mole or mass fractions are used
        if(useMoles)
        {
            std::cout<<"problem uses mole-fractions"<<std::endl;
        }
        else
        {
            std::cout<<"problem uses mass-fractions"<<std::endl;
        }
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
    { return name_; }

    /*!
     * \brief Returns the temperature [K]
     */
    Scalar temperature() const
    { return 293.15; }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param values Stores the value of the boundary type
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllDirichlet();
        return bcTypes;
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
         PrimaryVariables priVars;
         priVars.setState(wPhaseOnly);
         priVars[pressureIdx] = 1e5;
         priVars[switchIdx] = 1e-5 ;

         if (globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_)
             priVars[switchIdx] = 1e-3;

        return priVars;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^2*s) or kg/(m^2*s))
     */
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector values(0.0);
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
        return initial_(globalPos);
    }

    // \}

private:
    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * The internal method for the initial condition
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars.setState(wPhaseOnly);

        //mole-fraction formulation
        priVars[switchIdx] = 1e-5;
        priVars[pressureIdx] = 1e5;
        return priVars;
    }

    Scalar temperature_;
    static constexpr Scalar eps_ = 1e-6;

   std::string name_ ;


};

} //end namespace Dumux

#endif
