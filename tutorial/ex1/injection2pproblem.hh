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
 * \brief Non-isothermal gas injection problem where a gas (e.g. air) is injected into a fully
 *        water saturated medium. During buoyancy driven upward migration the gas
 *        passes a high temperature area.
 */

#ifndef DUMUX_INJECTION_PROBLEM_2P_HH
#define DUMUX_INJECTION_PROBLEM_2P_HH

#include <dumux/porousmediumflow/2p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/material/fluidsystems/h2on2.hh>

#include "injection2pspatialparams.hh"

namespace Dumux {

template <class TypeTag>
class InjectionProblem2P;

namespace Properties
{
NEW_TYPE_TAG(InjectionProblem2P, INHERITS_FROM(TwoP, InjectionSpatialParams));
NEW_TYPE_TAG(InjectionCCProblem2P, INHERITS_FROM(CCModel, InjectionProblem2P));

// Set the grid type
SET_TYPE_PROP(InjectionProblem2P, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(InjectionProblem2P, Problem, InjectionProblem2P<TypeTag>);

// Use the same fluid system as the 2p2c injection problem
SET_TYPE_PROP(InjectionProblem2P, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false>);

}

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
 *
 * To run the simulation execute the following line in shell:
 * <tt>./exercise1_2p -ParameterFile exercise1.input</tt>
 */
template<class TypeTag>
class InjectionProblem2P : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;


    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        contiNEqIdx = Indices::contiNEqIdx,

        // world dimension
        dimWorld = GridView::dimensionworld
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    InjectionProblem2P(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        // initialize the tables of the fluid system
        FluidSystem::init(/*tempMin=*/273.15,
                /*tempMax=*/423.15,
                /*numTemp=*/50,
                /*pMin=*/0.0,
                /*pMax=*/30e6,
                /*numP=*/300);

        // name of the problem and output file
        name_  = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        // depth of the aquifer, units: m
        aquiferDepth_  = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, AquiferDepth);
        // inflow rate of nitrogen water vapor mixture, units: kg/(s m^2)
        // the duration of the injection, units: second
        injectionDuration_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InjectionDuration);
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
    const std::string name() const
    { return "injection-2p"; }

    /*!
     * \brief Returns the temperature \f$ K \f$
     */
    Scalar temperature() const
    {
        return 273.15 + 30; // [K]
    }



    /*!
     * \brief Returns the source term
     *
     * \param values Stores the source values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param globalPos The global position
     */
    void sourceAtPos(PrimaryVariables &values,
                const GlobalPosition &globalPos) const
    {
        values = 0;
    }

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
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        if (globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initialAtPos(values, globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The globalPosition of the boundary interface
     */
     void neumannAtPos(PrimaryVariables &values,
                       const GlobalPosition &globalPos) const
    {
         // initialize values to zero, i.e. no-flow Neumann boundary conditions
         values = 0;
        values = 0;

        if (this->timeManager().time() + this->timeManager().timeStepSize() < injectionDuration_
            && globalPos[1] < 15 + eps_ && globalPos[1] > 7 - eps_ && globalPos[0] > 0.9*this->bBoxMax()[0])
        {
            // inject air. negative values mean injection
            values[Indices::contiNEqIdx] = -1e-4; // kg/(s*m^2)
            values[Indices::contiWEqIdx] = 0.0;
        }
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
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {

        // get the water density at atmospheric conditions
        const Scalar densityW = FluidSystem::H2O::liquidDensity(temperature(), 1.0e5);

        // assume an intially hydrostatic liquid pressure profile
        // note: we subtract rho_w*g*h because g is defined negative
        const Scalar pw = 1.0e5 - densityW*this->gravity()[dimWorld-1]*(aquiferDepth_ - globalPos[dimWorld-1]);

        values[Indices::pressureIdx] = pw;
        values[saturationIdx] = 0.0;
    }
    // \}

private:
    static constexpr Scalar eps_ = 1e-6;
    std::string name_; //! Problem name
    Scalar aquiferDepth_; //! Depth of the aquifer in m
    Scalar injectionDuration_; //! Duration of the injection in seconds
};
} //end namespace

#endif
