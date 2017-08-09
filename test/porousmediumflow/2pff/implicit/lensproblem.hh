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
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_LENSPROBLEM_FRACTIONAL_FLOW_HH
#define DUMUX_LENSPROBLEM_FRACTIONAL_FLOW_HH

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/air.hh>
#include <dumux/porousmediumflow/2pff/implicit/propertydefaults.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/implicit/cellcentered/tpfa/properties.hh>

#include "lensspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class LensProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
NEW_TYPE_TAG(LensProblem, INHERITS_FROM(TwoPFractionalFlow, LensSpatialParams));
NEW_TYPE_TAG(LensCCProblem, INHERITS_FROM(CCTpfaModel, LensProblem));

SET_TYPE_PROP(LensCCProblem, Grid, Dune::FoamGrid<1, 3>);

// Set the problem property
SET_TYPE_PROP(LensProblem, Problem, LensProblem<TypeTag>);

// Set the wetting phase
SET_PROP(LensProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(LensProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::GasPhase<Scalar, Air<Scalar> > type;
};

// Linear solver settings
SET_TYPE_PROP(LensCCProblem, LinearSolver, ILU0BiCGSTABBackend<TypeTag> );

NEW_PROP_TAG(BaseProblem);
SET_TYPE_PROP(LensCCProblem, BaseProblem, ImplicitPorousMediaProblem<TypeTag>);
}

/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_cc2p -parameterFile test_cc2p.input</tt>
 */
template <class TypeTag >
class LensProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    using ParentType = typename GET_PROP_TYPE(TypeTag, BaseProblem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using WettingPhase = typename GET_PROP_TYPE(TypeTag, WettingPhase);
    using NonwettingPhase = typename GET_PROP_TYPE(TypeTag, NonwettingPhase);

    enum {
        // equation indices
        contiNEqIdx = Indices::contiNEqIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,


        // world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    enum { adaptiveGrid = GET_PROP_VALUE(TypeTag, AdaptiveGrid) };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Sources = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    LensProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        temperature_ = 273.15 + 20; // -> 20°C
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
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
    {
        return name_;
    }

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * This problem assumes a uniform temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return temperature_; }

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
        BoundaryTypes values;
        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)) {
            values.setAllDirichlet();
        }
        else {
            values.setAllNeumann();
        }
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
        PrimaryVariables values(0.0);
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
    NeumannFluxes neumannAtPos(const GlobalPosition &globalPos) const
    {
        NeumannFluxes values(0.0);
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
        PrimaryVariables values(0.0);
        return values;
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

    Scalar temperature_;
    static constexpr Scalar eps_ = 1e-7;
    std::string name_;
};
} //end namespace

#endif
