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
 * \ingroup TwoPTwoCTests
 * \brief Problem where air is injected in a unsaturated porous medium. This test compares a mpnc problem with a 2p2c problem
 */
#ifndef DUMUX_TWOPTWOC_MPNC_PROBLEM_HH
#define DUMUX_TWOPTWOC_MPNC_PROBLEM_HH

#include <dune/common/parametertreeparser.hh>

#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/2p2c/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidstates/compositional.hh>

#include "2p2c_comparison_spatialparams.hh"
#include "vtkoutputfields.hh"

namespace Dumux
{
/*!
 * \ingroup TwoPTwoCTests
 * \briefProblem where air is injected in a unsaturated porous medium. This test compares a   mpnc problem with a 2p2c problem
 */
template <class TypeTag>
class TwoPTwoCComparisonProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoPTwoCComparisonTypeTag, INHERITS_FROM(TwoPTwoC, TwoPTwoCComparisonSpatialParams));
NEW_TYPE_TAG(TwoPTwoCComparisonBoxTypeTag, INHERITS_FROM(BoxModel, TwoPTwoCComparisonTypeTag));
NEW_TYPE_TAG(TwoPTwoCComparisonCCTypeTag, INHERITS_FROM(CCTpfaModel, TwoPTwoCComparisonTypeTag));

// Set the grid type
SET_TYPE_PROP(TwoPTwoCComparisonTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(TwoPTwoCComparisonTypeTag,
              Problem,
              TwoPTwoCComparisonProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(TwoPTwoCComparisonTypeTag,
              FluidSystem,
              FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), /*useComplexRelations=*/false>);

// decide which type to use for floating values (double / quad)
SET_TYPE_PROP(TwoPTwoCComparisonTypeTag, Scalar, double);
SET_INT_PROP(TwoPTwoCComparisonTypeTag, Formulation, TwoPTwoCFormulation::pnsw);

SET_BOOL_PROP(TwoPTwoCComparisonTypeTag, UseMoles, true);

SET_PROP(TwoPTwoCComparisonTypeTag, VtkOutputFields)
{
private:
   using FluidSystem =  typename GET_PROP_TYPE(TypeTag, FluidSystem);
   using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

public:
    using type = TwoPTwoCMPNCVtkOutputFields<FluidSystem, Indices>;
};
}


/*!
 * \ingroup TwoPTwoCTests
 * \brief Problem where air is injected in a unsaturated porous medium. This test compares a mpnc problem with a 2p2c problem
 *
 */
template <class TypeTag>
class TwoPTwoCComparisonProblem
    : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

    // world dimension
    enum {dimWorld = GridView::dimensionworld};
    enum {numPhases = GET_PROP_TYPE(TypeTag, ModelTraits)::numPhases()};
    enum {numComponents = GET_PROP_TYPE(TypeTag, ModelTraits)::numComponents()};
    enum {nPhaseIdx = FluidSystem::nPhaseIdx};
    enum {wPhaseIdx = FluidSystem::wPhaseIdx};
    enum {wCompIdx = FluidSystem::wCompIdx};
    enum {nCompIdx = FluidSystem::nCompIdx};
    enum
    {
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,
        contiH2OEqIdx = Indices::contiWEqIdx,
        contiN2EqIdx = Indices::contiNEqIdx
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    static constexpr bool isBox = GET_PROP_TYPE(TypeTag, FVGridGeometry)::discMethod == DiscretizationMethod::box;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    TwoPTwoCComparisonProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry)
    {
        temperature_ = 273.15 + 25; // -> 25°C

        // initialize the tables of the fluid system
        Scalar Tmin = temperature_ - 1.0;
        Scalar Tmax = temperature_ + 1.0;
        int nT = 3;

        Scalar pmin = 1.0e5 * 0.75;
        Scalar pmax = 2.0e5 * 1.25;
        int np = 1000;

        FluidSystem::init(Tmin, Tmax, nT, pmin, pmax, np);
        name_ = getParam<std::string>("Problem.Name");
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
    { return name_; }

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * \param globalPos The global position
     */
    Scalar temperature() const
    { return temperature_; }

    /*!
     * \name Boundary conditions
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if (onOutlet_(globalPos))
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
       return initial_(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann
     *        boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local index of the sub-control volume
     * \param boundaryFaceIdx The index of the boundary face
     *
     * Negative values mean influx.
     */
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);
        const auto& globalPos = scvf.ipGlobal();
        Scalar injectedAirMass = -1e-3;
        Scalar injectedAirMolarMass = injectedAirMass/FluidSystem::molarMass(nCompIdx);
        if (onInlet_(globalPos))
        {
          values[Indices::conti0EqIdx+1] = injectedAirMolarMass;
        }
        return values;
    }

    // \}

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        return initial_(globalPos);
    }

    // \}

private:
    // the internal method for the initial condition
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx] = 1e5;
        values[switchIdx] = 0.8;
        values.setState(Indices::bothPhases);

        return values;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        return x >= 60 - eps_ && y <= 10 + eps_;
    }

    bool onOutlet_(const GlobalPosition &globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        return x < eps_ && y <= 10 + eps_;
    }

    Scalar temperature_;
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
};
} //end namespace

#endif
