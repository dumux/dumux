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
 * \ingroup TwoPTests
 * \brief Non-isothermal gas injection problem where a gas (e.g. air) is injected into a fully
 *        water saturated medium. During buoyancy driven upward migration the gas
 *        passes a high temperature area.
 */

#ifndef DUMUX_INJECTION_PROBLEM_2PNI_HH
#define DUMUX_INJECTION_PROBLEM_2PNI_HH

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/components/n2.hh>

// use the spatial parameters as the injection problem of the 2p2c test program
#include <test/porousmediumflow/2p2c/implicit/injectionspatialparams.hh>

#ifndef GRIDTYPE // default to yasp grid if not provided by CMake
#define GRIDTYPE Dune::YaspGrid<2>
#endif

namespace Dumux {

//! Forward declaration of the problem class
template <class TypeTag> class InjectionProblem2PNI;

namespace Properties {
NEW_TYPE_TAG(InjectionProblem2PNI, INHERITS_FROM(TwoPNI, InjectionSpatialParams));
NEW_TYPE_TAG(InjectionBoxProblem2PNI, INHERITS_FROM(BoxModel, InjectionProblem2PNI));
NEW_TYPE_TAG(InjectionCCProblem2PNI, INHERITS_FROM(CCTpfaModel, InjectionProblem2PNI));

// Obtain grid type from COMPILE_DEFINITIONS
SET_TYPE_PROP(InjectionProblem2PNI, Grid, GRIDTYPE);

// Set the problem property
SET_TYPE_PROP(InjectionProblem2PNI, Problem, InjectionProblem2PNI<TypeTag>);

// Use the same fluid system as the 2p2c injection problem
SET_TYPE_PROP(InjectionProblem2PNI, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false>);
} // namespace Properties

/*!
 * \ingroup TwoPTests
 * \brief Non-isothermal gas injection problem where a gas (e.g. air) is injected into a fully
 *        water saturated medium. During buoyancy driven upward migration the gas
 *        passes a high temperature area.
 *
 * The domain is sized 60 m times 40 m. The rectangular area with the increased temperature (380 K)
 * starts at (20 m, 5 m) and ends at (30 m, 35 m)
 *
 * For the mass conservation equation neumann boundary conditions are used on
 * the top, on the bottom and on the right of the domain, while dirichlet conditions
 * apply on the left boundary.
 * For the energy conservation equation dirichlet boundary conditions are applied
 * on all boundaries.
 *
 * Gas is injected at the right boundary from 5 m to 15 m at a rate of
 * 0.001 kg/(s m), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * At the dirichlet boundaries a hydrostatic pressure, a gas saturation of zero and
 * a geothermal temperature gradient of 0.03 K/m are applied.
 *
 * This problem uses the \ref TwoPModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2pni -parameterFile test_box2pni.input</tt> or
 * <tt>./test_cc2pni -parameterFile test_cc2pni.input</tt>
 */
template<class TypeTag>
class InjectionProblem2PNI : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum
    {
        //! primary variable indices
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        temperatureIdx = Indices::temperatureIdx,

        //! equation indices
        contiNEqIdx = Indices::contiNEqIdx,
        energyEqIdx = Indices::energyEqIdx,

        //! phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        // world dimension
        dimWorld = GridView::dimensionworld
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Sources = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    InjectionProblem2PNI(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        maxDepth_ = 2700.0; // [m]

        // initialize the tables of the fluid system
        FluidSystem::init(/*tempMin=*/273.15,
                          /*tempMax=*/423.15,
                          /*numTemp=*/50,
                          /*pMin=*/0.0,
                          /*pMax=*/30e6,
                          /*numP=*/300);

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
     * \brief Returns the source term
     *
     * \param values Stores the source values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param globalPos The global position
     */
    Sources sourceAtPos(const GlobalPosition &globalPos) const
    {
        Sources values(0.0);
        values = 0;
        return values;
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
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] < eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        //! Use Dirichlet BCs everywhere for the temperature
        // values.setDirichlet(temperatureIdx);
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
        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 + (maxDepth_ - globalPos[dimWorld-1])*densityW*9.81;
        values[saturationIdx] = 0.0;
        values[temperatureIdx] = 283.0 + (maxDepth_ - globalPos[dimWorld-1])*0.03;
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param elemVolVars The element volume variables
     * \param scvf The sub-control volume face on which the BC is evaluated
     *
     * The \a values store the mass flux of each phase normal to the boundary.
     * Negative values indicate an inflow.
     */
    NeumannFluxes neumann(const Element &element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);
        const auto globalPos = scvf.ipGlobal();

        if (globalPos[1] < 13.75 + eps_ && globalPos[1] > 6.875 - eps_)
        {
            // inject air. negative values mean injection
            values[contiNEqIdx] = -1e-3; // kg/(s*m^2)

            // compute enthalpy flux associated with this injection [(J/(kg*s)]
            using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
            FluidState fs;

            const auto ic = initialAtPos(scvf.ipGlobal());
            fs.setPressure(wPhaseIdx, ic[pressureIdx]);
            fs.setPressure(nPhaseIdx, ic[pressureIdx]); // assume pressure equality here
            fs.setTemperature(ic[temperatureIdx]);

            // energy flux is mass flux times specific enthalpy
            values[energyEqIdx] = values[contiNEqIdx]*FluidSystem::enthalpy(fs, nPhaseIdx);
        }

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
        PrimaryVariables values;
        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 + (maxDepth_ - globalPos[1])*densityW*9.81;
        values[saturationIdx] = 0.0;
        values[temperatureIdx] = 283.0 + (maxDepth_ - globalPos[1])*0.03;

        if (globalPos[0] > 21.25 - eps_ && globalPos[0] < 28.75 + eps_ && globalPos[1] > 6.25 - eps_ && globalPos[1] < 33.75 + eps_)
            values[temperatureIdx] = 380;

        return values;
    }
    // \}

private:
    Scalar maxDepth_;
    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;
};

} //end namespace Dumux

#endif
