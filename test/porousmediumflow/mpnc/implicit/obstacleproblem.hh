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
 * \ingroup MPNCTests
 * \brief Problem where liquid water is injected -- by means of a
 *        Dirichlet condition on the lower right of the domain -- which has to go
 *        around an obstacle with \f$10^3\f$ lower permeability.
 */
#ifndef DUMUX_OBSTACLEPROBLEM_HH
#define DUMUX_OBSTACLEPROBLEM_HH

#include <dune/common/parametertreeparser.hh>

#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/mpnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>

#include "obstaclespatialparams.hh"

namespace Dumux
{
/*!
 * \ingroup MPNCTests
 * \brief Problem where liquid water is injected -- by means of a
 *        Dirichlet condition on the lower right of the domain -- which has to go
 *        around an obstacle with \f$10^3\f$ lower permeability.
 */
template <class TypeTag>
class ObstacleProblem;

namespace Properties
{
NEW_TYPE_TAG(ObstacleTypeTag, INHERITS_FROM(MPNC, ObstacleSpatialParams));
NEW_TYPE_TAG(ObstacleBoxTypeTag, INHERITS_FROM(BoxModel, ObstacleTypeTag));
NEW_TYPE_TAG(ObstacleCCTypeTag, INHERITS_FROM(CCTpfaModel, ObstacleTypeTag));

// Set the grid type
SET_TYPE_PROP(ObstacleTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ObstacleTypeTag,
              Problem,
              ObstacleProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(ObstacleTypeTag,
              FluidSystem,
              FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), /*useComplexRelations=*/false>);

// decide which type to use for floating values (double / quad)
SET_TYPE_PROP(ObstacleTypeTag, Scalar, double);

}


/*!
 * \ingroup MPNCModel
 * \ingroup ImplicitTestProblems
 * \brief Problem where liquid water is injected -- by means of a
 *        Dirichlet condition on the lower right of the domain -- which has to go
 *        around an obstacle with \f$10^3\f$ lower permeability.
 *
 * The domain is sized 60m times 40m and consists of two media, a
 * moderately permeable soil (\f$ K_0=10e-12 m^2\f$) and an obstacle
 * at \f$[10; 20]m \times [0; 35]m \f$ with a lower permeablility of
 * \f$ K_1=K_0/1000\f$.
 *
 * Initially the whole domain is filled with nitrogen, the temperature
 * is \f$25\symbol{23}C\f$ in the whole domain. The gas pressure in the
 * domain is 1 bar, except on the inlet (lower 10 meters of right hand
 * boundary) where the pressure is 2 bar.
 *
 * The boundary is no-flow except on the lower 10 meters of the left
 * and the right boundary which are Dirichlet conditions with the same
 * values as the initial condition.
 *
 * This problem uses the \ref MPNCModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_boxmpnc -parameterFile test_boxmpnc.input</tt> or
 * <tt>./test_ccmpnc -parameterFile test_ccmpnc.input</tt>
 */
template <class TypeTag>
class ObstacleProblem
    : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using ParameterCache = typename FluidSystem::ParameterCache;

    enum {dimWorld = GridView::dimensionworld};
    enum {numPhases = GET_PROP_VALUE(TypeTag, NumPhases)};
    enum {numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};
    enum {nPhaseIdx = FluidSystem::nPhaseIdx};
    enum {wPhaseIdx = FluidSystem::wPhaseIdx};
    enum {wCompIdx = FluidSystem::wCompIdx};
    enum {nCompIdx = FluidSystem::nCompIdx};
    enum {fug0Idx = Indices::fug0Idx};
    enum {s0Idx = Indices::s0Idx};
    enum {p0Idx = Indices::p0Idx};

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using PhaseVector = Dune::FieldVector<Scalar, numPhases>;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    ObstacleProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
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
     * \brief Write a restart file?
     */
    bool shouldWriteRestartFile() const
    {
        return ParentType::shouldWriteRestartFile();
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
     * \param globalPos The global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if (onInlet_(globalPos) || onOutlet_(globalPos))
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
    NumEqVector neumann(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf) const
    {
        return NumEqVector(0.0);
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all balance equations within a given
     *        sub-control-volume.
     *
     * \param values Stores the solution for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     *
     * Positive values mean that mass is created, negative ones mean that it vanishes.
     */
    //! \copydoc Dumux::ImplicitProblem::source()
    NumEqVector source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
       return NumEqVector(0.0);
    }

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
        FluidState fs;

        int refPhaseIdx;
        int otherPhaseIdx;

        // set the fluid temperatures
        fs.setTemperature(this->temperatureAtPos(globalPos));

        if (onInlet_(globalPos))
        {
            // only liquid on inlet
            refPhaseIdx = wPhaseIdx;
            otherPhaseIdx = nPhaseIdx;

            // set liquid saturation
            fs.setSaturation(wPhaseIdx, 1.0);

            // set pressure of the liquid phase
            fs.setPressure(wPhaseIdx, 2e5);

            // set the liquid composition to pure water
            fs.setMoleFraction(wPhaseIdx, nCompIdx, 0.0);
            fs.setMoleFraction(wPhaseIdx, wCompIdx, 1.0);
        }
        else {
            // elsewhere, only gas
            refPhaseIdx = nPhaseIdx;
            otherPhaseIdx = wPhaseIdx;

            // set gas saturation
            fs.setSaturation(nPhaseIdx, 1.0);

            // set pressure of the gas phase
            fs.setPressure(nPhaseIdx, 1e5);

            // set the gas composition to 99% nitrogen and 1% steam
            fs.setMoleFraction(nPhaseIdx, nCompIdx, 0.99);
            fs.setMoleFraction(nPhaseIdx, wCompIdx, 0.01);
        }

        // set the other saturation
        fs.setSaturation(otherPhaseIdx, 1.0 - fs.saturation(refPhaseIdx));

        // calulate the capillary pressure
        const auto& matParams =
            this->spatialParams().materialLawParamsAtPos(globalPos);
        PhaseVector pc;
        MaterialLaw::capillaryPressures(pc, matParams, fs);
        fs.setPressure(otherPhaseIdx,
                       fs.pressure(refPhaseIdx)
                       + (pc[otherPhaseIdx] - pc[refPhaseIdx]));

        // make the fluid state consistent with local thermodynamic
        // equilibrium
        using ComputeFromReferencePhase = ComputeFromReferencePhase<Scalar, FluidSystem>;

        ParameterCache paramCache;
        ComputeFromReferencePhase::solve(fs,
                                         paramCache,
                                         refPhaseIdx,
                                         /*setViscosity=*/false,
                                         /*setEnthalpy=*/false);

        ///////////
        // assign the primary variables
        ///////////

        // all N component fugacities
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
            values[fug0Idx + compIdx] = fs.fugacity(refPhaseIdx, compIdx);

        // first M - 1 saturations
        for (int phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
            values[s0Idx + phaseIdx] = fs.saturation(phaseIdx);

        // first pressure
        values[p0Idx] = fs.pressure(/*phaseIdx=*/0);
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
