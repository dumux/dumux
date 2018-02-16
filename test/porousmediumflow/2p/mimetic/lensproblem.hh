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

#ifndef DUMUX_MIMETIC_LENSPROBLEM_HH
#define DUMUX_MIMETIC_LENSPROBLEM_HH

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>
#include <dumux/porousmediumflow/2p/mimetic/model.hh>
//#include <dumux/porousmediumflow/2p/implicit/model.hh>
//#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

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
NEW_TYPE_TAG(LensProblem, INHERITS_FROM(TwoPMimetic, LensSpatialParams));
//NEW_TYPE_TAG(LensProblem, INHERITS_FROM(CCTpfaModel, TwoP, LensSpatialParams));

SET_TYPE_PROP(LensProblem, Grid, Dune::YaspGrid<2>);

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
    typedef FluidSystems::LiquidPhase<Scalar, DNAPL<Scalar> > type;
};

SET_TYPE_PROP(LensProblem, SpatialParams, LensSpatialParams<TypeTag> );

NEW_PROP_TAG(BaseProblem);
SET_TYPE_PROP(LensProblem, BaseProblem, ImplicitPorousMediaProblem<TypeTag>);

SET_TYPE_PROP(LensProblem, LinearSolver, SuperLUBackend<TypeTag> );

// Enable gravity
SET_BOOL_PROP(LensProblem, ProblemEnableGravity, true);

SET_BOOL_PROP(LensProblem, EnableGlobalFVGeometryCache, true);

SET_BOOL_PROP(LensProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(LensProblem, EnableGlobalVolumeVariablesCache, true);
}

/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeablility which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability. Note that
 * this problem is discretized using only two dimensions, so from the
 * point of view of the two-phase model, the depth of the domain
 * implicitly is 1 m everywhere.
 *
 * On the top and the bottom of the domain neumann boundary conditions
 * are used, while dirichlet conditions apply on the left and right
 * boundaries.
 *
 * DNAPL is injected at the top boundary from 3m to 4m at a rate of
 * 0.04 kg/(s m^2), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * The dirichlet boundaries on the left boundary is the hydrostatic
 * pressure scaled by a factor of 1.125, while on the right side it is
 * just the hydrostatic pressure. The DNAPL saturation on both sides
 * is zero.
 *
 * This problem uses the \ref TwoPModel.
 *
 * This problem should typically be simulated until \f$t_{\text{end}}
 * \approx 20\,000\;s\f$ is reached. A good choice for the initial time step
 * size is \f$t_{\text{inital}} = 250\;s\f$.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2p -parameterFile test_box2p.input</tt> or
 * <tt>./test_cc2p -parameterFile test_cc2p.input</tt>
 */
template <class TypeTag >
class LensProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

    enum {

        // primary variable indices
        pwIdx = Indices::pwIdx,
        snIdx = Indices::snIdx,
        facePressureWIdx = Indices::facePressureWIdx,
        facePressureNIdx = Indices::facePressureNIdx,

        // equation indices
        contiNEqIdx = Indices::contiNEqIdx,
        faceFluxBalanceNIdx = Indices::faceFluxBalanceNIdx,

        // phase indices
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,


        // world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    enum { adaptiveGrid = GET_PROP_VALUE(TypeTag, AdaptiveGrid) };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    LensProblem(TimeManager &timeManager,
                const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        temperature_ = 273.15 + 20; // -> 20°C

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);
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

    /*!
     * \brief Returns the source term
     *
     * \param values Stores the source values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param globalPos The global position
     */
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        return PrimaryVariables(0.0);
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
        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)) {
            values.setAllDirichlet();
        }
        else{
            values.setAllNeumann();
        }

        //values.setAllDirichlet();
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
        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;
        fluidState.setTemperature(temperature_);
        fluidState.setPressure(FluidSystem::wPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(FluidSystem::nPhaseIdx, /*pressure=*/1e5);

        Scalar densityW = FluidSystem::density(fluidState, FluidSystem::wPhaseIdx);

        Scalar height = this->bBoxMax()[1] - this->bBoxMin()[1];
        Scalar depth = this->bBoxMax()[1] - globalPos[1];
        Scalar alpha = 1 + 1.5/height;
        Scalar width = this->bBoxMax()[0] - this->bBoxMin()[0];
        Scalar factor = (width*alpha + (1.0 - alpha)*globalPos[0])/width;
        //factor = 1.0;

        // hydrostatic pressure scaled by alpha
        values[pwIdx] = 1e5 - factor*densityW*this->gravity()[1]*depth;
        values[snIdx] = 0.0;

        Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParamsAtPos(globalPos), 1.0 - values[snIdx]);
        values[facePressureWIdx] = values[pwIdx];
        values[facePressureNIdx] = pc + values[facePressureWIdx];

//        if (onUpperBoundary_(globalPos))
//        {
//            values[pwIdx] = 1e5 - densityW*this->gravity()[1]*depth;
//            values[snIdx] = 0.9;
//
//            Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParamsAtPos(globalPos), 1.0 - values[snIdx]);
//            values[facePressureWIdx] = values[pwIdx];
//            values[facePressureNIdx] = pc + values[facePressureWIdx];
//        }
//        else
//        {
//            values[pwIdx] = 1e5 - densityW*this->gravity()[1]*depth;
//            values[snIdx] = 0.1;
//
//            Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParamsAtPos(globalPos), 1.0 - values[snIdx]);
//            values[facePressureWIdx] = values[pwIdx];
//            values[facePressureNIdx] = pc + values[facePressureWIdx];
//        }

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
    PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        if (onInlet_(globalPos)) {
            values[contiNEqIdx] = -0.04; // kg / (m * s)
            values[faceFluxBalanceNIdx] = -0.04;
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
        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;
        fluidState.setTemperature(temperature_);
        fluidState.setPressure(FluidSystem::wPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(FluidSystem::nPhaseIdx, /*pressure=*/1e5);

        Scalar densityW = FluidSystem::density(fluidState, FluidSystem::wPhaseIdx);

        Scalar depth = this->bBoxMax()[1] - globalPos[1];

        // hydrostatic pressure
        values[pwIdx] = 1e5 - densityW*this->gravity()[1]*depth;
        values[snIdx] = 0.0;

        Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParamsAtPos(globalPos), 1.0-values[snIdx]);
        values[facePressureWIdx] = 0.0;
        values[facePressureNIdx] = 0.0;

        return values;
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }
    // \}

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& outputModule) const
    {
        auto& potW = outputModule.createScalarField("potW", 0);
        auto& potN = outputModule.createScalarField("potN", 0);

        typename GET_PROP_TYPE(TypeTag, FluidState) fluidState;
        fluidState.setTemperature(temperature_);
        fluidState.setPressure(FluidSystem::wPhaseIdx, /*pressure=*/1e5);
        fluidState.setPressure(FluidSystem::nPhaseIdx, /*pressure=*/1e5);

        Scalar densityW = FluidSystem::density(fluidState, FluidSystem::wPhaseIdx);
        Scalar densityN = FluidSystem::density(fluidState, FluidSystem::nPhaseIdx);

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->model().globalFvGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto ccDofIdx = scv.dofIndex();
                auto ccDofPosition = scv.dofPosition();

                auto elemVolVars = localView(this->model().curGlobalVolVars());
                elemVolVars.bind(element, fvGeometry, this->model().curSol());

                auto center = scv.center();
                Scalar depth = this->bBoxMax()[1] - center[1];
                potW[ccDofIdx] = elemVolVars[scv].pressure(wPhaseIdx) + densityW*this->gravity()[1]*depth;
                potN[ccDofIdx] = elemVolVars[scv].pressure(nPhaseIdx) + densityN*this->gravity()[1]*depth ;
            }
        }
    }

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

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = this->bBoxMax()[0] - this->bBoxMin()[0];
        Scalar lambda = (this->bBoxMax()[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda && lambda < 2.0/3.0;
    }

    Scalar temperature_;
    static constexpr Scalar eps_ = 3e-6;
    std::string name_;
};
} //end namespace

#endif
