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

#ifndef DUMUX_MIMETIC_CPTWOP_HH
#define DUMUX_MIMETIC_CPTWOP_HH

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/dnapl.hh>
#include <dumux/porousmediumflow/2p/mimetic/model.hh>
#include <dumux/discretization/staggered/mimetic/mimeticcpgeometryhelper.hh>
//#include <dumux/porousmediumflow/2p/implicit/model.hh>
//#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>

#include <dumux/io/cpgridcreator.hh>
#include <dumux/common/intersectionmapper.hh>

#include "cp2pspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class CPTwoPProblem;

//////////
// Specify the properties for the CPTwoP problem
//////////
namespace Properties
{
NEW_TYPE_TAG(CPTwoPProblem, INHERITS_FROM(TwoPMimetic, CPTwoPSpatialParams));
//NEW_TYPE_TAG(CPTwoPProblem, INHERITS_FROM(CCTpfaModel, TwoP, CPTwoPSpatialParams));

// Set the grid type
SET_TYPE_PROP(CPTwoPProblem, Grid, Dune::CpGrid);

// Set the grid creator
SET_TYPE_PROP(CPTwoPProblem, GridCreator, Dumux::CpGridCreator<TypeTag>);

// Set the problem property
SET_TYPE_PROP(CPTwoPProblem, Problem, CPTwoPProblem<TypeTag>);

// Set the wetting phase
SET_PROP(CPTwoPProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(CPTwoPProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, DNAPL<Scalar> > type;
};

SET_TYPE_PROP(CPTwoPProblem, SpatialParams, CPTwoPSpatialParams<TypeTag> );

NEW_PROP_TAG(BaseProblem);
SET_TYPE_PROP(CPTwoPProblem, BaseProblem, ImplicitPorousMediaProblem<TypeTag>);

// Enable gravity
SET_BOOL_PROP(CPTwoPProblem, ProblemEnableGravity, true);

SET_BOOL_PROP(CPTwoPProblem, EnableGlobalFVGeometryCache, true);

SET_BOOL_PROP(CPTwoPProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(CPTwoPProblem, EnableGlobalVolumeVariablesCache, true);

//SET_TYPE_PROP(CPTwoPProblem, LinearSolver, SuperLUBackend<TypeTag> );

SET_BOOL_PROP(CPTwoPProblem, VtkWriteFaceData, false);

// The geometry helper required for the stencils, etc.
//SET_PROP(CPTwoPProblem, StaggeredGeometryHelper)
//{
//private:
//    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
//public:
//    using type = MimeticCPGeometryHelper<GridView>;
//};

SET_TYPE_PROP(CPTwoPProblem, IntersectionMapper, Dumux::NonConformingGridIntersectionMapper<TypeTag>);
}

/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular CPTwoP
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
class CPTwoPProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
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
        contiWEqIdx = Indices::contiWEqIdx,
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
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
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
    CPTwoPProblem(TimeManager &timeManager,
                const GridView &gridView)
    : ParentType(timeManager, gridView), gravity_(0)
    {
        gravity_[dimWorld-1] = 9.81;
        temperature_ = 273.15 + 20; // -> 20°C

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);

        injectionElement_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Problem, InjectionElement);
        injectionRate_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InjectionRate);
        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLength);
        //this->timeManager().startNextEpisode(episodeLength_);

        useFixedTimeSteps_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Problem, UseFixedTimeSteps);

        timeStepSizes_ = std::vector<int>
            ({10000,15000,21250,30105,42645,60419,85591,121254,151568, 202091,
            269454,179636,254484,318106,344614,459486,612648,765809,957262,1276350,
            1382712,1728389,2016455,1833141,2291426,2864282,3102972,3102973,3361553,3641683,
            3641683,2801295,2334412,2334412,2723481,2723481,1937878,2700000,3000000,3400000,
            3800000,4290000,5110000,3600000,4200000,5000000,5700000,5800000,3400000,5700000,
            6300000,7000000,7000000,7000000,7000000,7000000,3000000,7000000,8000000,8000000,
            10000000,7000000,10000000,9000000,7000000,7000000,10000000,6000000,8000000,3000000,
            8000000,9000000,7000000,7000000,7000000,6000000,6000000,5000000,7000000,6000000,
            5000000,6000000,4000000,3000000,3000000,4000000,4000000,3000000,4000000,4000000,
            5000000,6000000,7000000,6000000,6000000,8000000,4000000,7000000,8000000,7000000,
            8000000,7000000,5000000,5000000,3000000,5000000,5000000,6000000,7000000,5000000,
            7000000,7000000,8000000,7000000,7000000,4000000,7000000,5000000,6000000,7000000,
            7000000,7000000,6000000,7000000,6000000,6000000,7000000,4000000,3000000,4000000});
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

    const GlobalPosition &gravity() const
    { return gravity_; }

    PrimaryVariables source(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume &scv) const
    {
        PrimaryVariables values(0.0);

        const GlobalPosition& globalPos = scv.center();

        Scalar deltaXMin = 1e16;
        Scalar deltaXMax = -1e16;
        Scalar deltaYMin = 1e16;
        Scalar deltaYMax = -1e16;
        Scalar deltaZMin = 1e16;
        Scalar deltaZMax = -1e16;

        GlobalPosition deltaXYZ(0);

        std::vector<GlobalPosition> corners;
        for(int i=0; i < element.geometry().corners(); i++)
        {
            auto corner = element.geometry().corner(i);
            deltaXMin = std::min(deltaXMin,corner[0]);
            deltaXMax = std::max(deltaXMax,corner[0]);
            deltaYMin = std::min(deltaYMin,corner[1]);
            deltaYMax = std::max(deltaYMax,corner[1]);
            deltaZMin = std::min(deltaZMin,corner[2]);
            deltaZMax = std::max(deltaZMax,corner[2]);
        }

        deltaXYZ[0] = deltaXMax - deltaXMin;
        deltaXYZ[1] = deltaYMax - deltaYMin;
        deltaXYZ[2] = deltaZMax - deltaZMin;

        Scalar P1x0 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, P1WellXCoord);
        Scalar P1y0 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, P1WellYCoord);
        Scalar P1deviationX = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, P1DeviationX);
        Scalar P1deviationY = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, P1DeviationY);
        Scalar P2x0 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, P2WellXCoord);
        Scalar P2y0 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, P2WellYCoord);
        Scalar P2deviationX = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, P2DeviationX);
        Scalar P2deviationY = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, P2DeviationY);
        Scalar I1x0 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, I1WellXCoord);
        Scalar I1y0 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, I1WellYCoord);
        Scalar I1deviationX = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, I1DeviationX);
        Scalar I1deviationY = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, I1DeviationY);
        Scalar I2x0 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, I2WellXCoord);
        Scalar I2y0 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, I2WellYCoord);
        Scalar I2deviationX = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, I2DeviationX);
        Scalar I2deviationY = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, I2DeviationY);

        Scalar height = deltaXYZ[2];
        Scalar pi = 4.0*atan(1.0);
        Scalar rw = 0.15;
        Scalar re = 0.14*std::sqrt(deltaXYZ[0]*deltaXYZ[0] + deltaXYZ[1]*deltaXYZ[1]);
        Scalar pbhI = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, BoreHolePressureI);
        Scalar pbhP = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, BoreHolePressureP);
        auto K = this->spatialParams().permeability(element, scv,
                this->model().elementSolution(element, this->model().curSol()));

        auto volVars = elemVolVars[scv];
        Scalar pw = volVars.pressure(wPhaseIdx);
        Scalar pn = volVars.pressure(nPhaseIdx);
        Scalar mobN = volVars.mobility(nPhaseIdx);
        Scalar densityW = volVars.density(wPhaseIdx);
        Scalar densityN = volVars.density(nPhaseIdx);
        Scalar viscosityW = volVars.fluidState().viscosity(wPhaseIdx);

        if(std::abs(globalPos[0]-I1x0) < I1deviationX && std::abs(globalPos[1]-I1y0) < I1deviationY)
            values[contiWEqIdx] = (2.0*pi*height*K[0][0]*densityW)/(std::log(re/rw)*viscosityW) * (pbhI - (pw -densityW*(this->gravity()*globalPos)));
        else if(std::abs(globalPos[0]-I2x0) < I2deviationX && std::abs(globalPos[1]-I2y0) < I2deviationY)
            values[contiWEqIdx] = (2.0*pi*height*K[0][0]*densityW)/(std::log(re/rw)*viscosityW) * (pbhI - (pw -densityW*(this->gravity()*globalPos)));
        else if(std::abs(globalPos[0]-P1x0) < P1deviationX && std::abs(globalPos[1]-P1y0) < P1deviationY)
            values[contiNEqIdx] = (2.0*pi*height*K[0][0]*densityN)/(std::log(re/rw))*mobN* (pbhP - (pn -densityN*(this->gravity()*globalPos)));
        else if(std::abs(globalPos[0]-P2x0) < P2deviationX && std::abs(globalPos[1]-P2y0) < P2deviationY)
             values[contiNEqIdx] = (2.0*pi*height*K[0][0]*densityN)/(std::log(re/rw))*mobN* (pbhP - (pn -densityN*(this->gravity()*globalPos)));

        values[contiWEqIdx] /= scv.volume();
        values[contiNEqIdx] /= scv.volume();

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

        if(globalPos[0] > 461000.0)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

//        if(globalPos[0] < this->bBoxMin()[0] + eps_ || globalPos[0] > this->bBoxMax()[0] - eps_)
//            values.setAllDirichlet();
//        else
//            values.setAllNeumann();

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

        values[pwIdx] = 1e5 + densityW*(this->gravity()*globalPos);
        values[snIdx] = 1.0;

        Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParamsAtPos(globalPos), 1.0 - values[snIdx]);
        values[facePressureWIdx] = values[pwIdx];
        values[facePressureNIdx] = pc + values[facePressureWIdx];

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
        values[pwIdx] = 1e5 + densityW*(this->gravity()*globalPos);
        values[snIdx] = 1.0;

        Scalar pc = MaterialLaw::pc(this->spatialParams().materialLawParamsAtPos(globalPos), 1.0 - values[snIdx]);
        values[facePressureWIdx] = 0.0;
        values[facePressureNIdx] = 0.0;

        return values;
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    /*!
     * \brief Called by TimeManager whenever a solution for a
     *        time step has been computed and the simulation time has
     *        been updated.
     *
     * \param dt The current time-step size
     */
    Scalar nextTimeStepSize(const Scalar dt)
    {
        if(useFixedTimeSteps_)
        {
            Scalar dtEpisodeEnd = episodeLength_ - fmod(this->timeManager().time(), episodeLength_);
            Scalar dtSuggested = timeStepSizes_[this->timeManager().timeStepIndex()];
            return std::min(dtEpisodeEnd,dtSuggested);
        }
        else
            return ParentType::nextTimeStepSize(dt);
    }


    bool shouldWriteOutput() const
    {
//        return this->timeManager().timeStepIndex() == 0 ||
//               this->timeManager().episodeWillBeOver() ||
//               this->timeManager().willBeFinished();
        return true;
    }

//    void episodeEnd()
//    {
//        this->timeManager().startNextEpisode(episodeLength_);
//    }
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
                potW[ccDofIdx] = elemVolVars[scv].pressure(wPhaseIdx) + densityW*this->gravity()[dimWorld-1]*depth;
                potN[ccDofIdx] = elemVolVars[scv].pressure(nPhaseIdx) + densityN*this->gravity()[dimWorld-1]*depth;
            }
        }
    }

private:
    Scalar temperature_;
    static constexpr Scalar eps_ = 3e-6;
    std::string name_;
    GlobalPosition gravity_;
    int injectionElement_;
    Scalar injectionRate_;
    Scalar episodeLength_;
    std::vector<int> timeStepSizes_;
    bool useFixedTimeSteps_;
};
} //end namespace

#endif
