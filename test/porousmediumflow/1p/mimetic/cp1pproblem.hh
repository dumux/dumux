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

#ifndef DUMUX_MIMETIC_CPONEP_HH
#define DUMUX_MIMETIC_CPONEP_HH

//#include <dumux/implicit/cellcentered/tpfa/properties.hh>
//#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include <dumux/porousmediumflow/1p/mimetic/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include <dumux/io/cpgridcreator.hh>
#include <dumux/common/intersectionmapper.hh>

#include "cp1pspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class CPOnePProblem;

//////////
// Specify the properties for the CPTwoP problem
//////////
namespace Properties
{
NEW_TYPE_TAG(CPOnePProblem, INHERITS_FROM(OnePMimetic, CPOnePSpatialParams));
//NEW_TYPE_TAG(CPOnePProblem, INHERITS_FROM(CCTpfaModel, OneP));

// Set the grid type
SET_TYPE_PROP(CPOnePProblem, Grid, Dune::CpGrid);

// Set the grid creator
SET_TYPE_PROP(CPOnePProblem, GridCreator, Dumux::CpGridCreator<TypeTag>);

// Set the problem property
SET_TYPE_PROP(CPOnePProblem, Problem, CPOnePProblem<TypeTag>);

SET_PROP(CPOnePProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, SimpleH2O<Scalar> > type;
};


SET_TYPE_PROP(CPOnePProblem, SpatialParams, CPOnePSpatialParams<TypeTag> );

NEW_PROP_TAG(BaseProblem);
SET_TYPE_PROP(CPOnePProblem, BaseProblem, ImplicitPorousMediaProblem<TypeTag>);

// Enable gravity
SET_BOOL_PROP(CPOnePProblem, ProblemEnableGravity, true);

SET_BOOL_PROP(CPOnePProblem, EnableGlobalFVGeometryCache, true);

SET_BOOL_PROP(CPOnePProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(CPOnePProblem, EnableGlobalVolumeVariablesCache, true);

SET_TYPE_PROP(CPOnePProblem, LinearSolver, SuperLUBackend<TypeTag> );

SET_TYPE_PROP(CPOnePProblem, IntersectionMapper, Dumux::NonConformingGridIntersectionMapper<TypeTag>);
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
class CPOnePProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum {
        // indices of the primary variables
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx,
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    CPOnePProblem(TimeManager &timeManager,
                const GridView &gridView)
    : ParentType(timeManager, gridView), gravity_(0)
    {
        temperature_ = 273.15 + 20; // -> 20°C

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);

        gravity_[dimWorld-1] = 9.81;

        injectionRate_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InjectionRate);
        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLength);
        this->timeManager().startNextEpisode(episodeLength_);
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
        //const GlobalPosition& deltaXYZ = scv.deltaXYZ;

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

        //Scalar height = deltaXYZ[2];
        //Scalar pi = 4.0*atan(1.0);
        //Scalar rw = 0.15;
        //Scalar re = 0.14*std::sqrt(deltaXYZ[0]*deltaXYZ[0] + deltaXYZ[1]*deltaXYZ[1]);
        Scalar pbhI = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, BoreHolePressureI);
        Scalar pbhP = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, BoreHolePressureP);
        auto K = this->spatialParams().permeability(element, scv,
                this->model().elementSolution(element, this->model().curSol()));

        Scalar pw =  elemVolVars[scv].pressure(conti0EqIdx);
        Scalar densityW =  elemVolVars[scv].density(conti0EqIdx);
        Scalar viscosityW =  elemVolVars[scv].fluidState().viscosity(conti0EqIdx);

//        if(std::abs(globalPos[0]-I1x0) < I1deviationX && std::abs(globalPos[1]-I1y0) < I1deviationY)
//            values[contiWEqIdx] = (2.0*pi*height*K[0][0]*densityW)/(std::log(re/rw)*viscosityW) * (pbhI - (pw -densityW*(this->gravity()*globalPos)));
//        else if(std::abs(globalPos[0]-I2x0) < I2deviationX && std::abs(globalPos[1]-I2y0) < I2deviationY)
//            values[contiWEqIdx] = (2.0*pi*height*K[0][0]*densityW)/(std::log(re/rw)*viscosityW) * (pbhI - (pw -densityW*(this->gravity()*globalPos)));
//        else if(std::abs(globalPos[0]-P1x0) < P1deviationX && std::abs(globalPos[1]-P1y0) < P1deviationY)
//            values[contiWEqIdx] = (2.0*pi*height*K[0][0]*densityW)/(std::log(re/rw)*viscosityW)* (pbhP - (pw -densityW*(this->gravity()*globalPos)));
//        else if(std::abs(globalPos[0]-P2x0) < P2deviationX && std::abs(globalPos[1]-P2y0) < P2deviationY)
//             values[contiWEqIdx] = (2.0*pi*height*K[0][0]*densityW)/(std::log(re/rw)*viscosityW)*(pbhP - (pw -densityW*(this->gravity()*globalPos)));

        if(std::abs(globalPos[0]-I1x0) < I1deviationX && std::abs(globalPos[1]-I1y0) < I1deviationY)
            values[conti0EqIdx] = densityW/viscosityW * (pbhI - (pw -densityW*(this->gravity()*globalPos)));
        else if(std::abs(globalPos[0]-I2x0) < I2deviationX && std::abs(globalPos[1]-I2y0) < I2deviationY)
            values[conti0EqIdx] = densityW/viscosityW * (pbhI - (pw -densityW*(this->gravity()*globalPos)));
        else if(std::abs(globalPos[0]-P1x0) < P1deviationX && std::abs(globalPos[1]-P1y0) < P1deviationY)
            values[conti0EqIdx] = densityW/viscosityW * (pbhP - (pw -densityW*(this->gravity()*globalPos)));
        else if(std::abs(globalPos[0]-P2x0) < P2deviationX && std::abs(globalPos[1]-P2y0) < P2deviationY)
             values[conti0EqIdx] = densityW/viscosityW * (pbhP - (pw -densityW*(this->gravity()*globalPos)));


        values[conti0EqIdx]/=scv.volume();

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
        Scalar densityW = 1000;
        PrimaryVariables values(1e5 + densityW*(this->gravity()*globalPos));

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
        Scalar densityW = 1000;
        PrimaryVariables values(1e5 + densityW*(this->gravity()*globalPos));

        return values;
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }
    // \}

private:
    Scalar temperature_;
    static constexpr Scalar eps_ = 3e-6;
    std::string name_;
    GlobalPosition gravity_;
    Scalar injectionRate_;
    Scalar episodeLength_;
};
} //end namespace

#endif
