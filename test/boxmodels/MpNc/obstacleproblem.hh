/*****************************************************************************
 *   Copyright (C) 2008-2010 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008-2009 by Bernd Flemisch                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/**
 * \file
 * \brief Problem where water and gas is injected by means of a
 *        dirchlet condition on the lower right of the domain which have to go
 *        around an obstacle with \f$10^3\f$ lower permeability.
 * \author Andreas Lauser, Klaus Mosthaf, Bernd Flemisch
 */
#ifndef DUMUX_OBSTACLEPROBLEM_HH
#define DUMUX_OBSTACLEPROBLEM_HH

#define USE_2P2C 0
#define OBSTACLE_FIND_CONVERGENCE_RADIUS 0

#include <dune/common/parametertreeparser.hh>

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#if USE_2P2C
#include <dumux/boxmodels/2p2c/2p2cmodel.hh>
#else
#include <dumux/boxmodels/MpNc/MpNcmodel.hh>
#endif

/*
#include <dumux/material/fluidsystems/simple_h2o_n2_system.hh>
#include <dumux/material/fluidsystems/simple_h2o_n2_tce_system.hh>
#include <dumux/material/fluidsystems/h2o_n2_system.hh>
#include <dumux/material/fluidsystems/h2o_n2_o2_system.hh>
#include <dumux/material/fluidsystems/h2o_h2_o2_system.hh>
#include <dumux/material/fluidsystems/h2o_h2_n2_o2_system.hh>
*/

#include <dumux/material/MpNcfluidsystems/h2on2fluidsystem.hh>

#include "obstaclespatialparameters.hh"

namespace Dumux
{

template <class TypeTag>
class ObstacleProblem;

namespace Properties
{
#if USE_2P2C
NEW_TYPE_TAG(ObstacleProblem, INHERITS_FROM(BoxTwoPTwoC, ObstacleSpatialParameters));
#else
NEW_TYPE_TAG(ObstacleProblem, INHERITS_FROM(BoxMPNC, ObstacleSpatialParameters));
#endif

// Set the grid type
SET_TYPE_PROP(ObstacleProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ObstacleProblem,
              Problem,
              Dumux::ObstacleProblem<TypeTag>);


// Set fluid configuration
SET_PROP(ObstacleProblem, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::H2ON2FluidSystem<Scalar> type;
};
//              Dumux::Simple_H2O_N2_TCE_System<TypeTag> );
//Dumux::H2O_N2_System<TypeTag> );

#if ! USE_2P2C
// Enable smooth upwinding?
SET_BOOL_PROP(ObstacleProblem, EnableSmoothUpwinding, false);

// Enable molecular diffusion of the components?
SET_BOOL_PROP(ObstacleProblem, EnableDiffusion, false);

// Use the chopped Newton method?
SET_BOOL_PROP(ObstacleProblem, NewtonEnableChop, true);
#endif

// Enable gravity
SET_BOOL_PROP(ObstacleProblem, EnableGravity, true);

// Write Newton convergence to disk?
SET_BOOL_PROP(ObstacleProblem, NewtonWriteConvergence, false);

// Use the line search strategy for the Newton update?
SET_BOOL_PROP(ObstacleProblem, NewtonUseLineSearch, false);

// Enable the re-use of the jacobian matrix whenever possible?
SET_BOOL_PROP(ObstacleProblem, EnableJacobianRecycling, true);

// Reassemble the jacobian matrix only where it changed?
SET_BOOL_PROP(ObstacleProblem, EnablePartialReassemble, true);

// use forward diffferences to approximate the partial derivatives
SET_INT_PROP(ObstacleProblem, NumericDifferenceMethod, +1);
}


/*!
 * \ingroup MPNCModel
 * \ingroup BoxTestProblems
 * \brief Problem where liquid water is injected by means of a
 *        dirchlet condition on the lower right of the domain which have to go
 *        around an obstacle with \f$10^3\f$ lower permeability.
 *
 * The domain is sized 60m times 40m and consists of two media, a
 * moderately permeable soil (\f$ K_0=10e-12 m^2\f$) and an obstacle
 * at \f$[10; 20]m \times [0; 35]m \f$ with a lower permeablility of
 * \f$ K_1=K_0/1000\f$.
 *
 * Initially the whole domain is filled by nitrogen, the temperature
 * is \f$20\symbol{23}C\f$ at the whole domain. The gas pressure at
 * the left side of the domain is 1 bar, at the right side it is 2 bar
 * with a linear gradient in between.
 *
 * The boundary is no-flow except on the lower 10 meters of the left
 * and the right boundary which is a Dirichlet condition with the same
 * values as the initial condition.
 */
template <class TypeTag>
class ObstacleProblem
#if USE_2P2C
    : public TwoPTwoCProblem<TypeTag>
#else
    : public MPNCProblem<TypeTag>
#endif
{
    typedef ObstacleProblem<TypeTag>             ThisType;
#if USE_2P2C
    typedef TwoPTwoCProblem<TypeTag> ParentType;
#else
    typedef MPNCProblem<TypeTag> ParentType;
#endif

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, PTAG(NumPhases)),
        numComponents = GET_PROP_VALUE(TypeTag, PTAG(NumComponents)),

        gPhaseIdx = FluidSystem::gPhaseIdx,
        //oPhaseIdx = FluidSystem::oPhaseIdx,
        lPhaseIdx = FluidSystem::lPhaseIdx,
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef Dune::FieldVector<typename GridView::Grid::ctype, dimWorld> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

#if USE_2P2C
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices)) Indices;

#else // ! USE_2P2C
    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MPNCIndices)) Indices;
    enum {
        fug0Idx = Indices::fug0Idx,
        S0Idx = Indices::S0Idx,
        p0Idx = Indices::p0Idx,
    };
#endif // USE_2P2C
public:
    ObstacleProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        temperature_ = 273.15 + 25; // -> 25°C

        // initialize the tables of the fluid system
        Scalar Tmin = temperature_ - 1.0;
        Scalar Tmax = temperature_ + 1.0;
        int nT = 10;
        Scalar pmin = 0.75 * 1e5;
        Scalar pmax = 1.25 * 2e5;
        int np = 1000;
        FluidSystem::init(Tmin, Tmax, nT, pmin, pmax, np);
    }

    /*!
     * \brief Called by Dumux::TimeManager in order to do a time
     *        integration on the model.
     */
    void timeIntegration()
    {
#if !OBSTACLE_FIND_CONVERGENCE_RADIUS
        if (GET_PROP_VALUE(TypeTag, PTAG(NewtonWriteConvergence)))
            this->newtonController().setMaxSteps(40);
        ParentType::timeIntegration();
#else
        std::cout << "Finding convergence radius\n";
        this->newtonController().setVerbose(false);
        this->newtonController().setMaxSteps(40);
        this->newtonController().setTargetSteps(20);
        Scalar successDt = 0.0;
        const int maxFails = 10;
        for (int i = 0; true; ++i) {
            std::cout << "Try dt of " << this->timeManager().timeStepSize() << "\n";
            std::cout.flush();

            if (i == maxFails && this->gridView().comm().rank() == 0)
                DUNE_THROW(Dune::MathError,
                           "Newton solver didn't converge after "
                           << maxFails
                           << " timestep divisions. dt="
                           << this->timeManager().timeStepSize());

            if (this->model().update(this->newtonMethod(),
                                     this->newtonController()))
            {
                // sucessfull update. remember time step size
                successDt = this->timeManager().timeStepSize();
                this->model().updateFailed();
                break;
            }

            if (i > 0 && this->gridView().comm().rank() == 0)
                std::cout << "Newton solver did not converge. Retrying with time step of "
                          << this->timeManager().timeStepSize() << "sec\n";

            // update failed
            Scalar dt = this->timeManager().timeStepSize();
            Scalar nextDt = dt / 2;
            this->timeManager().setTimeStepSize(nextDt);
            std::cout << "Failed for dt=" << dt << ". Try with half dt of " << nextDt << "\n";
        }

        // increase time step until the update fails
        while (true)
        {
            if (!this->model().update(this->newtonMethod(),
                                      this->newtonController()))
                break;

            // sucessfull update. increase time step size
            successDt = this->timeManager().timeStepSize();
            Scalar nextDt = successDt*1.25;
            this->timeManager().setTimeStepSize(nextDt);
            if (this->timeManager().timeStepSize() < nextDt) {
                std::cout << "End of simulation reached!\n";
                break;
            };
            std::cout << "Increase dt to " << nextDt << "\n";
            std::cout.flush();
            this->model().updateFailed();
        }

        // do a last update with the largest successful time step
        this->newtonController().setVerbose(true);
        this->timeManager().setTimeStepSize(successDt);
        std::cout << "Convergence radius is " << successDt << "\n";
        this->model().update(this->newtonMethod(), this->newtonController());
#endif
    }

    /*!
     * \brief Called directly after the time integration.
     */
    void postTimeStep()
    {
        // Calculate storage terms of the individual phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            PrimaryVariables phaseStorage;
            this->model().globalPhaseStorage(phaseStorage, phaseIdx);

            if (this->gridView().comm().rank() == 0) {
                std::cout
                    <<"Storage in "
                    << FluidSystem::phaseName(phaseIdx)
                    << "Phase: ["
                    << phaseStorage
                    << "]"
                    << "\n";
            }
        }

        // Calculate total storage terms
        PrimaryVariables storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout
                <<"Storage total: [" << storage << "]"
                << "\n";
        }
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    { return "obstacle"; }

    /*!
     * \brief Returns the temperature [K] within the domain.
     */
    Scalar temperature() const
    { return temperature_; };

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param vertex The vertex for which the boundary type is set
     */
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        if (onInlet_(globalPos) || onOutlet_(globalPos))
            values.setAllDirichlet();
        else
            values.setAllNeumann();
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex for which the boundary type is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        initial_(values, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 const Intersection &is,
                 int scvIdx,
                 int boundaryFaceIdx) const
    { values = 0; }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     */
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvElemGeom,
                int scvIdx) const
    {
        values = Scalar(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 int scvIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        initial_(values, globalPos);
        Valgrind::CheckDefined(values);
    }

    // \}

    /*!
     * \brief Write a restart file?
     */
    bool shouldWriteRestartFile() const
    {
        return ParentType::shouldWriteRestartFile();
    }


#if USE_2P2C
    /*!
     * \brief Return the initial phase state inside a control volume.
     */
   int initialPhasePresence(const Vertex &vert,
                             int &globalIdx,
                             const GlobalPosition &globalPos) const
    {
        if (onInlet_(globalPos))
            return Indices::lPhaseOnly;
        else
            return Indices::gPhaseOnly;
    };
#endif

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar pg;
        if (onInlet_(globalPos))
            pg = 2e5;
        else
            pg = 1e5;


#if USE_2P2C
        values[Indices::plIdx] = pg;

        if (onInlet_(globalPos))
            values[Indices::SgOrXIdx] = 0.0; // X^a_w
        else
            values[Indices::SgOrXIdx] = 0.0; // X^w_g
#else
        if (onInlet_(globalPos)) {
            // only liquid
            Scalar S[numPhases];
            Scalar p[numPhases];
            Scalar xl[numComponents];
            Scalar beta[numComponents];

            p[gPhaseIdx] = pg;
            p[lPhaseIdx] = pg;

            S[lPhaseIdx] = 1.0;
            S[gPhaseIdx] = 0.0;

            xl[FluidSystem::H2OIdx]     = 1.0;
            xl[FluidSystem::N2Idx]      = 0.0;
            beta[FluidSystem::H2OIdx]   = FluidSystem::H2O::vaporPressure(temperature_);
            beta[FluidSystem::N2Idx]    = Dumux::BinaryCoeff::H2O_N2::henry(temperature_);

            // assign the primary variables
            for (int i = 0; i < numComponents; ++i)
                values[fug0Idx + i] = xl[i]*beta[i];

            for (int i = 0; i < numPhases - 1; ++i)
                values[S0Idx + i] = S[i];

            values[p0Idx] = p[0];
        }
        else {
            // only gas
            Scalar S[numPhases];
            Scalar xg[numComponents];
            Scalar p[numPhases];

            S[lPhaseIdx] = 0.0;
            S[gPhaseIdx] = 1.0;

            p[lPhaseIdx] = pg;
            p[gPhaseIdx] = pg;


            xg[FluidSystem::H2OIdx] = 0.01;
            xg[FluidSystem::N2Idx] = 0.99;


            // assign the primary variables
            for (int i = 0; i < numComponents; ++i)
                values[fug0Idx + i] = xg[i]*pg;

            for (int i = 0; i < numPhases - 1; ++i)
                values[S0Idx + i] = S[i];

            values[p0Idx] = p[0];
        }
#endif // USE_2P2C
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        return x >= 60 - eps_ && y <= 10;
    };

    bool onOutlet_(const GlobalPosition &globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        return x < eps_ && y <= 10;
    };

    Scalar temperature_;
    static constexpr Scalar eps_ = 1e-6;
};
} //end namespace

#endif
