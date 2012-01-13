// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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

#include <dune/common/parametertreeparser.hh>

#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/boxmodels/MpNc/MpNcmodel.hh>

#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include "obstaclespatialparameters.hh"

namespace Dumux
{

template <class TypeTag>
class ObstacleProblem;

namespace Properties
{
NEW_TYPE_TAG(ObstacleProblem, INHERITS_FROM(BoxMPNC, ObstacleSpatialParameters));

// Set the grid type
SET_TYPE_PROP(ObstacleProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ObstacleProblem,
              Problem,
              Dumux::ObstacleProblem<TypeTag>);


// Set fluid configuration
SET_PROP(ObstacleProblem, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::FluidSystems::H2ON2<Scalar, /*useComplexRelations=*/false> type;
};
//              Dumux::Simple_H2O_N2_TCE_System<TypeTag> );
//Dumux::H2O_N2_System<TypeTag> );

// Enable smooth upwinding?
SET_BOOL_PROP(ObstacleProblem, EnableSmoothUpwinding, false);

// Enable molecular diffusion of the components?
SET_BOOL_PROP(ObstacleProblem, EnableDiffusion, false);

// Use the chopped Newton method?
SET_BOOL_PROP(ObstacleProblem, NewtonEnableChop, true);

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
    : public MPNCProblem<TypeTag>
{
    typedef MPNCProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents),

        gPhaseIdx = FluidSystem::gPhaseIdx,
        lPhaseIdx = FluidSystem::lPhaseIdx,

        H2OIdx = FluidSystem::H2OIdx,
        N2Idx = FluidSystem::N2Idx
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef Dune::FieldVector<typename GridView::Grid::ctype, dimWorld> GlobalPosition;

    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;

    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, MPNCIndices) Indices;
    enum {
        fug0Idx = Indices::fug0Idx,
        S0Idx = Indices::S0Idx,
        p0Idx = Indices::p0Idx
    };

public:
    ObstacleProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 1e-6;
        temperature_ = 273.15 + 25; // -> 25°C

        // initialize the tables of the fluid system
        Scalar Tmin = temperature_ - 1.0;
        Scalar Tmax = temperature_ + 1.0;
        int nT = 3;

        Scalar pmin = 1.0e5 * 0.75;
        Scalar pmax = 2.0e5 * 1.25;
        int np = 1000;

        FluidSystem::init(Tmin, Tmax, nT, pmin, pmax, np);
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
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
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

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> CompositionalFluidState;
        CompositionalFluidState fs;

        int refPhaseIdx;
        int otherPhaseIdx;

        // set the fluid temperatures
        fs.setTemperature(this->temperatureAtPos(globalPos));
            
        if (onInlet_(globalPos)) {
            // only liquid on inlet
            refPhaseIdx = lPhaseIdx;
            otherPhaseIdx = gPhaseIdx;
            
            // set liquid saturation
            fs.setSaturation(lPhaseIdx, 1.0);
            
            // set pressure of the liquid phase
            fs.setPressure(lPhaseIdx, 2e5);

            // set the liquid composition to pure water
            fs.setMoleFraction(lPhaseIdx, N2Idx, 0.0);
            fs.setMoleFraction(lPhaseIdx, H2OIdx, 1.0);
        }
        else {
            // elsewhere, only gas
            refPhaseIdx = gPhaseIdx;
            otherPhaseIdx = lPhaseIdx;

            // set gas saturation
            fs.setSaturation(gPhaseIdx, 1.0);
            
            // set pressure of the gas phase
            fs.setPressure(gPhaseIdx, 1e5);

            // set the gas composition to 99% nitrogen and 1% steam
            fs.setMoleFraction(gPhaseIdx, N2Idx, 0.99);
            fs.setMoleFraction(gPhaseIdx, H2OIdx, 0.01);
        };
        
        // set the other saturation
        fs.setSaturation(otherPhaseIdx, 1.0 - fs.saturation(refPhaseIdx));

        // calulate the capillary pressure
        const MaterialLawParams &matParams = 
            this->spatialParameters().materialLawParamsAtPos(globalPos);
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fs);
        fs.setPressure(otherPhaseIdx, 
                       fs.pressure(refPhaseIdx)
                       + (pC[otherPhaseIdx] - pC[refPhaseIdx]));
        
        // make the fluid state consistent with local thermodynamic
        // equilibrium
        typedef Dumux::ComputeFromReferencePhase<Scalar, FluidSystem> ComputeFromReferencePhase;

        typename FluidSystem::ParameterCache paramCache;
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
            values[S0Idx + phaseIdx] = fs.saturation(phaseIdx);

        // first pressure
        values[p0Idx] = fs.pressure(/*phaseIdx=*/0);
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
    Scalar eps_;
};
} //end namespace

#endif
