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
/**
 * \file
 * \brief Problem where liquid water is injected by means of a
 *        dirchlet condition on the left of the domain.
 *        Velocity according to Forchheimer.
 */
#ifndef DUMUX_FORCHEIMER1P_HH
#define DUMUX_FORCHEIMER1P_HH

#include <dune/common/parametertreeparser.hh>

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/implicit/mpnc/mpncmodel.hh>
#include <dumux/implicit/common/implicitporousmediaproblem.hh>

#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>
#include <dumux/material/constraintsolvers/computefromreferencephase.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

#include "forchheimerspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class Forchheimer1pProblem;

namespace Properties
{
NEW_TYPE_TAG(Forchheimer1pProblem, INHERITS_FROM(BoxMPNC, ForchheimerSpatialParams));

// Set the grid type
SET_TYPE_PROP(Forchheimer1pProblem, Grid, Dune::YaspGrid<1>);

// Set the problem property
SET_TYPE_PROP(Forchheimer1pProblem,
              Problem,
              Dumux::Forchheimer1pProblem<TypeTag>);


// Set fluid configuration
SET_PROP(Forchheimer1pProblem, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::FluidSystems::H2OAir<Scalar, Dumux::SimpleH2O<Scalar>, /*useComplexRelations=*/false> type;
};

// Enable molecular diffusion of the components?
SET_BOOL_PROP(Forchheimer1pProblem, EnableDiffusion, false);

// Enable gravity
SET_BOOL_PROP(Forchheimer1pProblem, ProblemEnableGravity, false);

// decide which type to use for floating values (double / quad)
SET_TYPE_PROP(Forchheimer1pProblem, Scalar, double);

// decide which to use for velocity calculation: Darcy / Forchheimer
SET_TYPE_PROP(Forchheimer1pProblem, BaseFluxVariables, ImplicitForchheimerFluxVariables<TypeTag>);

SET_BOOL_PROP(Forchheimer1pProblem, VtkAddVelocities, true);
}


/*!
 * \ingroup MPNCModel
 * \ingroup ImplicitTestProblems
 * \brief Problem where liquid water is injected by means of a
 *        dirichlet condition on the left of the domain.
 *        Velocity according to Forchheimer.
 *
 * The setup is for testing of the one phase Forchheimer relation.
 * The whole domain is filled with water at 1e5 Pa. On the left hand
 * side the pressure is raised to 42e5 Pa.
 *
 * The induced flow field is calculated by means of the Forchheimer relation.
 * This selection is chosen via the BaseFluxVariables property, which can also
 * be set to the Darcy relation.
 *
 * The setup is on purpose very simple allowing for the analytical calculation
 * of the correct velocity.
 * Just paste the following into e.g. matlab
 * (viscosity, forchheimer coefficient, density, permeability, pressure gradient)
 *  -> mu = 1e-03 ; cE =0.55; rho = 1000; K = 1e-12; a = mu  /( cE *rho* sqrt(K) );
 *  -> gradP = -41e5;  v = - a /2 + sqrt (a^2/4-sqrt(K)*gradP/(rho *cE))
 *
 * This problem uses the \ref MPNCModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_forchheimer1p -parameterFile test_forchheimer1p.input</tt>
 */
template <class TypeTag>
class Forchheimer1pProblem
    : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename FluidSystem::ParameterCache ParameterCache;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    // Grid and world dimension
    enum {dim = GridView::dimension};
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

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef Dune::FieldVector<typename GridView::Grid::ctype, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    Forchheimer1pProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        try
        {
            pMax_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.pMax);
            pMin_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.pMin);
            outputName_     = GET_RUNTIME_PARAM(TypeTag, std::string, Problem.outputName);
        }
        catch (Dumux::ParameterException &e) {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }
        catch (...) {
            std::cerr << "Unknown exception thrown!\n";
            exit(1);
        }

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
     * \brief User defined output after the time integration
     *
     * Will be called diretly after the time integration.
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
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    { return outputName_; }

    /*!
     * \brief Returns the temperature \f$ K \f$
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
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
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param vertex The vertex for which the boundary type is set
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        if (onLeftBoundary_(globalPos) or onRightBoundary_(globalPos))
            values.setAllDirichlet();
        else
            values.setAllNeumann();
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param vertex The vertex for which the boundary type is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
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
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &intersection,
                 const unsigned int scvIdx,
                 const unsigned int boundaryFaceIdx) const
    { values = 0.; }

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
     *         \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     *
     * Positive values mean that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const unsigned int scvIdx) const
    {
        values = Scalar(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values Stores the solution for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
        Valgrind::CheckDefined(values);
    }

    // \}

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        FluidState fs;

        int refPhaseIdx;
        int otherPhaseIdx;

        // set the fluid temperatures
        fs.setTemperature(this->temperatureAtPos(globalPos));


        if (onLeftBoundary_(globalPos)) {
            // set pressure of the liquid phase
            fs.setPressure(wPhaseIdx, pMax_);
        }
        else {
            const Scalar interpolation = (globalPos[0] - this->bBoxMin()[0] ) / (this->bBoxMax()[0] - this->bBoxMin()[0]);
            const Scalar interpolatedPressure = pMax_ - interpolation * (pMax_ - pMin_) ;

            // set pressure of the gas phase
            fs.setPressure(wPhaseIdx, interpolatedPressure);
        }

        const Scalar belowSolubilityMassFraction= 1e-5;
        // set the liquid composition to equilibrium concentration
        fs.setMoleFraction(wPhaseIdx, nCompIdx, belowSolubilityMassFraction);
        fs.setMoleFraction(wPhaseIdx, wCompIdx, 1.-belowSolubilityMassFraction);
        // only liquid on inlet
        refPhaseIdx = wPhaseIdx;
        otherPhaseIdx = nPhaseIdx;


        // set liquid saturation
        fs.setSaturation(wPhaseIdx, 1.);

        // set the other saturation
        fs.setSaturation(otherPhaseIdx, 1.0 - fs.saturation(refPhaseIdx));

        // calulate the capillary pressure
        const MaterialLawParams &matParams =
            this->spatialParams().materialLawParamsAtPos(globalPos);
        PhaseVector pc;
        MaterialLaw::capillaryPressures(pc, matParams, fs);
        fs.setPressure(otherPhaseIdx,
                       fs.pressure(refPhaseIdx)
                       + (pc[otherPhaseIdx] - pc[refPhaseIdx]));

        // make the fluid state consistent with local thermodynamic
        // equilibrium
        typedef Dumux::ComputeFromReferencePhase<Scalar, FluidSystem> ComputeFromReferencePhase;

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
    }

private:
    /*!
     * \brief Give back whether the testes position (input) is a specific region (left) in the domain
     */
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {       return globalPos[0] < this->bBoxMin()[0] + eps_;   }

    /*!
     * \brief Give back whether the testes position (input) is a specific region (right) in the domain
     */
    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {        return globalPos[0] > this->bBoxMax()[0] - eps_;    }

    /*!
     * \brief Give back whether the testes position (input) is a specific region (down, (gravityDir)) in the domain
     */
    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {        return globalPos[dimWorld-1] < this->bBoxMin()[dimWorld-1] + eps_;    }

    /*!
     * \brief Give back whether the testes position (input) is a specific region (up, (gravityDir)) in the domain
     */
    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {        return globalPos[dimWorld-1] > this->bBoxMax()[dimWorld-1] - eps_;    }

    Scalar temperature_;
    Scalar eps_;
    std::string outputName_;
    Scalar pMax_, pMin_ ;

};
} //end namespace

#endif
