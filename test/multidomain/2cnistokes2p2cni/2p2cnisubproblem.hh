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
 * \brief Non-isothermal two-phase two-component porous-medium subproblem
 *        with coupling at the top boundary.
 */
#ifndef DUMUX_2P2CNISUB_PROBLEM_HH
#define DUMUX_2P2CNISUB_PROBLEM_HH

#include <dune/common/float_cmp.hh>

#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/porousmediumflow/2p2c/implicit/model.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/multidomain/2cnistokes2p2cni/2p2cnicouplinglocalresidual.hh>
#include <dumux/multidomain/subdomainpropertydefaults.hh>
#include <dumux/multidomain/localoperator.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>

#include "2cnistokes2p2cnispatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class TwoPTwoCNISubProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoPTwoCNISubProblem,
    INHERITS_FROM(BoxTwoPTwoCNI, SubDomain, TwoCNIStokesTwoPTwoCNISpatialParams));

// Set the problem property
SET_TYPE_PROP(TwoPTwoCNISubProblem, Problem, TwoPTwoCNISubProblem<TTAG(TwoPTwoCNISubProblem)>);

// Use the 2p2cni local jacobian operator for the 2p2cniCoupling model
SET_TYPE_PROP(TwoPTwoCNISubProblem, LocalResidual, TwoPTwoCNICouplingLocalResidual<TypeTag>);

// choose pn and Sw as primary variables
SET_INT_PROP(TwoPTwoCNISubProblem, Formulation, TwoPTwoCFormulation::pnsw);

// the gas component balance (air) is replaced by the total mass balance
SET_INT_PROP(TwoPTwoCNISubProblem, ReplaceCompEqIdx, GET_PROP_TYPE(TypeTag, Indices)::contiNEqIdx);

// Used the fluid system from the coupled problem
SET_TYPE_PROP(TwoPTwoCNISubProblem,
              FluidSystem,
              typename GET_PROP_TYPE(typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag), FluidSystem));

// Somerton is used as model to compute the effective thermal heat conductivity
SET_TYPE_PROP(TwoPTwoCNISubProblem, ThermalConductivityModel,
              ThermalConductivitySomerton<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// use formulation based on mass fractions
SET_BOOL_PROP(TwoPTwoCNISubProblem, UseMoles, false);

// enable/disable velocity output
SET_BOOL_PROP(TwoPTwoCNISubProblem, VtkAddVelocity, true);

// Enable gravity
SET_BOOL_PROP(TwoPTwoCNISubProblem, ProblemEnableGravity, true);
}

/*!
 * \ingroup ImplicitTestProblems
 * \ingroup TwoPTwoCNIStokesTwoCNIModel
 * \brief Non-isothermal two-phase two-component porous-medium subproblem
 *        with coupling at the top boundary.
 *
 * The Darcy subdomain is sized 0.25m times 0.25m. All BCs for the balance
 * equations are set to Neumann no-flow, except for the top, where couplingInflow
 * conditions are applied.
 *
 * This sub problem uses the \ref TwoPTwoCNIModel. It is part of the 2p2cni model and
 * is combined with the stokes2cnisubproblem for the free flow domain.
 */
template <class TypeTag = TTAG(TwoPTwoCNISubProblem) >
class TwoPTwoCNISubProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef TwoPTwoCNISubProblem<TypeTag> ThisType;
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { // the equation indices
        contiTotalMassIdx = Indices::contiNEqIdx,
        contiWEqIdx = Indices::contiWEqIdx,
        energyEqIdx = Indices::energyEqIdx
    };
    enum { // the indices of the primary variables
        pNIdx = Indices::pressureIdx,
        sWIdx = Indices::switchIdx,
        temperatureIdx = Indices::temperatureIdx
    };
    enum { // the indices for the phase presence
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx
    };
    enum { // the indices for the phase presence
        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases
    };
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };
    enum { // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The sub-problem for the porous-medium subdomain
     *
     * \param timeManager The TimeManager which is used by the simulation
     * \param gridView The simulation's idea about physical space
     */
    TwoPTwoCNISubProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        Scalar noDarcyX = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyX);
        std::vector<Scalar> positions0 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<Scalar>, Grid, Positions0);
        std::vector<Scalar> positions1 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<Scalar>, Grid, Positions1);

        bBoxMin_[0] = std::max(positions0.front(),noDarcyX);
        bBoxMax_[0] = positions0.back();

        bBoxMin_[1] = positions1.front();
        bBoxMax_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePosY);

        runUpDistanceX_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, RunUpDistanceX); // first part of the interface without coupling
        initializationTime_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, InitTime);

        refTemperature_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, RefTemperature);
        refPressure_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, RefPressure);
        initialSw_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, InitialSw);

        freqMassOutput_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqMassOutput);

        storageLastTimestep_ = Scalar(0);
        lastMassOutputTime_ = Scalar(0);

        evaporationFile.open("evaporation.out");
        evaporationFile << "#Time[d]" << " "
                        << "WaterMass[kg]"
                        << std::endl;
        liveEvaporationRates_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Output, LiveEvaporationRates);

        outfile.open("storage.out");
        outfile << "Time;"
                << "TotalMassChange;"
                << "WaterMassChange;"
                << "IntEnergyChange;"
                << "WaterMass"
                << std::endl;
    }
   //! \brief The destructor
    ~TwoPTwoCNISubProblem()
    {
        evaporationFile.close();
        outfile.close();
    }

    // functions have to be overwritten, otherwise they remain uninitialized
    //! \copydoc Dumux::ImplicitProblem::bBoxMin()
    const GlobalPosition &bBoxMin() const
    { return bBoxMin_; }

    //! \copydoc Dumux::ImplicitProblem::bBoxMax()
    const GlobalPosition &bBoxMax() const
    { return bBoxMax_; }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string &name() const
    { return GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Output, NamePM); }

    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem.
     *
     * If you overload this method don't forget to call
     * ParentType::init()
     */
    void init()
    {
        ParentType::init();

        this->model().globalStorage(storageLastTimestep_);
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
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        Scalar time = this->timeManager().time();

        values.setAllNeumann();

        if (onLowerBoundary_(globalPos))
        {
            values.setDirichlet(temperatureIdx);
        }
        else if (onUpperBoundary_(globalPos)
                 && (globalPos[0] > runUpDistanceX_ - eps_)
                 && (time > initializationTime_))
        {
            values.setAllCouplingNeumann();
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * \param values The Neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^{\textrm{dim}-1} \cdot s )] \f$
     * \param globalPos The global position
     */
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        values = 0.;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Returns the source term
     *
     * \param values Stores the source values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param globalPos The global position
     */
    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
    {
        values = Scalar(0);
    }

    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        values = 0.;

        initial_(values, globalPos);
    }

    /*!
     * \brief Return the initial phase state inside a control volume.
     *
     * \param vertex The vertex
     * \param vIdxGlobal The global index of the vertex
     * \param globalPos The global position
     */
    int initialPhasePresence(const Vertex &vertex,
                             const int &vIdxGlobal,
                             const GlobalPosition &globalPos) const
    {
        return bothPhases;
    }

    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    {
        // Calculate masses
        PrimaryVariables storage;

        this->model().globalStorage(storage);
        const Scalar time = this->timeManager().time() +  this->timeManager().timeStepSize();

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0)
        {
            if (this->timeManager().timeStepIndex() % freqMassOutput_ == 0
                    || this->timeManager().episodeWillBeOver())
            {
                PrimaryVariables storageChange(0.);
                storageChange = storageLastTimestep_ - storage;

                assert( (Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(time - lastMassOutputTime_, 0.0, 1.0e-30)) );
                storageChange /= (time - lastMassOutputTime_);
                // 2d: interface length has to be accounted for
                // in order to obtain kg/m²s
                storageChange /= (bBoxMax_[0]-bBoxMin_[0]);

                std::cout << "Time: " << time
                          << " TotalMass: " << storage[contiTotalMassIdx]
                          << " WaterMass: " << storage[contiWEqIdx]
                          << " IntEnergy: " << storage[energyEqIdx]
                          << " WaterMassChange: " << storageChange[contiWEqIdx]
                          << std::endl;
                if (Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(this->timeManager().time(), 0.0, 1.0e-30))
                {
                    outfile << time << ";"
                            << storageChange[contiTotalMassIdx] << ";"
                            << storageChange[contiWEqIdx] << ";"
                            << storageChange[energyEqIdx] << ";"
                            << storage[contiWEqIdx]
                            << std::endl;

                    evaporationFile << time/86400.0 << " " << storageChange[contiWEqIdx]*86400.0 << std::endl;
                    gnuplot_.reset();
                    gnuplot_.setInteraction(liveEvaporationRates_);
                    gnuplot_.setXRange(0.0, time/86400.0);
                    gnuplot_.setYRange(0.0, 12.0);
                    gnuplot_.setXlabel("time [d]");
                    gnuplot_.setYlabel("evaporation rate [mm/d]");
                    gnuplot_.addFileToPlot("evaporation.out", "evaporation.out");
                    gnuplot_.plot("evaporation");
                }

                storageLastTimestep_ = storage;
                lastMassOutputTime_ = time;
            }
        }
    }

    // \}

private:
    /*!
     * \brief Internal method for the initial condition
     *        (reused for the dirichlet conditions!)
     */
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        values[pNIdx] = refPressure_
                        + 1000.*this->gravity()[1]*(globalPos[1]-bBoxMax_[1]);
        values[sWIdx] = initialSw_;
        values[temperatureIdx] = refTemperature_;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < bBoxMin_[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > bBoxMax_[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < bBoxMin_[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > bBoxMax_[1] - eps_; }

    static constexpr Scalar eps_ = 1e-8;
    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;

    int freqMassOutput_;

    PrimaryVariables storageLastTimestep_;
    Scalar lastMassOutputTime_;

    Scalar refTemperature_;
    Scalar refPressure_;
    Scalar initialSw_;

    Scalar runUpDistanceX_;
    Scalar initializationTime_;
    std::ofstream outfile;

    Dumux::GnuplotInterface<Scalar> gnuplot_;
    std::ofstream evaporationFile;
    bool liveEvaporationRates_;
};
} //end namespace Dumux

#endif
