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
 * \brief Isothermal two-phase two-component porous-medium subproblem
 *        with coupling at the top boundary.
 */
#ifndef DUMUX_2P2C_SUBPROBLEM_HH
#define DUMUX_2P2C_SUBPROBLEM_HH

#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/2p2c/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include "2cstokes2p2cspatialparams.hh"

#include <dumux/linear/seqsolverbackend.hh>

// coupling-specific includes
#include <dumux/multidomain/subproblemproperties.hh>

// TODO necessary?
//#include <dune/common/float_cmp.hh>
//
//#include <dumux/porousmediumflow/2p2c/implicit/indices.hh>
//#include <dumux/multidomain/2cstokes2p2c/2p2ccouplinglocalresidual.hh>
//#include <dumux/multidomain/subdomainpropertydefaults.hh>
//#include <dumux/multidomain/localoperator.hh>


namespace Dumux
{
template <class TypeTag>
class TwoPTwoCSubProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoPTwoCSubProblem, INHERITS_FROM(CCTpfaModel, TwoPTwoC, TwoCStokesTwoPTwoCSpatialParams));

// Set the problem property
SET_TYPE_PROP(TwoPTwoCSubProblem, Problem, TwoPTwoCSubProblem<TTAG(TwoPTwoCSubProblem)>);

//// TODO from stokesdarcy1p, here: set in global problem?
// Set the spatial parameters
SET_TYPE_PROP(TwoPTwoCSubProblem, SpatialParams, TwoCStokesTwoPTwoCSpatialParams<TypeTag>);

// Choose pn and Sw as primary variables
SET_INT_PROP(TwoPTwoCSubProblem, Formulation, TwoPTwoCFormulation::pnsw);

// The gas component balance (air) is replaced by the total mass balance
SET_INT_PROP(TwoPTwoCSubProblem, ReplaceCompEqIdx, GET_PROP_TYPE(TypeTag, Indices)::contiNEqIdx);

// Used the fluid system from the coupled problem
SET_TYPE_PROP(TwoPTwoCSubProblem, FluidSystem,
        typename GET_PROP_TYPE(typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag), FluidSystem));

// Use formulation based on mass fractions
SET_BOOL_PROP(TwoPTwoCSubProblem, UseMoles, false);

// Enable velocity output
SET_BOOL_PROP(TwoPTwoCSubProblem, VtkAddVelocity, true);

// Enable gravity
SET_BOOL_PROP(TwoPTwoCSubProblem, ProblemEnableGravity, false); // TODO ?

// TODO not needed in old multidomain, used in darcytestproblem 1p
// Set the grid type
#if ENABLE_3D
SET_TYPE_PROP(TwoPTwoCSubProblem, Grid, Dune::YaspGrid<3>);
#else
SET_TYPE_PROP(TwoPTwoCSubProblem, Grid, Dune::YaspGrid<2>);
#endif
//
//// Set the grid parameter group
SET_STRING_PROP(TwoPTwoCSubProblem, GridParameterGroup, "DarcyGrid");

SET_BOOL_PROP(TwoPTwoCSubProblem, EnableGlobalFVGeometryCache, true); // TODO default = false, but true needed for couplingmanager (access Darcy element)

NEW_PROP_TAG(GlobalProblemTypeTag);
NEW_PROP_TAG(CouplingManager);

SET_BOOL_PROP(TwoPTwoCSubProblem, EnableGlobalVolumeVariablesCache, true);

}

/*!
 * \ingroup ImplicitTestProblems
 * \ingroup TwoPTwoCStokesTwoCModel
 * \brief Isothermal two-phase two-component porous-medium subproblem
 *        with coupling at the top boundary.
 *
 * The Darcy subdomain is sized 0.25m times 0.25m. All BCs for the balance
 * equations are set to Neumann no-flow, except for the top, where couplingInflow
 * conditions are applied.
 *
 * This sub problem uses the \ref TwoPTwoCModel. It is part of the 2p2c model and
 * is combined with the stokes2csubproblem for the free flow domain.
 */
//template <class TypeTag = TTAG(TwoPTwoCSubProblem) > // TODO original
//class TwoPTwoCSubProblem : public ImplicitPorousMediaProblem<TypeTag>
template <class TypeTag> // TODO stokesdarcy1p
class TwoPTwoCSubProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables) ;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry) ;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;

    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    // copy some indices for convenience TODO from original 2cstokes2p2c
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq)};
    enum { // the equation indices
        contiTotalMassIdx = Indices::contiNEqIdx,
        contiWEqIdx = Indices::contiWEqIdx,
    };
    enum { // the indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx,
    };
    enum { // the indices for the phase presence
        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases
    };
    enum { // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using Intersection = typename GridView::Intersection;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using GlobalTypeTag = typename GET_PROP_TYPE(TypeTag, GlobalProblemTypeTag);
    using CouplingManager = typename GET_PROP_TYPE(GlobalTypeTag, CouplingManager);

    enum { dofCodim = 0 };

    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    /*!
     * \brief The sub-problem for the porous-medium subdomain
     *
     * \param timeManager The TimeManager which is used by the simulation
     * \param gridView The simulation's idea about physical space
     */
    TwoPTwoCSubProblem(TimeManager &timeManager, const GridView gridView)
        : ParentType(timeManager, gridView)
    {
        Scalar noDarcyX = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyX);
        std::vector<Scalar> positions0 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<Scalar>, Grid, Positions0);
        std::vector<Scalar> positions1 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<Scalar>, Grid, Positions1);

        using std::max;
        bBoxMin_[0] = max(positions0.front(),noDarcyX);
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

        outfile.open("storage.out");
        outfile << "Time;"
                << "TotalMassChange;"
                << "WaterMassChange;"
                << "WaterMass"
                << std::endl;
    }

    ~TwoPTwoCSubProblem()
    {
        outfile.close();
    }

    // functions have to be overwritten, otherwise they remain uninitialized
    //! \copydoc ImplicitProblem::bBoxMin()
    const GlobalPosition &bBoxMin() const
    { return bBoxMin_; }

    //! \copydoc ImplicitProblem::bBoxMax()
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
     * \brief Called by the TimeManager in order to
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

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * \param globalPos The global position
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return refTemperature_;
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
        Scalar time = this->timeManager().time();

        values.setAllNeumann();

        if (onUpperBoundary_(globalPos)
            && (globalPos[0] > runUpDistanceX_ - eps_)
            && (time > initializationTime_))
        {
                values.setAllCouplingNeumann();
        }
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        initial_(values, globalPos);

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * \param globalPos The global position
     */
    PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values = Scalar(0);

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Returns the source term
     *
     * \param globalPos The global position
     */
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        values = Scalar(0);

        return values;
    }

    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        initial_(values, globalPos);

        return values;
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
                    || this->timeManager().episodeWillBeFinished())
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
                          << " WaterMassChange: " << storageChange[contiWEqIdx]
                          << std::endl;
                if (Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(this->timeManager().time(), 0.0, 1.0e-30))
                    outfile << time << ";"
                            << storageChange[contiTotalMassIdx] << ";"
                            << storageChange[contiWEqIdx] << ";"
                            << storage[contiWEqIdx]
                            << std::endl;

                storageLastTimestep_ = storage;
                lastMassOutputTime_ = time;
            }
        }
    }

    // \}

    /*!
     * \brief Set the coupling manager
     * \param couplingManager The coupling manager
     *
     */
    void setCouplingManager(std::shared_ptr<CouplingManager> couplingManager)
    {
        couplingManager_ = couplingManager;
    }

    /*!
     * \brief Get the coupling manager
     *
     */
    CouplingManager& couplingManager() const
    { return *couplingManager_; }

    bool onCouplingInterface(const GlobalPosition &globalPos) const
    {return onLowerBoundary_(globalPos); }

private:
    /*!
     * \brief Internal method for the initial condition
     *        (reused for the dirichlet conditions!)
     */
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        values[pressureIdx] = refPressure_
                              + 1000.*this->gravity()[1]*(globalPos[1]-bBoxMax_[1]);
        values[switchIdx] = initialSw_;
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

    std::shared_ptr<CouplingManager> couplingManager_;

};
} //end namespace Dumux

#endif // DUMUX_2P2C_SUBPROBLEM_HH
