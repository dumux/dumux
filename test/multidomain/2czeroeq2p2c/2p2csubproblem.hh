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
 * \brief Isothermal two-phase two-component porous-medium subproblem
 *        with coupling at the top boundary.
 */
#ifndef DUMUX_TWOPTWOC_SUBPROBLEM_HH
#define DUMUX_TWOPTWOC_SUBPROBLEM_HH

#include <dumux/implicit/2p2c/2p2cindices.hh>
#include <dumux/implicit/common/implicitporousmediaproblem.hh>
#include <dumux/multidomain/common/subdomainpropertydefaults.hh>
#include <dumux/multidomain/common/multidomainlocaloperator.hh>
#include <dumux/multidomain/couplinglocalresiduals/2p2ccouplinglocalresidual.hh>

#include "2czeroeq2p2cspatialparameters.hh"

namespace Dumux
{
template <class TypeTag>
class TwoPTwoCSubProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoPTwoCSubProblem,
             INHERITS_FROM(BoxTwoPTwoC, SubDomain, TwoCZeroEqTwoPTwoCSpatialParams));

// Set the problem property
SET_TYPE_PROP(TwoPTwoCSubProblem, Problem, TwoPTwoCSubProblem<TTAG(TwoPTwoCSubProblem)>);

// Set the local residual extended for the coupling
SET_TYPE_PROP(TwoPTwoCSubProblem, LocalResidual, TwoPTwoCCouplingLocalResidual<TypeTag>);

// Set pn and Sw as primary variables
SET_INT_PROP(TwoPTwoCSubProblem, Formulation, TwoPTwoCFormulation::pnsw);

// Set the gas component balance (air) to be replaced by the total mass balance
SET_PROP(TwoPTwoCSubProblem, ReplaceCompEqIdx)
{
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    static const int value = Indices::contiNEqIdx;
};

// Used the fluid system from the coupled problem
SET_TYPE_PROP(TwoPTwoCSubProblem,
              FluidSystem,
              typename GET_PROP_TYPE(typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag), FluidSystem));

// Disable use of mole formulation
SET_BOOL_PROP(TwoPTwoCSubProblem, UseMoles, false);

// Enable velocity output
SET_BOOL_PROP(TwoPTwoCSubProblem, VtkAddVelocity, true);

// Enable gravity
SET_BOOL_PROP(TwoPTwoCSubProblem, ProblemEnableGravity, true);
}

/*!
 * \ingroup ImplicitTestProblems
 * \ingroup TwoPTwoCZeroEqTwoCModel
 * \brief Isothermal two-phase two-component porous-medium subproblem
 *        with coupling at the top boundary.
 *
 * The porous-medium subdomain is sized 0.25m times 0.25m. The boundary conditions
 * are Neumann no-flow everywhere, except at the top, where coupling conditions
 * are applied to all balance equations. They handle the exchange to the free-flow
 * subdomain.
 *
 * This subproblem uses the \ref TwoPTwoCModel. It is part of a multidomain model and
 * combined with the zeroeq2csubproblem for the free flow domain.
 */
template <class TypeTag = TTAG(TwoPTwoCSubProblem) >
class TwoPTwoCSubProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::Grid Grid;

    typedef TwoPTwoCSubProblem<TypeTag> ThisType;
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    // the type tag of the coupled problem
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag) CoupledTypeTag;

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

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

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
        Scalar xMin = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, XMin);

        bBoxMin_[0] = std::max(xMin,noDarcyX);
        bBoxMax_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, XMax);
        bBoxMin_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, YMin);
        bBoxMax_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePos);
        runUpDistanceX_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, RunUpDistanceX); // first part of the interface without coupling

        refTemperature_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, RefTemperaturePM);
        refPressure_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, RefPressurePM);
        refSw_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, RefSw);

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

    //! \brief The destructor
    ~TwoPTwoCSubProblem()
    {
        outfile.close();
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

    //! \copydoc ImplicitProblem::boundaryTypesAtPos()
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        values.setAllNeumann();

        if (onUpperBoundary_(globalPos))
        {
            if (globalPos[0] > runUpDistanceX_ - eps_)
                values.setAllCouplingInflow();
            else
                values.setAllNeumann();
        }
    }

    //! \copydoc ImplicitProblem::dirichletAtPos()
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    //! \copydoc ImplicitProblem::neumannAtPos()
    void neumannAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        values = 0.;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    //! \copydoc ImplicitProblem::sourceAtPos()
    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
    {
        values = Scalar(0);
    }

    //! \copydoc ImplicitProblem::initialAtPos()
    void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    // \}

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

                assert(time - lastMassOutputTime_ != 0);
                storageChange /= (time - lastMassOutputTime_);
                // 2d: interface length has to be accounted for
                // in order to obtain kg/m²s
                storageChange /= (bBoxMax_[0]-bBoxMin_[0]);

                std::cout << "Time: " << time
                          << " TotalMass: " << storage[contiTotalMassIdx]
                          << " WaterMass: " << storage[contiWEqIdx]
                          << " WaterMassChange: " << storageChange[contiWEqIdx]
                          << std::endl;
                if (this->timeManager().time() != 0.)
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

    /*!
     * \brief Determines if globalPos is a corner of the grid
     *
     * \param globalPos The global position
     */
    bool isCornerPoint(const GlobalPosition &globalPos)
    {
        return ((onLeftBoundary_(globalPos) && onLowerBoundary_(globalPos))
                || (onLeftBoundary_(globalPos) && onUpperBoundary_(globalPos))
                || (onRightBoundary_(globalPos) && onLowerBoundary_(globalPos))
                || (onRightBoundary_(globalPos) && onUpperBoundary_(globalPos)));
    }

    /*!
     * \brief Returns whether the position is an interface corner point
     *
     * This function is required in case of mortar coupling otherwise it should return false
     *
     * \param globalPos The global position
     */
    bool isInterfaceCornerPoint(const GlobalPosition &globalPos) const
    { return false; }

private:
    // Internal method for the initial condition (reused for the dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        values[pressureIdx] = refPressure_
                              + 1000. * this->gravity()[1] * (globalPos[1] - bBoxMax_[1]);
        values[switchIdx] = refSw_;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < bBoxMin_[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > bBoxMax_[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < bBoxMin_[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > bBoxMax_[1] - eps_; }

    bool onBoundary_(const GlobalPosition &globalPos) const
    {
        return (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)
                || onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos));
    }

    static constexpr Scalar eps_ = 1e-8;
    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;
    Scalar runUpDistanceX_;

    Scalar refTemperature_;
    Scalar refPressure_;
    Scalar refSw_;

    PrimaryVariables storageLastTimestep_;
    Scalar lastMassOutputTime_;
    int freqMassOutput_;
    std::ofstream outfile;
};
} //end namespace Dumux

#endif // DUMUX_TWOPTWOC_SUBPROBLEM_HH
