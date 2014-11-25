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

#include <dumux/implicit/2p2c/2p2cindices.hh>
#include <dumux/implicit/common/implicitporousmediaproblem.hh>
#include <dumux/multidomain/couplinglocalresiduals/2p2ccouplinglocalresidual.hh>
#include <dumux/multidomain/common/subdomainpropertydefaults.hh>
#include <dumux/multidomain/common/multidomainlocaloperator.hh>

#include "2cstokes2p2cspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class TwoPTwoCSubProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoPTwoCSubProblem,
             INHERITS_FROM(BoxTwoPTwoC, SubDomain, TwoCStokesTwoPTwoCSpatialParams));

// Set the problem property
SET_TYPE_PROP(TwoPTwoCSubProblem, Problem, TwoPTwoCSubProblem<TTAG(TwoPTwoCSubProblem)>);

// Set the local residual extended for the coupling
SET_TYPE_PROP(TwoPTwoCSubProblem, LocalResidual, TwoPTwoCCouplingLocalResidual<TypeTag>);

// choose pn and Sw as primary variables
SET_INT_PROP(TwoPTwoCSubProblem, Formulation, TwoPTwoCFormulation::pnsw);


// the gas component balance (air) is replaced by the total mass balance
SET_PROP(TwoPTwoCSubProblem, ReplaceCompEqIdx)
{
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    static const int value = Indices::contiNEqIdx;
};

// the fluidsystem is set in the coupled problem
SET_PROP(TwoPTwoCSubProblem, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag) CoupledTypeTag;
    typedef typename GET_PROP_TYPE(CoupledTypeTag, FluidSystem) FluidSystem;
public:
    typedef FluidSystem type;
};

// use formulation based on mass fractions
SET_BOOL_PROP(TwoPTwoCSubProblem, UseMoles, false);

// enable/disable velocity output
SET_BOOL_PROP(TwoPTwoCSubProblem, VtkAddVelocity, true);

// Enable gravity
SET_BOOL_PROP(TwoPTwoCSubProblem, ProblemEnableGravity, true);
}

/*!
 * \ingroup ImplicitTestProblems
 * \ingroup MultidomainProblems
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
        try
        {
            Scalar noDarcyX = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, NoDarcyX);
            Scalar xMin = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, XMin);

            bboxMin_[0] = std::max(xMin,noDarcyX);
            bboxMax_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, XMax);

            bboxMin_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, YMin);
            bboxMax_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePos);

            runUpDistanceX_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, RunUpDistanceX); // first part of the interface without coupling
            xMaterialInterface_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, MaterialInterfaceX);
            initializationTime_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, InitTime);

            refTemperature_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, RefTemperaturePM);
            refTemperatureFF_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefTemperature);
            refPressure_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, RefPressurePM);
            initialSw1_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, InitialSw1);
            initialSw2_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, InitialSw2);

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
        catch (Dumux::ParameterException &e) {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }
        catch (...) {
             std::cerr << "Unknown exception thrown!\n";
             exit(1);
        }
    }

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
    { return GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Vtk, NamePM); }

    /*!
     * \brief Called by the Dumux::TimeManager in order to
     *        initialize the problem.
     *
     * If you overload this method don't forget to call
     * ParentType::init()
     */
    void init()
    {
        // set the initial condition of the model
        ParentType::init();

        this->model().globalStorage(storageLastTimestep_);

        const int blModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, FreeFlow, UseBoundaryLayerModel);
        std::cout << "Using boundary layer model " << blModel << std::endl;
    }

    /*!
     * \brief Returns the temperature \f$ K \f$
     *
     * \param globalPos The global position
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return refTemperature_;
    };

    // \}

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given vertex
     *
     * \param values Stores the value of the boundary type
     * \param vertex The vertex
     */
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();
        Scalar time = this->timeManager().time();

        values.setAllNeumann();

        if (onUpperBoundary_(globalPos))
        {
            if (time > initializationTime_)
            {
                if (globalPos[0] > runUpDistanceX_ - eps_)
                    values.setAllCouplingInflow();
                else
                    values.setAllNeumann();
            }
            else
            {
//                values.setAllDirichlet();
                values.setAllNeumann();
            }
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The Dirichlet values for the primary variables
     * \param vertex The vertex representing the "half volume on the boundary"
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        initial_(values, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * \param values The Neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^{\textrm{dim}-1} \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param is The intersection between element and boundary
     * \param scvIdx The local subcontrolvolume index
     * \param boundaryFaceIdx The index of the boundary face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &is,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        values = Scalar(0);
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param values The source and sink values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^\textrm{dim} \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scvIdx The local subcontrolvolume index
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int scvIdx) const
    {
        values = Scalar(0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scvIdx The local subcontrolvolume index
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const int scvIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        values = Scalar(0);

        initial_(values, globalPos);
    }

    /*!
     * \brief Return the initial phase state inside a control volume.
     *
     * \param vertex The vertex
     * \param dofIdxGlobal The global index of the degree of freedom
     * \param globalPos The global position
     */
    int initialPhasePresence(const Vertex &vertex,
                             const int &dofIdxGlobal,
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
                storageChange /= (bboxMax_[0]-bboxMin_[0]);

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
        return ((onLeftBoundary_(globalPos) && onLowerBoundary_(globalPos)) ||
            (onLeftBoundary_(globalPos) && onUpperBoundary_(globalPos)) ||
            (onRightBoundary_(globalPos) && onLowerBoundary_(globalPos)) ||
            (onRightBoundary_(globalPos) && onUpperBoundary_(globalPos)));
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

    // \}
private:
    /*!
     * \brief Internal method for the initial condition
     *        (reused for the dirichlet conditions!)
     */
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        values[pressureIdx] = refPressure_
                + 1000.*this->gravity()[1]*(globalPos[1]-bboxMax_[1]);
        if (globalPos[0] < xMaterialInterface_)
            values[switchIdx] = initialSw1_;
        else
            values[switchIdx] = initialSw2_;

//        if (onUpperBoundary_(globalPos) && onLeftBoundary_(globalPos))
//        {
            //            values[switchIdx] = 0.99e-4;
//        }

        //        values[pressureIdx] = 1e5*(2.0 - globalPos[0]);
        //        if (onLeftBoundary_(globalPos))
        //            values[switchIdx] = 0.9e-4;
        //        else
        //            values[switchIdx] = 1e-4;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < bboxMin_[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > bboxMax_[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < bboxMin_[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > bboxMax_[1] - eps_; }

    bool onBoundary_(const GlobalPosition &globalPos) const
    {
        return (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)
                || onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos));
    }

    static constexpr Scalar eps_ = 1e-8;
    GlobalPosition bboxMin_;
    GlobalPosition bboxMax_;
    Scalar xMaterialInterface_;

    int freqMassOutput_;

    PrimaryVariables storageLastTimestep_;
    Scalar lastMassOutputTime_;

    Scalar refTemperature_;
    Scalar refTemperatureFF_;
    Scalar refPressure_;
    Scalar initialSw1_;
    Scalar initialSw2_;

    Scalar runUpDistanceX_;
    Scalar initializationTime_;
    std::ofstream outfile;    
};
} //end namespace Dumux

#endif // DUMUX_2P2C_SUBPROBLEM_HH
