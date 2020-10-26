// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup SequentialTwoPTwoCModel
 * \brief Base class for compositional pressure equations
 */
#ifndef DUMUX_FVPRESSURECOMPOSITIONAL_HH
#define DUMUX_FVPRESSURECOMPOSITIONAL_HH

#include <cmath>
#include <dune/common/float_cmp.hh>

// dumux environment
#include <dumux/common/math.hh>
#include <dumux/porousmediumflow/sequential/cellcentered/pressure.hh>
#include <dumux/material/constraintsolvers/compositionalflash.hh>
#include <dumux/porousmediumflow/2p2c/sequential/properties.hh>
#include <dumux/io/vtkmultiwriter.hh>

#include <dumux/common/deprecated.hh>

namespace Dumux {
/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief The finite volume model for the solution of the compositional pressure equation
 *
 * Provides the common ground to solve compositional pressure equations of the form
 * \f[
        c_{total}\frac{\partial p}{\partial t} + \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}}
        \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{\alpha} \bf{v}_{\alpha}\right)
         = \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}} q^{\kappa},
 * \f]
 * where \f$\bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 * \f$ c_{total} \f$ represents the total compressibility, for constant porosity this yields
 * \f$ - \frac{\partial V_{total}}{\partial p_{\alpha}} \f$,
 * \f$p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability,
 * \f$ \lambda_{\alpha} \f$ the phase mobility,
 * \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and
 * \f$ C^{\kappa} \f$ the total Component concentration.
 * See paper SPE 99619 or "Analysis of a Compositional Model for Fluid
 * Flow in Porous Media" by Chen, Qin and Ewing for derivation.
 *
 * Common functions such as output and the initialization procedure are provided here. Also,
 * private vector (the update estimate for the volume derivatives) are stored in this class,
 * as only derived classes (other compositional pressure models) need acess to it.
 * The partial derivatives of the actual fluid volume \f$ v_{total} \f$ are gained by using a secant method.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FVPressureCompositional
: public FVPressure<TypeTag>
{
    //the model implementation
    using Implementation = GetPropType<TypeTag, Properties::PressureModel>;
    using ParentType = FVPressure<TypeTag>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransportSolutionType = GetPropType<TypeTag, Properties::TransportSolutionType>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using CellData = GetPropType<TypeTag, Properties::CellData>;
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureN
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wCompIdx, nCompIdx = Indices::nCompIdx,
        contiWEqIdx = Indices::contiWEqIdx, contiNEqIdx = Indices::contiNEqIdx
    };
    enum
    {
        numPhases = getPropValue<TypeTag, Properties::NumPhases>(),
        numComponents = getPropValue<TypeTag, Properties::NumComponents>()
    };

    // using declarations to abbreviate a dune class...
    using Element = typename GridView::Traits::template Codim<0>::Entity;

    // convenience shortcuts for Vectors/Matrices
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using PhaseVector = Dune::FieldVector<Scalar, numPhases>;
    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;

public:
    //the variables object is initialized, non-compositional before and compositional after first pressure calculation
    void initialMaterialLaws(bool compositional);

    //initialization routine to prepare first timestep
    void initialize(bool solveTwice = false);

    /*!
     * \brief Compositional pressure solution routine: update estimate for secants, assemble, solve.
     *
     * An update estime (transport step acoording to old pressure field) determines changes in
     * mass, composition, which is used to calculate volume derivatives entering the pressure
     * equation, as well as an approximate guess for time step size for the storage terms in the
     * p.e.
     * Afterwards, the system is assembled and solved for pressure.
     */
    void update()
    {
        //pre-transport to estimate update vector
        Scalar dt_estimate = 0.;
        Dune::dinfo << "secant guess"<< std::endl;
        problem_.transportModel().update(problem_.timeManager().time(), dt_estimate, updateEstimate_, false);
        //last argument false in update() makes sure that this is estimate and no "real" transport step

        // if we just started a new episode, the TS size of the update Estimate is a better
        // estimate then the size of the last time step
        if(Dune::FloatCmp::eq<Scalar>(problem_.timeManager().time(), problem_.timeManager().episodeStartTime())
                && problem_.timeManager().episodeIndex() > 1)
            problem_.timeManager().setTimeStepSize(dt_estimate*getParam<Scalar>("Impet.CFLFactor"));

        updateEstimate_ *= problem_.timeManager().timeStepSize();

        problem_.pressureModel().assemble(false);           Dune::dinfo << "pressure calculation"<< std::endl;
        problem_.pressureModel().solve();

        return;
    }

    //constitutive functions are initialized and stored in the variables object
    void updateMaterialLaws(bool postTimeStep = false);

    //numerical volume derivatives wrt changes in mass, pressure
    void volumeDerivatives(const GlobalPosition&,const Element& ep);

    /*! \name general methods for output */
    //@{
    /*!
     * \brief Write data files
     *
     * Adds pressure-related quantities, including numerical things such as the volume Error
     * entering the pressure equation. Verobosity of the output can be triggered by
     * the property / parameter VtkOutputLevel, with 0 putting out only primary
     * variables and 4 being very verbose.
     * \tparam MultiWriter Class defining the output writer
     * \param writer  The output writer (usually a VTKMultiWriter object)
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        using ScalarSolutionType = typename GetProp<TypeTag, Properties::SolutionTypes>::ScalarSolution;
        int size = problem_.gridView().size(0);
        ScalarSolutionType *pressureW = writer.allocateManagedBuffer(size);
        ScalarSolutionType *pressureN = writer.allocateManagedBuffer(size);
        ScalarSolutionType *totalConcentration1 = writer.allocateManagedBuffer (size);
        ScalarSolutionType *totalConcentration2 = writer.allocateManagedBuffer (size);
        ScalarSolutionType *pc = writer.allocateManagedBuffer(size);
        ScalarSolutionType *saturationW = writer.allocateManagedBuffer(size);
        ScalarSolutionType *densityWetting = writer.allocateManagedBuffer(size);
        ScalarSolutionType *densityNonwetting = writer.allocateManagedBuffer(size);
        ScalarSolutionType *mobilityW = writer.allocateManagedBuffer (size);
        ScalarSolutionType *mobilityNW = writer.allocateManagedBuffer (size);
        ScalarSolutionType *massfraction1W = writer.allocateManagedBuffer (size);
        ScalarSolutionType *massfraction1NW = writer.allocateManagedBuffer (size);
        // numerical stuff
        ScalarSolutionType *volErr = writer.allocateManagedBuffer (size);


        for (int i = 0; i < size; i++)
        {
            // basic level 0 output
            CellData& cellData = problem_.variables().cellData(i);
            (*pressureW)[i] = cellData.pressure(wPhaseIdx);
            (*pressureN)[i] = cellData.pressure(nPhaseIdx);
            (*totalConcentration1)[i] = cellData.massConcentration(wCompIdx);
            (*totalConcentration2)[i] = cellData.massConcentration(nCompIdx);
            (*saturationW)[i] = cellData.saturation(wPhaseIdx);
            // output standard secondary variables
            if(problem_.vtkOutputLevel()>=1)
            {
                (*pc)[i] = cellData.capillaryPressure();
                (*densityWetting)[i] = cellData.density(wPhaseIdx);
                (*densityNonwetting)[i] = cellData.density(nPhaseIdx);
                (*mobilityW)[i] = cellData.mobility(wPhaseIdx);
                (*mobilityNW)[i] = cellData.mobility(nPhaseIdx);
                (*massfraction1W)[i] = cellData.massFraction(wPhaseIdx,wCompIdx);
                (*massfraction1NW)[i] = cellData.massFraction(nPhaseIdx,wCompIdx);
                (*volErr)[i] = cellData.volumeError();
            }
        }
        writer.attachCellData(*pressureW, "wetting pressure");
        writer.attachCellData(*pressureN, "nonwetting pressure");
        writer.attachCellData(*saturationW, "wetting saturation");
        writer.attachCellData(*totalConcentration1, "C^w from cellData");
        writer.attachCellData(*totalConcentration2, "C^n from cellData");
        if(problem_.vtkOutputLevel()>=1)
        {
            writer.attachCellData(*pc, "capillary pressure");
            writer.attachCellData(*densityWetting, "wetting density");
            writer.attachCellData(*densityNonwetting, "nonwetting density");
            writer.attachCellData(*mobilityW, "mobility w_phase");
            writer.attachCellData(*mobilityNW, "mobility nw_phase");
            std::ostringstream oss1, oss2;
            oss1 << "mass fraction " << FluidSystem::componentName(0) << " in " << FluidSystem::phaseName(0) << "-phase";
            writer.attachCellData(*massfraction1W, oss1.str());
            oss2 << "mass fraction " << FluidSystem::componentName(0) << " in " << FluidSystem::phaseName(1) << "-phase";
            writer.attachCellData(*massfraction1NW, oss2.str());
            writer.attachCellData(*volErr, "volume Error");
        }
        // verbose output with numerical details
        if(problem_.vtkOutputLevel()>=2)
        {
            // add debug stuff
            ScalarSolutionType *errorCorrPtr = writer.allocateManagedBuffer (size);
            ScalarSolutionType *dv_dpPtr = writer.allocateManagedBuffer (size);
            ScalarSolutionType *dV_dC1Ptr = writer.allocateManagedBuffer (size);
            ScalarSolutionType *dV_dC2Ptr = writer.allocateManagedBuffer (size);
            ScalarSolutionType *updEstimate1 = writer.allocateManagedBuffer (size);
            ScalarSolutionType *updEstimate2 = writer.allocateManagedBuffer (size);
            for (int i = 0; i < size; i++)
            {
                CellData& cellData = problem_.variables().cellData(i);
                (*errorCorrPtr)[i] = cellData.errorCorrection();
                (*dv_dpPtr)[i] = cellData.dv_dp();
                (*dV_dC1Ptr)[i] = cellData.dv(wCompIdx);
                (*dV_dC2Ptr)[i] = cellData.dv(nCompIdx);
                (*updEstimate1)[i] = updateEstimate_[0][i];
                (*updEstimate2)[i] = updateEstimate_[1][i];
            }
            writer.attachCellData(*errorCorrPtr, "Error Correction");
            writer.attachCellData(*dv_dpPtr, "dv_dp");
            writer.attachCellData(*dV_dC1Ptr, "dV_dC1");
            writer.attachCellData(*dV_dC2Ptr, "dV_dC2");
            writer.attachCellData(*updEstimate1, "updEstimate comp 1");
            writer.attachCellData(*updEstimate2, "updEstimate comp 2");
        }

        // very verbose output
        if(problem_.vtkOutputLevel()>=3)
        {
            ScalarSolutionType *pressurePV = writer.allocateManagedBuffer(size);
            ScalarSolutionType *viscosityWetting = writer.allocateManagedBuffer(size);
            ScalarSolutionType *viscosityNonwetting = writer.allocateManagedBuffer(size);
            // ScalarSolutionType *nun = writer.allocateManagedBuffer(size);
            // ScalarSolutionType *nuw = writer.allocateManagedBuffer(size);
            ScalarSolutionType *faceUpwindW = writer.allocateManagedBuffer(size);
            ScalarSolutionType *faceUpwindN = writer.allocateManagedBuffer(size);
            for (int i = 0; i < size; i++)
            {
                CellData& cellData = problem_.variables().cellData(i);
                (*viscosityWetting)[i] = cellData.viscosity(wPhaseIdx);
                (*viscosityNonwetting)[i] = cellData.viscosity(nPhaseIdx);
                // (*nun)[i] = cellData.phaseMassFraction(nPhaseIdx);
                // (*nuw)[i] = cellData.phaseMassFraction(wPhaseIdx);
                (*faceUpwindW)[i] = 0;
                (*faceUpwindN)[i] = 0;
                // run thorugh all local face idx and collect upwind information
                for(int fIdx = 0; fIdx<cellData.fluxData().size(); fIdx++)
                {
                    if(cellData.isUpwindCell(fIdx, contiWEqIdx))
                        (*faceUpwindW)[i] += pow(10,static_cast<double>(3-fIdx));
                    if(cellData.isUpwindCell(fIdx, contiNEqIdx))
                        (*faceUpwindN)[i] += pow(10,static_cast<double>(3-fIdx));
                }
            }
            //  writer.attachCellData(*nun, "phase mass fraction n-phase");
            //   writer.attachCellData(*nuw, "phase mass fraction w-phase");
            *pressurePV = this->pressure();
            writer.attachCellData(*faceUpwindW, "isUpwind w-phase");
            writer.attachCellData(*faceUpwindN, "isUpwind n-phase");
            writer.attachCellData(*pressurePV, "pressure (Primary Variable");
            writer.attachCellData(*viscosityWetting, "wetting viscosity");
            writer.attachCellData(*viscosityNonwetting, "nonwetting viscosity");
        }
        return;
    }

    /*!
     * \brief Write additional debug info in a special writer.
     *
     * To visualize the different steps through the initialization procedure,
     * we use very small pseudo time steps only for the writer!
     * This is only for debugging of the initialization procedure.
     * \param pseudoTS Time steps that only appear in the writer, not real.
     */
    void initializationOutput(double pseudoTS = 0.)
    {
        std::cout << "Writing debug for current time step\n";
        initializationOutputWriter_.beginWrite(problem_.timeManager().time() + pseudoTS);
        asImp_().addOutputVtkFields(initializationOutputWriter_);

#if DUNE_MINIMAL_DEBUG_LEVEL <= 2
        int size_ = problem_.gridView().size(0);
        // output porosity, permeability
        Dune::BlockVector<Dune::FieldVector<double,1> > *poroPtr = initializationOutputWriter_.allocateManagedBuffer (size_);
        Dune::BlockVector<Dune::FieldVector<double,1> > *permPtr = initializationOutputWriter_.allocateManagedBuffer (size_);

        Dune::BlockVector<Dune::FieldVector<double,1> > poro_(0.), perm_(0.);
        poro_.resize(size_); perm_.resize(size_);
        // iterate over all elements of domain
        for (const auto& element : elements(problem_.gridView()))
        {
            // get index
            int eIdxGlobal = problem_.variables().index(element);
            poro_[eIdxGlobal] = problem_.spatialParams().porosity(element);
            perm_[eIdxGlobal] = problem_.spatialParams().intrinsicPermeability(element)[0][0];
        }
        *poroPtr = poro_;
        *permPtr = perm_;
        initializationOutputWriter_.attachCellData(*poroPtr, "porosity");
        initializationOutputWriter_.attachCellData(*permPtr, "permeability");

        if(problem_.vtkOutputLevel()>=3)
        {
            Dune::BlockVector<Dune::FieldVector<double,1> > *permPtrY = initializationOutputWriter_.allocateManagedBuffer (size_);
            Dune::BlockVector<Dune::FieldVector<double,1> > *permPtrZ = initializationOutputWriter_.allocateManagedBuffer (size_);
            Dune::BlockVector<Dune::FieldVector<double,1> > permY_(0.), permZ_(0.);
            permY_.resize(size_); permZ_.resize(size_);
            // iterate over all elements of domain
            for (const auto& element : elements(problem_.gridView()))
            {
                // get index
                int eIdxGlobal = problem_.variables().index(element);
                if(dim >=2)
                    permY_[eIdxGlobal] = problem_.spatialParams().intrinsicPermeability(element)[1][1];
                if(dim >=3)
                    permZ_[eIdxGlobal] = problem_.spatialParams().intrinsicPermeability(element)[2][2];
            }
            if(dim >=2)
            {
                *permPtrY = permY_;
                initializationOutputWriter_.attachCellData(*permPtrY, "permeability Y");
            }
            if(dim >=3)
            {
                *permPtrZ = permZ_;
                initializationOutputWriter_.attachCellData(*permPtrZ, "permeability Z");
            }
        }
#endif

        initializationOutputWriter_.endWrite();
        return;
    }
    //@}

    /*!
     * \brief Constructs a FVPressureCompositional object
     * \param problem a problem class object
     */
    FVPressureCompositional(Problem& problem) : FVPressure<TypeTag>(problem),
        problem_(problem), initializationOutputWriter_(problem.gridView(),"initOutput2p2c"),
        maxError_(0.0), incp_(1.0e1)
    {
        updateEstimate_.resize(getPropValue<TypeTag, Properties::NumPhases>());
        for  (int i=0; i<getPropValue<TypeTag, Properties::NumPhases>(); i++)
            updateEstimate_[i].resize(problem.gridView().size(0));

        ErrorTermFactor_ = getParam<Scalar>("Impet.ErrorTermFactor");
        ErrorTermLowerBound_ = getParam<Scalar>("Impet.ErrorTermLowerBound", 0.2);
        ErrorTermUpperBound_ = getParam<Scalar>("Impet.ErrorTermUpperBound");

        if (pressureType == Indices::pressureGlobal)
        {
            DUNE_THROW(Dune::NotImplemented, "Global Pressure type not supported!");
        }
    }
protected:
    TransportSolutionType updateEstimate_; //!< Update estimate for changes in volume for the pressure equation
    Problem& problem_;

    //! output for the initialization procedure
    VtkMultiWriter<GridView> initializationOutputWriter_;

    Scalar maxError_; //!< Maximum volume error of all cells
    Scalar incp_; //!< Increment for the volume derivative w.r.t pressure
    Scalar ErrorTermFactor_; //!< Handling of error term: relaxation factor
    Scalar ErrorTermLowerBound_; //!< Handling of error term: lower bound for error dampening
    Scalar ErrorTermUpperBound_; //!< Handling of error term: upper bound for error dampening
    //! gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
    static constexpr int pressureType = getPropValue<TypeTag, Properties::PressureFormulation>();

private:
    Problem& problem()
    {
        return problem_;
    }
    const Problem& problem() const
    {
        return problem_;
    }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    {   return *static_cast<Implementation *>(this);}

    //! \copydoc IMPETProblem::asImp_()
    const Implementation &asImp_() const
    {   return *static_cast<const Implementation *>(this);}
};


/*!
 * \brief Initializes the simulation run
 *
 * Initializes the simulation to gain the initial pressure field.
 * Output throughout initialization procedure is only done in debug mode.
 *
 * \param solveTwice flag to determine possible iterations of the initialization process
 */
template<class TypeTag>
void FVPressureCompositional<TypeTag>::initialize(bool solveTwice)
{
    // inizialize matrix, RHS, and condition Matrix
    ParentType::initialize();

    // initialguess: set saturations, determine visco and mobility for initial pressure equation
    // at this moment, the pressure is unknown. Hence, dont regard compositional effects.
    Dune::dinfo << "first saturation guess"<<std::endl; //=J: initialGuess()
        problem_.pressureModel().initialMaterialLaws(false);
            #if DUNE_MINIMAL_DEBUG_LEVEL <= 3
                // this update produces fluxes on init pressures
                Scalar dummy;
                problem_.transportModel().update(0.,dummy, updateEstimate_, false);
                initializationOutput();
            #endif
    Dune::dinfo << "first pressure guess"<<std::endl;
        problem_.pressureModel().assemble(true);
        problem_.pressureModel().solve();
            #if DUNE_MINIMAL_DEBUG_LEVEL <= 3
                initializationOutput(1e-6);
            #endif
    // update the compositional variables (hence true)
    Dune::dinfo << "first guess for mass fractions"<<std::endl;//= J: transportInitial()
        problem_.pressureModel().initialMaterialLaws(true);

    // perform concentration update to determine secants
    Dune::dinfo << "secant guess"<< std::endl;
        Scalar dt_estimate = 0.;
        problem_.transportModel().update(0., dt_estimate, updateEstimate_, false);
        using std::min;
        dt_estimate = min ( problem_.timeManager().timeStepSize(), dt_estimate);
        //make sure the right time-step is used by all processes in the parallel case
        if (problem_.gridView().comm().size() > 1)
            dt_estimate = problem_.gridView().comm().min(dt_estimate);

        updateEstimate_ *= dt_estimate;
        //communicate in parallel case
//        problem_.variables().communicateUpdateEstimate();
            #if DUNE_MINIMAL_DEBUG_LEVEL <= 3
                initializationOutput(2e-6);
            #endif
    // pressure calculation
    Dune::dinfo << "second pressure guess"<<std::endl;
        problem_.pressureModel().assemble(false);
        problem_.pressureModel().solve();
            #if DUNE_MINIMAL_DEBUG_LEVEL <= 3
                initializationOutput(3e-6);
            #endif

        // update the compositional variables
        problem_.pressureModel().initialMaterialLaws(true);


    if (solveTwice)
    {
        Dune::BlockVector<Dune::FieldVector<Scalar, 1> > pressureOld(this->pressure());
        Dune::BlockVector<Dune::FieldVector<Scalar, 1> > pressureDiff;
        Scalar pressureNorm = 1.;   //dummy initialization to perform at least 1 iteration
        int numIter = 1;

        while (pressureNorm > 1e-5 && numIter < 10)
        {
            Scalar dt_dummy=0.;    // without this dummy, it will never converge!
            // update for secants
            problem_.transportModel().update(0., dt_dummy, updateEstimate_, false);   Dune::dinfo << "secant guess"<< std::endl;
            updateEstimate_ *= dt_estimate;

            //communicate in parallel case
//            problem_.variables().communicateUpdateEstimate();

            // pressure calculation
            problem_.pressureModel().assemble(false);             Dune::dinfo << "pressure guess number "<< numIter <<std::endl;
            problem_.pressureModel().solve();
            // update the compositional variables
            problem_.pressureModel().initialMaterialLaws(true);

            pressureDiff = pressureOld;
            pressureDiff -= this->pressure();
            pressureNorm = pressureDiff.infinity_norm();
            pressureOld = this->pressure();
            pressureNorm /= pressureOld.infinity_norm();

            numIter++;
        }
    }
    return;
}

/*!
 *  \brief initializes the fluid distribution and hereby the variables container
 *
 *  It differs from updateMaterialLaws() because there are two possible initial conditions:
 *  saturations and concentration.
 *  \param compositional flag that determines if compositional effects are regarded, i.e.
 *      a reasonable pressure field is known with which compositions can be calculated.
 */
template<class TypeTag>
void FVPressureCompositional<TypeTag>::initialMaterialLaws(bool compositional)
{
    // iterate through leaf grid an evaluate c0 at cell center
    for (const auto& element : elements(problem_.gridView()))
    {
        // get global coordinate of cell center
        GlobalPosition globalPos = element.geometry().center();

        // assign an Index for convenience
        int eIdxGlobal = problem_.variables().index(element);

        // get the temperature
        Scalar temperature_ = problem_.temperatureAtPos(globalPos);
        CellData& cellData = problem_.variables().cellData(eIdxGlobal);
        // acess the fluid state and prepare for manipulation
        FluidState& fluidState = cellData.manipulateFluidState();
        CompositionalFlash<Scalar, FluidSystem> flashSolver;

        // old material law interface is deprecated: Replace this by
        // const auto& fluidMatrixInteraction = spatialParams.fluidMatrixInteractionAtPos(element.geometry().center());
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem_.spatialParams(), element);

        // initial conditions
        PhaseVector pressure(0.);
        Scalar sat_0=0.;

        typename Indices::BoundaryFormulation icFormulation;
        problem_.initialFormulation(icFormulation, element);            // get type of initial condition

        if(!compositional) //means that we do the first approximate guess without compositions
        {
            // phase pressures are unknown, so start with an exemplary
            Scalar exemplaryPressure = problem_.referencePressure(element);
            pressure[wPhaseIdx] = pressure[nPhaseIdx] = this->pressure()[eIdxGlobal] = exemplaryPressure;
            if (icFormulation == Indices::saturation)  // saturation initial condition
            {
                sat_0 = problem_.initSat(element);
                flashSolver.saturationFlash2p2c(fluidState, sat_0, pressure, temperature_);
            }
            else if (icFormulation == Indices::concentration) // concentration initial condition
            {
                Scalar Z0 = problem_.initConcentration(element);
                flashSolver.concentrationFlash2p2c(fluidState, Z0, pressure, temperature_);
            }
        }
        else if(compositional)    //means we regard compositional effects since we know an estimate pressure field
        {
            if (icFormulation == Indices::saturation)  // saturation initial condition
            {
                //get saturation, determine pc
                sat_0 = problem_.initSat(element);
                Scalar pc=0.;
                if(getPropValue<TypeTag, Properties::EnableCapillarity>())
                {
                    pc = fluidMatrixInteraction.pc(sat_0);
                }
                else
                    pc = 0.;

                //determine phase pressures from primary pressure variable
                switch (pressureType)
                {
                    case pw:
                    {
                        pressure[wPhaseIdx] = this->pressure()[eIdxGlobal];
                        pressure[nPhaseIdx] = this->pressure()[eIdxGlobal] +pc;
                        break;
                    }
                    case pn:
                    {
                        pressure[wPhaseIdx] = this->pressure()[eIdxGlobal]-pc;
                        pressure[nPhaseIdx] = this->pressure()[eIdxGlobal];
                        break;
                    }
                }
                flashSolver.saturationFlash2p2c(fluidState, sat_0, pressure, temperature_);
            }
            else if (icFormulation == Indices::concentration) // concentration initial condition
            {
                Scalar Z0 = problem_.initConcentration(element);
                // If total concentrations are given at the boundary, saturation is unknown.
                // This may affect pc and hence p_alpha and hence again saturation -> iteration.

                // iterations in case of enabled capillary pressure
                if(getPropValue<TypeTag, Properties::EnableCapillarity>())
                {
                    //start with pc from last TS
                    Scalar pc(cellData.capillaryPressure());

                    int maxiter = 3;
                    //start iteration loop
                    for(int iter=0; iter < maxiter; iter++)
                    {
                        //determine phase pressures from primary pressure variable
                        switch (pressureType)
                        {
                            case pw:
                            {
                                pressure[wPhaseIdx] = this->pressure()[eIdxGlobal];
                                pressure[nPhaseIdx] = this->pressure()[eIdxGlobal] + pc;
                                break;
                            }
                            case pn:
                            {
                                pressure[wPhaseIdx] = this->pressure()[eIdxGlobal] - pc;
                                pressure[nPhaseIdx] = this->pressure()[eIdxGlobal];
                                break;
                            }
                        }

                        //store old pc
                        Scalar oldPc = pc;
                        //update with better pressures
                        flashSolver.concentrationFlash2p2c(fluidState, Z0, pressure, problem_.temperatureAtPos(globalPos));
                        pc = fluidMatrixInteraction.pc(fluidState.saturation(wPhaseIdx));
                        // TODO: get right criterion, do output for evaluation
                        //converge criterion
                        using std::abs;
                        if (abs(oldPc - pc) < 10.0)
                            iter = maxiter;

                        pc = fluidMatrixInteraction.pc(fluidState.saturation(wPhaseIdx));
                    }
                }
                else  // capillary pressure neglected
                {
                    pressure[wPhaseIdx] = pressure[nPhaseIdx]
                        = this->pressure()[eIdxGlobal];
                    flashSolver.concentrationFlash2p2c(fluidState, Z0, pressure, temperature_);
                }
            } //end conc initial condition
        } //end compositional

        cellData.calculateMassConcentration(problem_.spatialParams().porosity(element));

        problem_.transportModel().totalConcentration(wCompIdx,eIdxGlobal) = cellData.massConcentration(wCompIdx);
        problem_.transportModel().totalConcentration(nCompIdx,eIdxGlobal) = cellData.massConcentration(nCompIdx);

        // initialize mobilities
        cellData.setMobility(wPhaseIdx, fluidMatrixInteraction.krw(fluidState.saturation(wPhaseIdx))
                    / cellData.viscosity(wPhaseIdx));
        cellData.setMobility(nPhaseIdx, fluidMatrixInteraction.krn(fluidState.saturation(wPhaseIdx))
                    / cellData.viscosity(nPhaseIdx));

        // calculate perimeter used as weighting factor
        if(!compositional)
        {
            // run through all intersections with neighbors
            for (const auto& intersection : intersections(problem_.gridView(), element))
            {
                cellData.perimeter()
                        += intersection.geometry().volume();
            }
            cellData.globalIdx() = eIdxGlobal;

            // set dv to zero to prevent output errors
            cellData.dv_dp() = 0.;
            cellData.dv(wPhaseIdx) = 0.;
            cellData.dv(nPhaseIdx) = 0.;
        }
        // all
        cellData.reset();
    }
    return;
}

/*!
 * \brief Updates secondary variables
 *
 * A loop through all elements updates the secondary variables stored in the variableclass
 * by using the updated primary variables.
 * \param postTimeStep Flag indicating method is called from Problem::postTimeStep()
 */
template<class TypeTag>
void FVPressureCompositional<TypeTag>::updateMaterialLaws(bool postTimeStep)
{
    Scalar maxError = 0.;
    // iterate through leaf grid an evaluate c0 at cell center
    for (const auto& element : elements(problem().gridView()))
    {
        int eIdxGlobal = problem().variables().index(element);

        CellData& cellData = problem().variables().cellData(eIdxGlobal);

        asImp_().updateMaterialLawsInElement(element, postTimeStep);

        using std::max;
        maxError = max(maxError, fabs(cellData.volumeError()));
    }
    if (problem_.gridView().comm().size() > 1)
        maxError_ = problem_.gridView().comm().max(maxError_);

    maxError_ = maxError/problem().timeManager().timeStepSize();
    return;
}

/*!
 * \brief Partial derivatives of the volumes w.r.t. changes in total concentration and pressure
 *
 * This method calculates the volume derivatives via a secant method, where the
 * secants are gained in a pre-computational step via the transport equation and
 * the last TS size.
 * The partial derivatives w.r.t. mass are defined as
 * \f$ \frac{\partial v}{\partial C^{\kappa}} = \frac{\partial V}{\partial m^{\kappa}}\f$
 *
 * \param globalPos The global position of the current element
 * \param element The current element
 */
template<class TypeTag>
void FVPressureCompositional<TypeTag>::volumeDerivatives(const GlobalPosition& globalPos, const Element& element)
{
    // cell index
    int eIdxGlobal = problem_.variables().index(element);

    CellData& cellData = problem_.variables().cellData(eIdxGlobal);

    // get cell temperature
    Scalar temperature_ = cellData.temperature(wPhaseIdx);

    // initialize an Fluid state and a flash solver
    FluidState updFluidState;
    CompositionalFlash<Scalar, FluidSystem> flashSolver;

    /**********************************
     * a) get necessary variables
     **********************************/
    //determine phase pressures for flash calculation
    PhaseVector pressure(0.);
    for(int phaseIdx = 0; phaseIdx< numPhases; phaseIdx++)
            pressure[phaseIdx] = cellData.pressure(phaseIdx);

    // mass of components inside the cell
    ComponentVector mass(0.);
    for(int compIdx = 0; compIdx< numComponents; compIdx++)
        mass[compIdx] = cellData.massConcentration(compIdx);

    // actual fluid volume
    // see Fritz 2011 (Dissertation) eq.3.76
    Scalar specificVolume(0.); // = \sum_{\alpha} \nu_{\alpha} / \rho_{\alpha}
    for(int phaseIdx = 0; phaseIdx< numPhases; phaseIdx++)
        specificVolume += cellData.phaseMassFraction(phaseIdx) / cellData.density(phaseIdx);
    Scalar volalt = mass.one_norm() * specificVolume;
//    volalt = cellData.volumeError()+problem_.spatialParams().porosity(element);
        // = \sum_{\kappa} C^{\kappa} + \sum_{\alpha} \nu_{\alpha} / \rho_{\alpha}

    /**********************************
     * b) define increments
     **********************************/
    // increments for numerical derivatives
    ComponentVector massIncrement(0.);
    for(int compIdx = 0; compIdx< numComponents; compIdx++)
    {
        massIncrement[compIdx] = updateEstimate_[compIdx][eIdxGlobal];
        if(fabs(massIncrement[compIdx]) < 1e-8 * cellData.density(compIdx))
            massIncrement[compIdx] = 1e-8* cellData.density(compIdx);   // as phaseIdx = compIdx
    }
    Scalar& incp = incp_;

    /**********************************
     * c) Secant method for derivatives
     **********************************/

    // numerical derivative of fluid volume with respect to pressure
    PhaseVector p_(incp);
    p_ += pressure;
    Scalar Z0 = mass[0] / mass.one_norm();
    flashSolver.concentrationFlash2p2c(updFluidState, Z0, p_, temperature_);

    specificVolume=0.; // = \sum_{\alpha} \nu_{\alpha} / \rho_{\alpha}
    for(int phaseIdx = 0; phaseIdx< numPhases; phaseIdx++)
        specificVolume += updFluidState.phaseMassFraction(phaseIdx) / updFluidState.density(phaseIdx);
    Scalar dv_dp = ((mass.one_norm() * specificVolume) - volalt) /incp;

    if (dv_dp>0)
    {
        // dV_dp > 0 is unphysical: Try inverse increment for secant
        Dune::dinfo << "dv_dp larger 0 at Idx " << eIdxGlobal << " , try and invert secant"<< std::endl;

        p_ -= 2*incp;
        flashSolver.concentrationFlash2p2c(updFluidState, Z0, p_, temperature_);

        specificVolume=0.; // = \sum_{\alpha} \nu_{\alpha} / \rho_{\alpha}
        for(int phaseIdx = 0; phaseIdx< numPhases; phaseIdx++)
            specificVolume += updFluidState.phaseMassFraction(phaseIdx) / updFluidState.density(phaseIdx);
        dv_dp = ((mass.one_norm() * specificVolume) - volalt) /incp;

        // dV_dp > 0 is unphysical: Try inverse increment for secant
        if (dv_dp>0)
        {
            Dune::dwarn << "dv_dp still larger 0 after inverting secant at idx"<< eIdxGlobal<< std::endl;
            p_ += 2*incp;
            flashSolver.concentrationFlash2p2c(updFluidState, Z0, p_, temperature_);
            // neglect effects of phase split, only regard changes in phase densities
            specificVolume=0.; // = \sum_{\alpha} \nu_{\alpha} / \rho_{\alpha}
            for(int phaseIdx = 0; phaseIdx< numPhases; phaseIdx++)
                specificVolume += cellData.phaseMassFraction(phaseIdx) / updFluidState.density(phaseIdx);
            dv_dp = ((mass.one_norm() * specificVolume) - volalt) /incp;
            if (dv_dp>0)
            {
                std::cout << "dv_dp still larger 0 after both inverts at idx " << eIdxGlobal << std::endl;
                dv_dp = cellData.dv_dp();
            }
        }
    }
    cellData.dv_dp()=dv_dp;

    // numerical derivative of fluid volume with respect to mass of components
    for (int compIdx = 0; compIdx<numComponents; compIdx++)
    {
        mass[compIdx] +=  massIncrement[compIdx];
        Z0 = mass[0] / mass.one_norm();

        flashSolver.concentrationFlash2p2c(updFluidState, Z0, pressure, temperature_);

        specificVolume=0.; // = \sum_{\alpha} \nu_{\alpha} / \rho_{\alpha}
        for(int phaseIdx = 0; phaseIdx< numPhases; phaseIdx++)
            specificVolume += updFluidState.phaseMassFraction(phaseIdx) / updFluidState.density(phaseIdx);

        cellData.dv(compIdx) = ((mass.one_norm() * specificVolume) - volalt) / massIncrement[compIdx];
        mass[compIdx] -= massIncrement[compIdx];

        //check routines if derivatives are meaningful
        using std::isnan;
        using std::isinf;
        if (isnan(cellData.dv(compIdx)) || isinf(cellData.dv(compIdx)) )
        {
            DUNE_THROW(Dune::MathError, "NAN/inf of dV_dm. If that happens in first timestep, try smaller firstDt!");
        }
    }
    cellData.volumeDerivativesAvailable(true);
}

}//end namespace Dumux
#endif
