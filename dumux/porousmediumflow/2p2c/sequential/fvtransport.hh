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
 * \brief Finite volume discretization of the component transport equation
 */
#ifndef DUMUX_FVTRANSPORT2P2C_HH
#define DUMUX_FVTRANSPORT2P2C_HH

#include <cmath>
#include <unordered_map>

#include <dune/grid/common/gridenums.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/porousmediumflow/2p2c/sequential/properties.hh>
#include <dumux/material/constraintsolvers/compositionalflash.hh>
#include <dumux/common/math.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

#include <dumux/common/deprecated.hh>

namespace Dumux {
/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief Compositional transport step in a Finite Volume discretization.
 *
 *  The finite volume model for the solution of the transport equation for compositional
 *  two-phase flow.
 *  \f[
      \frac{\partial C^\kappa}{\partial t} = - \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha}
      \varrho_{alpha} \bf{v}_{\alpha}\right) + q^{\kappa},
 *  \f]
 *  where \f$ \bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 *  \f$ p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability, \f$ \lambda_{\alpha} \f$ the phase mobility,
 *  \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and \f$ C^{\kappa} \f$ the total Component concentration.
 *  The whole flux contribution for each cell is subdivided into a storage term, a flux term and a source term.
 *  Corresponding functions (<tt>getFlux()</tt> and <tt>getFluxOnBoundary()</tt>) are provided,
 *  internal sources are directly treated.
 *
 *  \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVTransport2P2C
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Implementation = GetPropType<TypeTag, Properties::TransportModel>;

    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using CellData = GetPropType<TypeTag, Properties::CellData>;

    using TransportSolutionType = GetPropType<TypeTag, Properties::TransportSolutionType>;

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
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx,
        contiWEqIdx=Indices::contiWEqIdx, contiNEqIdx=Indices::contiNEqIdx,
        NumPhases = getPropValue<TypeTag, Properties::NumPhases>(),
        NumComponents = getPropValue<TypeTag, Properties::NumComponents>()
    };

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;
    using PhaseVector = Dune::FieldVector<Scalar, NumPhases>;
    using ComponentVector = Dune::FieldVector<Scalar, NumComponents>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

public:
    //! @copydoc FVPressure::EntryType
    using EntryType = Dune::FieldVector<Scalar, 2>;
    using TimeStepFluxType = Dune::FieldVector<Scalar, 2>;

protected:
    struct LocalTimesteppingData
    {
        Dune::FieldVector<EntryType, 2*dim> faceFluxes;
        Dune::FieldVector<Scalar, 2*dim> faceTargetDt;
        Scalar dt;
        LocalTimesteppingData():faceFluxes(EntryType(0.0)), faceTargetDt(0.0), dt(0)
        {}
    };

    //! Acess function for the current problem
    Problem& problem()
    { return problem_; }

    void innerUpdate(TransportSolutionType& updateVec);

public:
    virtual void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet = false);

    void updateTransportedQuantity(TransportSolutionType& updateVector);

    void updateTransportedQuantity(TransportSolutionType& updateVector, Scalar dt);

    void updateConcentrations(TransportSolutionType& updateVector, Scalar dt);

    // Function which calculates the flux update
    void getFlux(ComponentVector& fluxEntries, EntryType& timestepFlux,
                 const Intersection& intersection, CellData& cellDataI);

    // Function which calculates the boundary flux update
    void getFluxOnBoundary(ComponentVector& fluxEntries, EntryType& timestepFlux,
                           const Intersection& intersection, const CellData& cellDataI);

    void evalBoundary(GlobalPosition globalPosFace,
                      const Intersection& intersection,
                      FluidState& BCfluidState,
                      PhaseVector& pressBound);


    /*!
     * \brief Set the initial values before the first pressure equation
     * This method is called before first pressure equation is solved from IMPET.
     */
    void initialize()
    {
        // resize update vector and set to zero
        int transportedQuantities = getPropValue<TypeTag, Properties::NumEq>() - 1; // NumEq - 1 pressure Eq
        totalConcentration_.resize(transportedQuantities);
        for (int eqNumber = 0; eqNumber < transportedQuantities; eqNumber++)
        {
            totalConcentration_[eqNumber].resize(problem().gridView().size(0));
            totalConcentration_[eqNumber] = 0;
        }
    }

    /*!
     * \brief Write transport variables into the output files
     * \param writer applied VTK-writer
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        if(problem().vtkOutputLevel()>3)
        {
            using ScalarSolutionType = typename GetProp<TypeTag, Properties::SolutionTypes>::ScalarSolution;
            int size = problem_.gridView().size(0);
            ScalarSolutionType *totalC1PV = writer.allocateManagedBuffer(size);
            ScalarSolutionType *totalC2PV = writer.allocateManagedBuffer(size);
            *totalC1PV = this->totalConcentration_[wCompIdx];
            *totalC2PV = this->totalConcentration_[nCompIdx];
            writer.attachCellData(*totalC1PV, "total Concentration w-Comp - vector");
            writer.attachCellData(*totalC2PV, "total Concentration n-Comp - vector");
        }
    }

    //! Function needed for restart option of the transport model: Write out
    void serializeEntity(std::ostream &outstream, const Element &element)
    {
        int eIdxGlobal = problem().variables().index(element);
        outstream << totalConcentration_[wCompIdx][eIdxGlobal]
                  << "  " << totalConcentration_[nCompIdx][eIdxGlobal];
    }
    //! Function needed for restart option of the transport model: Read in
    void deserializeEntity(std::istream &instream, const Element &element)
    {
        int eIdxGlobal = problem().variables().index(element);
        CellData& cellData = problem().variables().cellData(eIdxGlobal);
        instream >>  totalConcentration_[wCompIdx][eIdxGlobal]
                 >> totalConcentration_[nCompIdx][eIdxGlobal];
        cellData.setMassConcentration(wCompIdx, totalConcentration_[wCompIdx][eIdxGlobal]);
        cellData.setMassConcentration(nCompIdx, totalConcentration_[nCompIdx][eIdxGlobal]);
    }

    /*! \name Access functions for protected variables  */
    //@{
    /*!
     * \brief Return the vector of the transported quantity
     * For an immiscible IMPES scheme, this is the saturation. For compositional simulations, however,
     * the total concentration of all components is transported.
     * \param transportedQuantity Vector of both transported components
     */
    void getTransportedQuantity(TransportSolutionType& transportedQuantity)
    {
        // resize update vector and set to zero
        transportedQuantity.resize((getPropValue<TypeTag, Properties::NumEq>() - 1));
        for(int compIdx = 0; compIdx < (getPropValue<TypeTag, Properties::NumEq>() - 1); compIdx++)
            transportedQuantity[compIdx].resize(problem_.gridView().size(0));

        transportedQuantity = totalConcentration_;
    }

    /*! \name Access functions for protected variables  */
    //@{
    /*!
     * \brief Return the the total concentration stored in the transport vector
     * To get real cell values, do not acess this method, but rather
     * call the respective function in the cell data object.
     * \param compIdx The index of the component
     * \param eIdxGlobal The global index of the current cell.
     */
    Scalar& totalConcentration(int compIdx, int eIdxGlobal)
    {
        return totalConcentration_[compIdx][eIdxGlobal][0];
    }

    void getSource(Scalar& update, const Element& element, CellData& cellDataI)
    {}


    /*!
     * \brief Function to control the abort of the transport-sub-time-stepping
     * depending on a physical parameter range
     * \param entry Cell entries of the update vector
     */
    template<class DataEntry>
    bool inPhysicalRange(DataEntry& entry)
    {
        int numComp = getPropValue<TypeTag, Properties::NumEq>() - 1;
        for(int compIdx = 0; compIdx < numComp; compIdx++)
        {
            if (entry[compIdx] < -1.0e-6)
            {
                return false;
            }
        }
        return true;
    }

    //! Function to check if local time stepping is activated
    bool enableLocalTimeStepping()
    {
        return localTimeStepping_;
    }

    //@}
    /*!
     * \brief Constructs a FVTransport2P2C object
     * Currently, the compositional transport scheme can not be applied with a global pressure / total velocity
     * formulation.
     *
     * \param problem a problem class object
     */
    FVTransport2P2C(Problem& problem) :
        totalConcentration_(0.), problem_(problem),
        switchNormals(getParam<bool>("Impet.SwitchNormals")), accumulatedDt_(0),
        dtThreshold_(1e-6), subCFLFactor_(1.0)
    {
        restrictFluxInTransport_ = getParam<int>("Impet.RestrictFluxInTransport", 0);
        regulateBoundaryPermeability = getPropValue<TypeTag, Properties::RegulateBoundaryPermeability>();
        if(regulateBoundaryPermeability)
            minimalBoundaryPermeability = getParam<Scalar>("SpatialParams.MinBoundaryPermeability");

        Scalar cFLFactor = getParam<Scalar>("Impet.CFLFactor");
        using std::min;
        subCFLFactor_ = min(getParam<Scalar>("Impet.SubCFLFactor"), cFLFactor);
        verbosity_ = getParam<int>("TimeManager.SubTimestepVerbosity");

        localTimeStepping_ = subCFLFactor_/cFLFactor < 1.0 - dtThreshold_;

        if (localTimeStepping_)
            std::cout<<"max CFL-Number of "<<cFLFactor<<", max Sub-CFL-Number of "
                <<subCFLFactor_<<": Enable local time-stepping!" << std::endl;
    }

    virtual ~FVTransport2P2C()
    {
    }

protected:
    TransportSolutionType totalConcentration_; //!< private vector of transported primary variables
    Problem& problem_;
    bool impet_; //!< indicating if we are in an estimate (false) or real impet (true) step.
    int averagedFaces_; //!< number of faces were flux was restricted

    //! gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
    static const int pressureType = getPropValue<TypeTag, Properties::PressureFormulation>();
    //! Restriction of flux on new pressure field if direction reverses from the pressure equation
    int restrictFluxInTransport_;
    //! Enables regulation of permeability in the direction of a Dirichlet Boundary Condition
    bool regulateBoundaryPermeability;
    //! Minimal limit for the boundary permeability
    Scalar minimalBoundaryPermeability;
    bool switchNormals;
    //! Current time-interval in sub-time-stepping routine
    Scalar accumulatedDt_;
    //! Threshold for sub-time-stepping routine
    const Scalar dtThreshold_;
    //! Stores data for sub-time-stepping
    std::vector<LocalTimesteppingData> timeStepData_;

    void updatedTargetDt_(Scalar &dt);

    void resetTimeStepData_()
    {
        timeStepData_.clear();
        accumulatedDt_ = 0;
    }

    Scalar subCFLFactor_;
    bool localTimeStepping_;
    int verbosity_;


private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! @copydoc IMPETProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};


/*!
 *  \brief Calculate the update vector and determine timestep size
 *  This method calculates the update vector \f$ u \f$ of the discretized equation
 *  \f[
       C^{\kappa , new} = C^{\kappa , old} + u,
 *  \f]
 *  where \f$ u = \sum_{\gamma} \boldsymbol{v}_{\alpha} * \varrho_{\alpha} * X^{\kappa}_{\alpha} * \boldsymbol{n} * A_{\gamma} \f$,
 *  \f$ \boldsymbol{n} \f$ is the face normal and \f$ A_{\gamma} \f$ is the face area of face \f$ \gamma \f$.
 *
 *  In addition to the \a update vector, the recommended time step size \a dt is calculated
 *  employing a CFL condition.
 *
 *  \param t Current simulation time \f$\mathrm{[s]}\f$
 *  \param dt Time step size \f$\mathrm{[s]}\f$
 *  \param updateVec Update vector, or update estimate for secants, resp. Here in \f$\mathrm{[kg/m^3]}\f$
 *  \param impet Flag that determines if it is a real impet step or an update estimate for volume derivatives
 */
template<class TypeTag>
void FVTransport2P2C<TypeTag>::update(const Scalar t, Scalar& dt,
        TransportSolutionType& updateVec, bool impet)
{
    // initialize dt very large
    dt = 1E100;

    unsigned int size = problem_.gridView().size(0);
    if (localTimeStepping_)
    {
        if (this->timeStepData_.size() != size)
            this->timeStepData_.resize(size);
    }
    // store if we do update Estimate for flux functions
    impet_ = impet;
    averagedFaces_ = 0;

    // resize update vector and set to zero
    updateVec.resize(getPropValue<TypeTag, Properties::NumComponents>());
    updateVec[wCompIdx].resize(problem_.gridView().size(0));
    updateVec[nCompIdx].resize(problem_.gridView().size(0));
    updateVec[wCompIdx] = 0;
    updateVec[nCompIdx] = 0;

    // Cell which restricts time step size
    int restrictingCell = -1;

    ComponentVector entries(0.);
    EntryType timestepFlux(0.);
    // compute update vector
    for (const auto& element : elements(problem().gridView()))
    {
        // get cell infos
        int eIdxGlobalI = problem().variables().index(element);
        CellData& cellDataI = problem().variables().cellData(eIdxGlobalI);

        // some variables for time step calculation
        double sumfactorin = 0;
        double sumfactorout = 0;

        // run through all intersections with neighbors and boundary
        for (const auto& intersection : intersections(problem().gridView(), element))
        {
            int indexInInside = intersection.indexInInside();

            /****** interior face   *****************/
            if (intersection.neighbor())
                asImp_().getFlux(entries, timestepFlux, intersection, cellDataI);

            /******  Boundary Face   *****************/
            if (intersection.boundary())
                asImp_().getFluxOnBoundary(entries, timestepFlux, intersection, cellDataI);


            if (localTimeStepping_)
            {
                LocalTimesteppingData& localData = timeStepData_[eIdxGlobalI];
                if (localData.faceTargetDt[indexInInside] < accumulatedDt_ + dtThreshold_)
                {
                    localData.faceFluxes[indexInInside] = entries;
                }
            }
            else
            {
            // add to update vector
                updateVec[wCompIdx][eIdxGlobalI] += entries[wCompIdx];
                updateVec[nCompIdx][eIdxGlobalI] += entries[nCompIdx];
            }

            // for time step calculation
            sumfactorin += timestepFlux[0];
            sumfactorout += timestepFlux[1];

        }// end all intersections

        if (localTimeStepping_)
        {
            LocalTimesteppingData& localData = timeStepData_[eIdxGlobalI];
            for (int i=0; i < 2*dim; i++)
            {
                updateVec[wCompIdx][eIdxGlobalI] += localData.faceFluxes[i][wCompIdx];
                updateVec[nCompIdx][eIdxGlobalI] += localData.faceFluxes[i][nCompIdx];
            }
        }

        /***********     Handle source term     ***************/
        PrimaryVariables q(NAN);
        problem().source(q, element);
        updateVec[wCompIdx][eIdxGlobalI] += q[contiWEqIdx];
        updateVec[nCompIdx][eIdxGlobalI] += q[contiNEqIdx];

        // account for porosity in fluxes for time-step
        using std::max;
        sumfactorin = max(sumfactorin,sumfactorout)
                        / problem().spatialParams().porosity(element);

        //calculate time step
        if (localTimeStepping_)
        {
            timeStepData_[eIdxGlobalI].dt = 1./sumfactorin;
            if ( 1./sumfactorin < dt)
            {
                dt = 1./sumfactorin;
                restrictingCell= eIdxGlobalI;
            }
        }
        else
        {
            if ( 1./sumfactorin < dt)
            {
                dt = 1./sumfactorin;
                restrictingCell= eIdxGlobalI;
            }
        }
    } // end grid traversal

#if HAVE_MPI
    // communicate updated values
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using ElementMapper = typename SolutionTypes::ElementMapper;
    using DataHandle = VectorCommDataHandleEqual<ElementMapper, Dune::BlockVector<Dune::FieldVector<Scalar, 1> >, 0/*elementCodim*/>;
    for (int i = 0; i < updateVec.size(); i++)
    {
        DataHandle dataHandle(problem_.variables().elementMapper(), updateVec[i]);
        problem_.gridView().template communicate<DataHandle>(dataHandle,
                                                            Dune::InteriorBorder_All_Interface,
                                                            Dune::ForwardCommunication);
    }

    if (localTimeStepping_)
    {
        using TimeDataHandle = VectorCommDataHandleEqual<ElementMapper, std::vector<LocalTimesteppingData>, 0/*elementCodim*/>;

        TimeDataHandle timeDataHandle(problem_.elementMapper(), timeStepData_);
        problem_.gridView().template communicate<TimeDataHandle>(timeDataHandle,
                                                         Dune::InteriorBorder_All_Interface,
                                                         Dune::ForwardCommunication);
    }

    dt = problem_.gridView().comm().min(dt);
#endif

    if(impet)
    {
        Dune::dinfo << "Timestep restricted by CellIdx " << restrictingCell << " leads to dt = "
                <<dt * getParam<Scalar>("Impet.CFLFactor")<< std::endl;
        if (averagedFaces_ != 0)
            Dune::dinfo  << " Averageing done for " << averagedFaces_ << " faces. "<< std::endl;
    }
}
/*!
 * \brief Updates the transported quantity once an update is calculated.
 *
 * This method updates both, the internal transport solution vector and the entries in the cellData.
 * \param updateVector Update vector, or update estimate for secants, resp. Here in \f$\mathrm{[kg/m^3]}\f$
 */
template<class TypeTag>
void FVTransport2P2C<TypeTag>::updateTransportedQuantity(TransportSolutionType& updateVector)
{
    if (this->enableLocalTimeStepping())
        this->innerUpdate(updateVector);
    else
        updateConcentrations(updateVector, problem().timeManager().timeStepSize());
}

/*!
 * \brief Updates the transported quantity once an update is calculated.
 *
 * This method updates both, the internal transport solution vector and the entries in the cellData.
 * \param updateVector Update vector, or update estimate for secants, resp. Here in \f$\mathrm{[kg/m^3]}\f$
 * \param dt Time step size \f$\mathrm{[s]}\f$
 */
template<class TypeTag>
void FVTransport2P2C<TypeTag>::updateTransportedQuantity(TransportSolutionType& updateVector, Scalar dt)
{
    updateConcentrations(updateVector, dt);
}

/*!
 * \brief Updates the concentrations once an update is calculated.
 *
 * This method updates both, the internal transport solution vector and the entries in the cellData.
 * \param updateVector Update vector, or update estimate for secants, resp. Here in \f$\mathrm{[kg/m^3]}\f$
 * \param dt Time step size \f$\mathrm{[s]}\f$
 */
template<class TypeTag>
void FVTransport2P2C<TypeTag>::updateConcentrations(TransportSolutionType& updateVector, Scalar dt)
{
    // loop thorugh all elements
    for (int i = 0; i< problem().gridView().size(0); i++)
    {
        CellData& cellDataI = problem().variables().cellData(i);
        for(int compIdx = 0; compIdx < getPropValue<TypeTag, Properties::NumComponents>(); compIdx++)
        {
            totalConcentration_[compIdx][i] += (updateVector[compIdx][i]*=dt);
            cellDataI.setMassConcentration(compIdx, totalConcentration_[compIdx][i]);
        }
    }
}

/*!
 * \brief  Get flux at an interface between two cells
 * The flux through \f$ \gamma \f$  is calculated according to the underlying pressure field,
 * calculated by the pressure model.
 *  \f[ - A_{\gamma} \mathbf{n}^T_{\gamma} \mathbf{K}  \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha}
     \mathbf{d}_{ij}  \left( \frac{p_{\alpha,j}^t - p^{t}_{\alpha,i}}{\Delta x} + \varrho_{\alpha} \mathbf{g}^T \mathbf{d}_{ij} \right)
                \sum_{\kappa} X^{\kappa}_{\alpha} \f]
 * Here, \f$ \mathbf{d}_{ij} \f$ is the normalized vector connecting the cell centers, and \f$ \mathbf{n}_{\gamma} \f$
 * represents the normal of the face \f$ \gamma \f$. Due to the nature of the primay Variable, the (volume-)specific
 * total mass concentration, this represents a mass flux per cell volume.
 *
 * \param fluxEntries The flux entries, mass influx from cell \f$j\f$ to \f$i\f$.
 * \param timestepFlux flow velocities for timestep estimation
 * \param intersection The intersection
 * \param cellDataI The cell data for cell \f$i\f$
 */
template<class TypeTag>
void FVTransport2P2C<TypeTag>::getFlux(ComponentVector& fluxEntries,
                                        Dune::FieldVector<Scalar, 2>& timestepFlux,
                                        const Intersection& intersection,
                                        CellData& cellDataI)
{
    fluxEntries = 0.;
    timestepFlux = 0.;
    // cell information
    auto elementI = intersection.inside();
    int eIdxGlobalI = problem().variables().index(elementI);

    // get position
    const GlobalPosition globalPos = elementI.geometry().center();
    const GlobalPosition& gravity_ = problem().gravity();
    // cell volume, assume linear map here
    Scalar volume = elementI.geometry().volume();

    // get values of cell I
    Scalar pressI = problem().pressureModel().pressure(eIdxGlobalI);
    Scalar pcI = cellDataI.capillaryPressure();
    DimMatrix K_I(problem().spatialParams().intrinsicPermeability(elementI));

    // old material law interface is deprecated: Replace this by
    // const auto& fluidMatrixInteraction = problem().spatialParams.fluidMatrixInteractionAtPos(elementI.geometry().center());
    // after the release of 3.3, when the deprecated interface is no longer supported
    const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem().spatialParams(), elementI);

    PhaseVector SmobI(0.);
    using std::max;
    SmobI[wPhaseIdx] = max((cellDataI.saturation(wPhaseIdx)
                            - fluidMatrixInteraction.pcSwCurve().effToAbsParams().swr())
                            , 1e-2);
    SmobI[nPhaseIdx] = max((cellDataI.saturation(nPhaseIdx)
                                - fluidMatrixInteraction.pcSwCurve().effToAbsParams().snr())
                            , 1e-2);

    Scalar densityWI (0.), densityNWI(0.);
    densityWI= cellDataI.density(wPhaseIdx);
    densityNWI = cellDataI.density(nPhaseIdx);

    // face properties
    GlobalPosition unitOuterNormal = intersection.centerUnitOuterNormal();
    if (switchNormals)
        unitOuterNormal *= -1.0;
    Scalar faceArea = intersection.geometry().volume();

//    // local interface index
//    const int indexInInside = intersection.indexInInside();
//    const int indexInOutside = intersection.indexInOutside();

    // create vector for timestep and for update
    Dune::FieldVector<Scalar, 2> factor (0.);
    Dune::FieldVector<Scalar, 2> updFactor (0.);

    PhaseVector potential(0.);

    // access neighbor
    auto neighbor = intersection.outside();
    int eIdxGlobalJ = problem().variables().index(neighbor);
    CellData& cellDataJ = problem().variables().cellData(eIdxGlobalJ);

    // neighbor cell center in global coordinates
    const GlobalPosition& globalPosNeighbor = neighbor.geometry().center();

    // distance vector between barycenters
    GlobalPosition distVec = globalPosNeighbor - globalPos;
    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    GlobalPosition unitDistVec(distVec);
    unitDistVec /= dist;

    // phase densities in neighbor
    Scalar densityWJ (0.), densityNWJ(0.);
    densityWJ = cellDataJ.density(wPhaseIdx);
    densityNWJ = cellDataJ.density(nPhaseIdx);

    // average phase densities with central weighting
    double densityW_mean = (densityWI + densityWJ) * 0.5;
    double densityNW_mean = (densityNWI + densityNWJ) * 0.5;

    double pressJ = problem().pressureModel().pressure(eIdxGlobalJ);
    Scalar pcJ = cellDataJ.capillaryPressure();

    // compute mean permeability
    DimMatrix meanK_(0.);
    harmonicMeanMatrix(meanK_,
            K_I,
            problem().spatialParams().intrinsicPermeability(neighbor));
    Dune::FieldVector<Scalar,dim> K(0);
    meanK_.umv(unitDistVec,K);

    // determine potentials for upwind
    switch (pressureType)
    {
    case pw:
    {
        potential[wPhaseIdx] = (pressI - pressJ) / (dist);
        potential[nPhaseIdx] = (pressI - pressJ + pcI - pcJ) / (dist);
        break;
    }
    case pn:
    {
        potential[wPhaseIdx] = (pressI - pressJ - pcI + pcJ) / (dist);
        potential[nPhaseIdx] = (pressI - pressJ) / (dist);
        break;
    }
    }
    // add gravity term
    potential[nPhaseIdx] +=  (gravity_ * unitDistVec) * densityNW_mean;
    potential[wPhaseIdx] +=  (gravity_ * unitDistVec) * densityW_mean;

    potential[wPhaseIdx] *= fabs(K * unitOuterNormal);
    potential[nPhaseIdx] *= fabs(K * unitOuterNormal);

    // determine upwinding direction, perform upwinding if possible
    Dune::FieldVector<bool, NumPhases> doUpwinding(true);
    PhaseVector lambda(0.);
    for(int phaseIdx = 0; phaseIdx < NumPhases; phaseIdx++)
    {
        int contiEqIdx = 0;
        if(phaseIdx == wPhaseIdx)
            contiEqIdx = contiWEqIdx;
        else
            contiEqIdx = contiNEqIdx;

        if(!impet_ || restrictFluxInTransport_==0) // perform a strict uwpind scheme
        {
            if (potential[phaseIdx] > 0.)
            {
                lambda[phaseIdx] = cellDataI.mobility(phaseIdx);
                cellDataI.setUpwindCell(intersection.indexInInside(), contiEqIdx, true);

            }
            else if(potential[phaseIdx] < 0.)
            {
                lambda[phaseIdx] = cellDataJ.mobility(phaseIdx);
                cellDataI.setUpwindCell(intersection.indexInInside(), contiEqIdx, false);
            }
            else
            {
                doUpwinding[phaseIdx] = false;
                cellDataI.setUpwindCell(intersection.indexInInside(), contiEqIdx, false);
                cellDataJ.setUpwindCell(intersection.indexInOutside(), contiEqIdx, false);
            }
        }
        else // Transport after PE with check on flow direction
        {
            bool wasUpwindCell = cellDataI.isUpwindCell(intersection.indexInInside(), contiEqIdx);

            if (potential[phaseIdx] > 0. && wasUpwindCell)
                lambda[phaseIdx] = cellDataI.mobility(phaseIdx);
            else if (potential[phaseIdx] < 0. && !wasUpwindCell)
                lambda[phaseIdx] = cellDataJ.mobility(phaseIdx);
            // potential direction does not coincide with that of P.E.
            else if(restrictFluxInTransport_ == 2)   // perform central averaging for all direction changes
                doUpwinding[phaseIdx] = false;
            else    // i.e. restrictFluxInTransport == 1
            {
               //check if harmonic weighting is necessary
                if (potential[phaseIdx] > 0. && (Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(cellDataJ.mobility(phaseIdx), 0.0, 1.0e-30)   // check if outflow induce neglected (i.e. mob=0) phase flux
                       || (cellDataI.wasRefined() && cellDataJ.wasRefined() && elementI.father() == neighbor.father())))
                    lambda[phaseIdx] = cellDataI.mobility(phaseIdx);
                else if (potential[phaseIdx] < 0. && (Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(cellDataI.mobility(phaseIdx), 0.0, 1.0e-30) // check if inflow induce neglected phase flux
                        || (cellDataI.wasRefined() && cellDataJ.wasRefined() && elementI.father() == neighbor.father())))
                    lambda[phaseIdx] = cellDataJ.mobility(phaseIdx);
                else
                    doUpwinding[phaseIdx] = false;
            }

        }

        // do not perform upwinding if so desired
        if(!doUpwinding[phaseIdx])
        {
            //a) no flux if there wont be any flux regardless how to average/upwind
            if(cellDataI.mobility(phaseIdx)+cellDataJ.mobility(phaseIdx)==0.)
            {
                potential[phaseIdx] = 0;
                continue;
            }

            //b) perform harmonic averaging
            fluxEntries[wCompIdx] -= potential[phaseIdx] * faceArea / volume
                    * harmonicMean(cellDataI.massFraction(phaseIdx, wCompIdx) * cellDataI.mobility(phaseIdx) * cellDataI.density(phaseIdx),
                            cellDataJ.massFraction(phaseIdx, wCompIdx) * cellDataJ.mobility(phaseIdx) * cellDataJ.density(phaseIdx));
            fluxEntries[nCompIdx] -= potential[phaseIdx] * faceArea / volume
                    * harmonicMean(cellDataI.massFraction(phaseIdx, nCompIdx) * cellDataI.mobility(phaseIdx) * cellDataI.density(phaseIdx),
                            cellDataJ.massFraction(phaseIdx, nCompIdx) * cellDataJ.mobility(phaseIdx) * cellDataJ.density(phaseIdx));
            // c) timestep control
            // for timestep control : influx
            using std::max;
            timestepFlux[0] += max(0.,
                    - potential[phaseIdx] * faceArea / volume
                      * harmonicMean(cellDataI.mobility(phaseIdx),cellDataJ.mobility(phaseIdx)));
            // outflux
            timestepFlux[1] += max(0.,
                    potential[phaseIdx] * faceArea / volume
                    * harmonicMean(cellDataI.mobility(phaseIdx),cellDataJ.mobility(phaseIdx))/SmobI[phaseIdx]);

            //d) output
            if(!(cellDataI.wasRefined() && cellDataJ.wasRefined() && elementI.father() == neighbor.father())
                    && eIdxGlobalI > eIdxGlobalJ) //(only for one side)
            {
                averagedFaces_++;
                #if DUNE_MINIMAL_DEBUG_LEVEL < 3
                // verbose (only for one side)
                if(eIdxGlobalI > eIdxGlobalJ)
                    Dune::dinfo << "harmonicMean flux of phase" << phaseIdx <<" used from cell" << eIdxGlobalI<< " into " << eIdxGlobalJ
                    << " ; TE upwind I = "<< cellDataI.isUpwindCell(intersection.indexInInside(), contiEqIdx)
                    << " but pot = "<< potential[phaseIdx] <<  std::endl;
                #endif
            }

            //e) stop further standard calculations
            potential[phaseIdx] = 0;
        }
    }

    // calculate and standardized velocity
    using std::max;
    double velocityJIw = max((-lambda[wPhaseIdx] * potential[wPhaseIdx]) * faceArea / volume, 0.0);
    double velocityIJw = max(( lambda[wPhaseIdx] * potential[wPhaseIdx]) * faceArea / volume, 0.0);
    double velocityJIn = max((-lambda[nPhaseIdx] * potential[nPhaseIdx]) * faceArea / volume, 0.0);
    double velocityIJn = max(( lambda[nPhaseIdx] * potential[nPhaseIdx]) * faceArea / volume, 0.0);

    // for timestep control : influx
    timestepFlux[0] += velocityJIw + velocityJIn;

    double foutw = velocityIJw/SmobI[wPhaseIdx];
    double foutn = velocityIJn/SmobI[nPhaseIdx];
    using std::isnan;
    using std::isinf;
    if (isnan(foutw) || isinf(foutw) || foutw < 0) foutw = 0;
    if (isnan(foutn) || isinf(foutn) || foutn < 0) foutn = 0;
    timestepFlux[1] += foutw + foutn;

    fluxEntries[wCompIdx] +=
        velocityJIw * cellDataJ.massFraction(wPhaseIdx, wCompIdx) * densityWJ
        - velocityIJw * cellDataI.massFraction(wPhaseIdx, wCompIdx) * densityWI
        + velocityJIn * cellDataJ.massFraction(nPhaseIdx, wCompIdx) * densityNWJ
        - velocityIJn * cellDataI.massFraction(nPhaseIdx, wCompIdx) * densityNWI;
    fluxEntries[nCompIdx] +=
        velocityJIw * cellDataJ.massFraction(wPhaseIdx, nCompIdx) * densityWJ
        - velocityIJw * cellDataI.massFraction(wPhaseIdx, nCompIdx) * densityWI
        + velocityJIn * cellDataJ.massFraction(nPhaseIdx, nCompIdx) * densityNWJ
        - velocityIJn * cellDataI.massFraction(nPhaseIdx, nCompIdx) * densityNWI;


//    //add dispersion
//    for (int phaseIdx = 0; phaseIdx < NumPhases; ++phaseIdx)
//    {
//
//        //calculate porous media diffusion coefficient
//        // for phases which exist in both finite volumes
//        if (cellDataI.saturation(phaseIdx) <= 0 ||
//                cellDataJ.saturation(phaseIdx) <= 0)
//        {
//            continue;
//        }
//
//        // calculate tortuosity at the nodes i and j needed
//        // for porous media diffusion coefficient
//        Scalar poroI = problem().spatialParams().porosity(elementI);
//        Scalar poroJ = problem().spatialParams().porosity(neighbor);
//        Scalar tauI =
//            1.0/(poroI * poroI) *
//            pow(poroI * cellDataI.saturation(phaseIdx), 7.0/3);
//        Scalar tauJ =
//            1.0/(poroJ * poroJ) *
//            pow(poroJ * cellDataJ.saturation(phaseIdx), 7.0/3);
//        // Diffusion coefficient in the porous medium
//        // -> harmonic mean
//        Scalar porousDiffCoeff
//                = harmonicMean(poroI * cellDataI.saturation(phaseIdx) * tauI
//                        * FluidSystem::binaryDiffusionCoefficient(cellDataI.fluidState(), phaseIdx, wCompIdx, nCompIdx),
//                        poroJ * cellDataJ.saturation(phaseIdx) * tauJ
//                        * FluidSystem::binaryDiffusionCoefficient(cellDataJ.fluidState(), phaseIdx, wCompIdx, nCompIdx));
//        Scalar averagedMolarDensity = (cellDataI.density(phaseIdx) / cellDataI.fluidState().averageMolarMass(phaseIdx)
//                +cellDataJ.density(phaseIdx) / cellDataJ.fluidState().averageMolarMass(phaseIdx))*0.5;
//
//        for (int compIdx = 0; compIdx < NumComponents; ++compIdx)
//        {
//
//            Scalar gradx = (cellDataJ.moleFraction(phaseIdx, compIdx)-cellDataI.moleFraction(phaseIdx, compIdx) )/ dist;
//            Scalar gradX = (cellDataJ.massFraction(phaseIdx, compIdx)-cellDataI.massFraction(phaseIdx, compIdx) )/ dist;
//
//            fluxEntries[compIdx] += gradx * averagedMolarDensity * porousDiffCoeff *  faceArea / volume;
//        }
//    }

    return;
}

/*!
 * \brief Get flux on Boundary
 *
 * The flux through \f$ \gamma \f$  is calculated according to the underlying pressure field,
 * calculated by the pressure model.
 *  \f[ - A_{\gamma}  \mathbf{n}^T_{\gamma} \mathbf{K} \mathbf{d}_{i-Boundary}
      \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha} \sum_{\kappa} X^{\kappa}_{\alpha}
    \left( \frac{p_{\alpha,Boundary}^t - p^{t}_{\alpha,i}}{\Delta x} + \varrho_{\alpha}\mathbf{g}^T \mathbf{d}_{i-Boundary}\right) \f]
 * Here, \f$ \mathbf{u} \f$ is the normalized vector connecting the cell centers, and \f$ \mathbf{n}_{\gamma_{ij}} \f$
 * represents the normal of the face \f$ \gamma{ij} \f$.
 * \param fluxEntries The flux entries, mass influx from cell \f$j\f$ to \f$i\f$.
 * \param timestepFlux flow velocities for timestep estimation
 * \param intersection The intersection
 * \param cellDataI The cell data for cell \f$i\f$
 */
template<class TypeTag>
void FVTransport2P2C<TypeTag>::getFluxOnBoundary(ComponentVector& fluxEntries,
                                                 Dune::FieldVector<Scalar, 2>& timestepFlux,
                                                 const Intersection& intersection,
                                                 const CellData& cellDataI)
{
    using std::max;
    // cell information
    auto elementI = intersection.inside();
    int eIdxGlobalI = problem().variables().index(elementI);

    // get position
    const GlobalPosition globalPos = elementI.geometry().center();

    // cell volume, assume linear map here
    Scalar volume = elementI.geometry().volume();
    const GlobalPosition& gravity_ = problem().gravity();
    // get values of cell I
    Scalar pressI = problem().pressureModel().pressure(eIdxGlobalI);
    Scalar pcI = cellDataI.capillaryPressure();
    DimMatrix K_I(problem().spatialParams().intrinsicPermeability(elementI));

    if(regulateBoundaryPermeability)
    {
        int axis = intersection.indexInInside() / 2;
        if(K_I[axis][axis] < minimalBoundaryPermeability)
            K_I[axis][axis] = minimalBoundaryPermeability;
    }

    // old material law interface is deprecated: Replace this by
    // const auto& fluidMatrixInteraction = problem().spatialParams.fluidMatrixInteractionAtPos(elementI.geometry().center());
    // after the release of 3.3, when the deprecated interface is no longer supported
    const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem().spatialParams(), elementI);

    Scalar SwmobI = max((cellDataI.saturation(wPhaseIdx)
                            - fluidMatrixInteraction.pcSwCurve().effToAbsParams().swr())
                            , 1e-2);
    Scalar SnmobI = max((cellDataI.saturation(nPhaseIdx)
                                - fluidMatrixInteraction.pcSwCurve().effToAbsParams().snr())
                            , 1e-2);

    Scalar densityWI (0.), densityNWI(0.);
    densityWI= cellDataI.density(wPhaseIdx);
    densityNWI = cellDataI.density(nPhaseIdx);

    // face properties
    GlobalPosition unitOuterNormal = intersection.centerUnitOuterNormal();
    if (switchNormals)
        unitOuterNormal *= -1.0;
    Scalar faceArea = intersection.geometry().volume();

    // create vector for timestep and for update
    Dune::FieldVector<Scalar, 2> factor (0.);
    Dune::FieldVector<Scalar, 2> updFactor (0.);

    PhaseVector potential(0.);
    // center of face in global coordinates
    const GlobalPosition& globalPosFace = intersection.geometry().center();

    // distance vector between barycenters
    GlobalPosition distVec = globalPosFace - globalPos;
    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    GlobalPosition unitDistVec(distVec);
    unitDistVec /= dist;

    //instantiate a fluid state
    FluidState BCfluidState;

    //get boundary type
    BoundaryTypes bcTypes;
    problem().boundaryTypes(bcTypes, intersection);

    /**********         Dirichlet Boundary        *************/
    if (bcTypes.isDirichlet(contiWEqIdx)) // if contiWEq is Dirichlet, so is contiNEq
    {
        //get dirichlet pressure boundary condition
        PhaseVector pressBound(0.);
        Scalar pcBound (0.);

        // read boundary values
        this->evalBoundary(globalPosFace,
                           intersection,
                           BCfluidState,
                           pressBound);

        // determine fluid properties at the boundary
        Scalar densityWBound = BCfluidState.density(wPhaseIdx);
        Scalar densityNWBound = BCfluidState.density(nPhaseIdx);
        Scalar viscosityWBound = FluidSystem::viscosity(BCfluidState, wPhaseIdx);
        Scalar viscosityNWBound = FluidSystem::viscosity(BCfluidState, nPhaseIdx);
        if(getPropValue<TypeTag, Properties::EnableCapillarity>())
            pcBound = (BCfluidState.pressure(nPhaseIdx) - BCfluidState.pressure(wPhaseIdx));
        // average
        double densityW_mean = (densityWI + densityWBound) / 2;
        double densityNW_mean = (densityNWI + densityNWBound) / 2;

        // prepare K
        Dune::FieldVector<Scalar,dim> K(0);
        K_I.umv(unitDistVec,K);

        //calculate potential gradient
        switch (pressureType)
        {
            case pw:
            {
                potential[wPhaseIdx] = (pressI - pressBound[wPhaseIdx]) / dist;
                potential[nPhaseIdx] = (pressI + pcI - pressBound[wPhaseIdx] - pcBound)
                        / dist;
                break;
            }
            case pn:
            {
                potential[wPhaseIdx] = (pressI - pcI - pressBound[nPhaseIdx] + pcBound)
                        / dist;
                potential[nPhaseIdx] = (pressI - pressBound[nPhaseIdx]) / dist;
                break;
            }
        }
        potential[wPhaseIdx] += (gravity_ * unitDistVec) * densityW_mean;
        potential[nPhaseIdx] += (gravity_ * unitDistVec) * densityNW_mean;

        potential[wPhaseIdx] *= fabs(K * unitOuterNormal);
        potential[nPhaseIdx] *= fabs(K * unitOuterNormal);

        // old material law interface is deprecated: Replace this by
        // const auto& fluidMatrixInteraction = problem().spatialParams.fluidMatrixInteractionAtPos(elementI.geometry().center());
        // after the release of 3.3, when the deprecated interface is no longer supported
        const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem().spatialParams(), elementI);

        // do upwinding for lambdas
        PhaseVector lambda(0.);
        if (potential[wPhaseIdx] >= 0.)
            lambda[wPhaseIdx] = cellDataI.mobility(wPhaseIdx);
        else
            {
            if(getPropValue<TypeTag, Properties::BoundaryMobility>()==Indices::satDependent)
                lambda[wPhaseIdx] = BCfluidState.saturation(wPhaseIdx) / viscosityWBound;
            else
                lambda[wPhaseIdx] = fluidMatrixInteraction.krw(BCfluidState.saturation(wPhaseIdx)) / viscosityWBound;
            }
        if (potential[nPhaseIdx] >= 0.)
            lambda[nPhaseIdx] = cellDataI.mobility(nPhaseIdx);
        else
            {
            if(getPropValue<TypeTag, Properties::BoundaryMobility>()==Indices::satDependent)
                lambda[nPhaseIdx] = BCfluidState.saturation(nPhaseIdx) / viscosityNWBound;
            else
                lambda[nPhaseIdx] = fluidMatrixInteraction.krn(BCfluidState.saturation(wPhaseIdx)) / viscosityNWBound;
            }
        // calculate and standardized velocity

        double velocityJIw = max((-lambda[wPhaseIdx] * potential[wPhaseIdx]) * faceArea / volume, 0.0);
        double velocityIJw = max(( lambda[wPhaseIdx] * potential[wPhaseIdx]) * faceArea / volume, 0.0);
        double velocityJIn = max((-lambda[nPhaseIdx] * potential[nPhaseIdx]) * faceArea / volume, 0.0);
        double velocityIJn = max(( lambda[nPhaseIdx] * potential[nPhaseIdx]) * faceArea / volume, 0.0);

        // for timestep control
        timestepFlux[0] = velocityJIw + velocityJIn;

        double foutw = velocityIJw/SwmobI;
        double foutn = velocityIJn/SnmobI;
        using std::isnan;
        using std::isinf;
        if (isnan(foutw) || isinf(foutw) || foutw < 0) foutw = 0;
        if (isnan(foutn) || isinf(foutn) || foutn < 0) foutn = 0;
        timestepFlux[1] = foutw + foutn;

        fluxEntries[wCompIdx] =
            + velocityJIw * BCfluidState.massFraction(wPhaseIdx, wCompIdx) * densityWBound
            - velocityIJw * cellDataI.massFraction(wPhaseIdx, wCompIdx) * densityWI
            + velocityJIn * BCfluidState.massFraction(nPhaseIdx, wCompIdx) * densityNWBound
            - velocityIJn * cellDataI.massFraction(nPhaseIdx, wCompIdx) * densityNWI ;
        fluxEntries[nCompIdx] =
            velocityJIw * BCfluidState.massFraction(wPhaseIdx, nCompIdx) * densityWBound
            - velocityIJw * cellDataI.massFraction(wPhaseIdx, nCompIdx) * densityWI
            + velocityJIn * BCfluidState.massFraction(nPhaseIdx, nCompIdx) * densityNWBound
            - velocityIJn * cellDataI.massFraction(nPhaseIdx, nCompIdx) * densityNWI ;
    }//end dirichlet boundary
    else if (bcTypes.isNeumann(contiWEqIdx))
    {
        // Convention: outflow => positive sign : has to be subtracted from update vec
        PrimaryVariables J(NAN);
        problem().neumann(J, intersection);
        fluxEntries[wCompIdx] = - J[contiWEqIdx] * faceArea / volume;
        fluxEntries[nCompIdx] = - J[contiNEqIdx] * faceArea / volume;

        // for timestep control
        #define cflIgnoresNeumann
        #ifdef cflIgnoresNeumann
        timestepFlux[0] = 0;
        timestepFlux[1] = 0;
        #else
        double inflow = updFactor[wCompIdx] / densityW + updFactor[nCompIdx] / densityNW;
        if (inflow>0)
            {
            timestepFlux[0] = updFactor[wCompIdx] / densityW
                        + updFactor[nCompIdx] / densityNW;    // =factor in
            timestepFlux[1] = -(updFactor[wCompIdx] / densityW /SwmobI
                        + updFactor[nCompIdx] / densityNW / SnmobI);    // =factor out
            }
        else
        {
            timestepFlux[0] = -(updFactor[wCompIdx] / densityW
                            + updFactor[nCompIdx] / densityNW);    // =factor in
            timestepFlux[1] = updFactor[wCompIdx] / densityW /SwmobI
                        + updFactor[nCompIdx] / densityNW / SnmobI;    // =factor out
        }
        #endif
    }//end neumann boundary
    return;
}


/*!
 * \brief Evaluate the boundary conditions
 *
 *  As the transport primary variable in this formulation is the total component
 *  concentration, \f$ C^{\kappa} \f$ it seems natural that the boundary values
 *  are also total concentrations. However, as for the initial conditions, it is
 *  possible to define boundaries by means of a saturation. This choice determines
 *  which version of flash calculation is necessary to get to the composition at
 *  the boundary.
 *  \param globalPosFace The global position of the face of the current boundary cell
 *  \param intersection The current intersection
 *  \param BCfluidState FluidState object that is used for the boundary
 *  \param pressBound Boundary values of phase pressures
 */
template<class TypeTag>
void FVTransport2P2C<TypeTag>::evalBoundary(GlobalPosition globalPosFace,
                                            const Intersection& intersection,
                                            FluidState &BCfluidState,
                                            PhaseVector &pressBound)
{
    // prepare a flash solver
    CompositionalFlash<Scalar, FluidSystem> flashSolver;

    auto element = intersection.inside();
    // read boundary values
    PrimaryVariables primaryVariablesOnBoundary(0.);
    problem().dirichlet(primaryVariablesOnBoundary, intersection);

    // old material law interface is deprecated: Replace this by
    // const auto& fluidMatrixInteraction = problem().spatialParams.fluidMatrixInteractionAtPos(element.geometry().center());
    // after the release of 3.3, when the deprecated interface is no longer supported
    const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem().spatialParams(), element);

    // read boundary type
    typename Indices::BoundaryFormulation bcType;
    problem().boundaryFormulation(bcType, intersection);
    if (bcType == Indices::saturation)
    {
        Scalar satBound = primaryVariablesOnBoundary[contiWEqIdx];

        if(getPropValue<TypeTag, Properties::EnableCapillarity>())
        {
            Scalar pcBound = fluidMatrixInteraction.pc(satBound);
            switch (pressureType)
            {
            case pw:
            {
                pressBound[wPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx];
                pressBound[nPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx] + pcBound;
                break;
            }
            case pn:
            {
                pressBound[wPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx] - pcBound;
                pressBound[nPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx];
                break;
            }
            }
        }
        else // capillarity neglected
            pressBound[wPhaseIdx] = pressBound[nPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx];

        flashSolver.saturationFlash2p2c(BCfluidState, satBound, pressBound,
                problem().temperatureAtPos(globalPosFace));
    }
    else if (bcType == Indices::concentration)
    {
        // saturation and hence pc and hence corresponding pressure unknown
        pressBound[wPhaseIdx] = pressBound[nPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx];
        Scalar Z0Bound = primaryVariablesOnBoundary[contiWEqIdx];
        flashSolver.concentrationFlash2p2c(BCfluidState, Z0Bound, pressBound,
            problem().temperatureAtPos(globalPosFace));

        if(getPropValue<TypeTag, Properties::EnableCapillarity>())
        {
            Scalar pcBound = fluidMatrixInteraction.pc(BCfluidState.saturation(wPhaseIdx));
            int maxiter = 3;
            //start iteration loop
            for(int iter=0; iter < maxiter; iter++)
            {
                //prepare pressures to enter flash calculation
                switch (pressureType)
                {
                case pw:
                {
                    pressBound[wPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx];
                    pressBound[nPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx]
                                                                                        + pcBound;
                    break;
                }
                case pn:
                {
                    pressBound[wPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx]
                                                                                         - pcBound;
                    pressBound[nPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx];
                    break;
                }
                }

                //store old pc
                Scalar oldPc = pcBound;
                //update with better pressures
                flashSolver.concentrationFlash2p2c(BCfluidState, Z0Bound, pressBound,
                        problem().temperatureAtPos(globalPosFace));
                pcBound = fluidMatrixInteraction.pc(BCfluidState.saturation(wPhaseIdx));
                // TODO: get right criterion, do output for evaluation
                //converge criterion
                using std::abs;
                if (abs(oldPc - pcBound) < 10.0)
                    iter = maxiter;
            }
        }
    }
    else
        DUNE_THROW(Dune::NotImplemented, "Boundary Formulation neither Concentration nor Saturation??");
}

template<class TypeTag>
void FVTransport2P2C<TypeTag>::updatedTargetDt_(Scalar &dt)
{
    dt = std::numeric_limits<Scalar>::max();

    // update target time-step-sizes
    for (const auto& element : elements(problem_.gridView()))
    {
#if HAVE_MPI
        if (element.partitionType() != Dune::InteriorEntity)
        {
            continue;
        }
#endif

        // cell index
        int eIdxGlobalI = problem_.variables().index(element);

        LocalTimesteppingData& localDataI = timeStepData_[eIdxGlobalI];


        using FaceDt = std::unordered_map<int, Scalar>;
        FaceDt faceDt;

        // run through all intersections with neighbors and boundary
        for (const auto& intersection : intersections(problem_.gridView(), element))
        {
            int indexInInside = intersection.indexInInside();
            using std::min;
            if (intersection.neighbor())
            {
                auto neighbor = intersection.outside();
                int eIdxGlobalJ = problem_.variables().index(neighbor);

                int levelI = element.level();
                int levelJ = neighbor.level();


                if (eIdxGlobalI < eIdxGlobalJ && levelI <= levelJ)
                {
                    LocalTimesteppingData& localDataJ = timeStepData_[eIdxGlobalJ];

                    int indexInOutside = intersection.indexInOutside();

                    if (localDataI.faceTargetDt[indexInInside] < accumulatedDt_ + dtThreshold_
                        || localDataJ.faceTargetDt[indexInOutside] < accumulatedDt_ + dtThreshold_)
                    {
                        Scalar timeStep  = min(localDataI.dt, localDataJ.dt);

                        if (levelI < levelJ)
                        {
                            typename FaceDt::iterator it = faceDt.find(indexInInside);
                            if (it != faceDt.end())
                            {
                                it->second = min(it->second, timeStep);
                            }
                            else
                            {
                                faceDt.insert(std::make_pair(indexInInside, timeStep));
                            }
                        }
                        else
                        {
                            localDataI.faceTargetDt[indexInInside] += subCFLFactor_ * timeStep;
                            localDataJ.faceTargetDt[indexInOutside] += subCFLFactor_ * timeStep;
                        }

                        dt = min(dt, timeStep);
                    }
                }
            }
            else if (intersection.boundary())
            {
                if (localDataI.faceTargetDt[indexInInside] < accumulatedDt_ + dtThreshold_)
                {
                    localDataI.faceTargetDt[indexInInside] += subCFLFactor_ * localDataI.dt;
                    dt = min(dt, subCFLFactor_ * localDataI.dt);
                }
            }
        }
        if (faceDt.size() > 0)
        {
            typename FaceDt::iterator it = faceDt.begin();
            for (;it != faceDt.end();++it)
            {
                localDataI.faceTargetDt[it->first] += subCFLFactor_ * it->second;
            }

            for (const auto& intersection : intersections(problem_.gridView(), element))
            {
                if (intersection.neighbor())
                {
                    int indexInInside = intersection.indexInInside();

                    it = faceDt.find(indexInInside);
                    if (it != faceDt.end())
                    {
                        auto neighbor = intersection.outside();
                        int eIdxGlobalJ = problem_.variables().index(neighbor);

                        LocalTimesteppingData& localDataJ = timeStepData_[eIdxGlobalJ];

                        int indexInOutside = intersection.indexInOutside();

                        localDataJ.faceTargetDt[indexInOutside] += subCFLFactor_ * it->second;
                    }
                }
            }
        }
    }

#if HAVE_MPI
    // communicate updated values
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using ElementMapper = typename SolutionTypes::ElementMapper;
    using TimeDataHandle = VectorCommDataHandleEqual<ElementMapper, std::vector<LocalTimesteppingData>, 0/*elementCodim*/>;

    TimeDataHandle timeDataHandle(problem_.elementMapper(), timeStepData_);
    problem_.gridView().template communicate<TimeDataHandle>(timeDataHandle,
                                                         Dune::InteriorBorder_All_Interface,
                                                         Dune::ForwardCommunication);

    dt = problem_.gridView().comm().min(dt);
#endif
}

template<class TypeTag>
void FVTransport2P2C<TypeTag>::innerUpdate(TransportSolutionType& updateVec)
{
    if (localTimeStepping_)
    {
        Scalar realDt = problem_.timeManager().timeStepSize();

        Scalar subDt = realDt;

        updatedTargetDt_(subDt);

        Scalar accumulatedDtOld = accumulatedDt_;
        accumulatedDt_ += subDt;

        Scalar t = problem_.timeManager().time();

        if (accumulatedDt_ < realDt)
        {
            using std::min;
            while(true)
            {
                Scalar dtCorrection = min(0.0, realDt - accumulatedDt_);
                subDt += dtCorrection;

                if (verbosity_ > 0)
                    std::cout<<"    Sub-time-step size: "<<subDt<< std::endl;

                bool stopTimeStep = false;
                int size = problem_.gridView().size(0);
                for (int i = 0; i < size; i++)
                {
                    EntryType newVal(0);
                    int transportedQuantities = getPropValue<TypeTag, Properties::NumEq>() - 1; // NumEq - 1 pressure Eq
                    for (int eqNumber = 0; eqNumber < transportedQuantities; eqNumber++)
                    {
                        newVal[eqNumber] = totalConcentration_[eqNumber][i];
                        newVal[eqNumber] += updateVec[eqNumber][i] * subDt;
                    }
                    if (!asImp_().inPhysicalRange(newVal))
                    {
                        stopTimeStep = true;

                        break;
                    }
                }

#if HAVE_MPI
                int rank = 0;
                if (stopTimeStep)
                    rank = problem_.gridView().comm().rank();

                rank = problem_.gridView().comm().max(rank);
                problem_.gridView().comm().broadcast(&stopTimeStep,1,rank);
#endif


                if (stopTimeStep && accumulatedDtOld > dtThreshold_)
                {
                    problem_.timeManager().setTimeStepSize(accumulatedDtOld);
                    break;
                }
                else
                {
                    asImp_().updateTransportedQuantity(updateVec, subDt);
                }


                if (accumulatedDt_ >= realDt)
                {
                    break;
                }

                problem_.pressureModel().updateMaterialLaws();
                problem_.model().updateTransport(t, subDt, updateVec);

                updatedTargetDt_(subDt);

                accumulatedDtOld = accumulatedDt_;
                accumulatedDt_ += subDt;
            }
        }
        else
        {
            asImp_().updateTransportedQuantity(updateVec, realDt);
        }

        resetTimeStepData_();
    }
}

} // end namespace Dumux
#endif
