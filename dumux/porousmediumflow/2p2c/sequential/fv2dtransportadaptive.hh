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
 * \brief Finite volume discretization of the component transport equation.
 */
#ifndef DUMUX_FV2DTRANSPORT2P2C_ADAPTIVE_HH
#define DUMUX_FV2DTRANSPORT2P2C_ADAPTIVE_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

#include "adaptiveproperties.hh"
#include "fvtransport.hh"

#include <dumux/common/deprecated.hh>

namespace Dumux {
/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief Compositional Transport step in a Finite Volume discretization for a adaptive 2D-grid
 *
 *  The finite volume model for the solution of the transport equation for compositional
 *  two-phase flow.
 *  \f[
      \frac{\partial C^\kappa}{\partial t} = - \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha}
      \varrho_{alpha} \bf{v}_{\alpha}\right) + q^{\kappa},
 *  \f]
 *  where \f$ \bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 *  \f$ p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability,
 *  \f$ \lambda_{\alpha} \f$ the phase mobility,
 *  \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and
 *  \f$ C^{\kappa} \f$ the total Component concentration.
 *
 *  \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FV2dTransport2P2CAdaptive : public FVTransport2P2C<TypeTag>
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using CellData = GetPropType<TypeTag, Properties::CellData>;

    using TransportSolutionType = GetPropType<TypeTag, Properties::TransportSolutionType>;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld,
        NumPhases = getPropValue<TypeTag, Properties::NumPhases>()
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
        contiWEqIdx=Indices::contiWEqIdx, contiNEqIdx=Indices::contiNEqIdx
    };

    using IntersectionIterator = typename GridView::IntersectionIterator;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;
    using TransmissivityMatrix = Dune::FieldVector<Scalar,dim+1>;
    using PhaseVector = Dune::FieldVector<Scalar, NumPhases>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    //! Acess function for the current problem
    Problem& problem()
    { return problem_; }

    using LocalTimesteppingData = typename FVTransport2P2C<TypeTag>::LocalTimesteppingData;

public:
    virtual void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impes = false);

    void getMpfaFlux(Dune::FieldVector<Scalar, 2>&, Dune::FieldVector<Scalar, 2>&,
            const IntersectionIterator&, CellData&);

    /*!
     * \brief Constructs a FV2dTransport2P2CAdaptive object
     *
     * The compositional transport scheme can not be applied with a global pressure / total velocity
     * formulation. This is a 2D-specific implementation! In case of 3d, use the class
     * FV3dTransport2P2CAdaptive
     *
     * \param problem a problem class object
     */
    FV2dTransport2P2CAdaptive(Problem& problem) : FVTransport2P2C<TypeTag>(problem),
        problem_(problem)
    {
        enableMPFA = getParam<bool>("GridAdapt.EnableMultiPointFluxApproximation");
    }

    virtual ~FV2dTransport2P2CAdaptive()
    {
    }

protected:
    Problem& problem_;

    bool enableMPFA; //!> Specifies if the MPFA is used on hanging nodes
    //! gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
    static const int pressureType = getPropValue<TypeTag, Properties::PressureFormulation>();
};

/*!
 * \brief Calculate the update vector and determine timestep size
 *
 * This method calculates the update vector \f$ u \f$ of the discretized equation
 * \f[
      C^{\kappa , new} = C^{\kappa , old} + u,
 * \f]
 * where \f$ u = \sum_{\gamma} \boldsymbol{v}_{\alpha} * \varrho_{\alpha} * X^{\kappa}_{\alpha} * \boldsymbol{n} * A_{\gamma} \f$,
 * \f$ \boldsymbol{n} \f$ is the face normal and \f$ A_{\gamma} \f$ is the face area of face \f$ \gamma \f$.
 * In addition to the \a update vector, the recommended time step size \a dt is calculated
 * employing a CFL condition. This method uses a standard \a Tpfa method for regular fluxes,
 * and a \a MPFA can be used near hanging nodes.
 * The lengths of the vectors are resized to agree with the current grid resolution.
 *
 * \param t Current simulation time \f$\mathrm{[s]}\f$
 * \param[out] dt Time step size \f$\mathrm{[s]}\f$
 * \param[out] updateVec Update vector, or update estimate for secants, resp. Here in \f$\mathrm{[kg/m^3]}\f$
 * \param impet Flag that determines if it is a real impet step or an update estimate for volume derivatives
 */
template<class TypeTag>
void FV2dTransport2P2CAdaptive<TypeTag>::update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet)
{
    this->impet_ = impet;
    // initialize dt very large
    dt = 1E100;
    this->averagedFaces_ = 0;

    // resize update vector and set to zero
    unsigned int size_ = problem_.gridView().size(0);
    updateVec.resize(getPropValue<TypeTag, Properties::NumComponents>());
    updateVec[wCompIdx].resize(size_);
    updateVec[nCompIdx].resize(size_);
    updateVec[wCompIdx] = 0;
    updateVec[nCompIdx] = 0;

    if (this->enableLocalTimeStepping())
    {
        if (this->timeStepData_.size() != size_)
            this->timeStepData_.resize(size_);
    }

    //also resize PV vector if necessary
    if(this->totalConcentration_.size() != size_)
    {
        this->totalConcentration_[wCompIdx].resize(size_);
        this->totalConcentration_[nCompIdx].resize(size_);
        // copy data //TODO: remove this, remove PM Transport Vector!!
        // loop through all elements
        for (int i = 0; i< problem().gridView().size(0); i++)
        {
            CellData& cellDataI = problem().variables().cellData(i);
            for(int compIdx = 0; compIdx < getPropValue<TypeTag, Properties::NumComponents>(); compIdx++)
            {
                this->totalConcentration_[compIdx][i]
                        = cellDataI.totalConcentration(compIdx);
            }
        }
    }

    // Cell which restricts time step size
    int restrictingCell = -1;

    Dune::FieldVector<Scalar, 2> entries(0.), timestepFlux(0.);
    // compute update vector
    for (const auto& element : elements(problem().gridView()))
    {
        // get cell infos
        int globalIdxI = problem().variables().index(element);
        CellData& cellDataI = problem().variables().cellData(globalIdxI);

        // some variables for time step calculation
        double sumfactorin = 0;
        double sumfactorout = 0;

        if (this->enableLocalTimeStepping())
        {
            LocalTimesteppingData& localData = this->timeStepData_[globalIdxI];
            for (int i = 0; i < 2*dim; i++)
            {
                if (localData.faceTargetDt[i] < this->accumulatedDt_ + this->dtThreshold_)
                {
                    localData.faceFluxes[i] = 0.0;
                }
            }
        }

        // run through all intersections with neighbors and boundary
        const auto isEndIt = problem().gridView().iend(element);
        for (auto isIt = problem().gridView().ibegin(element); isIt != isEndIt; ++isIt)
        {
            const auto& intersection = *isIt;

            int indexInInside = intersection.indexInInside();

            // handle interior face
            if (intersection.neighbor())
            {
                if (enableMPFA && (intersection.outside().level() != element.level()))
                    getMpfaFlux(entries, timestepFlux, isIt, cellDataI);
                else
                    this->getFlux(entries, timestepFlux, intersection, cellDataI);
            }

            //     Boundary Face
            if (intersection.boundary())
            {
                this->getFluxOnBoundary(entries, timestepFlux, intersection, cellDataI);
            }

            if (this->enableLocalTimeStepping())
            {
                LocalTimesteppingData& localData = this->timeStepData_[globalIdxI];

                if (localData.faceTargetDt[indexInInside] < this->accumulatedDt_ + this->dtThreshold_)
                {
                    localData.faceFluxes[indexInInside] += entries;
                }
            }
            else
            {
            // add to update vector
            updateVec[wCompIdx][globalIdxI] += entries[wCompIdx];
            updateVec[nCompIdx][globalIdxI] += entries[nCompIdx];
            }

            // for time step calculation
            sumfactorin += timestepFlux[0];
            sumfactorout += timestepFlux[1];

        }// end all intersections

        if (this->enableLocalTimeStepping())
        {
            LocalTimesteppingData& localData = this->timeStepData_[globalIdxI];
            for (int i=0; i < 2*dim; i++)
            {
                updateVec[wCompIdx][globalIdxI] += localData.faceFluxes[i][wCompIdx];
                updateVec[nCompIdx][globalIdxI] += localData.faceFluxes[i][nCompIdx];
            }
        }

        /***********     Handle source term     ***************/
        PrimaryVariables q(NAN);
        problem().source(q, element);
        updateVec[wCompIdx][globalIdxI] += q[contiWEqIdx];
        updateVec[nCompIdx][globalIdxI] += q[contiNEqIdx];

        // account for porosity
        using std::max;
        sumfactorin = max(sumfactorin,sumfactorout)
                        / problem().spatialParams().porosity(element);

        //calculate time step
        if (this->enableLocalTimeStepping())
        {
            this->timeStepData_[globalIdxI].dt = 1./sumfactorin;
            if ( 1./sumfactorin < dt)
            {
                dt = 1./sumfactorin;
                restrictingCell= globalIdxI;
            }
        }
        else
        {
            if ( 1./sumfactorin < dt)
            {
                dt = 1./sumfactorin;
                restrictingCell= globalIdxI;
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
    dt = problem_.gridView().comm().min(dt);
#endif

    if(impet)
    {
        Dune::dinfo << "Timestep restricted by CellIdx " << restrictingCell << " leads to dt = "
                <<dt * getParam<Scalar>("Impet.CFLFactor")<< std::endl;
        if(this->averagedFaces_ != 0)
            Dune::dinfo  << " Averageing done for " << this->averagedFaces_ << " faces. "<< std::endl;
    }
    return;
}

/*!
 * \brief Compute flux over an irregular interface using a \a mpfa method
 *
 * A mpfa l-method is applied to calculate fluxes near hanging nodes, using:
 * \f[
      - \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha}
        \left( \sum_k \tau_{2k} p^t_{\alpha,k} + \varrho_{\alpha} \sum_k \tau_{2k} \mathbf{g}^T \mathbf{x}_{k} \right)
                \sum_{\kappa} X^{\kappa}_{\alpha}
    \f]
 *
 * We provide two options: Calculating the flux expressed by twice the flux
 * through the one unique interaction region on the hanging node if one
 * halfedge is stored (eg on boundaries). Or using the second interaction
 * region covering neighboring cells.
 *
 * \param fluxEntries The flux entries, mass influx from cell \f$j\f$ to \f$i\f$.
 * \param timestepFlux flow velocities for timestep estimation
 * \param intersectionIterator Iterator of the intersection between cell I and J
 * \param cellDataI The cell data for cell \f$i\f$
 */
template<class TypeTag>
void FV2dTransport2P2CAdaptive<TypeTag>::getMpfaFlux(Dune::FieldVector<Scalar, 2>& fluxEntries, Dune::FieldVector<Scalar, 2>& timestepFlux,
                                        const IntersectionIterator& intersectionIterator, CellData& cellDataI)
{
    fluxEntries = 0.;
    timestepFlux = 0.;
    // cell information
    auto elementI = intersectionIterator->inside();
    int globalIdxI = problem().variables().index(elementI);

    // get position
    const GlobalPosition globalPos = elementI.geometry().center();
    const GlobalPosition& gravity_ = problem().gravity();
    // cell volume, assume linear map here
    Scalar volume = elementI.geometry().volume();

    // get values of cell I
    Scalar pressI = problem().pressureModel().pressure(globalIdxI);
    Scalar pcI = cellDataI.capillaryPressure();

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

    PhaseVector potential(0.);

    // access neighbor
    auto neighbor = intersectionIterator->outside();
    int globalIdxJ = problem().variables().index(neighbor);
    CellData& cellDataJ = problem().variables().cellData(globalIdxJ);

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

    double pressJ = problem().pressureModel().pressure(globalIdxJ);
    Scalar pcJ = cellDataJ.capillaryPressure();

    // determine potentials for upwind
        /** get geometrical Info, transmissibility matrix */
        GlobalPosition globalPos3(0.);
        int globalIdx3=-1;
        TransmissivityMatrix T(0.);
        IntersectionIterator additionalIsIt = intersectionIterator;
        TransmissivityMatrix additionalT(0.);

        int halfedgesStored
            = problem().variables().getMpfaData(*intersectionIterator, additionalIsIt, T, additionalT,
                                                globalPos3, globalIdx3);
        if (halfedgesStored == 0)
            halfedgesStored = problem().pressureModel().computeTransmissibilities(intersectionIterator,additionalIsIt,
                                                                                  T,additionalT, globalPos3, globalIdx3 );

        // acess cell 3 and prepare mpfa
        Scalar press3 = problem().pressureModel().pressure(globalIdx3);
        CellData& cellData3 = problem().variables().cellData(globalIdx3);
        Scalar pc3 = cellData3.capillaryPressure();
        Scalar temp1 = globalPos * gravity_;
        Scalar temp2 = globalPosNeighbor * gravity_;
        Scalar temp3 = globalPos3 * gravity_;
        if(pressureType==pw)
        {
            potential[wPhaseIdx] += (pressI-temp1*densityW_mean) * T[2]
                            +(pressJ-temp2*densityW_mean) * T[0]
                            +(press3- temp3*densityW_mean) * T[1];
            potential[nPhaseIdx] += (pressI+pcI-temp1*densityNW_mean) * T[2]
                            +(pressJ+pcJ-temp2*densityNW_mean) * T[0]
                            +(press3+pc3- temp3*densityNW_mean) * T[1];
            // second half edge, if there is one
            if(halfedgesStored == 2)
            {
                int AdditionalIdx = problem().variables().index(additionalIsIt->outside());
                CellData& cellDataAdditional = problem().variables().cellData(AdditionalIdx);
                potential[wPhaseIdx] += (pressI-temp1*densityW_mean) * additionalT[2]
                                +(pressJ-temp2*densityW_mean) * additionalT[0]
                                +(problem().pressureModel().pressure(AdditionalIdx)
                                    -(additionalIsIt->outside().geometry().center()*gravity_*densityW_mean)
                                  ) * additionalT[1];
                potential[nPhaseIdx] += (pressI+pcI-temp1*densityNW_mean) * additionalT[2]
                                +(pressJ+pcJ-temp2*densityNW_mean) * additionalT[0]
                                +(problem().pressureModel().pressure(AdditionalIdx)
                                    + cellDataAdditional.capillaryPressure()
                                    -(additionalIsIt->outside().geometry().center()*gravity_*densityNW_mean)
                                  ) * additionalT[1];
            }
        }
        else if(pressureType==pn)
        {
            potential[wPhaseIdx] += (pressI-pcI-temp1*densityW_mean) * T[2]
                          + (pressJ-pcJ-temp2*densityW_mean) * T[0]
                          + (press3-pc3- temp3*densityW_mean) * T[1];
            potential[nPhaseIdx] += (pressI-temp1*densityNW_mean) * T[2]
                          + (pressJ-temp2*densityNW_mean) * T[0]
                          + (press3-temp3*densityNW_mean) * T[1];
            // second half edge, if there is one
            if(halfedgesStored == 2)
            {
                int AdditionalIdx = problem().variables().index(additionalIsIt->outside());
                CellData& cellDataAdditional = problem().variables().cellData(AdditionalIdx);

                potential[wPhaseIdx] += (pressI-pcI-temp1*densityW_mean) * additionalT[2]
                                +(pressJ-pcJ-temp2*densityW_mean) * additionalT[0]
                                +(problem().pressureModel().pressure(AdditionalIdx)
                                   - cellDataAdditional.capillaryPressure()
                                   -(additionalIsIt->outside().geometry().center()*gravity_*densityW_mean)
                                  ) * additionalT[1];
                potential[nPhaseIdx] += (pressI-temp1*densityNW_mean) * additionalT[2]
                                +(pressJ-temp2*densityNW_mean) * additionalT[0]
                                +(problem().pressureModel().pressure(AdditionalIdx)
                                    -(additionalIsIt->outside().geometry().center()*gravity_*densityNW_mean))
                                    * additionalT[1];
            }
        }
    //end of mpfa specific stuff

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

        if(!this->impet_ or !this->restrictFluxInTransport_) // perform a strict uwpind scheme
        {
            if(potential[phaseIdx] > 0.)
            {
                lambda[phaseIdx] = cellDataI.mobility(phaseIdx);
                if(elementI.level()>neighbor.level())
                    cellDataI.setUpwindCell(intersectionIterator->indexInInside(), contiEqIdx, true);
            }
            else if(potential[phaseIdx] < 0.)
            {
                lambda[phaseIdx] = cellDataJ.mobility(phaseIdx);
                if(elementI.level()>neighbor.level())
                    cellDataI.setUpwindCell(intersectionIterator->indexInInside(), contiEqIdx, false);
            }
            else
            {
                doUpwinding[phaseIdx] = false;
                if(elementI.level()>neighbor.level())
                cellDataI.setUpwindCell(intersectionIterator->indexInInside(), contiEqIdx, false);
                else
                cellDataJ.setUpwindCell(intersectionIterator->indexInOutside(), contiEqIdx, false);
            }
        }
        else // Transport after PE with check on flow direction
        {
            bool cellIwasUpwindCell;
            //get the information from smaller (higher level) cell, as its IS is unique
            if(elementI.level()>neighbor.level())
                cellIwasUpwindCell = cellDataI.isUpwindCell(intersectionIterator->indexInInside(), contiEqIdx);
            else // reverse neighbors information gathered
                cellIwasUpwindCell = !cellDataJ.isUpwindCell(intersectionIterator->indexInOutside(), contiEqIdx);

            if (potential[phaseIdx] > 0. && cellIwasUpwindCell)
                lambda[phaseIdx] = cellDataI.mobility(phaseIdx);
            else if (potential[phaseIdx] < 0. && !cellIwasUpwindCell)
                lambda[phaseIdx] = cellDataJ.mobility(phaseIdx);
            // potential direction does not coincide with that of P.E.
            else if(this->restrictFluxInTransport_ == 2)   // perform central averaging for all direction changes
                doUpwinding[phaseIdx] = false;
            else    // i.e. restrictFluxInTransport == 1
            {
                //check if harmonic weighting is necessary
                // check if outflow induce neglected (i.e. mob=0) phase flux
                if (potential[phaseIdx] > 0. && Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(cellDataJ.mobility(phaseIdx), 0.0, 1.0e-30))
                    lambda[phaseIdx] = cellDataI.mobility(phaseIdx);
                // check if inflow induce neglected phase flux
                else if (potential[phaseIdx] < 0. && Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(cellDataI.mobility(phaseIdx), 0.0, 1.0e-30))
                    lambda[phaseIdx] = cellDataJ.mobility(phaseIdx);
                else
                    doUpwinding[phaseIdx] = false;
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
                fluxEntries[wCompIdx] -= potential[phaseIdx] / volume
                        * harmonicMean(cellDataI.massFraction(phaseIdx, wCompIdx) * cellDataI.mobility(phaseIdx) * cellDataI.density(phaseIdx),
                                cellDataJ.massFraction(phaseIdx, wCompIdx) * cellDataJ.mobility(phaseIdx) * cellDataJ.density(phaseIdx));
                fluxEntries[nCompIdx] -= potential[phaseIdx] / volume
                        * harmonicMean(cellDataI.massFraction(phaseIdx, nCompIdx) * cellDataI.mobility(phaseIdx) * cellDataI.density(phaseIdx),
                                cellDataJ.massFraction(phaseIdx, nCompIdx) * cellDataJ.mobility(phaseIdx) * cellDataJ.density(phaseIdx));
                // c) timestep control
                // for timestep control : influx
                using std::max;
                timestepFlux[0] += max(0.,
                        - potential[phaseIdx] / volume
                          * harmonicMean(cellDataI.mobility(phaseIdx),cellDataJ.mobility(phaseIdx)));
                // outflux
                timestepFlux[1] += max(0.,
                        potential[phaseIdx]  / volume
                        * harmonicMean(cellDataI.mobility(phaseIdx),cellDataJ.mobility(phaseIdx))/SmobI[phaseIdx]);

                //d) output (only for one side)
                this->averagedFaces_++;
                #if DUNE_MINIMAL_DEBUG_LEVEL < 3
                // verbose (only for one side)
                if(globalIdxI > globalIdxJ)
                    Dune::dinfo << "harmonicMean flux of phase" << phaseIdx <<" used from cell" << globalIdxI<< " into " << globalIdxJ
                    << " ; TE upwind I = "<< cellDataI.isUpwindCell(intersectionIterator->indexInInside(), contiEqIdx)
                    << " but pot = "<< potential[phaseIdx] <<  " \n";
                #endif

                //e) stop further standard calculations
                potential[phaseIdx] = 0;
            }
        }
    }

    // calculate and standardized velocity
    using std::max;
    double velocityJIw = max((-lambda[wPhaseIdx] * potential[wPhaseIdx]) / volume, 0.0);
    double velocityIJw = max(( lambda[wPhaseIdx] * potential[wPhaseIdx]) / volume, 0.0);
    double velocityJIn = max((-lambda[nPhaseIdx] * potential[nPhaseIdx]) / volume, 0.0);
    double velocityIJn = max(( lambda[nPhaseIdx] * potential[nPhaseIdx]) / volume, 0.0);

    // for timestep control
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
    return;
}

} // end namespace Dumux
#endif
