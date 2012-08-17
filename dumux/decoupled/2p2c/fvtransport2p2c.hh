// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff, Benjamin Faigle                     *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
#ifndef DUMUX_FVTRANSPORT2P2C_HH
#define DUMUX_FVTRANSPORT2P2C_HH

#include <dune/grid/common/gridenums.hh>
#include <dumux/decoupled/2p2c/2p2cproperties.hh>
#include <dumux/common/math.hh>

/**
 * @file
 * @brief  Finite Volume discretization of the component transport equation
 * @author Markus Wolff, Jochen Fritz, Benjamin Faigle
 */

namespace Dumux
{
//! Miscible Transport step in a Finite Volume discretization
/*! \ingroup multiphase
 *  The finite volume model for the solution of the transport equation for compositional
 *  two-phase flow.
 *  \f[
      \frac{\partial C^\kappa}{\partial t} = - \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{alpha} \bf{v}_{\alpha}\right) + q^{\kappa},
 *  \f]
 *  where \f$ \bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 *  \f$ p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability, \f$ \lambda_{\alpha} \f$ the phase mobility,
 *  \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and \f$ C^{\kappa} \f$ the total Component concentration.
 *  The whole flux contribution for each cell is subdivided into a storage term, a flux term and a source term.
 *  Corresponding functions (<tt>getFlux()</tt> and <tt>getFluxOnBoundary()</tt>) are provided,
 *  internal sources are directly treated.
 *
 *
 *  \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVTransport2P2C
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename SpatialParams::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    typedef typename GET_PROP_TYPE(TypeTag, TransportSolutionType) TransportSolutionType;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld,
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNW,
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx,
        contiWEqIdx=Indices::contiWEqIdx, contiNEqIdx=Indices::contiNEqIdx,
        NumPhases = GET_PROP_VALUE(TypeTag, NumPhases)
    };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> DimMatrix;
    typedef Dune::FieldVector<Scalar, NumPhases> PhaseVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    //! Acess function for the current problem
    Problem& problem()
    {return problem_;};

public:
    virtual void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet = false);

    void updateTransportedQuantity(TransportSolutionType& updateVector);

    // Function which calculates the flux update
    void getFlux(Dune::FieldVector<Scalar, 2>&, Dune::FieldVector<Scalar, 2>&,
            const Intersection&, CellData&);

    // Function which calculates the boundary flux update
    void getFluxOnBoundary(Dune::FieldVector<Scalar, 2>&, Dune::FieldVector<Scalar, 2>&,
                            const Intersection&, const CellData&);

    void evalBoundary(GlobalPosition,const Intersection&,FluidState &, PhaseVector &);

    //! Set the initial values before the first pressure equation
    /*!
     * This method is called before first pressure equation is solved from Dumux::IMPET.
     */
    void initialize()
    {
        totalConcentration_[wCompIdx].resize(problem_.gridView().size(0));
        totalConcentration_[nCompIdx].resize(problem_.gridView().size(0));
    };

    //! \brief Write data files
     /*  \param writer applied VTK-writer */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        typedef typename GET_PROP(TypeTag, SolutionTypes)::ScalarSolution ScalarSolutionType;
        int size = problem_.gridView().size(0);
        ScalarSolutionType *totalC1PV = writer.allocateManagedBuffer(size);
        ScalarSolutionType *totalC2PV = writer.allocateManagedBuffer(size);
        *totalC1PV = this->totalConcentration_[wCompIdx];
        *totalC2PV = this->totalConcentration_[nCompIdx];
        writer.attachCellData(*totalC1PV, "total Concentration w-Comp");
        writer.attachCellData(*totalC2PV, "total Concentration n-Comp");
    }

    //! Function needed for restart option of the transport model: Write out
    void serializeEntity(std::ostream &outstream, const Element &element)
    {
        int globalIdx = problem().variables().index(element);
        outstream << totalConcentration_[wCompIdx][globalIdx]
                  << "  " << totalConcentration_[nCompIdx][globalIdx];
    }
    //! Function needed for restart option of the transport model: Read in
    void deserializeEntity(std::istream &instream, const Element &element)
    {
        int globalIdx = problem().variables().index(element);
        instream >>  totalConcentration_[wCompIdx][globalIdx]
                 >> totalConcentration_[nCompIdx][globalIdx];
    }

    /*! \name Access functions for protected variables  */
    //@{
    //! Return the vector of the transported quantity
    /*! For an immiscible IMPES scheme, this is the saturation. For Miscible simulations, however,
     *  the total concentration of all components is transported.
     */
    TransportSolutionType& transportedQuantity() DUNE_DEPRECATED
    {
        return totalConcentration_;
    }
    //! \copydoc transportedQuantity()
    void getTransportedQuantity(TransportSolutionType& transportedQuantity)
    {
        // resize update vector and set to zero
        transportedQuantity.resize(GET_PROP_VALUE(TypeTag, NumComponents));
        transportedQuantity[wCompIdx].resize(problem_.gridView().size(0));
        transportedQuantity[nCompIdx].resize(problem_.gridView().size(0));

        transportedQuantity = totalConcentration_;
    }
    //! \copydoc transportedQuantity()
    Scalar& totalConcentration(int compIdx, int globalIdx)
    {
        return totalConcentration_[compIdx][globalIdx][0];
    }
    //@}
    //! Constructs a FVTransport2P2C object
    /*!
     * Currently, the miscible transport scheme can not be applied with a global pressure / total velocity
     * formulation.
     *
     * \param problem a problem class object
     */
    FVTransport2P2C(Problem& problem) :
        totalConcentration_(0.),problem_(problem), switchNormals(false)
    {
        totalConcentration_.resize(GET_PROP_VALUE(TypeTag, NumComponents));
        totalConcentration_[wCompIdx].resize(problem.gridView().size(0));
        totalConcentration_[nCompIdx].resize(problem.gridView().size(0));

        restrictFluxInTransport_ = GET_PARAM_FROM_GROUP(TypeTag,int, Impet, RestrictFluxInTransport);
    }

    virtual ~FVTransport2P2C()
    {
    }

protected:
    TransportSolutionType totalConcentration_;
    Problem& problem_;
    bool impet_;
    int averagedFaces_;

    static const int pressureType = GET_PROP_VALUE(TypeTag, PressureFormulation); //!< gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
    int restrictFluxInTransport_;
    bool switchNormals;
};

//! \brief Calculate the update vector and determine timestep size
/*!
 *  This method calculates the update vector \f$ u \f$ of the discretized equation
 *  \f[
       C^{\kappa , new} = C^{\kappa , old} + u,
 *  \f]
 *  where \f$ u = \sum_{element faces} \boldsymbol{v}_{\alpha} * \varrho_{\alpha} * X^{\kappa}_{\alpha} * \boldsymbol{n} * A_{element face} \f$,
 *  \f$ \boldsymbol{n} \f$ is the face normal and \f$ A_{element face} \f$ is the face area.
 *
 *  In addition to the \a update vector, the recommended time step size \a dt is calculated
 *  employing a CFL condition.
 *  This method = old concentrationUpdate()
 *
 *  \param t Current simulation time \f$\mathrm{[s]}\f$
 *  \param[out] dt Time step size \f$\mathrm{[s]}\f$
 *  \param[out] updateVec Update vector, or update estimate for secants, resp. Here in \f$\mathrm{[kg/m^3]}\f$
 *  \param impet Flag that determines if it is a real impet step or an update estimate for volume derivatives
 */
template<class TypeTag>
void FVTransport2P2C<TypeTag>::update(const Scalar t, Scalar& dt,
		TransportSolutionType& updateVec, bool impet)
{
    // initialize dt very large
    dt = 1E100;
    // store if we do update Estimate for flux functions
    impet_ = impet;
    averagedFaces_ = 0;

    // resize update vector and set to zero
    updateVec.resize(GET_PROP_VALUE(TypeTag, NumComponents));
    updateVec[wCompIdx].resize(problem_.gridView().size(0));
    updateVec[nCompIdx].resize(problem_.gridView().size(0));
    updateVec[wCompIdx] = 0;
    updateVec[nCompIdx] = 0;

    // Cell which restricts time step size
    int restrictingCell = -1;

    PhaseVector entries(0.), timestepFlux(0.);
    // compute update vector
    ElementIterator eItEnd = problem().gridView().template end<0> ();
    for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // get cell infos
        int globalIdxI = problem().variables().index(*eIt);
        CellData& cellDataI = problem().variables().cellData(globalIdxI);

        // some variables for time step calculation
        double sumfactorin = 0;
        double sumfactorout = 0;

        // run through all intersections with neighbors and boundary
        IntersectionIterator isItEnd = problem().gridView().iend(*eIt);
        for (IntersectionIterator isIt = problem().gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
        {

            /****** interior face   *****************/
            if (isIt->neighbor())
                getFlux(entries, timestepFlux, *isIt, cellDataI);

            /******  Boundary Face   *****************/
            if (isIt->boundary())
                getFluxOnBoundary(entries, timestepFlux, *isIt, cellDataI);

            // add to update vector
            updateVec[wCompIdx][globalIdxI] += entries[wCompIdx];
            updateVec[nCompIdx][globalIdxI] += entries[nCompIdx];

            // for time step calculation
            sumfactorin += timestepFlux[0];
            sumfactorout += timestepFlux[1];

        }// end all intersections

        /***********     Handle source term     ***************/
        PrimaryVariables q(NAN);
        problem().source(q, *eIt);
        updateVec[wCompIdx][globalIdxI] += q[contiWEqIdx];
        updateVec[nCompIdx][globalIdxI] += q[contiNEqIdx];

        // account for porosity in fluxes for time-step
        sumfactorin = std::max(sumfactorin,sumfactorout)
                        / problem().spatialParams().porosity(*eIt);

        if ( 1./sumfactorin < dt)
        {
            dt = 1./sumfactorin;
            restrictingCell= globalIdxI;
        }
    } // end grid traversal
    if(impet)
    {
        Dune::dinfo << "Timestep restricted by CellIdx " << restrictingCell << " leads to dt = "
                <<dt * GET_PARAM_FROM_GROUP(TypeTag, Scalar, Impet, CFLFactor)<< std::endl;
    	if(averagedFaces_ != 0)
            Dune::dinfo  << " Averageing done for " << averagedFaces_ << " faces. "<< std::endl;
    }
    return;
}
/*	Updates the transported quantity once an update is calculated.
 *  This method updates both, the internal transport solution vector and the entries in the cellData.
 *  \param updateVec Update vector, or update estimate for secants, resp. Here in \f$\mathrm{[kg/m^3]}\f$
 *
 */
template<class TypeTag>
void FVTransport2P2C<TypeTag>::updateTransportedQuantity(TransportSolutionType& updateVector)
{
    Scalar dt = problem().timeManager().timeStepSize();
    // loop thorugh all elements
    for (int i = 0; i< problem().gridView().size(0); i++)
    {
        CellData& cellDataI = problem().variables().cellData(i);
        for(int compIdx = 0; compIdx < GET_PROP_VALUE(TypeTag, NumComponents); compIdx++)
        {
            totalConcentration_[compIdx][i] += (updateVector[compIdx][i]*=dt);
            cellDataI.setMassConcentration(compIdx, totalConcentration_[compIdx][i]);
        }
    }
}

//! Get flux at an interface between two cells
/** The flux thorugh \f$ \gamma{ij} \f$  is calculated according to the underlying pressure field,
 * calculated by the pressure model.
 *  \f[ - A_{\gamma_{ij}}  \cdot \mathbf{u} \cdot (\mathbf{n}_{\gamma_{ij}} \cdot \mathbf{u})\cdot \mathbf{K}
      \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha} \sum_{\kappa} X^{\kappa}_{\alpha}
    \left( \frac{p_{\alpha,j}^t - p^{t}_{\alpha,i}}{\Delta x} + \varrho_{\alpha} \mathbf{g}\right) \f]
 * Here, \f$ \mathbf{u} \f$ is the normalized vector connecting the cell centers, and \f$ \mathbf{n}_{\gamma_{ij}} \f$
 * represents the normal of the face \f$ \gamma{ij} \f$. Due to the nature of the Primay Variable, the (volume-)specific
 * total mass concentration, this represents a mass flux per cell volume.
 * \param fluxEntries The flux entries, mass influx from cell \f$j\f$ to \f$i\f$.
 * \param timestepFlux flow velocities for timestep estimation
 * \param intersection The intersection
 * \param cellDataI The cell data for cell \f$i\f$
 */
template<class TypeTag>
void FVTransport2P2C<TypeTag>::getFlux(Dune::FieldVector<Scalar, 2>& fluxEntries,
                                        Dune::FieldVector<Scalar, 2>& timestepFlux,
                                        const Intersection& intersection,
                                        CellData& cellDataI)
{
    fluxEntries = 0.;
    timestepFlux = 0.;
    // cell information
    ElementPointer elementPtrI= intersection.inside();
    int globalIdxI = problem().variables().index(*elementPtrI);

    // get position
    const GlobalPosition globalPos = elementPtrI->geometry().center();
    const GlobalPosition& gravity_ = problem().gravity();
    // cell volume, assume linear map here
    Scalar volume = elementPtrI->geometry().volume();

    // get values of cell I
    Scalar pressI = problem().pressureModel().pressure(globalIdxI);
    Scalar pcI = cellDataI.capillaryPressure();
    DimMatrix K_I(problem().spatialParams().intrinsicPermeability(*elementPtrI));

    PhaseVector SmobI(0.);
    SmobI[wPhaseIdx] = std::max((cellDataI.saturation(wPhaseIdx)
                            - problem().spatialParams().materialLawParams(*elementPtrI).Swr())
                            , 1e-2);
    SmobI[nPhaseIdx] = std::max((cellDataI.saturation(nPhaseIdx)
                                - problem().spatialParams().materialLawParams(*elementPtrI).Snr())
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
    ElementPointer neighborPtr = intersection.outside();
    int globalIdxJ = problem().variables().index(*neighborPtr);
    CellData& cellDataJ = problem().variables().cellData(globalIdxJ);

    // neighbor cell center in global coordinates
    const GlobalPosition& globalPosNeighbor = neighborPtr->geometry().center();

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

    // compute mean permeability
    DimMatrix meanK_(0.);
    Dumux::harmonicMeanMatrix(meanK_,
            K_I,
            problem().spatialParams().intrinsicPermeability(*neighborPtr));
    Dune::FieldVector<Scalar,dim> K(0);
    meanK_.umv(unitDistVec,K);

    // determine potentials for upwind
    switch (pressureType)
    {
    case pw:
    {
        potential[wPhaseIdx] = (K * unitOuterNormal) * (pressI - pressJ) / (dist);
        potential[nPhaseIdx] = (K * unitOuterNormal) * (pressI - pressJ + pcI - pcJ) / (dist);
        break;
    }
    case pn:
    {
        potential[wPhaseIdx] = (K * unitOuterNormal) * (pressI - pressJ - pcI + pcJ) / (dist);
        potential[nPhaseIdx] = (K * unitOuterNormal) * (pressI - pressJ) / (dist);
        break;
    }
    }
    // add gravity term
    potential[nPhaseIdx] +=  (K * gravity_)  * (unitOuterNormal * unitDistVec) * densityNW_mean;
    potential[wPhaseIdx] +=  (K * gravity_)  * (unitOuterNormal * unitDistVec) * densityW_mean;

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

        if(!impet_ or restrictFluxInTransport_==0) // perform a strict uwpind scheme
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
            else if(restrictFluxInTransport_ == 2)   // perform central averageing for all direction changes
                doUpwinding[phaseIdx] = false;
            else    // i.e. restrictFluxInTransport == 1
            {
               //check if harmonic weithing is necessary
                if (!wasUpwindCell && (cellDataJ.mobility(phaseIdx) != 0.   // check if outflow induce neglected (i.e. mob=0) phase flux
                       or (cellDataI.wasRefined() && cellDataJ.wasRefined() && elementPtrI->father() == neighborPtr->father())))
                    lambda[phaseIdx] = cellDataI.mobility(phaseIdx);
                else if (wasUpwindCell && (cellDataI.mobility(phaseIdx) != 0. // check if inflow induce neglected phase flux
                        or (cellDataI.wasRefined() && cellDataJ.wasRefined() && elementPtrI->father() == neighborPtr->father())))
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

            //b) perform harmonic averageing
            fluxEntries[wCompIdx] -= potential[phaseIdx] * faceArea / volume
                    * harmonicMean(cellDataI.massFraction(phaseIdx, wCompIdx) * cellDataI.mobility(phaseIdx) * cellDataI.density(phaseIdx),
                            cellDataJ.massFraction(phaseIdx, wCompIdx) * cellDataJ.mobility(phaseIdx) * cellDataJ.density(phaseIdx));
            fluxEntries[nCompIdx] -= potential[phaseIdx] * faceArea / volume
                    * harmonicMean(cellDataI.massFraction(phaseIdx, nCompIdx) * cellDataI.mobility(phaseIdx) * cellDataI.density(phaseIdx),
                            cellDataJ.massFraction(phaseIdx, nCompIdx) * cellDataJ.mobility(phaseIdx) * cellDataJ.density(phaseIdx));
            // c) timestep control
            // for timestep control : influx
            timestepFlux[0] += std::max(0.,
                    - potential[phaseIdx] * faceArea / volume
                      * harmonicMean(cellDataI.mobility(phaseIdx),cellDataJ.mobility(phaseIdx)));
            // outflux
            timestepFlux[1] += std::max(0.,
                    potential[phaseIdx] * faceArea / volume
                    * harmonicMean(cellDataI.mobility(phaseIdx),cellDataJ.mobility(phaseIdx))/SmobI[phaseIdx]);

            //d) output (only for one side)
            averagedFaces_++;
            #if DUNE_MINIMAL_DEBUG_LEVEL < 3
            // verbose (only for one side)
            if(globalIdxI > globalIdxJ)
                Dune::dinfo << "harmonicMean flux of phase" << phaseIdx <<" used from cell" << globalIdxI<< " into " << globalIdxJ
                << " ; TE upwind I = "<< cellDataI.isUpwindCell(intersection.indexInInside(), contiEqIdx) << " but pot = "<< potential[phaseIdx] <<  " \n";
            #endif

            //e) stop further standard calculations
            potential[phaseIdx] = 0;
        }
    }

    // calculate and standardized velocity
    double velocityJIw = std::max((-lambda[wPhaseIdx] * potential[wPhaseIdx]) * faceArea / volume, 0.0);
    double velocityIJw = std::max(( lambda[wPhaseIdx] * potential[wPhaseIdx]) * faceArea / volume, 0.0);
    double velocityJIn = std::max((-lambda[nPhaseIdx] * potential[nPhaseIdx]) * faceArea / volume, 0.0);
    double velocityIJn = std::max(( lambda[nPhaseIdx] * potential[nPhaseIdx]) * faceArea / volume, 0.0);

    // for timestep control : influx
    timestepFlux[0] += velocityJIw + velocityJIn;

    double foutw = velocityIJw/SmobI[wPhaseIdx];
    double foutn = velocityIJn/SmobI[nPhaseIdx];
    if (std::isnan(foutw) || std::isinf(foutw) || foutw < 0) foutw = 0;
    if (std::isnan(foutn) || std::isinf(foutn) || foutn < 0) foutn = 0;
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
//! Get flux on Boundary
/** The flux thorugh \f$ \gamma{ij} \f$  is calculated according to the underlying pressure field,
 * calculated by the pressure model.
 *  \f[ - A_{\gamma_{ij}}  \cdot \mathbf{u} \cdot (\mathbf{n}_{\gamma_{ij}} \cdot \mathbf{u})\cdot \mathbf{K}
      \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha} \sum_{\kappa} X^{\kappa}_{\alpha}
    \left( \frac{p_{\alpha,j}^t - p^{t}_{\alpha,i}}{\Delta x} + \varrho_{\alpha} \mathbf{g}\right) \f]
 * Here, \f$ \mathbf{u} \f$ is the normalized vector connecting the cell centers, and \f$ \mathbf{n}_{\gamma_{ij}} \f$
 * represents the normal of the face \f$ \gamma{ij} \f$.
 * \param fluxEntries The flux entries, mass influx from cell \f$j\f$ to \f$i\f$.
 * \param timestepFlux flow velocities for timestep estimation
 * \param intersection The intersection
 * \param cellDataI The cell data for cell \f$i\f$
 */
template<class TypeTag>
void FVTransport2P2C<TypeTag>::getFluxOnBoundary(Dune::FieldVector<Scalar, 2>& fluxEntries,
                                                    Dune::FieldVector<Scalar, 2>& timestepFlux,
                                                    const Intersection& intersection,
                                                    const CellData& cellDataI)
{
    // cell information
    ElementPointer elementPtrI= intersection.inside();
    int globalIdxI = problem().variables().index(*elementPtrI);

    // get position
    const GlobalPosition globalPos = elementPtrI->geometry().center();

    // cell volume, assume linear map here
    Scalar volume = elementPtrI->geometry().volume();
    const GlobalPosition& gravity_ = problem().gravity();
    // get values of cell I
    Scalar pressI = problem().pressureModel().pressure(globalIdxI);
    Scalar pcI = cellDataI.capillaryPressure();
    DimMatrix K_I(problem().spatialParams().intrinsicPermeability(*elementPtrI));

    Scalar SwmobI = std::max((cellDataI.saturation(wPhaseIdx)
                            - problem().spatialParams().materialLawParams(*elementPtrI).Swr())
                            , 1e-2);
    Scalar SnmobI = std::max((cellDataI.saturation(nPhaseIdx)
                                - problem().spatialParams().materialLawParams(*elementPtrI).Snr())
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
        if(GET_PROP_VALUE(TypeTag, EnableCapillarity))
            pcBound = BCfluidState.capillaryPressure();
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
                potential[wPhaseIdx] = (K * unitOuterNormal) *
                        (pressI - pressBound[wPhaseIdx]) / dist;
                potential[nPhaseIdx] = (K * unitOuterNormal) *
                        (pressI + pcI - pressBound[wPhaseIdx] - pcBound)
                        / dist;
                break;
            }
            case pn:
            {
                potential[wPhaseIdx] = (K * unitOuterNormal) *
                        (pressI - pcI - pressBound[nPhaseIdx] + pcBound)
                        / dist;
                potential[nPhaseIdx] = (K * unitOuterNormal) *
                        (pressI - pressBound[nPhaseIdx]) / dist;
                break;
            }
        }
        potential[wPhaseIdx] += (K * gravity_)  * (unitOuterNormal * unitDistVec) * densityW_mean;;
        potential[nPhaseIdx] += (K * gravity_)  * (unitOuterNormal * unitDistVec) * densityNW_mean;;

        // do upwinding for lambdas
        PhaseVector lambda(0.);
        if (potential[wPhaseIdx] >= 0.)
            lambda[wPhaseIdx] = cellDataI.mobility(wPhaseIdx);
        else
            {
            if(GET_PROP_VALUE(TypeTag, BoundaryMobility)==Indices::satDependent)
                lambda[wPhaseIdx] = BCfluidState.saturation(wPhaseIdx) / viscosityWBound;
            else
                lambda[wPhaseIdx] = MaterialLaw::krw(
                        problem().spatialParams().materialLawParams(*elementPtrI), BCfluidState.saturation(wPhaseIdx))
                        / viscosityWBound;
            }
        if (potential[nPhaseIdx] >= 0.)
            lambda[nPhaseIdx] = cellDataI.mobility(nPhaseIdx);
        else
            {
            if(GET_PROP_VALUE(TypeTag, BoundaryMobility)==Indices::satDependent)
                lambda[nPhaseIdx] = BCfluidState.saturation(nPhaseIdx) / viscosityNWBound;
            else
                lambda[nPhaseIdx] = MaterialLaw::krn(
                        problem().spatialParams().materialLawParams(*elementPtrI), BCfluidState.saturation(wPhaseIdx))
                        / viscosityNWBound;
            }
        // calculate and standardized velocity
        double velocityJIw = std::max((-lambda[wPhaseIdx] * potential[wPhaseIdx]) * faceArea / volume, 0.0);
        double velocityIJw = std::max(( lambda[wPhaseIdx] * potential[wPhaseIdx]) * faceArea / volume, 0.0);
        double velocityJIn = std::max((-lambda[nPhaseIdx] * potential[nPhaseIdx]) * faceArea / volume, 0.0);
        double velocityIJn = std::max(( lambda[nPhaseIdx] * potential[nPhaseIdx]) * faceArea / volume, 0.0);

        // for timestep control
        timestepFlux[0] = velocityJIw + velocityJIn;

        double foutw = velocityIJw/SwmobI;
        double foutn = velocityIJn/SnmobI;
        if (std::isnan(foutw) || std::isinf(foutw) || foutw < 0) foutw = 0;
        if (std::isnan(foutn) || std::isinf(foutn) || foutn < 0) foutn = 0;
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

//! evaluate the boundary conditions
/*!
 *  As the transport primary variable in this formulation is the total component
 *  concentration, \f$ C^{\kappa} \f$ it seems natural that the boundary values
 *  are also total concentrations. However, as for the initial conditions, it is
 *  possible to define boundaries by means of a saturation. This choice determines
 *  which version of flash calculation is necessary to get to the composition at
 *  the boundary.
 *  \param globalPosFace Face of the current boundary
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
    CompositionalFlash<TypeTag> flashSolver;

    const ElementPointer eIt= intersection.inside();
    // read boundary values
    PrimaryVariables primaryVariablesOnBoundary(0.);
    problem().dirichlet(primaryVariablesOnBoundary, intersection);

    // read boundary type
    typename Indices::BoundaryFormulation bcType;
    problem().boundaryFormulation(bcType, intersection);
    if (bcType == Indices::saturation)
    {
        Scalar satBound = primaryVariablesOnBoundary[contiWEqIdx];
        if(GET_PROP_VALUE(TypeTag, EnableCapillarity))
        {
            Scalar pcBound = MaterialLaw::pC(problem().spatialParams().materialLawParams(*eIt),
                    satBound);
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
                problem().spatialParams().porosity(*eIt), problem().temperatureAtPos(globalPosFace));
    }
    else if (bcType == Indices::concentration)
    {
        // saturation and hence pc and hence corresponding pressure unknown
        pressBound[wPhaseIdx] = pressBound[nPhaseIdx] = primaryVariablesOnBoundary[Indices::pressureEqIdx];
        Scalar Z1Bound = primaryVariablesOnBoundary[contiWEqIdx];
        flashSolver.concentrationFlash2p2c(BCfluidState, Z1Bound, pressBound,
        	problem().spatialParams().porosity(*eIt), problem().temperatureAtPos(globalPosFace));

        if(GET_PROP_VALUE(TypeTag, EnableCapillarity))
        {
            Scalar pcBound = MaterialLaw::pC(problem().spatialParams().materialLawParams(*eIt),
                    BCfluidState.saturation(wPhaseIdx));
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
                flashSolver.concentrationFlash2p2c(BCfluidState, Z1Bound, pressBound,
                        problem().spatialParams().porosity(*eIt), problem().temperatureAtPos(globalPosFace));
                pcBound = MaterialLaw::pC(problem().spatialParams().materialLawParams(*eIt),
                        BCfluidState.saturation(wPhaseIdx));
                // TODO: get right criterion, do output for evaluation
                //converge criterion
                if (abs(oldPc-pcBound)<10)
                    iter = maxiter;
            }
        }
    }
    else
        DUNE_THROW(Dune::NotImplemented, "Boundary Formulation neither Concentration nor Saturation??");
}

}
#endif
