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
 * \brief The problem class for the coupling of a non-isothermal two-component Stokes
 *        and a non-isothermal two-phase two-component Darcy model.
 *
 * The problem class for the coupling of a non-isothermal two-component Stokes (stokes2cn)
 * and a non-isothermal two-phase two-component Darcy model (2p2cni).
 * It uses the 2p2cniCoupling model and the Stokes2cnicoupling model and provides
 * the problem specifications for common parameters of the two submodels.
 * The initial and boundary conditions of the submodels are specified in the two subproblems,
 * 2p2cnisubproblem.hh and stokes2cnisubproblem.hh, which are accessible via the coupled problem.
 */

#ifndef DUMUX_2CNISTOKES2P2CNIPROBLEM_HH
#define DUMUX_2CNISTOKES2P2CNIPROBLEM_HH

#include <dune/common/float_cmp.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/multidomaingrid.hh>
#include <dune/grid/io/file/dgfparser.hh>

#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>
#include <dumux/multidomain/common/multidomainproblem.hh>
#include <dumux/multidomain/2cnistokes2p2cni/2cnistokes2p2cnilocaloperator.hh>
#include <dumux/multidomain/2cnistokes2p2cni/2cnistokes2p2cnipropertydefaults.hh>

#include <dumux/linear/seqsolverbackend.hh>
#ifdef HAVE_PARDISO
#include <dumux/linear/pardisobackend.hh>
#endif // HAVE_PARDISO

#include "2cnistokes2p2cnispatialparams.hh"
#include "stokes2cnisubproblem.hh"
#include "2p2cnisubproblem.hh"

namespace Dumux
{
template <class TypeTag>
class TwoCNIStokesTwoPTwoCNIProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoCNIStokesTwoPTwoCNIProblem, INHERITS_FROM(TwoCNIStokesTwoPTwoCNI));

// Set the grid type
#ifdef HAVE_UG
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, Grid, Dune::UGGrid<2>);
#elif HAVE_ALUGRID || HAVE_DUNE_ALUGRID
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);
#else
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, Grid, Dune::YaspGrid<2>);
#endif

// Set the global problem
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, Problem, TwoCNIStokesTwoPTwoCNIProblem<TypeTag>);

// Set the local coupling operator
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, MultiDomainCouplingLocalOperator, 
              Dumux::TwoCNIStokesTwoPTwoCNILocalOperator<TypeTag>);

// Set the two sub-problems of the global problem
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, SubDomain1TypeTag, TTAG(Stokes2cniSubProblem));
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, SubDomain2TypeTag, TTAG(TwoPTwoCNISubProblem));

// Set the global problem in the context of the two sub-problems
SET_TYPE_PROP(Stokes2cniSubProblem, MultiDomainTypeTag, TTAG(TwoCNIStokesTwoPTwoCNIProblem));
SET_TYPE_PROP(TwoPTwoCNISubProblem, MultiDomainTypeTag, TTAG(TwoCNIStokesTwoPTwoCNIProblem));

// Set the other sub-problem for each of the two sub-problems
SET_TYPE_PROP(Stokes2cniSubProblem, OtherSubDomainTypeTag, TTAG(TwoPTwoCNISubProblem));
SET_TYPE_PROP(TwoPTwoCNISubProblem, OtherSubDomainTypeTag, TTAG(Stokes2cniSubProblem));

// Set the spatial parameters used for the problems
SET_TYPE_PROP(Stokes2cniSubProblem, SpatialParams, Dumux::TwoCNIStokesTwoPTwoCNISpatialParams<TypeTag>);
SET_TYPE_PROP(TwoPTwoCNISubProblem, SpatialParams, Dumux::TwoCNIStokesTwoPTwoCNISpatialParams<TypeTag>);

// Set the fluid system
SET_PROP(TwoCNIStokesTwoPTwoCNIProblem, FluidSystem)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dumux::FluidSystems::H2OAir<Scalar,
        Dumux::H2O<Scalar>,
        /*useComplexrelations=*/true> type;
};

#ifdef HAVE_PARDISO
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, LinearSolver, PardisoBackend<TypeTag>);
#else
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, LinearSolver, SuperLUBackend<TypeTag>);
#endif // HAVE_PARDISO
}

/*!
 * \ingroup ImplicitTestProblems
 * \ingroup TwoPTwoCNIStokesTwoCNIModel
 * \brief The problem class for the coupling of a non-isothermal two-component Stokes
 *        and a non-isothermal two-phase two-component Darcy model.
 *
 * The problem class for the coupling of a non-isothermal two-component Stokes (stokes2cn)
 * and a non-isothermal two-phase two-component Darcy model (2p2cni).
 * It uses the 2p2cniCoupling model and the Stokes2cnicoupling model and provides
 * the problem specifications for common parameters of the two submodels.
 * The initial and boundary conditions of the submodels are specified in the two subproblems,
 * 2p2cnisubproblem.hh and stokes2cnisubproblem.hh, which are accessible via the coupled problem.
 */
template <class TypeTag = TTAG(TwoCNIStokesTwoPTwoCNIProblem) >
class TwoCNIStokesTwoPTwoCNIProblem : public MultiDomainProblem<TypeTag>
{
    typedef TwoCNIStokesTwoPTwoCNIProblem<TypeTag> ThisType;
    typedef MultiDomainProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) Stokes2cniTypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) TwoPTwoCNITypeTag;

    typedef typename GET_PROP_TYPE(Stokes2cniTypeTag, Problem) Stokes2cniSubProblem;
    typedef typename GET_PROP_TYPE(TwoPTwoCNITypeTag, Problem) TwoPTwoCNISubProblem;

    typedef typename GET_PROP_TYPE(Stokes2cniTypeTag, GridView) Stokes2cniGridView;
    typedef typename GET_PROP_TYPE(TwoPTwoCNITypeTag, GridView) TwoPTwoCNIGridView;

    typedef typename GET_PROP_TYPE(Stokes2cniTypeTag, PrimaryVariables) Stokes2cniPrimaryVariables;
    typedef typename GET_PROP_TYPE(TwoPTwoCNITypeTag, PrimaryVariables) TwoPTwoCNIPrimaryVariables;

    typedef typename GET_PROP_TYPE(Stokes2cniTypeTag, ElementSolutionVector) ElementSolutionVector1;
    typedef typename GET_PROP_TYPE(TwoPTwoCNITypeTag, ElementSolutionVector) ElementSolutionVector2;

    typedef typename GET_PROP_TYPE(Stokes2cniTypeTag, ElementVolumeVariables) ElementVolumeVariables1;
    typedef typename GET_PROP_TYPE(TwoPTwoCNITypeTag, ElementVolumeVariables) ElementVolumeVariables2;

    typedef typename GET_PROP_TYPE(Stokes2cniTypeTag, FluxVariables) BoundaryVariables1;

    typedef typename GET_PROP_TYPE(Stokes2cniTypeTag, FVElementGeometry) FVElementGeometry1;
    typedef typename GET_PROP_TYPE(TwoPTwoCNITypeTag, FVElementGeometry) FVElementGeometry2;

    typedef typename GET_PROP_TYPE(TypeTag, Grid) HostGrid;
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGrid) MDGrid;
    typedef typename MDGrid::LeafGridView MDGridView;
    typedef typename MDGrid::SubDomainGrid SDGrid;

    typedef typename MDGrid::Traits::template Codim<0>::Entity MDElement;
    typedef typename MDGrid::Traits::template Codim<0>::EntityPointer MDElementPointer;
    typedef typename Stokes2cniGridView::template Codim<0>::Entity SDElement1;
    typedef typename TwoPTwoCNIGridView::template Codim<0>::Entity SDElement2;
    typedef typename SDGrid::Traits::template Codim<0>::EntityPointer SDElementPointer;

    typedef typename GET_PROP_TYPE(Stokes2cniTypeTag, Indices) Stokes2cniIndices;
    typedef typename GET_PROP_TYPE(TwoPTwoCNITypeTag, Indices) TwoPTwoCNIIndices;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { dim = Stokes2cniGridView::dimension };
    enum { // indices in the Stokes domain
        momentumXIdx1 = Stokes2cniIndices::momentumXIdx, //!< Index of the x-component of the momentum balance
        momentumYIdx1 = Stokes2cniIndices::momentumYIdx, //!< Index of the y-component of the momentum balance
        momentumZIdx1 = Stokes2cniIndices::momentumZIdx, //!< Index of the z-component of the momentum balance
        massBalanceIdx1 = Stokes2cniIndices::massBalanceIdx, //!< Index of the mass balance
        transportEqIdx1 = Stokes2cniIndices::transportEqIdx, //!< Index of the transport equation
        energyEqIdx1 = Stokes2cniIndices::energyEqIdx //!< Index of the transport equation
    };
    enum { // indices of the PVs in the Darcy domain
        massBalanceIdx2 = TwoPTwoCNIIndices::pressureIdx,
        switchIdx2 = TwoPTwoCNIIndices::switchIdx,
        temperatureIdx2 = TwoPTwoCNIIndices::temperatureIdx
    };
    enum { // indices of the balance equations
        contiTotalMassIdx2 = TwoPTwoCNIIndices::contiNEqIdx,
        contiWEqIdx2 = TwoPTwoCNIIndices::contiWEqIdx,
        energyEqIdx2 = TwoPTwoCNIIndices::energyEqIdx
    };
    enum { transportCompIdx1 = Stokes2cniIndices::transportCompIdx };
    enum {
        wCompIdx2 = TwoPTwoCNIIndices::wCompIdx,
        nCompIdx2 = TwoPTwoCNIIndices::nCompIdx
    };
    enum { phaseIdx = GET_PROP_VALUE(Stokes2cniTypeTag, PhaseIdx) };
    enum {
        numEq1 = GET_PROP_VALUE(Stokes2cniTypeTag, NumEq),
        numEq2 = GET_PROP_VALUE(TwoPTwoCNITypeTag, NumEq)
    };

    typedef Dune::FieldVector<Scalar, dim> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> FieldVector;

    typedef typename MDGrid::template Codim<0>::LeafIterator ElementIterator;
    typedef typename MDGrid::LeafSubDomainInterfaceIterator SDInterfaceIterator;

public:
    /*!
     * \brief The problem for the coupling of Stokes and Darcy flow
     *
     * \param mdGrid The multidomain grid
     * \param timeManager The time manager
     */
    TwoCNIStokesTwoPTwoCNIProblem(MDGrid &mdGrid,
    							  TimeManager &timeManager)
        : ParentType(mdGrid, timeManager)
    {
        interfacePosY_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePosY);
        noDarcyX_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyX);
        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLength);
        initializationTime_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, InitTime);
        dtInit_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, DtInitial);

        // define output options
        freqRestart_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqRestart);
        freqOutput_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqOutput);
        freqFluxOutput_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqFluxOutput);
        freqVaporFluxOutput_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqVaporFluxOutput);

        stokes2cni_ = this->sdID1();
        twoPtwoCNI_ = this->sdID2();

        initializeGrid();

        // initialize the tables of the fluid system
        FluidSystem::init(/*tempMin=*/273.15, /*tempMax=*/373.15, /*numTemp=*/200,
                          /*pMin=*/1e3, /*pMax=*/2e5, /*numP=*/200);

        if (initializationTime_ > 0.0)
            this->timeManager().startNextEpisode(initializationTime_);
        else
            this->timeManager().startNextEpisode(episodeLength_);
    }

    /*!
     * \brief The destructor
     */
    ~TwoCNIStokesTwoPTwoCNIProblem()
    {
        fluxFile_.close();
    };

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

        std::cout << "Writing flux data at interface\n";
        if (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(this->timeManager().time(), 0.0, 1.0e-30))
        {
            fluxFile_.open("fluxes.out");
            fluxFile_ << "Time;"
                      << "TotalWaterVaporFluxFF;" // from residuum
                      << "AdvWaterVaporFluxFF;" // from gradients (imprecise)
                      << "DiffWaterVaporFluxFF;" // from gradients (imprecise)
                      << "EnergyFluxFF;"
                      << "TotalWaterComponentFluxPM;"
                      << "WaterFluxLiquidPhasePM;"
                      << "WaterFluxGasPhasePM;"
                      << "EnergyFluxPM"
                      << std::endl;
            counter_ = 1;
        }
        else
            fluxFile_.open("fluxes.out", std::ios_base::app);
    }

    /*!
     * \brief Initialization of the grids
     *
     * This function splits the multidomain grid in the two
     * individual subdomain grids and takes care of parallelization.
     */
    void initializeGrid()
    {
        MDGrid& mdGrid = this->mdGrid();
        mdGrid.startSubDomainMarking();

        // subdivide grid in two subdomains
        ElementIterator eEndIt = mdGrid.template leafend<0>();
        for (ElementIterator eIt = mdGrid.template leafbegin<0>();
             eIt != eEndIt; ++eIt)
        {
            // this is required for parallelization
            // checks if element is within a partition
            if (eIt->partitionType() != Dune::InteriorEntity)
                continue;

            GlobalPosition globalPos = eIt->geometry().center();

            if (globalPos[1] > interfacePosY_)
                mdGrid.addToSubDomain(stokes2cni_,*eIt);
            else
                if(globalPos[0] > noDarcyX_)
                    mdGrid.addToSubDomain(twoPtwoCNI_,*eIt);
        }
        mdGrid.preUpdateSubDomains();
        mdGrid.updateSubDomains();
        mdGrid.postUpdateSubDomains();

        gridinfo(this->sdGrid1());
        gridinfo(this->sdGrid2());
    }

    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    {
        // call the postTimeStep function of the subproblems
        this->sdProblem1().postTimeStep();
        this->sdProblem2().postTimeStep();

        if (shouldWriteFluxFile() || shouldWriteVaporFlux())
        {
            counter_ = this->sdProblem1().currentVTKFileNumber() + 1;

            calculateFirstInterfaceFluxes();
            calculateSecondInterfaceFluxes();
        }
    }

    /*!
     * \brief Called when the end of an simulation episode is reached.
     *
     * Typically a new episode should be started in this method.
     */
    void episodeEnd()
    {
        this->timeManager().startNextEpisode(episodeLength_);
        if (this->timeManager().time() <= initializationTime_ + dtInit_)
        {
            std::cout << "setting timeStepSize to " << dtInit_ << std::endl;
            this->timeManager().setTimeStepSize(dtInit_);
        }
    }

    /*!
     * \brief Calculates fluxes and coupling terms at the interface
     *        for the Stokes model.
     * 
     * Flux output files are created and the summarized flux is
     * written to a file.
     */
    void calculateFirstInterfaceFluxes()
    {
        const MDGrid& mdGrid = this->mdGrid();
        ElementVolumeVariables1 elemVolVarsPrev1, elemVolVarsCur1;
        Scalar totalWaterVaporFlux = 0.;
        Scalar advectiveWaterVaporFlux = 0.;
        Scalar diffusiveWaterVaporFlux = 0.;
        Scalar energyFlux = 0.;

        // count number of elements to determine number of interface nodes
        int numElements = 0;
        const SDInterfaceIterator endIfIt = mdGrid.leafSubDomainInterfaceEnd(stokes2cni_, twoPtwoCNI_);
        for (SDInterfaceIterator ifIt =
                 mdGrid.leafSubDomainInterfaceBegin(stokes2cni_, twoPtwoCNI_); ifIt != endIfIt; ++ifIt)
            numElements++;

        const int numInterfaceVertices = numElements + 1;
        std::vector<InterfaceFluxes<numEq1> > outputVector(numInterfaceVertices); // vector for the output of the fluxes
        FVElementGeometry1 fvGeometry1;
        int interfaceVertIdx = -1;

        // loop over the element faces on the interface
        for (SDInterfaceIterator ifIt =
                 mdGrid.leafSubDomainInterfaceBegin(stokes2cni_, twoPtwoCNI_); ifIt != endIfIt; ++ifIt)
        {
            const int firstFaceIdx = ifIt->indexInFirstCell();
            const MDElementPointer mdElementPointer1 = ifIt->firstCell(); // ATTENTION!!!
            const MDElement& mdElement1 = *mdElementPointer1;     // Entity pointer has to be copied before.
            const SDElementPointer sdElementPointer1 = this->sdElementPointer1(mdElement1);
            const SDElement1& sdElement1 = *sdElementPointer1;
            fvGeometry1.update(this->sdGridView1(), sdElement1);

            const Dune::ReferenceElement<typename MDGrid::ctype,dim>& referenceElement1 =
                Dune::ReferenceElements<typename MDGrid::ctype,dim>::general(mdElement1.type());
            const int numVerticesOfFace = referenceElement1.size(firstFaceIdx, 1, dim);

            // evaluate residual of the sub model without boundary conditions (stabilization is removed)
            // the element volume variables are updated here
            this->localResidual1().evalNoBoundary(sdElement1, fvGeometry1,
                                                  elemVolVarsPrev1, elemVolVarsCur1);

            for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++)
            {
                const int vertInElem1 = referenceElement1.subEntity(firstFaceIdx, 1, nodeInFace, dim);
                const FieldVector& vertexGlobal = mdElement1.geometry().corner(vertInElem1);
                const unsigned firstGlobalIdx = this->mdVertexMapper().map(stokes2cni_, mdElement1, vertInElem1, dim);
                const ElementSolutionVector1& firstVertexResidual = this->localResidual1().residual();

                // loop over all interface vertices to check if vertex id is already in stack
                bool existing = false;
                for (int interfaceVertex=0; interfaceVertex < numInterfaceVertices; ++interfaceVertex)
                {
                    if (firstGlobalIdx == outputVector[interfaceVertex].vIdxGlobal)
                    {
                        existing = true;
                        interfaceVertIdx = interfaceVertex;
                        break;
                    }
                }

                if (!existing)
                    interfaceVertIdx++;

                if (shouldWriteFluxFile()) // compute only if required
                {
                    outputVector[interfaceVertIdx].interfaceVertex = interfaceVertIdx;
                    outputVector[interfaceVertIdx].vIdxGlobal = firstGlobalIdx;
                    outputVector[interfaceVertIdx].xCoord = vertexGlobal[0];
                    outputVector[interfaceVertIdx].yCoord = vertexGlobal[1];
                    outputVector[interfaceVertIdx].count += 1;
                    for (int eqIdx=0; eqIdx < numEq1; ++eqIdx)
                        outputVector[interfaceVertIdx].residual[eqIdx] +=
                            firstVertexResidual[vertInElem1][eqIdx];
                }

                // compute summarized fluxes for output
                if (shouldWriteVaporFlux())
                {
                    int boundaryFaceIdx =
                        fvGeometry1.boundaryFaceIndex(firstFaceIdx, nodeInFace);

                    const BoundaryVariables1 boundaryVars1(this->sdProblem1(),
                                                           sdElement1,
                                                           fvGeometry1,
                                                           boundaryFaceIdx,
                                                           elemVolVarsCur1,
                                                           /*onBoundary=*/true);

                    totalWaterVaporFlux += firstVertexResidual[vertInElem1][transportEqIdx1];
                    advectiveWaterVaporFlux += computeAdvectiveVaporFluxes1(elemVolVarsCur1, boundaryVars1, vertInElem1);
                    diffusiveWaterVaporFlux += computeDiffusiveVaporFluxes1(elemVolVarsCur1, boundaryVars1, vertInElem1);
                    energyFlux += firstVertexResidual[vertInElem1][energyEqIdx1];
                }
            }
        } // end loop over element faces on interface

        if (shouldWriteFluxFile())
        {
            std::cout << "Writing flux file\n";
            char outputname[20];
            sprintf(outputname, "%s%05d%s","fluxesFF_", counter_,".out");
            std::ofstream outfile(outputname, std::ios_base::out);
            outfile << "XCoordFF;"
                    << "TotalMassFluxFF;"
                    << "TotalComponentMassFluxFF;"
                    << "TotalEnergyFluxFF"
                    << std::endl;
            for (int interfaceVertIdx=0; interfaceVertIdx < numInterfaceVertices; interfaceVertIdx++)
            {
                if (outputVector[interfaceVertIdx].count > 2)
                    std::cerr << "too often at one node!!";

                if (outputVector[interfaceVertIdx].count==2)
                    outfile << outputVector[interfaceVertIdx].xCoord << ";"
                            << outputVector[interfaceVertIdx].residual[massBalanceIdx1] << ";" // total mass flux
                            << outputVector[interfaceVertIdx].residual[transportEqIdx1] << ";" // total flux of component
                            << outputVector[interfaceVertIdx].residual[energyEqIdx1] // total flux of heat
                            << std::endl;
            }
            outfile.close();
        }
        if (shouldWriteVaporFlux())
            fluxFile_ << this->timeManager().time() + this->timeManager().timeStepSize() << ";"
                      << totalWaterVaporFlux << ";"
                      << advectiveWaterVaporFlux << ";"
                      << diffusiveWaterVaporFlux << ";"
                      << energyFlux << ";";
    }

    /*!
     * \brief Calculates fluxes and coupling terms at the interface
     *        for the Darcy model.
     *
     * Flux output files are created and the summarized flux is written
     * to a file.
     */
    void calculateSecondInterfaceFluxes()
    {
        const MDGrid& mdGrid = this->mdGrid();
        ElementVolumeVariables2 elemVolVarsPrev2, elemVolVarsCur2;

        Scalar totalWaterComponentFlux = 0.;
        Scalar energyFlux = 0.;
        Scalar waterFluxGasPhase = 0.;

        // count number of elements to determine number of interface nodes
        int numElements = 0;
        const SDInterfaceIterator endIfIt = mdGrid.leafSubDomainInterfaceEnd(stokes2cni_, twoPtwoCNI_);
        for (SDInterfaceIterator ifIt =
                 mdGrid.leafSubDomainInterfaceBegin(stokes2cni_, twoPtwoCNI_); ifIt != endIfIt; ++ifIt)
            numElements++;

        const int numInterfaceVertices = numElements + 1;
        std::vector<InterfaceFluxes<numEq2> > outputVector(numInterfaceVertices); // vector for the output of the fluxes
        FVElementGeometry2 fvGeometry2;
        int interfaceVertIdx = -1;

        // loop over the element faces on the interface
        for (SDInterfaceIterator ifIt =
                 mdGrid.leafSubDomainInterfaceBegin(stokes2cni_, twoPtwoCNI_); ifIt != endIfIt; ++ifIt)
        {
            const int secondFaceIdx = ifIt->indexInSecondCell();
            const MDElementPointer mdElementPointer2 = ifIt->secondCell(); // ATTENTION!!!
            const MDElement& mdElement2 = *mdElementPointer2;     // Entity pointer has to be copied before.
            const SDElementPointer sdElementPointer2 = this->sdElementPointer2(mdElement2);
            const SDElement2& sdElement2 = *sdElementPointer2;
            fvGeometry2.update(this->sdGridView2(), sdElement2);

            const Dune::ReferenceElement<typename MDGrid::ctype,dim>& referenceElement2 =
                Dune::ReferenceElements<typename MDGrid::ctype,dim>::general(mdElement2.type());
            const int numVerticesOfFace = referenceElement2.size(secondFaceIdx, 1, dim);

            // evaluate residual of the sub model without boundary conditions
            this->localResidual2().evalNoBoundary(sdElement2, fvGeometry2,
                                                  elemVolVarsPrev2, elemVolVarsCur2);
            // evaluate the vapor fluxes within each phase
            this->localResidual2().evalPhaseFluxes();

            for (int nodeInFace = 0; nodeInFace < numVerticesOfFace; nodeInFace++)
            {
                const int vertInElem2 = referenceElement2.subEntity(secondFaceIdx, 1, nodeInFace, dim);
                const FieldVector& vertexGlobal = mdElement2.geometry().corner(vertInElem2);
                const unsigned secondGlobalIdx = this->mdVertexMapper().map(twoPtwoCNI_, mdElement2, vertInElem2, dim);
                const ElementSolutionVector2& secondVertexResidual = this->localResidual2().residual();

                bool existing = false;
                // loop over all interface vertices to check if vertex id is already in stack
                for (int interfaceVertex=0; interfaceVertex < numInterfaceVertices; ++interfaceVertex)
                {
                    if (secondGlobalIdx == outputVector[interfaceVertex].vIdxGlobal)
                    {
                        existing = true;
                        interfaceVertIdx = interfaceVertex;
                        break;
                    }
                }

                if (!existing)
                    interfaceVertIdx++;

                if (shouldWriteFluxFile())
                {
                    outputVector[interfaceVertIdx].interfaceVertex = interfaceVertIdx;
                    outputVector[interfaceVertIdx].vIdxGlobal = secondGlobalIdx;
                    outputVector[interfaceVertIdx].xCoord = vertexGlobal[0];
                    outputVector[interfaceVertIdx].yCoord = vertexGlobal[1];
                    for (int eqIdx=0; eqIdx < numEq2; ++eqIdx)
                        outputVector[interfaceVertIdx].residual[eqIdx] += secondVertexResidual[vertInElem2][eqIdx];
                    outputVector[interfaceVertIdx].count += 1;
                }
                if (shouldWriteVaporFlux())
                {
                    if (!existing) // add phase storage only once per vertex
                        waterFluxGasPhase +=
                            this->localResidual2().evalPhaseStorageDerivative(vertInElem2);

                    totalWaterComponentFlux += secondVertexResidual[vertInElem2][contiWEqIdx2];
                    waterFluxGasPhase += this->localResidual2().elementFluxes(vertInElem2);
                    energyFlux += secondVertexResidual[vertInElem2][energyEqIdx2];
                }
            }
        }

        if (shouldWriteFluxFile())
        {
            char outputname[20];
            sprintf(outputname, "%s%05d%s","fluxesPM_", counter_,".out");
            std::ofstream outfile(outputname, std::ios_base::out);
            outfile << "XCoordPM;"
                    << "TotalMassFluxPM;"
                    << "TotalComponentMassFluxPM;"
                    << "TotalEnergyFluxPM"
                    << std::endl;

            for (int interfaceVertIdx=0; interfaceVertIdx < numInterfaceVertices; interfaceVertIdx++)
            {
                if (outputVector[interfaceVertIdx].count > 2)
                    std::cerr << "too often at one node!!";

                if (outputVector[interfaceVertIdx].count==2)
                    outfile << outputVector[interfaceVertIdx].xCoord << ";"
                            << outputVector[interfaceVertIdx].residual[contiTotalMassIdx2] << ";" // total mass flux
                            << outputVector[interfaceVertIdx].residual[contiWEqIdx2] << ";" // total flux of component
                            << outputVector[interfaceVertIdx].residual[energyEqIdx2] // total heat flux
                            << std::endl;
            }
            outfile.close();
        }
        if (shouldWriteVaporFlux())
        {
            Scalar waterFluxLiquidPhase = totalWaterComponentFlux - waterFluxGasPhase;
            fluxFile_ << totalWaterComponentFlux << ";"
                      << waterFluxLiquidPhase << ";"
                      << waterFluxGasPhase << ";"
                      << energyFlux
                      << std::endl;
        }
    }

    /*!
     * \brief Returns the advective vapor fluxes
     *
     * The phaseIdx and transportCompIdx1 are predefined
     *
     * \todo boundaryVars1 violates naming convention
     *
     * \param elemVolVars1 All volume variables for the element
     * \param boundaryVars1 Flux variables
     * \param vertInElem1 Vertex index for the inside element
     */
    Scalar computeAdvectiveVaporFluxes1(const ElementVolumeVariables1& elemVolVars1,
                                        const BoundaryVariables1& boundaryVars1,
                                        int vertInElem1)
    {
        Scalar advFlux = elemVolVars1[vertInElem1].density() *
            elemVolVars1[vertInElem1].massFraction(transportCompIdx1) *
            boundaryVars1.normalVelocity();
        return advFlux;
    }

    /*!
     * \brief Returns the diffusive vapor fluxes
     *
     * The transportCompIdx1 is predefined
     *
     * \todo boundaryVars1 violates naming convention
     *
     * \param elemVolVars1 All volume variables for the element
     * \param boundaryVars1 Flux variables
     * \param vertInElem1 Vertex index for the inside elements
     */
    Scalar computeDiffusiveVaporFluxes1(const ElementVolumeVariables1& elemVolVars1,
                                        const BoundaryVariables1& boundaryVars1,
                                        int vertInElem1)
    {
        Scalar diffFlux = (boundaryVars1.moleFractionGrad(transportCompIdx1) *
                          boundaryVars1.face().normal) *
                          boundaryVars1.diffusionCoeff(transportCompIdx1) *
                          boundaryVars1.molarDensity() *
                          FluidSystem::molarMass(transportCompIdx1);
        return diffFlux;
    }

    /*!
     * \brief Returns true if a restart file should be written to
     *        disk.
     *
     * The default behavior is to write one restart file every 5 time
     * steps. This file is intended to be overwritten by the
     * implementation.
     */
    bool shouldWriteRestartFile() const
    {
        if ((this->timeManager().timeStepIndex() > 0 &&
             (this->timeManager().timeStepIndex() % freqRestart_ == 0))
            // also write a restart file at the end of each episode
            || this->timeManager().episodeWillBeOver())
            return true;
        return false;
    }

    /*!
     * \brief Returns true if the current solution should be written to
     *        disk (i.e. as a VTK file)
     *
     * The default behavior is to write out the solution for
     * every time step. This function is intended to be overwritten by the
     * implementation.
     */
    bool shouldWriteOutput() const
    {
        if (this->timeManager().timeStepIndex() % freqOutput_ == 0
            || this->timeManager().episodeWillBeOver())
            return true;
        return false;
    }

    /*!
     * \brief Returns true if a file with the fluxes across the
     *        free-flow -- porous-medium interface should be
     *        written to disk.
     */
    bool shouldWriteFluxFile() const
    {
        if (this->timeManager().timeStepIndex() % freqFluxOutput_ == 0
            || this->timeManager().episodeWillBeOver())
            return true;
        return false;
    }

    /*!
     * \brief Returns true if the summarized vapor fluxes
     *        across the free-flow -- porous-medium interface,
     *        representing the evaporation rate (related to the
     *        interface area), should be written.
     */
    bool shouldWriteVaporFlux() const
    {
        if (this->timeManager().timeStepIndex() % freqVaporFluxOutput_ == 0
            || this->timeManager().episodeWillBeOver())
            return true;
        return false;
    }

    /*!
     * \brief Returns a pointer to the Stokes problem
     */
    Stokes2cniSubProblem& stokes2cniProblem()
    { return this->sdProblem1(); }
    const Stokes2cniSubProblem& stokes2cniProblem() const
    { return this->sdProblem1(); }

    /*!
     * \brief Returns a pointer to the Darcy problem
     */
    TwoPTwoCNISubProblem& twoPtwoCNIProblem()
    { return this->sdProblem2(); }
    const TwoPTwoCNISubProblem& twoPtwoCNIProblem() const
    { return this->sdProblem2(); }

private:
    typename MDGrid::SubDomainIndex stokes2cni_;
    typename MDGrid::SubDomainIndex twoPtwoCNI_;

    unsigned counter_;
    unsigned freqRestart_;
    unsigned freqOutput_;
    unsigned freqFluxOutput_;
    unsigned freqVaporFluxOutput_;

    Scalar interfacePosY_;
    Scalar noDarcyX_;
    Scalar episodeLength_;
    Scalar initializationTime_;
    Scalar dtInit_;

    template <int numEq>
    struct InterfaceFluxes
    {
        unsigned count;
        unsigned interfaceVertex;
        unsigned vIdxGlobal;
        Scalar xCoord;
        Scalar yCoord;
        Dune::FieldVector<Scalar, numEq> residual;

        InterfaceFluxes()
        {
            count = 0;
            interfaceVertex = 0;
            vIdxGlobal = 0;
            xCoord = 0.0;
            yCoord = 0.0;
            residual = 0.0;
        }
    };
    std::ofstream fluxFile_;
};

} // namespace Dumux

#endif // DUMUX_2CNISTOKES2P2CNIPROBLEM_HH
