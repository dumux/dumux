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
 * \brief A test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 */
#ifndef DUMUX_EL2P_DECOUPLED_PROBLEM_SIMPLE_HH
#define DUMUX_EL2P_DECOUPLED_PROBLEM_SIMPLE_HH

#include <valgrind/callgrind.h>

#include <dumux/geomechanics/el2p/decoupled/simplecoupledproblem.hh> // copied from dumux-devel

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/common/timer.hh>

#include "decoupledelasticproblem_simple.hh"
#include "decoupled2pproblem_simple.hh"
#include <vector>
#include "time.h"

#include <dumux/implicit/adaptive/gridadaptproperties.hh>

namespace Dumux
{
template <class TypeTag>
class El2p_DecoupledProblem_Simple;

namespace Properties
{
NEW_TYPE_TAG(El2p_DecoupledProblem_Simple, INHERITS_FROM(SimpleCoupled, GridAdapt));

// Set the problem property
SET_TYPE_PROP(El2p_DecoupledProblem_Simple, Problem,
              Dumux::El2p_DecoupledProblem_Simple<TypeTag>);

SET_INT_PROP(El2p_DecoupledProblem_Simple, ImplicitMaxTimeStepDivisions, 5);

// Set the two sub-problems of the global problem
//SET_TYPE_PROP(El2p_DecoupledProblem_Simple, SubProblem1TypeTag, TTAG(BioTranspProblem));
//SET_TYPE_PROP(El2p_DecoupledProblem_Simple, SubProblem2TypeTag, TTAG(BioChemProblem));
SET_TYPE_PROP(El2p_DecoupledProblem_Simple, SubProblem1TypeTag, TTAG(TwoP_TestProblem_Simple));
SET_TYPE_PROP(El2p_DecoupledProblem_Simple, SubProblem2TypeTag, TTAG(El2P_TestProblem_Simple));

SET_PROP(TwoP_TestProblem_Simple, ParameterTree)
{private:
    typedef typename GET_PROP(TTAG(El2p_DecoupledProblem_Simple), ParameterTree) ParameterTree;
public:
    typedef typename ParameterTree::type type;

    static type &tree()
    { return ParameterTree::tree(); }

    static type &compileTimeParams()
    { return ParameterTree::compileTimeParams(); }

    static type &runTimeParams()
    { return ParameterTree::runTimeParams(); }

    static type &deprecatedRunTimeParams()
    { return ParameterTree::deprecatedRunTimeParams(); }

    static type &unusedNewRunTimeParams()
    { return ParameterTree::unusedNewRunTimeParams(); }

};

SET_PROP(El2P_TestProblem_Simple, ParameterTree)
{private:
    typedef typename GET_PROP(TTAG(El2p_DecoupledProblem_Simple), ParameterTree) ParameterTree;
public:
    typedef typename ParameterTree::type type;

    static type &tree()
    { return ParameterTree::tree(); }

    static type &compileTimeParams()
    { return ParameterTree::compileTimeParams(); }

    static type &runTimeParams()
    { return ParameterTree::runTimeParams(); }

    static type &deprecatedRunTimeParams()
    { return ParameterTree::deprecatedRunTimeParams(); }

    static type &unusedNewRunTimeParams()
    { return ParameterTree::unusedNewRunTimeParams(); }

};
}

template <class TypeTag>
class El2p_DecoupledProblem_Simple : public SimpleCoupledProblem<TypeTag>
{
    typedef SimpleCoupledProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // obtain the type tags of the subproblems
    typedef typename GET_PROP_TYPE(TypeTag, SubProblem1TypeTag) SubTypeTag1;
    typedef typename GET_PROP_TYPE(TypeTag, SubProblem2TypeTag) SubTypeTag2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, Problem) TwoPSpatialParams_Simple;
    typedef typename GET_PROP_TYPE(SubTypeTag2, Problem) El2P_TestProblem_Simple;

    typedef typename GET_PROP_TYPE(SubTypeTag1, FVElementGeometry) FVElementGeometryTranspProblem;
    typedef typename GET_PROP_TYPE(SubTypeTag1, ElementVolumeVariables) ElementVolumeVariablesTranspProblem;

    typedef typename GET_PROP_TYPE(SubTypeTag2, FVElementGeometry) FVElementGeometryMechanicsProblem;
    typedef typename GET_PROP_TYPE(SubTypeTag2, ElementVolumeVariables) ElementVolumeVariablesMechanicsProblem;

    typedef typename GET_PROP_TYPE(SubTypeTag1, PrimaryVariables) PrimaryVariablesTranspProblem;
    typedef typename GET_PROP_TYPE(SubTypeTag2, PrimaryVariables) PrimaryVariablesMechanicsProblem;

    typedef typename GET_PROP_TYPE(SubTypeTag1, Indices) Indices;
    enum {
        nPhaseIdx = Indices::nPhaseIdx,
        wPhaseIdx = Indices::wPhaseIdx
    };
    enum {
        // indices of the primary variables
            pressureIdx = Indices::pwIdx,
            saturationIdx = Indices::snIdx
    };

    typedef typename GET_PROP_TYPE(SubTypeTag1, GridView) GridView;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(SubTypeTag1, SolutionVector) SolutionVectorTranspProblem;
    typedef typename GET_PROP_TYPE(SubTypeTag2, SolutionVector) SolutionVectorMechanicsProblem;

    typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
    typedef Dune::BlockVector<Dune::FieldVector<double, dim> > VectorField;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dim> LocalPosition;

    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

    typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> LocalFiniteElementCache;
    typedef typename LocalFiniteElementCache::FiniteElementType LocalFiniteElement;

    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;

public:
    El2p_DecoupledProblem_Simple(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView),
    gridView_(gridView)
    {

        try
        {
            name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                                 std::string,
                                                 Problem,
                                                 Name);

            maxTimeStepSize_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, MaxTimeStepSize);

            previousTimeStepSize_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TimeManager.DtInitial);

            numIterations_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TimeManager.NumIterationsInit);

            episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLengthInit);

//             maxCouplingError_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
//                                                             Scalar,
//                                                             Problem,
//                                                             MaxCouplingError);
//             minPvValue_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
//                                                             Scalar,
//                                                             Problem,
//                                                             MinPvValue);
//             timeIntegrationIdx_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
//                                                                     Scalar,
//                                                                     Problem,
//                                                                     TimeIntegrationIdx);
//             injectionParameters_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Injection, InjectionParamFile);
//             dtmin_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, DtMin);
            }

    catch (Dumux::ParameterException &e) {
        std::cerr << e << ". Abort!\n";
        exit(1) ;
    }

        tInitEnd_ = GET_RUNTIME_PARAM(TypeTag, Scalar,TimeManager.TInitEnd);
        ParentType::timeManager_.startNextEpisode(tInitEnd_);
        // transfer the episode index to spatial parameters
        // (during intialization episode hydraulic different parameters might be applied)
        MechanicsProblem().spatialParams().setEpisode(ParentType::timeManager_.episodeIndex());

        output_ = true;
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    {
        return name_.c_str();
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    void init()
    {
        ParentType::init();

        // create the required scalar and vector fields
        unsigned numVert = gridView_.size(dim);
        unsigned numElements = gridView_.size(0);

//         std::cout << "numVert = " << numVert << std::endl;
//         std::cout << "numElements = " << numElements << std::endl;

        ScalarField rank;
        rank.resize(numElements);

        lastSol_.resize(numVert);

        std::vector<Scalar> zerosNumVert(numVert, 0.0);
        std::vector<int> zerosNumVertInt(numVert, 0);
        std::vector<Scalar> zerosNumElements(numElements, 0.0);

        pwVector_.resize(numElements);
        pnVector_.resize(numElements);
        pcVector_.resize(numElements);
        SwVector_.resize(numElements);
        SnVector_.resize(numElements);
        rhonVector_.resize(numElements);
        rhowVector_.resize(numElements);

        numScvOfNode_ = zerosNumVert;
        BNodeWiseAveraged_ = zerosNumVert;

        effPorosityVector_.resize(numElements);
        effPermeabilityVector_.resize(numElements);
        volumetricStrain_.resize(numElements);
        deltaVolumetricStrainOldIteration_.resize(numElements);
//         dUVector_.resize(numElements);

        FVElementGeometryTranspProblem fvGeometryTranspProblem;

        for (const auto& element : elements(gridView_))
        {
            if(element.partitionType() == Dune::InteriorEntity)
            {
                int eIdx = TranspProblem().model().elementMapper().index(element);

                fvGeometryTranspProblem.update(gridView_, element);

                int numScv = fvGeometryTranspProblem.numScv;

                effPorosityVector_[eIdx].resize(numScv);
                effPermeabilityVector_[eIdx].resize(numScv);

                volumetricStrain_[eIdx].resize(numScv);
                deltaVolumetricStrainOldIteration_[eIdx].resize(numScv);;

                pwVector_[eIdx].resize(numScv);
                pnVector_[eIdx].resize(numScv);
                pcVector_[eIdx].resize(numScv);
                SwVector_[eIdx].resize(numScv);
                SnVector_[eIdx].resize(numScv);
                rhonVector_[eIdx].resize(numScv);
                rhowVector_[eIdx].resize(numScv);

                for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
                {
                    volumetricStrain_[eIdx][scvIdx] = 0.0;
                    deltaVolumetricStrainOldIteration_[eIdx][scvIdx] = 0.0;
//                         DimVector zeroes(0.0);
//                         dUVector_[eIdx][scvIdx] = zeroes;
                }
            }
        }

        MechanicsProblem().setpw() = pwVector_;
        MechanicsProblem().setpn() = pnVector_;
        MechanicsProblem().setpc() = pcVector_;
        MechanicsProblem().setSw() = SwVector_;
        MechanicsProblem().setSn() = SnVector_;
        MechanicsProblem().setRhon() = rhonVector_;
        MechanicsProblem().setRhow() = rhowVector_;

        TranspProblem().setpWOldIteration_() = pwVector_;
        TranspProblem().setpNOldIteration_() = pnVector_;
        TranspProblem().setsWOldIteration_() = SwVector_;
        TranspProblem().setsNOldIteration_() = SnVector_;

        TranspProblem().setEffPorosity() = effPorosityVector_;
        TranspProblem().setEffPermeability() = effPermeabilityVector_;

        TranspProblem().setDeltaVolumetricStrainOldIteration() = deltaVolumetricStrainOldIteration_;
//         TranspProblem().setdU() = dUVector_;

        MechanicsProblem().setEffPorosity() = effPorosityVector_;
        MechanicsProblem().setEffPermeability() = effPermeabilityVector_;

//         this->setOutput(true);

//         for (int vIdxGlobal = 0; vIdxGlobal < numVert; ++vIdxGlobal)
//         {
//             lastSol_[vIdxGlobal] = 0.0;
//         }

/*
        ParentType::init();
        updateParamsMassBalance();

        std::cout << "************************************" << std::endl;*/

        std::cout << "initialisation done" << std::endl;
    }

    /*!
     * \brief Called by Dumux::TimeManager in order to do a time
     *        integration on the model.
     */
    void timeIntegration()
    {

        CALLGRIND_START_INSTRUMENTATION;
        Dune::Timer Timer;

//         const int maxFails =
//                 GET_PARAM_FROM_GROUP(TypeTag, int, Implicit, MaxTimeStepDivisions);
        const int maxFails = 5;
        for (int i = 0; i < maxFails; ++i) {
            try {

                MechanicsProblem().spatialParams().setEpisode(ParentType::timeManager_.episodeIndex());

//                 if (ParentType::timeManager_.time() > 172800 - eps_)
//                 {
//                     Scalar oldTimeStep = ParentType::timeManager_.timeStepSize();
//                     Scalar newTimeStepSize = TranspProblem().newtonController().suggestTimeStepSize(oldTimeStep);
//
//                     std::cout << "newTimeStepSize " << newTimeStepSize << std::endl;
//                     ParentType::timeManager_.setTimeStepSize(newTimeStepSize);
//
//                     TranspProblem().setEpisodeLength(newTimeStepSize);
//                     MechanicsProblem().setEpisodeLength(newTimeStepSize);
//                 }

                Scalar time = ParentType::timeManager_.time() + ParentType::timeManager_.timeStepSize();

                if((time - tInitEnd_) > eps_)
                {
                    numIterations_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TimeManager.NumIterations);
                }

                for (int iterator = 1; iterator < (numIterations_ + 1); iterator++)
                {
                    if( ((time - tInitEnd_) < eps_) && (iterator == 1) )
                    {
                        initializePoroPerm();
                    }

                    TranspProblem().setIteration() = iterator;
                    runTransport();

//                     std::cout << "norm of curSol is " << TranspProblem().model().curSol().two_norm() << std::endl;
                    SolutionVectorTranspProblem b(TranspProblem().model().curSol());
                    TranspProblem().setEvalGlobalResidual(true);
                    res_ = TranspProblem().model().globalResidual(b);
                    std::cout << "global residual = " << res_ << std::endl;
                    TranspProblem().setEvalGlobalResidual(false);

                    // set timestep size to 0 as Timemanager will do time_ += dt and
                    // this leads to wrong time in the output
                    // the preTimestep method takes care of setting the appropriate timestepsize
                    TranspProblem().timeManager().setTimeStepSize(0.0);


                    std::cout << "Ran transport" << std::endl;

                    updateParamsMomentumBalance();

                    std::cout << "updateParamsMomentumBalance" << std::endl;

                    // When the pore compressibility is used,there is no need to run the mechanics
                    if (!(GET_RUNTIME_PARAM(TypeTag, bool, PoreCompressibility.UsePoreCompressibility)))
                    {
                        runMechanics();

                        // et timestep size to 0 - see above
                        MechanicsProblem().timeManager().setTimeStepSize(0.0);

                        std::cout << "Ran mechanics" << std::endl;

                        updateParamsMassBalance();
                    }

                    std::cout
                    << "----------------------------- \n"
                    << "Iteration " << iterator << " done \n"
                    << "----------------------------- \n"
                    << std::endl;


//                     std::ofstream file2;
//                     std::string outname2 = "bOutput.txt";
//                     file2.open(outname2, std::ios::out);
//                     if (file2.fail())
//                         throw std::ios_base::failure(std::strerror(errno));
//         //             for (int i = 0; i < b.base().size(); i++)
//         //                 file2 << b.base()[i] << std::endl;
//                     printvector(file2, b, "b", "row", 5, 1, 5);

                    // After last iteration of the initialization
                    if( (std::abs(time - tInitEnd_) < eps_) && (iterator == numIterations_) )
                    {
                        TranspProblem().setCoupled(true);
                        MechanicsProblem().setCoupled(true);
                        TranspProblem().setInitializationRun(false);

//                         TranspProblem().setOutput(true);
//                         MechanicsProblem().setOutput(true);

                        // set old iteration values to initial values
                        TranspProblem().setpWOldIteration_() = pwVector_;
                        TranspProblem().setpNOldIteration_() = pnVector_;
                        TranspProblem().setsWOldIteration_() = SwVector_;
                        TranspProblem().setsNOldIteration_() = SnVector_;

//                         for (int scvIdx = 0; scvIdx < pwVector_[12873].size(); ++scvIdx)
//                         {
//                             std::cout << "pwVector_[12873][" << scvIdx << "] is " << pwVector_[12873][scvIdx] << std::endl;
//                         }
                    }

                    unsigned numVert = gridView_.size(0);

                    std::vector<Scalar> relDiffPw;
                    relDiffPw.resize(numVert);

                    for (int vIdxGlobal = 0; vIdxGlobal < numVert; ++vIdxGlobal)
                    {
                        relDiffPw[vIdxGlobal] = std::abs(TranspProblem().model().curSol()[vIdxGlobal][0] - lastSol_[vIdxGlobal][0]) /TranspProblem().model().curSol()[vIdxGlobal][0];
                    }

                    auto biggest = std::max_element(std::begin(relDiffPw), std::end(relDiffPw));

                    std::cout << "Max relDiffPw is " << *biggest
                              << " at position " << std::distance(std::begin(relDiffPw), biggest) << std::endl;
//                     std::cout << "Max relDiffPw is " << relDiffPw[2925] << " at position 2925" << std::endl;
//                     std::cout << "pwCurrent     is " << TranspProblem().model().curSol()[2925][0] << " at position 2925" << std::endl;
//                     std::cout << "pwLastIter    is " <<  lastSol_[2925][0] << " at position 2925" << std::endl;

                    // finish timestep if initialization is done and change between iteration sufficiently small
                    Scalar maxRelDiffPw = GET_RUNTIME_PARAM(TypeTag, Scalar, TimeManager.MaxRelDiffPw);
                    Scalar maxAbsoluteGlobalResidual = GET_RUNTIME_PARAM(TypeTag, Scalar, Newton.MaxAbsoluteGlobalResidual);
//                     if( ( ((time - tInitEnd_) > eps_) && (res_ < maxAbsoluteGlobalResidual) && iterator > 1)
                    if((time > tInitEnd_ + eps_) && (*biggest < maxRelDiffPw) || (iterator == numIterations_))
                    {

                        // switch off all output for intr. measurement
//                         TranspProblem().setOutput(false);
//                         MechanicsProblem().setOutput(false);
//                         this->setOutput(false);

                        CALLGRIND_STOP_INSTRUMENTATION;
                        CALLGRIND_DUMP_STATS;
                        std::cout << "DuneTimer:" << Timer.elapsed() << std::endl;

//                         || (iterator == numIterations_) )
                        break;
//                     TranspProblem().setEvalOriginalRhs(true);
                    }

                    lastSol_ = TranspProblem().model().curSol();

                    // If iterations continue, start the Newton with the previous solution, not the new one
                    TranspProblem().model().curSol() = TranspProblem().model().prevSol();
                    MechanicsProblem().model().curSol() = MechanicsProblem().model().prevSol();
//
                }

                std::cout
                << "*************************************** \n"
                << "Timestep for t = " << time << " done \n"
                << "*************************************** \n"
                << std::endl;

                // after iterations are finished, update prevSol for the next timestep
                TranspProblem().model().prevSol() = TranspProblem().model().curSol();
                MechanicsProblem().model().prevSol() = MechanicsProblem().model().curSol();

                TranspProblem().setEffPorosityOldTimestep() = effPorosityVector_;

                return;
            }
            catch (Dune::MathError &e) {
                std::cerr << "Caught exception: " << e << std::endl;

                Scalar dt = ParentType::timeManager_.timeStepSize();
                Scalar nextDt = dt / 2;
                ParentType::timeManager_.setTimeStepSize(nextDt);
//                 uCur_ = uPrev_;
                // reset function to avoid new solutions in case of non-convergences
//                 resetTranspSol();
//                 resetChemSol();
//                 reduceChemTimestep();

                // update failed
                std::cout << "Macro time manager didn't find a solution with dt="<<dt<<" seconds. Retrying with time step of "
                          << nextDt << " seconds\n";
            }
        }

//         std::cout<<"\n the transport solution is: \n"<< TranspProblem().model().curSol() <<std::endl;
//         std::cout<<"\n the chemistry solution is: \n"<< MechanicsProblem().model().curSol() <<std::endl;

        DUNE_THROW(Dune::MathError,
                   "Macro time manager didn't find a solution after "
                   << maxFails
                   << " time-step divisions. dt="
                   << ParentType::timeManager_.timeStepSize());
    }

    /*
     *\brief Initialise effPorosity, effPermeability and volumetric strain vectors with the intrinsic values
     */
    void initializePoroPerm()
    {
        FVElementGeometryTranspProblem fvGeometryTranspProblem;
        FVElementGeometryMechanicsProblem fvGeometryMechanicsProblem;

        unsigned numVert = gridView_.size(dim);

        for (const auto& element : elements(gridView_))
        {
            if(element.partitionType() == Dune::InteriorEntity)
            {
                fvGeometryTranspProblem.update(gridView_, element);
                fvGeometryMechanicsProblem.update(gridView_, element);

                int numScv = fvGeometryTranspProblem.numScv;
                int eIdx = MechanicsProblem().model().elementMapper().index(element);

                for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
                {
                    LocalPosition cellCenter = element.geometry().center();

                    effPorosityVector_[eIdx][scvIdx] = MechanicsProblem().spatialParams().porosityAtPos(cellCenter);

                    effPermeabilityVector_[eIdx][scvIdx] = MechanicsProblem().spatialParams().intrinsicPermeabilityAtPos(cellCenter)[0][0];

                    int dofIdxGlobal = MechanicsProblem().model().vertexMapper().subIndex(element, scvIdx, dim);

//                     if (eIdx == 211)
//                     {
//                         std::cout << "B = " << MechanicsProblem().spatialParams().lameParams(element, fvGeometryMechanicsProblem, scvIdx)[1] << std::endl;
//                         std::cout << "BNode[" << dofIdxGlobal << "] before = " << BNodeWiseAveraged_[dofIdxGlobal] << std::endl;
//                     }

//                     BNodeWiseAveraged_[dofIdxGlobal] += (1.0/ MechanicsProblem().spatialParams().lameParams(element, fvGeometryMechanicsProblem, scvIdx)[1]);
                    BNodeWiseAveraged_[dofIdxGlobal] += MechanicsProblem().spatialParams().lameParams(element, fvGeometryMechanicsProblem, scvIdx)[1];
                    numScvOfNode_[dofIdxGlobal] += 1;

                    if (eIdx == 211)
                    {
//                         std::cout << "B[" << dofIdxGlobal << "] after  = " << MechanicsProblem().spatialParams().lameParams(element, fvGeometryMechanicsProblem, scvIdx)[1] << std::endl;
//                         std::cout << "BNode[" << dofIdxGlobal << "] after  = " << BNodeWiseAveraged_[dofIdxGlobal] << std::endl;
                    }
                }
            }
        }


        for (int vIdx = 0; vIdx < numVert; ++vIdx)
        {
//             BNodeWiseAveraged_[vIdx] = numScvOfNode_[vIdx]*1.0/BNodeWiseAveraged_[vIdx];
             BNodeWiseAveraged_[vIdx] = BNodeWiseAveraged_[vIdx]/numScvOfNode_[vIdx];
//             std::cout << "BNode[" << vIdx << "] final  = " << BNodeWiseAveraged_[vIdx] << std::endl;
        }

        TranspProblem().setEffPorosity() = effPorosityVector_;
        TranspProblem().setEffPermeability() = effPermeabilityVector_;

        MechanicsProblem().setEffPorosity() = effPorosityVector_;
        MechanicsProblem().setEffPermeability() = effPermeabilityVector_;
    }

    /*
     *\brief Calculate the effective pressure from pw, pn, sw and sn.
     */
    void updateParamsMomentumBalance()
    {
        unsigned numElements = gridView_.size(0);

        ScalarField rank;
        rank.resize(numElements);

        FVElementGeometryTranspProblem fvGeometryTranspProblem;
        ElementVolumeVariablesTranspProblem elemVolVarsTranspProblem;

        for (const auto& element : elements(gridView_))
        {
            if(element.partitionType() == Dune::InteriorEntity)
            {
                fvGeometryTranspProblem.update(gridView_, element);
                elemVolVarsTranspProblem.update(TranspProblem(), element, fvGeometryTranspProblem, false);

                int numScv = fvGeometryTranspProblem.numScv;

                typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> LocalFiniteElementCache;
                typedef typename LocalFiniteElementCache::FiniteElementType LocalFiniteElement;

                int eIdx = MechanicsProblem().model().elementMapper().index(element);

                for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
                {

                    rank[eIdx] = gridView_.comm().rank();

//                     Scalar pw = elemVolVarsTranspProblem[scvIdx].pressure(wPhaseIdx);
//                     Scalar pn = elemVolVarsTranspProblem[scvIdx].pressure(nPhaseIdx);
//                     Scalar pc = elemVolVarsTranspProblem[scvIdx].capillaryPressure();
//                     Scalar sw = elemVolVarsTranspProblem[scvIdx].saturation(wPhaseIdx);
//                     Scalar sn = elemVolVarsTranspProblem[scvIdx].saturation(nPhaseIdx);

                    pwVector_[eIdx][scvIdx] = elemVolVarsTranspProblem[scvIdx].pressure(wPhaseIdx);
                    pnVector_[eIdx][scvIdx] = elemVolVarsTranspProblem[scvIdx].pressure(nPhaseIdx);
                    pcVector_[eIdx][scvIdx] = elemVolVarsTranspProblem[scvIdx].capillaryPressure();
                    SwVector_[eIdx][scvIdx] = elemVolVarsTranspProblem[scvIdx].saturation(wPhaseIdx);
                    SnVector_[eIdx][scvIdx] = elemVolVarsTranspProblem[scvIdx].saturation(nPhaseIdx);
                    rhowVector_[eIdx][scvIdx] = elemVolVarsTranspProblem[scvIdx].density(wPhaseIdx);
                    rhonVector_[eIdx][scvIdx] = elemVolVarsTranspProblem[scvIdx].density(nPhaseIdx);

                    effPorosityVector_[eIdx][scvIdx] = elemVolVarsTranspProblem[scvIdx].effPorosity();

//                     if (eIdx == 210)
//                     {
//                         std::cout << "pwVector_[ " << eIdx << "][" << scvIdx << "] is " << pwVector_[eIdx][scvIdx] << std::endl;
//                     }
                }
            }

        }

        MechanicsProblem().setpw() = pwVector_;
        MechanicsProblem().setpn() = pnVector_;
        MechanicsProblem().setpc() = pcVector_;
        MechanicsProblem().setSw() = SwVector_;
        MechanicsProblem().setSn() = SnVector_;
        MechanicsProblem().setRhon() = rhonVector_;
        MechanicsProblem().setRhow() = rhowVector_;

        // pwVector_ etc. is only updated for the fixed stress scheme
        // When the pore compressibility is used, pwVector_ stores the initial values
        if (!(GET_RUNTIME_PARAM(TypeTag, bool, PoreCompressibility.UsePoreCompressibility)))
        {
//             TranspProblem().setpWOldIteration_() = pwVector_;
//             TranspProblem().setpNOldIteration_() = pnVector_;
//             TranspProblem().setsWOldIteration_() = SwVector_;
//             TranspProblem().setsNOldIteration_() = SnVector_;
        }

//         // check of get function
//         for (const auto& element : elements(gridView_))
//         {
//             if(element.partitionType() == Dune::InteriorEntity)
//             {
//                 std::cout << "Checking the getEffPressureOldIteration_ function" << std::endl;
//
//                 int eIdx = TranspProblem().model().elementMapper().index(element);
//
//                 Scalar effPressureOld = TranspProblem().getEffPressureOldIteration_(element);
//
//                 std::cout << "effPressureOld[ " << eIdx << "] is " << effPressureOld << std::endl;
//
//             }
//         }
    }

    /*
     *\brief Calculate solidity according to the solid phase in MechanicsProblem.
     */
    void updateParamsMassBalance()
    {
        // create the required scalar and vector fields
//         unsigned numVert = gridView_.size(dim);
//         unsigned numElements = gridView_.size(0);

        FVElementGeometryTranspProblem fvGeometryTranspProblem;
        FVElementGeometryMechanicsProblem fvGeometryMechanicsProblem;

//         ElementVolumeVariablesMechanicsProblem elemVolVarsMechanicsProblem;


//         ScalarField rank;
//         rank.resize(numElements);
//
//         std::vector<Scalar> zerosNumVert(numVert, 0.0);
//         std::vector<int> zerosNumVertInt(numVert, 0);
//
//         numScvOfNode_ = zerosNumVert;
//
//         std::ifstream inputfile("xOutput.txt");
//         if (inputfile.fail())
//             throw std::ios_base::failure(std::strerror(errno));
//
//         std::cout << "Read xOutput.txt" << std::endl;


//         std::string line;
//         std::vector<std::string> xOutput;
//         while (std::getline(inputfile, line))
//         {
//             xOutput.push_back(line);
//         }

        for (const auto& element : elements(gridView_))
        {
            if(element.partitionType() == Dune::InteriorEntity)
            {
//                 elemVolVarsMechanicsProblem.update(MechanicsProblem(), element, fvGeometryMechanicsProblem, false);

                fvGeometryMechanicsProblem.update(gridView_, element);
                fvGeometryTranspProblem.update(gridView_, element);

                int eIdx = MechanicsProblem().model().elementMapper().index(element);
                int numScv = fvGeometryMechanicsProblem.numScv;

                typedef typename GET_PROP_TYPE(SubTypeTag2, GridFunctionSpace) GridFunctionSpace;
                typedef Dune::PDELab::LocalFunctionSpace<GridFunctionSpace> LocalFunctionSpace;

                const GridFunctionSpace& gridFunctionSpace = MechanicsProblem().model().jacobianAssembler().gridFunctionSpace();

                const typename GridFunctionSpace::Ordering& ordering = gridFunctionSpace.ordering();
                LocalFunctionSpace localFunctionSpace(gridFunctionSpace);
                localFunctionSpace.bind(element);

                // type of function space for solid displacement vector
                // number of degrees of freedom for each displacement value (here number of element vertices)
                const unsigned int dispSize = localFunctionSpace.child(0).size();
                // type of function space of solid displacement value (one for each coordinate direction)
                typedef typename LocalFunctionSpace::template Child<0>::Type ScalarDispLFS;
                typedef typename ScalarDispLFS::Traits::FiniteElementType::
                        Traits::LocalBasisType::Traits::RangeFieldType RF;
                typedef typename ScalarDispLFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RT_V;
                typedef typename ScalarDispLFS::Traits::FiniteElementType::
                        Traits::LocalBasisType::Traits::JacobianType JacobianType_V;

                typedef typename GridView::template Codim<0>::Entity::Geometry::JacobianInverseTransposed JacobianInverseTransposed;

                for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
                {
                    LocalPosition scvCenter = fvGeometryMechanicsProblem.subContVol[scvIdx].localCenter;

                    // copy values of previous solution into prevSolutionValues Vector
                    std::vector<Scalar> prevSolutionValues(localFunctionSpace.size());
                    // copy values of current solution into curSolutionValues Vector
                    std::vector<Scalar> curSolutionValues(localFunctionSpace.size());
//                     std::vector<Scalar> fullyCoupledSolutionValues(localFunctionSpace.size());


                    for (typename LocalFunctionSpace::Traits::IndexContainer::size_type k=0; k<localFunctionSpace.size(); ++k)
                    {
                        const typename GridFunctionSpace::Ordering::Traits::DOFIndex& di = localFunctionSpace.dofIndex(k);
                        typename GridFunctionSpace::Ordering::Traits::ContainerIndex ci;
                        ordering.mapIndex(di.view(),ci);
                        prevSolutionValues[k] = MechanicsProblem().model().prevSol()[ci];
                        curSolutionValues[k] = MechanicsProblem().model().curSol()[ci];
//                         fullyCoupledSolutionValues[k] = std::stod(xOutput[xOutput.size()/2 + ci[0]]);
//                         std::cout << "fullyCoupledSolutionValues[" << k << "] = " << fullyCoupledSolutionValues[k] << ", calculated from xOutput[" << ci[0] << "]" << std::endl;
                    }

                    // evaluate gradient of displacement shape functions at the center of
                    // the sub control volume in the reference element
                    std::vector<JacobianType_V> vRefShapeGradient(dispSize);
                    localFunctionSpace.child(0).finiteElement().localBasis().evaluateJacobian(scvCenter, vRefShapeGradient);

                    // transform gradient to element in global coordinates
                    const JacobianInverseTransposed jacInvT = element.geometry().jacobianInverseTransposed(scvCenter);
                    std::vector<Dune::FieldVector<RF,dim> > vShapeGradient(dispSize);

                    // loop over element vertices
                    for (size_t i = 0; i < dispSize; i++)
                    {
                        vShapeGradient[i] = 0.0;
                        jacInvT.umv(vRefShapeGradient[i][0],vShapeGradient[i]);
                    }

                    DimMatrix gradU(0.0);
                    DimMatrix gradUFullyCoupled(0.0);
                    for(int coordDir = 0; coordDir < dim; ++coordDir) {
                        const ScalarDispLFS& uLFS = localFunctionSpace.child(coordDir);
                        // compute gradient of u
                        for (size_t i = 0; i < dispSize; i++)
                        {
                            gradU[coordDir].axpy(curSolutionValues[uLFS.localIndex(i)],vShapeGradient[i]);
//                             gradUFullyCoupled[coordDir].axpy(fullyCoupledSolutionValues[uLFS.localIndex(i)],vShapeGradient[i]);
                        }
//                         gradU[coordDir].axpy(elemVolVarsMechanicsProblem[i].displacement(coordDir), vShapeGradient[i]);
                    }

                    volumetricStrain_[eIdx][scvIdx] = 0.0;

                    for(int coordDir = 0; coordDir < dim; ++coordDir)
                    {
                        volumetricStrain_[eIdx][scvIdx] += gradU[coordDir][coordDir];
                    }

                    #if PROBLEM_IS_CC==1
                    deltaVolumetricStrainOldIteration_[eIdx][0] += volumetricStrain_[eIdx][scvIdx];
                    #else
                    deltaVolumetricStrainOldIteration_[eIdx][scvIdx] = volumetricStrain_[eIdx][scvIdx];
                    #endif

                    //Calculation of the solid displacement gradients.
                    for(int coordDir = 0; coordDir < dim; ++coordDir)
                    {
                        // get displacement function space for coordinate direction coordDir
                        const ScalarDispLFS & scalarDispLFS = localFunctionSpace.child(coordDir);
                        std::vector<RT_V> vShape(dispSize);
                        // evaluate shape functions of all element vertices for current integration point and write it into vector vShape
                        scalarDispLFS.finiteElement().localBasis().evaluateFunction(scvCenter, vShape);

//                         for (size_t i = 0; i < dispSize; i++)
//                         {
//                             dUVector_[eIdx][scvIdx][coordDir] +=( curSolutionValues[scalarDispLFS.localIndex(i)]*vShape[i]- prevSolutionValues[scalarDispLFS.localIndex(i)])*vShape[i];
//                         }
                    }
                }

                #if PROBLEM_IS_CC==1
                deltaVolumetricStrainOldIteration_[eIdx][0] = deltaVolumetricStrainOldIteration_[eIdx][0] / numScv;
                #endif
            }
        }

        if (TranspProblem().coupled() == true)
        {
//             TranspProblem().setdU() = dUVector_;
//             TranspProblem().setDeltaVolumetricStrainOldIteration() = deltaVolumetricStrainOldIteration_;
        }

//         // check of get function
//         for (const auto& element : elements(gridView_))
//         {
//             if(element.partitionType() == Dune::InteriorEntity)
//             {
//                 int eIdx = MechanicsProblem().model().elementMapper().index(element);
//                 int numScv = element.subEntities(dim);
//
//                 for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
//                 {
//                     DimVector dUVector = TranspProblem().getdU(element, fvGeometryTranspProblem, scvIdx);
//                     Scalar volumetricStrainOldIteration = TranspProblem().getDeltaVolumetricStrainOldIteration(element, fvGeometryTranspProblem, scvIdx);

//                     if (std::abs(volumetricStrainOldIteration) > 1e-6)
//                     {
//                         std::cout << "dUVector[ " << eIdx << "][" << scvIdx << "] is " << dUVector[0] << " " << dUVector[1] << std::endl;
//                         std::cout << "volumetricStrainOldIteration[ " << eIdx << "][" << scvIdx << "] is " << volumetricStrainOldIteration << std::endl;
//                     }
//                     if (eIdx == 12873)
//                     {
//                         Scalar effPorosityVectorOldIteration = MechanicsProblem().getEffPorosityOldIteration(element, fvGeometryMechanicsProblem, scvIdx);
//                         std::cout << "effPorosityVectorOldIteration[ " << eIdx << "][" << scvIdx << "] is " << effPorosityVectorOldIteration << std::endl;
//                     }
//                 }
//             }
//         }
    }

    void runTransport()
    {
        runSubProblem1();
    }
    void runMechanics()
    {
        runSubProblem2();
    }

    Scalar getSubProblem1TimeStepSize(Scalar timeStepSize)
    {
        return std::max<Scalar>(MechanicsProblem().getPreviousTimeStepSize(), timeStepSize);
    }

    void runSubProblem1()   //Transport
    {
        ParentType::subTimeManager1_.setTime(ParentType::timeManager_.time());
        ParentType::subTimeManager1_.setEndTime(ParentType::timeManager_.time() + ParentType::timeManager_.timeStepSize());
        ParentType::subTimeManager1_.setTimeStepSize(ParentType::subTimeManager1_.previousTimeStepSize());
        ParentType::subTimeManager1_.run();
    }

    void runSubProblem2()   //Mechanics
    {
        ParentType::subTimeManager2_.setTime(ParentType::timeManager_.time());
        ParentType::subTimeManager2_.setEndTime(ParentType::timeManager_.time() + ParentType::timeManager_.timeStepSize());
        ParentType::subTimeManager2_.setTimeStepSize(ParentType::subTimeManager2_.previousTimeStepSize());
        ParentType::subTimeManager2_.run();
    }

    Scalar maxTimeStepSize() const
    {
        return maxTimeStepSize_;
    }

    // allows to change the output_ variable which defines if output is written
    void setOutput(bool output)
    {
        output_ = output;
    }

    // returns true if the current solution should be written to
    // disk (i.e. as a VTK file)
    // during initialization no output is written
    // during actual simulation output is written initially and
    // at episode/simulation end
    bool shouldWriteOutput()
    {
        return output_;
    }

    // returns true if the current solution should be written to
    // disk (i.e. as a drs file)
    bool shouldWriteRestartFile() const
    {
        return output_;
    }

        /*!
     * \brief Define length of next episode
     */
    void episodeEnd()
    {
//         Scalar oldTimeStep = ParentType::timeManager_.timeStepSize();
        //calls suggestTimeStepSize function, which returns the new suggested TimeStepSize for no failure
        //and the failureTimeStepSize for anyFailure = true
//         double newTimeStepSize = TranspProblem().newtonController().suggestTimeStepSize(oldTimeStep);
//         std::cout << "el2p: newTimeStepSize is " << newTimeStepSize << "\n";
        Scalar time = ParentType::timeManager_.time();

//         std::vector<double> timeAntonio = {0.864000E+01, 0.259200E+02, 0.604800E+02, 0.129600E+03, 0.267840E+03, 0.544320E+03, 0.109728E+04, 0.220320E+04, 0.441504E+04, 0.883872E+04, 0.176861E+05, 0.353808E+05, 0.707702E+05, 0.141549E+06, 0.283107E+06, 0.424665E+06, 0.566222E+06, 0.849338E+06, 0.129600E+07, 0.218932E+07, 0.259200E+07, 0.339735E+07, 0.432000E+07, 0.616530E+07, 0.777600E+07, 0.109974E+08, 0.155520E+08, 0.246612E+08, 0.311040E+08, 0.432000E+08, 0.604800E+08, 0.777600E+08, 0.103680E+09, 0.129600E+09, 0.155520E+09};

//         std::vector<double> timeAntonio = {172800, 275200, 396562.9629629629, 540400.548696845, 710874.7243814459, 912918.1918594914, 1152377.116277916, 1436180.285959011, 1591220.906432943, 1760617.13991372, 1992752.719128119, 2246382.333454961, 2546980.394879368, 2875411.610139367, 3295073.718527144, 3792451.032171917, 4427988.710718015, 5063526.389264113, 5875602.311850795, 6606470.642178807, 7472684.9596046, 8419104.306421671, 8726400, 9672819.346817071, 10619238.69363414, 11471016.1057695, 12480530.07570771, 13306496.0511117, 14208940.35757161, 15194944.32203707, 16082347.89005599, 16445376.62242736, 16909246.66934634, 17257149.20453557, 17366400, 17714302.53518923, 18094418.26808117, 18509729.90216681, 18755840.50014348, 19093103.17144485, 19492821.89298721, 19929551.60726498, 20406719.25804995, 20836170.14375642, 21345148.97126038, 21995510.80640434, 22706091.32998755, 23482466.34649513, 24402614.51420782, 25407961.58633835, 26006400, 27011747.07213053, 28110181.836125, 29310323.52271155, 30510465.2092981, 31821731.12612415, 32538070.83994579, 33519721.55888655, 34592265.86291441, 34646400, 35718944.30402787, 36890798.2658361, 38171157.22410805, 39688619.6931711, 40930179.89513178, 42401658.65301111, 43286400, 44757878.75787932, 46082209.63997071, 47897033.44135521, 49530374.86260127, 51314951.60062935, 51926400, 53710976.73802809, 54946452.9412783, 56525116.97876469, 57945914.61250244, 59498267.58269739, 60566400, 62118752.97019495, 64102315.09877738, 66269540.38741375, 69038772.70067132, 69206400, 71975632.31325758, 75001275.02589086, 76946331.05544081, 77846400, 79791456.02954996, 82276805.40064158, 85452529.59703642, 86486400, 88074262.09819742, 90250221.26980129, 92829135.84355403, 95126400, 97705314.57375275, 101000594.3068812, 103766400, 107061679.7331285, 110662077.9600652, 112406400, 116006798.2269367, 120607307.0724669, 121046400, 125646908.8455302, 127826097.2460445, 129686400, 131865588.4005143, 134650106.9122826, 138326400, 142002693.0877174, 146700178.6998007, 146966400, 151663885.6120833, 155520000};
//
//         if(ParentType::timeManager_.episodeIndex()  == 1)
//         {
//             episodeLength_ = timeAntonio[0] - time;
//         }
//         else if(ParentType::timeManager_.episodeIndex()  <= timeAntonio.size())
//         {
//             episodeLength_ = timeAntonio[ParentType::timeManager_.episodeIndex()-1] - time;
//         }
//         else
//         {
//             episodeLength_ = 8640000;
//         }

//         if(ParentType::timeManager_.episodeIndex()  <= 2)
//         {
//             episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLengthMainSimulation);
//         }
//
// //         if( ParentType::timeManager_.episodeIndex() >= 2 && ParentType::timeManager_.episodeIndex() < 9)
// //         {
// //             episodeLength_ = episodeLength_ * 2.0;
// //         }
// //
// //         if( ParentType::timeManager_.episodeIndex() == 12)
// //             episodeLength_ = 327600;

//         if( ParentType::timeManager_.episodeIndex() >= 3 && ParentType::timeManager_.episodeIndex() < 20)
//         {
//             episodeLength_ = episodeLength_ * 2.0;
//         }
//
//         if( ParentType::timeManager_.episodeIndex() == 5)
//             episodeLength_ = 604800;
//         if( ParentType::timeManager_.episodeIndex() == 7)
//             episodeLength_ = 86400;
//         if( ParentType::timeManager_.episodeIndex() == 11)
//             episodeLength_ = 518400;
//         if( ParentType::timeManager_.episodeIndex() == 14)
//             episodeLength_ = 345600;
//         if( ParentType::timeManager_.episodeIndex() == 18)
//             episodeLength_ = 2937600;
//         if( ParentType::timeManager_.episodeIndex() == 20)
//             episodeLength_ = 9676800;
//         if( ParentType::timeManager_.episodeIndex() == 21)
//             episodeLength_ = 12096000;
//         if( ParentType::timeManager_.episodeIndex() == 22)
//             episodeLength_ = 17280000;
//         if( ParentType::timeManager_.episodeIndex() == 24)
//             episodeLength_ = 25920000;


//         if( ParentType::timeManager_.episodeIndex() >= 2 && ParentType::timeManager_.episodeIndex() < 13)
//         {
//             episodeLength_ = episodeLength_ * 2.0;
//         }
//
//         if( ParentType::timeManager_.episodeIndex() == 13)
//             episodeLength_ = 1026000;
//
//         if( ParentType::timeManager_.episodeIndex() == 14)
//         {
//             episodeLength_ = 15768000;
//         }

//         if( ParentType::timeManager_.episodeIndex() == 10)
//             episodeLength_ = 389.0;
//         if( ParentType::timeManager_.episodeIndex() == 12)
//             episodeLength_ = 122.0;
//         if( ParentType::timeManager_.episodeIndex() == 16)
//             episodeLength_ = 92.0;
//         if( ParentType::timeManager_.episodeIndex() == 21)
//             episodeLength_ = 840;
//         if( ParentType::timeManager_.episodeIndex() == 23)
//             episodeLength_ = 1820;
//         if( ParentType::timeManager_.episodeIndex() == 24)
//             episodeLength_ = 100;

//         if( ((ParentType::timeManager_.time() - 47304000) > eps_) )
// //         {
//             episodeLength_ = 15768000;
//         }

        std::cout << "el2p: episodeLength_ " <<  episodeLength_ << "\n";
        ParentType::timeManager_.startNextEpisode(episodeLength_);
        ParentType::timeManager_.setTimeStepSize(episodeLength_);
        std::cout << "el2p: TimeStepSize_ " <<  ParentType::timeManager_.timeStepSize() << "\n";

        TranspProblem().setEpisodeLength(episodeLength_);
        MechanicsProblem().setEpisodeLength(episodeLength_);

    }

    TwoPSpatialParams_Simple& TranspProblem()
    { return this->subProblem1(); }
    const TwoPSpatialParams_Simple& TranspProblem() const
    { return this->subProblem1(); }

    El2P_TestProblem_Simple& MechanicsProblem()
    { return this->subProblem2(); }
    const El2P_TestProblem_Simple& MechanicsProblem() const
    { return this->subProblem2(); }

private:
    std::string name_;

//     SolutionVectorMechanicsProblem uCur_;
//     SolutionVectorMechanicsProblem uPrev_;

    Scalar timeStepSizeFactor_;

    Scalar couplingError_ = 1e-10;
    Scalar maxCouplingError_ = 1e-4;
    Scalar minPvValue_ = 1e-10;
    Scalar dtmin_ = 20;
    Scalar maxTimeStepSize_;
    static constexpr Scalar eps_ = 1e-10;

    bool output_;

    Scalar episodeLength_;

    GridView gridView_;

//     std::vector<Dune::FieldVector<Scalar, (1 << dim)>> pwVector_, pnVector_, pcVector_, SwVector_, SnVector_, rhonVector_, rhowVector_;
//     std::vector<Dune::FieldVector<Scalar, (1 << dim)>> effPorosityVector_, effPermeabilityVector_;
//
// //     std::vector<std::vector<DimVector>> dUVector_;
//
//     std::vector<Dune::FieldVector<Scalar, (1 << dim)>> volumetricStrain_;
//     std::vector<Dune::FieldVector<Scalar, (1 << dim)>> deltaVolumetricStrainOldIteration_;

    std::vector<std::vector<Scalar>> pwVector_, pnVector_, pcVector_, SwVector_, SnVector_, rhonVector_, rhowVector_;
    std::vector<std::vector<Scalar>> effPorosityVector_, effPermeabilityVector_;

//     std::vector<std::vector<DimVector>> dUVector_;

    std::vector<std::vector<Scalar>> volumetricStrain_;
    std::vector<std::vector<Scalar>> deltaVolumetricStrainOldIteration_;

    std::vector<Scalar> epiEnd_;
    std::string injectionParameters_;

    Scalar tInitEnd_;
    int numIterations_;

    std::vector<Scalar> BNodeWiseAveraged_;
    std::vector<Scalar> numScvOfNode_;

    Scalar previousTimeStepSize_;

    SolutionVectorTranspProblem lastSol_;

    Scalar res_;

};

} //end namespace

#endif
