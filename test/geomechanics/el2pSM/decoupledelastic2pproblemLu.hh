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
#ifndef DUMUX_EL2P_2PDECOUPLED_PROBLEM_HH
#define DUMUX_EL2P_2PDECOUPLED_PROBLEM_HH

#include <dumux/geomechanics/el2p/decoupled/simplecoupledproblem.hh> // copied from dumux-devel

#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>

#include "decoupledelasticproblemLu.hh"
#include "decoupled2pproblemLu.hh"
#include <vector>
#include "time.h"

#include <dumux/implicit/adaptive/gridadaptproperties.hh>

namespace Dumux
{
template <class TypeTag>
class El2p_2pDecoupledProblem;

namespace Properties
{
NEW_TYPE_TAG(El2p_2pDecoupledProblem, INHERITS_FROM(SimpleCoupled, GridAdapt));

// Set the problem property
SET_TYPE_PROP(El2p_2pDecoupledProblem, Problem,
              Dumux::El2p_2pDecoupledProblem<TypeTag>);

SET_INT_PROP(El2p_2pDecoupledProblem, ImplicitMaxTimeStepDivisions, 5);

// Set the two sub-problems of the global problem
//SET_TYPE_PROP(El2p_2pDecoupledProblem, SubProblem1TypeTag, TTAG(BioTranspProblem));
//SET_TYPE_PROP(El2p_2pDecoupledProblem, SubProblem2TypeTag, TTAG(BioChemProblem));
SET_TYPE_PROP(El2p_2pDecoupledProblem, SubProblem1TypeTag, TTAG(TwoP_TestProblem));
SET_TYPE_PROP(El2p_2pDecoupledProblem, SubProblem2TypeTag, TTAG(El2P_TestProblem));

SET_PROP(TwoP_TestProblem, ParameterTree)
{private:
    typedef typename GET_PROP(TTAG(El2p_2pDecoupledProblem), ParameterTree) ParameterTree;
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

SET_PROP(El2P_TestProblem, ParameterTree)
{private:
    typedef typename GET_PROP(TTAG(El2p_2pDecoupledProblem), ParameterTree) ParameterTree;
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
class El2p_2pDecoupledProblem : public SimpleCoupledProblem<TypeTag>
{
    typedef SimpleCoupledProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // obtain the type tags of the subproblems
    typedef typename GET_PROP_TYPE(TypeTag, SubProblem1TypeTag) SubTypeTag1;
    typedef typename GET_PROP_TYPE(TypeTag, SubProblem2TypeTag) SubTypeTag2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, Problem) TwoP_TestProblem;
    typedef typename GET_PROP_TYPE(SubTypeTag2, Problem) El2P_TestProblem;

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
    El2p_2pDecoupledProblem(TimeManager &timeManager, const GridView &gridView)
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

        this->setOutput(true);

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
//         const int maxFails =
//                 GET_PARAM_FROM_GROUP(TypeTag, int, Implicit, MaxTimeStepDivisions);
        const int maxFails = 5;
        for (int i = 0; i < maxFails; ++i) {
            try {

                MechanicsProblem().spatialParams().setEpisode(ParentType::timeManager_.episodeIndex());

//                 if (ParentType::timeManager_.time() > eps_)
//                 {
//                     Scalar oldTimeStep = ParentType::timeManager_.timeStepSize();
//                     Scalar newTimeStepSize = TranspProblem().newtonController().suggestTimeStepSize(oldTimeStep);
//                     std::cout << "newTimeStepSize " << newTimeStepSize << std::endl;
//                     ParentType::timeManager_.setTimeStepSize(newTimeStepSize);
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
                    if( (std::abs(time - tInitEnd_)  < eps_) && (iterator == numIterations_) )
                    {
                        TranspProblem().setCoupled(true);
                        MechanicsProblem().setCoupled(true);
                        TranspProblem().setInitializationRun(false);

                        TranspProblem().setOutput(true);
                        MechanicsProblem().setOutput(true);

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
//                         || (iterator == numIterations_) )
                        break;
//                     TranspProblem().setEvalOriginalRhs(true);

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
                int numScv = fvGeometryTranspProblem.numScv;

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

                    deltaVolumetricStrainOldIteration_[eIdx][scvIdx] = volumetricStrain_[eIdx][scvIdx];

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
//         //calls suggestTimeStepSize function, which returns the new suggested TimeStepSize for no failure
//         //and the failureTimeStepSize for anyFailure = true
//         double newTimeStepSize = TranspProblem().newtonController().suggestTimeStepSize(oldTimeStep);
//         std::cout << "el2p: newTimeStepSize is " << newTimeStepSize << "\n";
//         Scalar time = ParentType::timeManager_.time() + ParentType::timeManager_.timeStepSize();
        Scalar time = ParentType::timeManager_.time();
//         std::cout << "el2p: time= " <<  time << "\n";

//         if(ParentType::timeManager_.episodeIndex()  <= 2)
//         {
//             episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLengthMainSimulation);
//         }

        std::vector<double> timeOutput = {0, 0.01, 0.0125, 0.016157, 0.022185, 0.03212, 0.048495, 0.075482, 0.11996, 0.19327, 0.31409, 0.51322, 0.84142, 1.3823, 1.8281, 2.5215, 3.6001, 5.278, 7.888, 12.19, 18.881, 29.29, 44.518, 68.207, 105.05, 162.37, 251.54, 381.98, 584.89, 881.74, 1288.5, 1586.1, 2049, 2387.6, 2635.2, 2827.9, 3127.5, 3593.6, 4318.7, 5446.6, 7096.7, 9510.7, 13042, 17555, 23739, 32213, 44611, 62748, 74335, 86400, 98465, 107290, 120200, 129640, 142590, 152050, 165900, 172800, 186650, 196780, 211600, 221760, 236610, 247480, 259200, 270920, 288070, 300610, 318960, 332280, 345600, 358920, 378410, 392660, 413510, 432000, 452850, 479500, 497760, 518400, 539040, 569250, 604800, 640350, 664710, 691200, 717690, 756440, 777600, 816350, 840170, 864000, 887830, 920480, 950400, 983050, 1030800, 1036800, 1084600, 1123200, 1171000, 1209600, 1257400, 1296000, 1343800, 1382400, 1430200, 1468800, 1516600, 1555200, 1603000, 1641600, 1689400, 1728000, 1775800, 1814400, 1862200, 1900800, 1948600, 1987200, 2035000, 2073600, 2121400, 2160000, 2207800, 2246400, 2294200, 2332800, 2380600, 2419200, 2467000, 2505600, 2553400, 2592000, 2639800, 2678400, 2726200, 2764800, 2812600, 2831900, 2841500, 2851200, 2860900, 2875900, 2897900, 2930000, 2937600, 2969800, 3016800, 3024000, 3026900, 3031500, 3038600, 3043500, 3050600, 3061700, 3070400, 3083800, 3104700, 3110400, 3131300, 3161800, 3196800, 3231800, 3282900, 3283200, 3308800, 3346200, 3369600, 3407000, 3456000, 3505000, 3523700, 3533000, 3537700, 3542400, 3547100, 3550700, 3553700, 3558700, 3566300, 3578300, 3596900, 3624100, 3628800, 3656000, 3675900, 3705000, 3715200, 3744300, 3786900, 3801600, 3844200, 3888000, 3931800, 3974400, 4018200, 4060800, 4104600, 4147200, 4191000, 4233600, 4277400, 4320000, 4363800, 4385100, 4395700, 4406400, 4417100, 4425300, 4438200, 4458300, 4487600, 4492800, 4522100, 4565000, 4579200, 4622100, 4665600, 4687300, 4695300, 4707700, 4717300, 4724800, 4736400, 4744200, 4752000, 4759800, 4765800, 4775300, 4789900, 4812700, 4838400, 4864100, 4901600, 4924800, 4962300, 5011200, 5060100, 5097600, 5146500, 5184000, 5232900, 5270400, 5319300, 5328600, 5332300, 5335300, 5337600, 5341300, 5344300, 5346700, 5350800, 5353800, 5356800, 5359800, 5362100, 5365900, 5368900, 5371400, 5375400, 5378600, 5381000, 5384800, 5388000, 5392900, 5400500, 5412300, 5430600, 5443200, 5461600, 5490200, 5529600, 5569000, 5616000, 5663000, 5702400, 5749400, 5788800, 5835800, 5875200, 5922200, 5961600, 6008600, 6048000, 6095000, 6134400, 6181400, 6220800, 6267800, 6307200, 6354200, 6393600, 6440600, 6480000, 6527000, 6566400, 6613400, 6652800, 6699800, 6739200, 6786200, 6825600, 6872600, 6892300, 6912000, 6931700, 6960600, 6998400, 7036200, 7084800, 7109100, 7128000, 7141800, 7163300, 7171200, 7192700, 7226100, 7257600, 7291000, 7339900, 7344000, 7392900, 7430400, 7479300, 7498000, 7507400, 7516800, 7521500, 7528800, 7540100, 7548900, 7562700, 7584000, 7603200, 7624500, 7655800, 7689600, 7723400, 7772900, 7776000, 7825500, 7862400, 7911900, 7948800, 7998300, 8035200, 8084700, 8121600, 8171100, 8208000, 8257500, 8294400, 8343900, 8380800, 8430300, 8467200, 8516700, 8553600, 8603100, 8640000, 8689500, 8726400, 8775900, 8812800, 8837500, 8876000, 8899200, 8937700, 8985600, 9033500, 9072000, 9119900, 9158400, 9206300, 9215900, 9223400, 9235800, 9244800, 9257100, 9276300, 9306200, 9331200, 9361000, 9404700, 9417600, 9461300, 9504000, 9547700, 9590400, 9634100, 9676800, 9720500, 9763200, 9806900, 9849600, 9893300, 9936000, 9979700, 10022000, 10066000, 10109000, 10152000, 10195000, 10239000, 10282000, 10325000, 10368000, 10412000, 10454000, 10498000, 10541000, 10584000, 10627000, 10671000, 10714000, 10757000, 10800000, 10844000, 10886000, 10930000, 10973000, 11016000, 11059000, 11103000, 11146000, 11189000, 11232000, 11276000, 11318000, 11362000, 11405000, 11448000, 11491000, 11535000, 11578000, 11621000, 11664000, 11708000, 11750000, 11794000, 11837000, 11880000, 11923000, 11967000, 12010000, 12053000, 12096000, 12140000, 12182000, 12226000, 12269000, 12312000, 12355000, 12399000, 12442000, 12485000, 12528000, 12572000, 12614000, 12658000, 12701000, 12744000, 12787000, 12831000, 12874000, 12917000, 12960000};//output time is taken from the output of fully coupled model to get rid of the error due to time step difference


         if(ParentType::timeManager_.episodeIndex()  <= timeOutput.size())
         {
             episodeLength_ = timeOutput[ParentType::timeManager_.episodeIndex() + 1 ] - time;
         }
         else
         {
             episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLengthMainSimulation);
         }

//         std::cout << "el2p: episodeLength_ " <<  episodeLength_ << "\n";
//         std::cout << "el2p: episodeIndex " <<  ParentType::timeManager_.episodeIndex() << "\n";

        ParentType::timeManager_.startNextEpisode(episodeLength_);
        ParentType::timeManager_.setTimeStepSize(episodeLength_);
        std::cout << "el2p: TimeStepSize_ " <<  ParentType::timeManager_.timeStepSize() << "\n";

        TranspProblem().setEpisodeLength(episodeLength_);
        MechanicsProblem().setEpisodeLength(episodeLength_);

        time = ParentType::timeManager_.time() + ParentType::timeManager_.timeStepSize();
        std::cout << "el2p: time= " <<  time << "\n";

    }

    TwoP_TestProblem& TranspProblem()
    { return this->subProblem1(); }
    const TwoP_TestProblem& TranspProblem() const
    { return this->subProblem1(); }

    El2P_TestProblem& MechanicsProblem()
    { return this->subProblem2(); }
    const El2P_TestProblem& MechanicsProblem() const
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
