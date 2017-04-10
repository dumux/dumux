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

#include "decoupledelasticproblem.hh"
#include "decoupled2pproblem.hh"
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

        output_ = false;
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
        // create the required scalar and vector fields
        unsigned numVert = gridView_.size(dim);
        unsigned numElements = gridView_.size(0);

        ScalarField rank;
        rank.resize(numElements);

        std::vector<Scalar> zerosNumVert(numVert, 0.0);
        std::vector<int> zerosNumVertInt(numVert, 0);
        std::vector<Scalar> zerosNumElements(numElements, 0.0);

        pwVector_ = zerosNumVert;
        pnVector_ = zerosNumVert;
        pcVector_ = zerosNumVert;
        SwVector_ = zerosNumVert;
        SnVector_ = zerosNumVert;
        rhonVector_ = zerosNumVert;
        rhowVector_ = zerosNumVert;

        factorP_ = zerosNumVert;
        localFactorP_ = zerosNumVert;
        numScvOfNode_ = zerosNumVert;

        effPressureVector_ = zerosNumElements;

        effPorosityVector_.resize(numElements);
        effPorosityVectorOldIteration_.resize(numElements);
        effPorosityVectorFirstIteration_.resize(numElements);
        volumetricStrainOldIteration_.resize(numElements);
        effPorosityVectorOldTimestep_.resize(numElements);
        lowerLimitEffPorosity_.resize(numElements);
        upperLimitEffPorosity_.resize(numElements);
        effPermeabilityVector_.resize(numElements);
        dUVector_.resize(numElements);

        MechanicsProblem().setEffPressure() = effPressureVector_;
        MechanicsProblem().setpw() = pwVector_;
        MechanicsProblem().setpn() = pnVector_;
        MechanicsProblem().setpc() = pcVector_;
        MechanicsProblem().setSw() = SwVector_;
        MechanicsProblem().setSn() = SnVector_;
        MechanicsProblem().setRhon() = rhonVector_;
        MechanicsProblem().setRhow() = rhowVector_;

        for (const auto& element : elements(gridView_))
        {
            if(element.partitionType() == Dune::InteriorEntity)
            {
                for (int eIdx = 0; eIdx < numElements; ++eIdx)
                {
                    int numScv = element.subEntities(dim);

                    effPorosityVector_[eIdx].resize(numScv);
                    effPorosityVectorOldIteration_[eIdx].resize(numScv);
                    effPorosityVectorFirstIteration_[eIdx].resize(numScv);
                    volumetricStrainOldIteration_[eIdx].resize(numScv);
                    effPorosityVectorOldTimestep_[eIdx].resize(numScv);
                    lowerLimitEffPorosity_[eIdx].resize(numScv);
                    upperLimitEffPorosity_[eIdx].resize(numScv);
                    effPermeabilityVector_[eIdx].resize(numScv);
                    dUVector_[eIdx].resize(numScv);
                }
            }
        }

        TranspProblem().setEffPorosity() = effPorosityVector_;
        TranspProblem().setEffPorosityOldTimestep() = effPorosityVectorOldTimestep_;
        TranspProblem().setEffPorosityOldIteration() = effPorosityVectorOldIteration_;
        TranspProblem().setEffPermeability() = effPermeabilityVector_;

        MechanicsProblem().setEffPorosity() = effPorosityVector_;
        MechanicsProblem().setEffPorosityOldIteration() = effPorosityVectorOldIteration_;
        MechanicsProblem().setEffPorosityOldTimestep() = effPorosityVectorOldTimestep_;
        MechanicsProblem().setEffPermeability() = effPermeabilityVector_;

        MechanicsProblem().setFactorP() = factorP_;

        this->setOutput(true);

        ParentType::init();
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

                if(ParentType::timeManager_.time() + ParentType::timeManager_.timeStepSize() > tInitEnd_ + eps_)
                {
                    numIterations_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TimeManager.NumIterations);
                }

                for (int iterator = 1; iterator < (numIterations_ + 1); iterator++)
                {

                    if( (ParentType::timeManager_.time() + ParentType::timeManager_.timeStepSize() < tInitEnd_ + eps_) && (iterator == 1) )
                    {
                        initializePoroPerm();
                    }

//                     if (iterator == (numIterations_ - 1))
//                     {
//                         MechanicsProblem().setOutput(true);
//                         std::cout << "Set MechanicsProblem output to true." << std::endl;
//                     }

                    runTransport();

                    std::cout << "Ran transport" << std::endl;

                    updateParamsMomentumBalance();

                    std::cout << "updateParamsMomentumBalance" << std::endl;

                    runMechanics();

                    std::cout << "Ran mechanics" << std::endl;

                    updateParamsMassBalance(iterator);

                    std::cout
                    << "----------------------------- \n"
                    << "Iteration " << iterator << " done \n"
                    << "----------------------------- \n"
                    << std::endl;
                }

                Scalar time = ParentType::timeManager_.time() + ParentType::timeManager_.timeStepSize();

                std::cout
                << "*************************************** \n"
                << "Timestep for t = " << time << " done \n"
                << "*************************************** \n"
                << std::endl;

                // after iterations are finished, update prevSol for the next timestep

                TranspProblem().model().prevSol() = TranspProblem().model().curSol();
                MechanicsProblem().model().prevSol() = MechanicsProblem().model().curSol();

                effPorosityVectorOldTimestep_ = effPorosityVector_;
                TranspProblem().setEffPorosityOldTimestep() = effPorosityVectorOldTimestep_;

                // check of get function
                FVElementGeometryMechanicsProblem fvGeometryMechanicsProblem;

                for (const auto& element : elements(gridView_))
                {
                    if(element.partitionType() == Dune::InteriorEntity)
                    {
                        fvGeometryMechanicsProblem.update(gridView_, element);

                        int eIdx = MechanicsProblem().model().elementMapper().index(element);
                        int numScv = element.subEntities(dim);

                        for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
                        {
            //                             Scalar effPorosityVectorOld = TranspProblem().getEffPorosityOld(element, fvGeometryTranspProblem, scvIdx);
            //                             Scalar effPorosityVector = TranspProblem().getEffPorosity(element, fvGeometryTranspProblem, scvIdx);
            //
            //                             if (std::abs(effPorosityVector - effPorosityVectorOld) > 1e-6)
            //                             {
            //                                 std::cout << "effPorosityVector_[ " << eIdx << "][" << scvIdx << "] is " << effPorosityVector << std::endl;
            //                                 std::cout << "effPorosityVectorOld_[ " << eIdx << "][" << scvIdx << "] is " << effPorosityVectorOld << std::endl;
            //                             }
                            if (eIdx == 12873)
                            {
                                Scalar effPorosity = MechanicsProblem().getEffPorosity(element, fvGeometryMechanicsProblem, scvIdx);
                                std::cout << "effPorosity from getFunction[ " << eIdx << "][" << scvIdx << "] is " << effPorosity << std::endl;
                            }
                        }
                    }
                }
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

        for (const auto& element : elements(gridView_))
        {
            if(element.partitionType() == Dune::InteriorEntity)
            {
                int numScv = element.subEntities(dim);
                int eIdx = MechanicsProblem().model().elementMapper().index(element);

                fvGeometryTranspProblem.update(gridView_, element);

                for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
                {
                    LocalPosition cellCenter = element.geometry().center();

                    effPorosityVector_[eIdx][scvIdx] = MechanicsProblem().spatialParams().porosity(cellCenter);
                    effPorosityVectorOldTimestep_[eIdx][scvIdx] = MechanicsProblem().spatialParams().porosity(cellCenter);
                    effPorosityVectorOldIteration_[eIdx][scvIdx] = MechanicsProblem().spatialParams().porosity(cellCenter);
                    effPorosityVectorFirstIteration_[eIdx][scvIdx] = MechanicsProblem().spatialParams().porosity(cellCenter);
                    lowerLimitEffPorosity_[eIdx][scvIdx] = MechanicsProblem().spatialParams().porosity(cellCenter);
                    upperLimitEffPorosity_[eIdx][scvIdx] = MechanicsProblem().spatialParams().porosity(cellCenter);

                    effPermeabilityVector_[eIdx][scvIdx] = MechanicsProblem().spatialParams().intrinsicPermeabilityAtPos(cellCenter)[0][0];

                    volumetricStrainOldIteration_[eIdx][scvIdx] = 0.0;
                }
            }
        }

        TranspProblem().setEffPorosity() = effPorosityVector_;
        TranspProblem().setEffPorosityOldTimestep() = effPorosityVectorOldTimestep_;
        TranspProblem().setEffPermeability() = effPermeabilityVector_;

        MechanicsProblem().setEffPorosity() = effPorosityVector_;
        MechanicsProblem().setEffPorosityOldIteration() = effPorosityVectorOldIteration_;
        MechanicsProblem().setEffPermeability() = effPermeabilityVector_;

        TranspProblem().setVolumetricStrainOldIteration() = volumetricStrainOldIteration_;

        MechanicsProblem().setFactorP() = factorP_;
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
                int numScv = element.subEntities(dim);

                fvGeometryTranspProblem.update(gridView_, element);
                elemVolVarsTranspProblem.update(TranspProblem(), element, fvGeometryTranspProblem, false);


                int eIdx = MechanicsProblem().model().elementMapper().index(element);

                effPressureVector_[eIdx] = 0.0;

                for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
                {

                    rank[eIdx] = gridView_.comm().rank();


                    int dofIdxGlobal = MechanicsProblem().model().vertexMapper().subIndex(element, scvIdx, dim);

                    Scalar pw = elemVolVarsTranspProblem[scvIdx].pressure(wPhaseIdx);
                    Scalar pn = elemVolVarsTranspProblem[scvIdx].pressure(nPhaseIdx);
                    Scalar pc = elemVolVarsTranspProblem[scvIdx].capillaryPressure();
                    Scalar sw = elemVolVarsTranspProblem[scvIdx].saturation(wPhaseIdx);
                    Scalar sn = elemVolVarsTranspProblem[scvIdx].saturation(nPhaseIdx);

                    pwVector_[dofIdxGlobal] = pw;
                    pnVector_[dofIdxGlobal] = pn;
                    pcVector_[dofIdxGlobal] = pc;
                    SwVector_[dofIdxGlobal] = sw;
                    SnVector_[dofIdxGlobal] = sn;

                    effPressureVector_[eIdx] += (pn * sn
                                                    + pw * sw)
                                                    / numScv;

                    rhonVector_[dofIdxGlobal] = elemVolVarsTranspProblem[scvIdx].density(nPhaseIdx);
                    rhowVector_[dofIdxGlobal] = elemVolVarsTranspProblem[scvIdx].density(wPhaseIdx);
                }
            }

        }

        MechanicsProblem().setEffPressure() = effPressureVector_;
        MechanicsProblem().setpw() = pwVector_;
        MechanicsProblem().setpn() = pnVector_;
        MechanicsProblem().setpn() = pcVector_;
        MechanicsProblem().setSw() = SwVector_;
        MechanicsProblem().setSn() = SnVector_;
        MechanicsProblem().setRhon() = rhonVector_;
        MechanicsProblem().setRhow() = rhowVector_;
    }

    /*
     *\brief Calculate solidity according to the solid phase in MechanicsProblem.
     */
    void updateParamsMassBalance(int iteration)
    {
        // create the required scalar and vector fields
        unsigned numVert = gridView_.size(dim);
        unsigned numElements = gridView_.size(0);

        FVElementGeometryMechanicsProblem fvGeometryMechanicsProblem;
        ElementVolumeVariablesMechanicsProblem elemVolVarsMechanicsProblem;

        FVElementGeometryTranspProblem fvGeometryTranspProblem;
        ElementVolumeVariablesTranspProblem elemVolVarsTranspProblem;

        ScalarField rank;
        rank.resize(numElements);

        std::vector<Scalar> zerosNumVert(numVert, 0.0);
        std::vector<int> zerosNumVertInt(numVert, 0);

        factorP_ = zerosNumVert;
        localFactorP_ = zerosNumVert;
        numScvOfNode_ = zerosNumVert;

        for (const auto& element : elements(gridView_))
        {
            if(element.partitionType() == Dune::InteriorEntity)
            {
                int eIdx = MechanicsProblem().model().elementMapper().index(element);
                int numScv = element.subEntities(dim);

                fvGeometryMechanicsProblem.update(gridView_, element);
                elemVolVarsMechanicsProblem.update(MechanicsProblem(), element, fvGeometryMechanicsProblem, false);

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
                    for (typename LocalFunctionSpace::Traits::IndexContainer::size_type k=0; k<localFunctionSpace.size(); ++k)
                    {
                        const typename GridFunctionSpace::Ordering::Traits::DOFIndex& di = localFunctionSpace.dofIndex(k);
                        typename GridFunctionSpace::Ordering::Traits::ContainerIndex ci;
                        ordering.mapIndex(di.view(),ci);
                        prevSolutionValues[k] = MechanicsProblem().model().prevSol()[ci];
                        curSolutionValues[k] = MechanicsProblem().model().curSol()[ci];
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
                    for(int coordDir = 0; coordDir < dim; ++coordDir) {
                        const ScalarDispLFS& uLFS = localFunctionSpace.child(coordDir);
                        // compute gradient of u
                        for (size_t i = 0; i < dispSize; i++)
                            gradU[coordDir].axpy(curSolutionValues[uLFS.localIndex(i)],vShapeGradient[i]);
                    }

//                     if (gradU[0][0] > eps_)
//                     {
//                         std::cout << "gradU is " << gradU << std::endl;
//                     }

                    // calculate gradients
//                     DimVector tmp(0.0);
//                     for (int idx = 0;
//                             idx < fvGeometryMechanicsProblem.numScv;
//                             idx++) // loop over adjacent vertices
//                     {
//                         // FE gradient at vertex idx
//                         const DimVector &centerGrad = fvGeometryMechanicsProblem.subContVol[scvIdx].grad[idx];
//
//                         if (centerGrad[0] > eps_)
//                         {
//                             std::cout << "centerGrad is " << centerGrad << std::endl;
//                         }
//
//
//                         // the displacement vector gradient
//                         for (int coordIdx = 0; coordIdx < dim; ++coordIdx) {
//                             tmp = centerGrad;
//                             tmp *= elemVolVarsMechanicsProblem[idx].displacement(coordIdx);
//                             gradU[coordIdx] += tmp;
//                         }
//                     }

                    volumetricStrainOldIteration_[eIdx][scvIdx] = 0.0;

                    for(int coordDir = 0; coordDir < dim; ++coordDir)
                    {
                        volumetricStrainOldIteration_[eIdx][scvIdx] += gradU[coordDir][coordDir];
//                         volumetricStrain += gradU[coordDir][coordDir];
                    }


                    Scalar intrinsicPorosity = MechanicsProblem().spatialParams().porosity(scvCenter);

                    if (eIdx == 12873)
                    {
//                         std::cout << "volumetricStrain of " << eIdx << "][" << scvIdx << "] is " << volumetricStrain << std::endl;

                        std::cout << "effPorosityVector_ " << eIdx << "][" << scvIdx << "] before is " << effPorosityVector_[eIdx][scvIdx] << std::endl;
                        Scalar effPorosityOldTimestep = TranspProblem().getEffPorosityOldTimestep(element, fvGeometryTranspProblem, scvIdx);
                        std::cout << "effPorosityVectorOldTimestep_ " << eIdx << "][" << scvIdx << "] before is " << effPorosityOldTimestep << std::endl;
                    }

                    effPorosityVector_[eIdx][scvIdx] = 1 - (1 - intrinsicPorosity)*exp(-volumetricStrainOldIteration_[eIdx][scvIdx]);

                    Scalar factor = pow( (effPorosityVector_[eIdx][scvIdx]/ intrinsicPorosity), 15);

                    //intrinsicPermeability is a DimMatrix, but we assume an isotropic K here and just used one entry to calculate Keff
                    effPermeabilityVector_[eIdx][scvIdx] = MechanicsProblem().spatialParams().intrinsicPermeability(element, fvGeometryMechanicsProblem, scvIdx)[0][0] * factor;

                    //Calculation of the solid displacement gradients.
                    for(int coordDir = 0; coordDir < dim; ++coordDir)
                    {
                        // get displacement function space for coordinate direction coordDir
                        const ScalarDispLFS & scalarDispLFS = localFunctionSpace.child(coordDir);
                        std::vector<RT_V> vShape(dispSize);
                        // evaluate shape functions of all element vertices for current integration point and write it into vector vShape
                        scalarDispLFS.finiteElement().localBasis().evaluateFunction(scvCenter, vShape);

                        dUVector_[eIdx][scvIdx] = 0.0;

                        for (size_t i = 0; i < dispSize; i++){
                            dUVector_[eIdx][scvIdx][coordDir] +=( curSolutionValues[scalarDispLFS.localIndex(i)]*vShape[i]- prevSolutionValues[scalarDispLFS.localIndex(i)])*vShape[i];                       }

                        dUVector_[eIdx][scvIdx][coordDir] = curSolutionValues[coordDir] - prevSolutionValues[coordDir];
                    }

                    if(ParentType::timeManager_.time() + ParentType::timeManager_.timeStepSize() > tInitEnd_ + eps_)
                    {

                        if (eIdx == 12873)
                        {
                        std::cout << "LowerLimitEffPorosity[" << eIdx << "][" << scvIdx << "] is "
                        << lowerLimitEffPorosity_[eIdx][scvIdx] << "\n"
                        << "UpperLimitEffPorosity[" << eIdx << "][" << scvIdx << "] is "
                        << upperLimitEffPorosity_[eIdx][scvIdx] << "\n"
                        << "effPorosityVector_[" << eIdx << "][" << scvIdx << "] is "
                        << effPorosityVector_[eIdx][scvIdx] << std::endl;
                        }

                        Scalar effPorosityChange = effPorosityVector_[eIdx][scvIdx] - effPorosityVectorOldIteration_[eIdx][scvIdx];

                        Scalar dampingFactor = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Damping, DampingFactor);

                        // damp porosity update
                        effPorosityVector_[eIdx][scvIdx] = effPorosityVectorOldIteration_[eIdx][scvIdx] + dampingFactor * effPorosityChange;

                        if (eIdx == 12873)
                        {
                        std::cout << "damped change [" << eIdx << "][" << scvIdx << "] is " << dampingFactor * effPorosityChange << "\n"
                        << "effPorosityVector_[" << eIdx << "][" << scvIdx << "] is "
                        << effPorosityVector_[eIdx][scvIdx] << std::endl;
                        }

//                                     std::cout << "effPorosityVector_ damped [" << eIdx << "][" << scvIdx << "] is " << effPorosityVector_[eIdx][scvIdx] << std::endl;


                        // if damped porosity is still above/below limits
                        if ( (effPorosityVector_[eIdx][scvIdx] > upperLimitEffPorosity_[eIdx][scvIdx] + eps_) || (effPorosityVector_[eIdx][scvIdx] < lowerLimitEffPorosity_[eIdx][scvIdx] - eps_) )
                        {
                            if (eIdx == 12873)
                            {
//                                     DUNE_THROW(Dune::MathError,
//                                     "Porosity update is diverging! \n"
//                                     << "LowerLimitEffPorosity[" << eIdx << "][" << scvIdx << "] is "
//                                     << lowerLimitEffPorosity_[eIdx][scvIdx] << "\n"
//                                     << "UpperLimitEffPorosity[" << eIdx << "][" << scvIdx << "] is "
//                                     << upperLimitEffPorosity_[eIdx][scvIdx] << "\n"
//                                     << "effPorosityVector_[" << eIdx << "][" << scvIdx << "] is "
//                                     << effPorosityVector_[eIdx][scvIdx]);
                            std::cout << "Not confined within limits!" << std::endl;
                            std::cout << "effPorosityVector_[" << eIdx << "][" << scvIdx << "] is " << effPorosityVector_[eIdx][scvIdx] << std::endl;
                            }
                        }
//                             }
//                         }
                    }

                    effPorosityVectorOldIteration_[eIdx][scvIdx] = effPorosityVector_[eIdx][scvIdx];

//                     if (eIdx == 12873)
//                     {
//                         std::cout << "effPorosityVectorOldIteration_[" << eIdx << "][" << scvIdx << "]:" << std::endl;
//                         std::cout << effPorosityVectorOldIteration_[eIdx][scvIdx] << std::endl;
//                     }

                }

//                 if(ParentType::timeManager_.time() + ParentType::timeManager_.timeStepSize() > tInitEnd_ - eps_)
//                 {
// //                     std::cout << "Updating pressure solution" << std::endl;
//                     // declare and initialize start and end of vertex iterator
//                     for (int scvIdx = 0; scvIdx < numScv; ++scvIdx)
//                     {
//                         int dofIdxGlobal = MechanicsProblem().model().vertexMapper().subIndex(element, scvIdx, dim);
//
//                         TranspProblem().model().curSol()[dofIdxGlobal][pressureIdx] *=  factorP_[dofIdxGlobal];
//                     }
//                 }
            }
        }


        TranspProblem().setEffPorosity() = effPorosityVector_;
        TranspProblem().setEffPorosityOldTimestep() = effPorosityVectorOldTimestep_;
        TranspProblem().setEffPermeability() = effPermeabilityVector_;
        TranspProblem().setdU() = dUVector_;
        TranspProblem().setEffPorosityOldIteration() = volumetricStrainOldIteration_;

        MechanicsProblem().setEffPorosity() = effPorosityVector_;
        MechanicsProblem().setEffPorosityOldIteration() = effPorosityVectorOldIteration_;
        MechanicsProblem().setEffPermeability() = effPermeabilityVector_;

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
// //                     Scalar effPorosityVectorOld = TranspProblem().getEffPorosityOld(element, fvGeometryTranspProblem, scvIdx);
// //                     Scalar effPorosityVector = TranspProblem().getEffPorosity(element, fvGeometryTranspProblem, scvIdx);
// //
// //                     if (std::abs(effPorosityVector - effPorosityVectorOld) > 1e-6)
// //                     {
// //                         std::cout << "effPorosityVector_[ " << eIdx << "][" << scvIdx << "] is " << effPorosityVector << std::endl;
// //                         std::cout << "effPorosityVectorOld_[ " << eIdx << "][" << scvIdx << "] is " << effPorosityVectorOld << std::endl;
// //                     }
//                     if (eIdx == 12873)
//                     {
//                         Scalar effPorosityVectorOldIteration = MechanicsProblem().getEffPorosityOldIteration(element, fvGeometryMechanicsProblem, scvIdx);
//                         std::cout << "effPorosityVectorOldIteration[ " << eIdx << "][" << scvIdx << "] is " << effPorosityVectorOldIteration << std::endl;
//                     }
//                 }
//             }
//         }

        MechanicsProblem().setFactorP() = factorP_;
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
        Scalar oldTimeStep = ParentType::timeManager_.timeStepSize();
        //calls suggestTimeStepSize function, which returns the new suggested TimeStepSize for no failure
        //and the failureTimeStepSize for anyFailure = true
        double newTimeStepSize = TranspProblem().newtonController().suggestTimeStepSize(oldTimeStep);
        std::cout << "newTimeStepSize is " << newTimeStepSize << "\n";
        if(ParentType::timeManager_.time() + ParentType::timeManager_.timeStepSize() > tInitEnd_ + eps_)
        {
            MechanicsProblem().setCoupled(true);
            std::cout << "Coupled set to true \n";
            TranspProblem().setInitializationRun(false);

            episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLengthMainSimulation);
            std::cout << "episodeLength_ is " << episodeLength_ << "\n";
        }
        else
        {
            episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLengthMainSimulation);
            std::cout << "episodeLength_ is " << episodeLength_ << "\n";
        }
//         if(ParentType::timeManager_.time() > 12000 - eps_)
//         {
//             episodeLength_ = 100.0;
//         }

        ParentType::timeManager_.startNextEpisode(episodeLength_);
        ParentType::timeManager_.setTimeStepSize(newTimeStepSize);
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

    std::vector<Scalar> effPressureVector_;
    std::vector<Scalar> pwVector_, pnVector_, pcVector_, SwVector_, SnVector_, rhonVector_, rhowVector_;
    std::vector<std::vector<Scalar>> effPorosityVector_, effPermeabilityVector_;
    std::vector<std::vector<Scalar>> effPorosityVectorOldIteration_;
    std::vector<std::vector<Scalar>> effPorosityVectorFirstIteration_;
    std::vector<std::vector<Scalar>> volumetricStrainOldIteration_;
    std::vector<std::vector<Scalar>> effPorosityVectorOldTimestep_;
    std::vector<std::vector<Scalar>> lowerLimitEffPorosity_, upperLimitEffPorosity_;

    std::vector<std::vector<DimVector>> dUVector_;

    std::vector<Scalar> epiEnd_;
    std::string injectionParameters_;

    Scalar tInitEnd_;
    int numIterations_;

    std::vector<Scalar> factorP_;
    std::vector<Scalar> localFactorP_;
    std::vector<Scalar> numScvOfNode_;

    Scalar previousTimeStepSize_;

};

} //end namespace

#endif
