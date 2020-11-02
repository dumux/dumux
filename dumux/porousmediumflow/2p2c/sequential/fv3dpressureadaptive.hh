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
 * \brief Finite volume diffusion model.
 */
#ifndef DUMUX_FV3DPRESSURE2P2C_ADAPTIVE_HH
#define DUMUX_FV3DPRESSURE2P2C_ADAPTIVE_HH

// dune environent:
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

// dumux environment
#include <dumux/porousmediumflow/2p2c/sequential/fvpressuremultiphysics.hh>
#include <dumux/common/math.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/porousmediumflow/2p2c/sequential/adaptiveproperties.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

// include pressure model from Markus
#include <dumux/porousmediumflow/sequential/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/2p2c/sequential/fvmpfal3dinteractionvolumecontaineradaptive.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/lmethod/3dtransmissibilitycalculator.hh>

namespace Dumux {
namespace Properties {
template<class TypeTag>
struct MPFAInteractionVolume<TypeTag, TTag::SequentialTwoPTwoCAdaptive> { using type = FvMpfaL3dInteractionVolumeAdaptive<TypeTag>; };
template<class TypeTag>
struct MPFAInteractionVolumeContainer<TypeTag, TTag::SequentialTwoPTwoCAdaptive> { using type = FvMpfaL3d2P2CInteractionVolumeContainerAdaptive<TypeTag>; };
} // end namespace Properties

/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief The finite volume model for the solution of the compositional pressure equation
 *
 * Provides a Finite Volume implementation for the pressure equation of a compressible
 * system with two components. An IMPES-like method is used for the sequential
 * solution of the problem.  Diffusion is neglected, capillarity can be regarded.
 * Isothermal conditions and local thermodynamic
 * equilibrium are assumed.  Gravity is included.
 * \f[
        c_{total}\frac{\partial p}{\partial t} + \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}} \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{\alpha} \bf{v}_{\alpha}\right)
         = \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}} q^{\kappa},
 * \f]
 * where \f$\bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 * \f$ c_{total} \f$ represents the total compressibility, for constant porosity this yields \f$ - \frac{\partial V_{total}}{\partial p_{\alpha}} \f$,
 * \f$p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability, \f$ \lambda_{\alpha} \f$ the phase mobility,
 * \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and \f$ C^{\kappa} \f$ the total Component concentration.
 * See paper SPE 99619 or "Analysis of a Compositional Model for Fluid
 * Flow in Porous Media" by Chen, Qin and Ewing for derivation.
 *
 * The pressure base class FVPressure assembles the matrix and right-hand-side vector and solves for the pressure vector,
 * whereas this class provides the actual entries for the matrix and RHS vector.
 * The partial derivatives of the actual fluid volume \f$ v_{total} \f$ are gained by using a secant method.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class FV3dPressure2P2CAdaptive
: public FVPressure2P2CMultiPhysics<TypeTag>
{
    //the model implementation
    using Implementation = GetPropType<TypeTag, Properties::PressureModel>;
    using ParentType = FVPressure2P2CMultiPhysics<TypeTag>;
    using BaseType = FVPressure<TypeTag>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using CellData = GetPropType<TypeTag, Properties::CellData>;
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld,
        NumPhases = getPropValue<TypeTag, Properties::NumPhases>(), NumComponents = getPropValue<TypeTag, Properties::NumComponents>()
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureN,
        pGlobal = Indices::pressureGlobal,
        Sw = Indices::saturationW,
        Sn = Indices::saturationN
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx,
        contiWEqIdx = Indices::contiWEqIdx, contiNEqIdx = Indices::contiNEqIdx
    };
    enum
    {
        rhs = BaseType::rhs, matrix = BaseType::matrix,
    };

    // using declarations to abbreviate several dune classes...
    using Vertex = typename GridView::Traits::template Codim<dim>::Entity;
    using Element = typename GridView::Traits::template Codim<0>::Entity;

    using Grid = typename GridView::Grid;
    using Intersection = typename GridView::Intersection;
    using IntersectionIterator = typename GridView::IntersectionIterator;

    // convenience shortcuts for Vectors/Matrices
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using TransmissivityMatrix = Dune::FieldVector<Scalar,dim+1>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;
    using PhaseVector = Dune::FieldVector<Scalar, getPropValue<TypeTag, Properties::NumPhases>()>;
    using ComponentVector = Dune::FieldVector<Scalar, getPropValue<TypeTag, Properties::NumComponents>()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    // the typenames used for the stiffness matrix and solution vector
    using Matrix = GetPropType<TypeTag, Properties::PressureCoefficientMatrix>;
    using RHSVector = GetPropType<TypeTag, Properties::PressureRHSVector>;

    // Dumux MPFA types
    using InteractionVolumeContainer = GetPropType<TypeTag, Properties::MPFAInteractionVolumeContainer>;
    using InteractionVolume = typename InteractionVolumeContainer::InteractionVolume;

protected:
    //! \cond \private
    Problem& problem()
    {
        return this->problem_;
    }
    const Problem& problem() const
    {
        return this->problem_;
    }
    //! \endcond

public:
    void update()
    {
        //! update of interaction volumes
        //! \todo maybe only do this if the grid changed
        if(enableMPFA && maxInteractionVolumes>1)
        {
            if(!interactionVolumesContainer_)
                interactionVolumesContainer_ =
                        new InteractionVolumeContainer(problem());

            interactionVolumesContainer_->update();
        }
        asImp_().initializeMatrix();
        ParentType::update();
    }
    //initializes the matrix to store the system of equations
    void initializeMatrix();
    void initializeMatrixRowSize();
    void initializeMatrixIndices();

    void initialize(bool solveTwice = false){
        if(enableMPFA && maxInteractionVolumes>1)
        {
            if(!interactionVolumesContainer_)
                interactionVolumesContainer_ =
                        new InteractionVolumeContainer(problem());

            interactionVolumesContainer_->update();
        }
        ParentType::initialize(solveTwice);
    }

    //function which assembles the system of equations to be solved
    void assemble(bool first);

    void getMpfaFlux(const IntersectionIterator&, const CellData&);

    void get1pMpfaFlux(const IntersectionIterator&, const CellData&);

    //constitutive functions are initialized and stored in the variables object
    void updateMaterialLaws(bool fromPostTimestep = false);

    // mpfa transmissibilities
    int computeTransmissibilities(const IntersectionIterator&,
                                    TransmissivityMatrix&,
                                    GlobalPosition&,
                                    int&,
                                    GlobalPosition&,
                                    int&);

    //! Adapt primary variables vector after adapting the grid
    void adaptPressure()
    {
        int gridSize = problem().gridView().size(0);
        this->pressure().resize(gridSize);

        for (const auto& element : elements(problem().gridView()))
        {
            // get the global index of the cell
            int eIdxGlobalI = problem().variables().index(element);

            // assemble interior element contributions
            if (element.partitionType() == Dune::InteriorEntity)
            {
                this->pressure()[eIdxGlobalI]
                      = problem().variables().cellData(eIdxGlobalI).pressure(this->pressureType);
            }
        }
#if HAVE_MPI
    // communicate updated values
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using ElementMapper = typename SolutionTypes::ElementMapper;
    using PressureSolution = GetPropType<TypeTag, Properties::PressureSolutionVector>;
    using DataHandle = VectorCommDataHandleEqual<ElementMapper, PressureSolution, 0/*elementCodim*/>;

        DataHandle dataHandle(problem().variables().elementMapper(), this->pressure());
        problem().gridView().template communicate<DataHandle>(dataHandle,
                                                            Dune::InteriorBorder_All_Interface,
                                                            Dune::ForwardCommunication);
#endif
    }

    /*!
     * \brief Constructs a FVPressure2P2C object
     * \param problem a problem class object
     */
    FV3dPressure2P2CAdaptive(Problem& problem) : FVPressure2P2CMultiPhysics<TypeTag>(problem),
            enableMPFA(false),
            interactionVolumesContainer_(0), mpfal3DTransmissibilityCalculator_(problem)
    {
        enableVolumeIntegral_ = this->enableVolumeIntegral;
        enableMPFA = getParam<bool>("GridAdapt.EnableMultiPointFluxApproximation");

        maxInteractionVolumes = getParam<int>("GridAdapt.MaxInteractionVolumes");
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    {   return *static_cast<Implementation *>(this);}

    //! \copydoc IMPETProblem::asImp_()
    const Implementation &asImp_() const
    {   return *static_cast<const Implementation *>(this);}

    int searchCommonVertex_(const Intersection& is, Vertex& vertex)
    {
        /******* get corner of interest ************/
        // search through corners of large cell with isIt
        int localIdxLarge = 0;
        for(localIdxLarge = 0; localIdxLarge < is.inside().subEntities(dim); ++localIdxLarge)
        {
            auto vLarge = is.inside().template subEntity<dim>(localIdxLarge);

            // search through corners of small cell with isIt
            for(int verticeSmall = 0; verticeSmall<is.outside().subEntities(dim); ++verticeSmall)
            {
                auto vSmall = is.outside().template subEntity<dim>(verticeSmall);

                if(problem().variables().index(vSmall) == problem().variables().index(vLarge) )
                {
                    vertex = vSmall;
                    return localIdxLarge;
                }
            }
        }
        return -1;
    }

protected:
    int transmissibilityAdapter_(const IntersectionIterator& isIt,
                                    InteractionVolume& interactionVolume,
                                    const int& subVolumeFaceIdx,
                                    bool properFluxDirection,
                                    Element& additional2,
                                    Element& additional3,
                                    TransmissivityMatrix& additionalT);

    std::map<int, std::vector<int> > irregularCellMap_; //!< Container to store all cell's Indice with a hanging node
    bool enableVolumeIntegral_; //!< Calculates the volume integral (on by default)
    bool enableMPFA; //!< Enables mpfa on hanging nodes (on by default)
    int maxInteractionVolumes; //!< Maximum number of interaction volumes considered (4 by default)

    //! A pointer to the adaptive interaction volumes container
    InteractionVolumeContainer* interactionVolumesContainer_;
    //! The common implementation to calculate the Transmissibility with the mpfa-L-method
    FvMpfaL3dTransmissibilityCalculator<TypeTag> mpfal3DTransmissibilityCalculator_;
};

//! \copydoc FV2dPressure2P2CAdaptive::initializeMatrix()
template<class TypeTag>
void FV3dPressure2P2CAdaptive<TypeTag>::initializeMatrix()
{
    int gridSize_ = problem().gridView().size(0);
    // update RHS vector, matrix
    this->A_.setSize (gridSize_,gridSize_); //
    this->f_.resize(gridSize_);
    irregularCellMap_.clear();

    this->initializeMatrixRowSize();
    this->A_.endrowsizes();
    this->initializeMatrixIndices();
    this->A_.endindices();
}

//!Initialize the row sizes of the sparse global matrix
template<class TypeTag>
void FV3dPressure2P2CAdaptive<TypeTag>::initializeMatrixRowSize()
{
    // determine matrix row sizes
    for (const auto& element : elements(problem().gridView()))
    {
        // cell index
        int eIdxGlobalI = problem().variables().index(element);
        CellData& cellDataI = problem().variables().cellData(eIdxGlobalI);

        // initialize row size
        int rowSize = 1;

        // set perimeter to zero
        cellDataI.perimeter() = 0;

        // prepare storage for all found 3rd cells
        std::vector<int> foundAdditionals;

        int numberOfIntersections = 0;
        // run through all intersections with neighbors
        for (const auto& intersection : intersections(problem().gridView(), element))
        {
            cellDataI.perimeter() += intersection.geometry().volume();
            numberOfIntersections++;
            if (intersection.neighbor())
            {
                rowSize++;

                // special treatment for hanging nodes in the mpfa case
                if (enableMPFA && (element.level() < intersection.outside().level()))
                {
                    // each cell might have 4 hanging nodes with 4 irregular neighbors each
                    // get global ID of Interface from larger cell
                    int intersectionID = problem().grid().localIdSet().subId(element,
                            intersection.indexInInside(), 1);
                    //index outside
                    int eIdxGlobalJ = problem().variables().index(intersection.outside());

                    // add Entry of current neighbor cell to the IS seen from large cell
                    irregularCellMap_[intersectionID].push_back(eIdxGlobalJ);
                }
            }
        }
        cellDataI.fluxData().resize(numberOfIntersections);
        this->A_.incrementrowsize(eIdxGlobalI, rowSize);
    } //end first loop (that already reserved enough space for the MPFA connections on hanging nodes)

    // second loop to determine which matrix entries will be occupied
    if(enableMPFA && maxInteractionVolumes>1)
    {
        //prepare map for additional cell-connection through mpfa
        std::multimap<int, int> addionalRelations;
        using IntPair = std::pair<int,int>;
        std::pair<std::multimap<int,int>::iterator,std::multimap<int,int>::iterator> range;
        std::multimap<int,int>::iterator rangeIt;

        // second loop for further sub-faces
        for (const auto& element : elements(problem().gridView()))
        {
            // cell index
            int eIdxGlobalI = problem().variables().index(element);
            // run through all intersections with neighbors
            const auto isEndIt = problem().gridView().iend(element);
            for (auto isIt = problem().gridView().ibegin(element); isIt != isEndIt; ++isIt)
            {
                const auto& intersection = *isIt;

                if (intersection.neighbor())
                {
                    //index outside
                    int eIdxGlobalJ = problem().variables().index(intersection.outside());

                    // if mpfa is used, more entries might be needed if all interactionRegions are regarded
                    if (intersection.outside().level() > element.level()) //look from larger cell
                    {
                        // Prepare MPFA
                        /* get geometric Info, transmissibility matrix */
                        GlobalPosition globalPos3(0.);
                        int eIdxGlobal3=-1;
                        GlobalPosition globalPos4(0.);
                        int eIdxGlobal4=-1;
                        TransmissivityMatrix T(0.);

                        int interactionRegions
                            = problem().variables().getMpfaData3D(intersection, T,
                                                    globalPos3, eIdxGlobal3, globalPos4, eIdxGlobal4  );
                        if (interactionRegions == 0)
                            interactionRegions = problem().pressureModel().computeTransmissibilities(isIt,T,
                                    globalPos3, eIdxGlobal3,  globalPos4, eIdxGlobal4  );
                        if(!interactionRegions)
                            Dune::dgrave << "something went wrong getting mpfa data on cell " << eIdxGlobalI << std::endl;
                        if (interactionRegions == 1) // no second subface
                            continue;

                        //loop over all found interaction regions
                        for (int cocumber=1; cocumber<interactionRegions; cocumber++ )
                        {
                            problem().variables().getMpfaData3D(intersection, T,
                                   globalPos3, eIdxGlobal3, globalPos4, eIdxGlobal4, cocumber);
                            // indices
                            int additionalIdx2 = eIdxGlobal3;
                            int additionalIdx3 = eIdxGlobal4;

                            bool addIndex = true;

                            // check if additionals are "normal" neighbors of I
                            bool additional2isNeighbor(false), additional3isNeighbor(false);
                            // run through all intersections with neighbors if eIt
                            for (const auto& checkIntersection
                                 : intersections(problem().gridView(), element))
                            {
                                if (checkIntersection.neighbor())
                                {
                                    if(additionalIdx2==problem().variables().index(checkIntersection.outside()))
                                        additional2isNeighbor = true;
                                    if(additionalIdx3 == problem().variables().index(checkIntersection.outside()))
                                        additional3isNeighbor = true;
                                }

                            }
                            // if not "normal" neighbor, increase row size
                            if(!additional2isNeighbor)
                            {
                                // check if relation not already added
                                IntPair intPair(eIdxGlobalI,additionalIdx2);
                                if(eIdxGlobalI > additionalIdx2)
                                {
                                    using std::swap;
                                    swap(intPair.first, intPair.second);
                                }
                                range = addionalRelations.equal_range(intPair.first);
                                for (rangeIt=range.first; range.first!=range.second
                                                          && rangeIt!=range.second; ++rangeIt)
                                    if((*rangeIt).second == intPair.second)
                                        addIndex = false;
                                if(addIndex)
                                {
                                    this->A_.incrementrowsize(eIdxGlobalI);
                                    // add space for additional itsself
                                    this->A_.incrementrowsize(additionalIdx2);
                                    // mark relation as added
                                    addionalRelations.insert(intPair);
                                }
                            }

                            // if not "normal" neighbor, increase row size
                            if(!additional3isNeighbor)
                            {
                                addIndex = true;
                                // check if relation not already added
                                IntPair intPair(eIdxGlobalI,additionalIdx3);
                                if(eIdxGlobalI > additionalIdx3)
                                {
                                    using std::swap;
                                    swap(intPair.first, intPair.second);
                                }
                                range = addionalRelations.equal_range(intPair.first);
                                for (rangeIt=range.first; range.first!=range.second
                                                          && rangeIt!=range.second; ++rangeIt)
                                    if((*rangeIt).second == intPair.second)
                                        addIndex = false;
                                if(addIndex)
                                {
                                    this->A_.incrementrowsize(eIdxGlobalI);
                                    // add space for additional itsself
                                    this->A_.incrementrowsize(additionalIdx3);
                                    // mark relation as added
                                    addionalRelations.insert(intPair);
                                }
                            }

                            //reset bools for J
                            additional2isNeighbor = additional3isNeighbor = false;
                            // run through all intersections with neighbors of J
                            for (const auto& checkIntersection
                                 : intersections(problem().gridView(), intersection.outside()))
                            {
                                if (checkIntersection.neighbor())
                                {
                                    if(additionalIdx2 == problem().variables().index(checkIntersection.outside()))
                                        additional2isNeighbor = true;
                                    if(additionalIdx3 == problem().variables().index(checkIntersection.outside()))
                                        additional3isNeighbor = true;
                                }
                            }

                            // if not "normal" neighbor, increase row size
                            if(!additional2isNeighbor)
                            {
                                addIndex = true;
                                // check if relation not already added
                                IntPair intPair(eIdxGlobalJ,additionalIdx2);
                                if(eIdxGlobalJ > additionalIdx2)
                                {
                                    using std::swap;
                                    swap(intPair.first, intPair.second);
                                }
                                range = addionalRelations.equal_range(intPair.first);
                                for (rangeIt=range.first; range.first!=range.second
                                                          && rangeIt!=range.second; ++rangeIt)
                                    if((*rangeIt).second == intPair.second)
                                        addIndex = false;
                                if(addIndex)
                                {
                                    this->A_.incrementrowsize(eIdxGlobalJ);
                                    // add space for additional itsself
                                    this->A_.incrementrowsize(additionalIdx2);
                                    // mark relation as added
                                    addionalRelations.insert(intPair);
                                }
                            }

                            // if not "normal" neighbor, increase row size
                            if(!additional3isNeighbor)
                            {
                                addIndex = true;
                                // check if relation not already added
                                IntPair intPair(eIdxGlobalJ,additionalIdx3);
                                if(eIdxGlobalJ > additionalIdx3)
                                {
                                    using std::swap;
                                    swap(intPair.first, intPair.second);
                                }
                                range = addionalRelations.equal_range(intPair.first);
                                for (rangeIt=range.first; range.first!=range.second
                                                          && rangeIt!=range.second; ++rangeIt)
                                    if((*rangeIt).second == intPair.second)
                                        addIndex = false;
                                if(addIndex)
                                {
                                    this->A_.incrementrowsize(eIdxGlobalJ);
                                    // add space for additional itsself
                                    this->A_.incrementrowsize(additionalIdx3);
                                    // mark relation as added
                                    addionalRelations.insert(intPair);
                                }
                            }
                        }
                    }
                }
            }
        } //end second loop
    }
}

//!Determine position of matrix entries
/* Method adds TPFA and MPFA matrix entries
 */
template<class TypeTag>
void FV3dPressure2P2CAdaptive<TypeTag>::initializeMatrixIndices()
{
    // determine position of matrix entries
    for (const auto& element : elements(problem().gridView()))
    {
        // cell index
        int eIdxGlobalI = problem().variables().index(element);

        // add diagonal index
        this->A_.addindex(eIdxGlobalI, eIdxGlobalI);

        // run through all intersections with neighbors
        const auto isEndIt = problem().gridView().iend(element);
        for (auto isIt = problem().gridView().ibegin(element); isIt != isEndIt; ++isIt)
        {
            const auto& intersection = *isIt;

            if (intersection.neighbor())
            {
                // access neighbor
                int eIdxGlobalJ = problem().variables().index(intersection.outside());

                // add off diagonal index
                this->A_.addindex(eIdxGlobalI, eIdxGlobalJ);

                // special treatment for hanging nodes in the mpfa case
                if (enableMPFA && (element.level() < intersection.outside().level()))
                {
                    // prepare stuff to enter transmissibility calculation
                    GlobalPosition globalPos3(0.);
                    int eIdxGlobal3=-1;
                    GlobalPosition globalPos4(0.);
                    int eIdxGlobal4=-1;
                    TransmissivityMatrix T(0.);
                    TransmissivityMatrix additionalT(0.);

                    int interactionRegions
                        = problem().variables().getMpfaData3D(intersection, T, globalPos3, eIdxGlobal3, globalPos4, eIdxGlobal4  );
                    if (interactionRegions == 0)
                        interactionRegions = problem().pressureModel().computeTransmissibilities(isIt,T,
                                globalPos3, eIdxGlobal3,  globalPos4, eIdxGlobal4  );

                    for (int cocumber=1; cocumber<interactionRegions; cocumber++ )
                    {
                        problem().variables().getMpfaData3D(intersection, T,
                               globalPos3, eIdxGlobal3, globalPos4, eIdxGlobal4, cocumber);

                        // add off diagonal index in both directions!!
                        this->A_.addindex(eIdxGlobalI, eIdxGlobal3);
                        this->A_.addindex(eIdxGlobal3, eIdxGlobalI);
                        this->A_.addindex(eIdxGlobalI, eIdxGlobal4);
                        this->A_.addindex(eIdxGlobal4, eIdxGlobalI);
                        this->A_.addindex(eIdxGlobalJ, eIdxGlobal3);
                        this->A_.addindex(eIdxGlobal3, eIdxGlobalJ);
                        this->A_.addindex(eIdxGlobalJ, eIdxGlobal4);
                        this->A_.addindex(eIdxGlobal4, eIdxGlobalJ);
                    }
                }
            }
        }
    }
}

//! \copydoc FV2dPressure2P2CAdaptive::assemble()
template<class TypeTag>
void FV3dPressure2P2CAdaptive<TypeTag>::assemble(bool first)
{
    if(first)
    {
        BaseType::assemble(true);
        return;
    }

    // initialization: set matrix A_ to zero
    this->A_ = 0;
    this->f_ = 0;

    for (const auto& element : elements(problem().gridView()))
    {
        // get the global index of the cell
        int eIdxGlobalI = problem().variables().index(element);

        // assemble interior element contributions
        if (element.partitionType() == Dune::InteriorEntity)
        {
                // get the cell data
            CellData& cellDataI = problem().variables().cellData(eIdxGlobalI);

            Dune::FieldVector<Scalar, 2> entries(0.);

            /*****  source term ***********/
#ifndef noMultiphysics
            if(cellDataI.subdomain() != 2)
                problem().pressureModel().get1pSource(entries,element, cellDataI);
            else
#endif
                problem().pressureModel().getSource(entries,element, cellDataI, first);

            this->f_[eIdxGlobalI] += entries[rhs];

            /*****  flux term ***********/
            // iterate over all faces of the cell
            const auto isEndIt = problem().gridView().iend(element);
            for (auto isIt = problem().gridView().ibegin(element); isIt != isEndIt; ++isIt)
            {
                const auto& intersection = *isIt;

                /************* handle interior face *****************/
                if (intersection.neighbor())
                {
                    auto neighbor = intersection.outside();
                    int eIdxGlobalJ = problem().variables().index(neighbor);
                    //check for hanging nodes
                    //take a hanging node never from the element with smaller level!
                    bool haveSameLevel = (element.level() == neighbor.level());
                    // calculate only from one side, but add matrix entries for both sides
                    // the last condition is needed to properly assemble in the presence
                    // of ghost elements
                    if (getPropValue<TypeTag, Properties::VisitFacesOnlyOnce>()
                        && (eIdxGlobalI > eIdxGlobalJ) && haveSameLevel
                        && neighbor.partitionType() == Dune::InteriorEntity)
                        continue;

                    entries = 0;
                    //check for hanging nodes
                    if(!haveSameLevel && enableMPFA)
                    {
                        if (cellDataI.subdomain() != 2
                            || problem().variables().cellData(eIdxGlobalJ).subdomain() != 2) // cell in the 1p domain
                        {
                            asImp_().get1pMpfaFlux(isIt, cellDataI);
                        }
                        else
                        {
                            asImp_().getMpfaFlux(isIt, cellDataI);
                        }
                    }
                    else
                    {
                        CellData cellDataJ = problem().variables().cellData(eIdxGlobalJ);
                        if (cellDataI.subdomain() != 2
                            || problem().variables().cellData(eIdxGlobalJ).subdomain() != 2) // cell in the 1p domain
                        {
                            asImp_().get1pFlux(entries, intersection, cellDataI);
                        }
                        else
                        {
                            asImp_().getFlux(entries, intersection, cellDataI, first);
                        }

                        //set right hand side
                        this->f_[eIdxGlobalI] -= entries[rhs];

                        // set diagonal entry
                        this->A_[eIdxGlobalI][eIdxGlobalI] += entries[matrix];

                            // set off-diagonal entry
                        this->A_[eIdxGlobalI][eIdxGlobalJ] -= entries[matrix];

                        // The second condition is needed to not spoil the ghost element entries
                        if (getPropValue<TypeTag, Properties::VisitFacesOnlyOnce>()
                            && neighbor.partitionType() == Dune::InteriorEntity)
                        {
                            this->f_[eIdxGlobalJ] += entries[rhs];
                            this->A_[eIdxGlobalJ][eIdxGlobalJ] += entries[matrix];
                            this->A_[eIdxGlobalJ][eIdxGlobalI] -= entries[matrix];
                        }
                    }
                } // end neighbor

                /************* boundary face ************************/
                else
                {
                    entries = 0;
                    if (cellDataI.subdomain() != 2) //the current cell in the 1p domain
                        asImp_().get1pFluxOnBoundary(entries, intersection, cellDataI);
                    else
                        asImp_().getFluxOnBoundary(entries, intersection, cellDataI, first);

                    //set right hand side
                    this->f_[eIdxGlobalI] += entries[rhs];
                    // set diagonal entry
                    this->A_[eIdxGlobalI][eIdxGlobalI] += entries[matrix];
                }
            } //end interfaces loop

            /*****  storage term ***********/
            if (cellDataI.subdomain() != 2) //the current cell in the 1p domain
                asImp_().get1pStorage(entries, element, cellDataI);
            else
                asImp_().getStorage(entries, element, cellDataI, first);

            this->f_[eIdxGlobalI] += entries[rhs];
            // set diagonal entry
            this->A_[eIdxGlobalI][eIdxGlobalI] += entries[matrix];
        }
        // assemble overlap and ghost element contributions
        else
        {
            this->A_[eIdxGlobalI] = 0.0;
            this->A_[eIdxGlobalI][eIdxGlobalI] = 1.0;
            this->f_[eIdxGlobalI] = this->pressure()[eIdxGlobalI];
        }
    } // end grid traversal
}

//! Compute flux through an irregular interface using a \a mpfa method
/*! A mpfa l-method is applied to calculate fluxes near hanging nodes, using:
 * \f[
      - \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha}
        \left( \sum_k \tau_{2k} p^t_{\alpha,k} + \varrho_{\alpha} \sum_k \tau_{2k} \mathbf{g}^T \mathbf{x}_{k} \right)
                \sum_{\kappa} X^{\kappa}_{\alpha} \frac{\partial v_{t}}{\partial C^{\kappa}}
      + \frac{ V_i}{U_i} \sum_{\alpha} \varrho_{\alpha} \lambda_{\alpha}
       \left( \sum_k \tau_{2k} p^t_{\alpha,k} + \varrho_{\alpha} \sum_k \tau_{2k} \mathbf{g}^T \mathbf{x}_{k} \right)
          \sum_{\kappa} X^{\kappa}_{\alpha} \frac{\frac{\partial v_{t,j}}{\partial C^{\kappa}_j}
          -\frac{\partial v_{t,i}}{\partial C^{\kappa}_i}}{\Delta x}
    \f]
 *
 * Depending on the applied interaction volumes (e.g. on the boundary only
 * a single interaction volume might be applied), the stencil also spans
 * over neighboring cells. The contribution in other cells than I or J make
 * it necessary that the matrix and rhs entries are filled up within this function.
 * \param isIt Iterator of the intersection between cell I and J
 * \param cellDataI Data of cell I
 */
template<class TypeTag>
void FV3dPressure2P2CAdaptive<TypeTag>::getMpfaFlux(const IntersectionIterator& isIt,
                                                    const CellData& cellDataI)
{
    const auto& intersection = *isIt;

    // acess Cell I
    auto elementI = intersection.inside();
    int eIdxGlobalI = problem().variables().index(elementI);

    // get global coordinate of cell center
    const GlobalPosition& globalPos = elementI.geometry().center();

    // cell volume & perimeter, assume linear map here
    Scalar volume = elementI.geometry().volume();
    Scalar perimeter = cellDataI.perimeter();

    // get absolute permeability
    DimMatrix permeabilityI(problem().spatialParams().intrinsicPermeability(elementI));

    // access neighbor
    auto neighbor = intersection.outside();
    int eIdxGlobalJ = problem().variables().index(neighbor);
    CellData& cellDataJ = problem().variables().cellData(eIdxGlobalJ);

    // gemotry info of neighbor
    const GlobalPosition& globalPosNeighbor = neighbor.geometry().center();

    // distance vector between barycenters
    GlobalPosition distVec = globalPosNeighbor - globalPos;

    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    GlobalPosition unitDistVec(distVec);
    unitDistVec /= dist;

    DimMatrix permeabilityJ
        = problem().spatialParams().intrinsicPermeability(neighbor);

    // compute vectorized permeabilities
    DimMatrix meanPermeability(0);
    harmonicMeanMatrix(meanPermeability, permeabilityI, permeabilityJ);

    Dune::FieldVector<Scalar, dim> permeability(0);
    meanPermeability.mv(unitDistVec, permeability);

    // get average density for gravity flux
    PhaseVector rhoMean(0.);
    for (int phaseIdx=0; phaseIdx<NumPhases; phaseIdx++)
        rhoMean[phaseIdx] =0.5 * (cellDataI.density(phaseIdx)+cellDataJ.density(phaseIdx));

    // reset potential gradients
    Dune::FieldVector<Scalar,2> potential(0.);

    // determine volume derivatives in neighbor
    if (!cellDataJ.hasVolumeDerivatives())
        asImp_().volumeDerivatives(globalPosNeighbor, neighbor);

    ComponentVector dv_dC(0.), graddv_dC(0.);
    for (int compIdx = 0; compIdx < NumComponents; ++compIdx)
    {
        dv_dC[compIdx]= (cellDataJ.dv(compIdx) // dV/dm1= dv/dC^1
                + cellDataI.dv(compIdx)) * 0.5;
        graddv_dC[compIdx] = (cellDataJ.dv(compIdx)
                                - cellDataI.dv(compIdx)) / dist;
    }
//                    potentialW = problem().variables().potentialWetting(eIdxGlobalI, isIndex);
//                    potentialNW = problem().variables().potentialNonwetting(eIdxGlobalI, isIndex);
//
//                    densityW = (potentialW > 0.) ? densityWI : densityWJ;
//                    densityNW = (potentialNW > 0.) ? densityNWI : densityNWJ;
//
//                    densityW = (potentialW == 0.) ? rhoMeanW : densityW;
//                    densityNW = (potentialNW == 0.) ? rhoMeanNW : densityNW;
    //jochen: central weighting for gravity term

        // Prepare MPFA
        /* get geometrical Info, transmissibility matrix */
        GlobalPosition globalPos3(0.);
        int eIdxGlobal3=-1;
        GlobalPosition globalPos4(0.);
        int eIdxGlobal4=-1;
        TransmissivityMatrix T(0.);

            // prepare second interaction region
            GlobalPosition globalPosAdditional3(0.);
            int eIdxGlobalAdditional3=-1;
            GlobalPosition globalPosAdditional4(0.);
            int eIdxGlobalAdditional4=-1;

            TransmissivityMatrix additionalT(0.);

        int interactionRegions
            = problem().variables().getMpfaData3D(intersection, T,
                                    globalPos3, eIdxGlobal3, globalPos4, eIdxGlobal4  );
        if (interactionRegions == 0)
            interactionRegions = problem().pressureModel().computeTransmissibilities(isIt,T,
                    globalPos3, eIdxGlobal3,  globalPos4, eIdxGlobal4  );
        if(!interactionRegions)
            Dune::dgrave << "something went wrong getting mpfa data on cell " << eIdxGlobalI << std::endl;

        // shortcurts mpfa case
        CellData& cellData3 = problem().variables().cellData(eIdxGlobal3);
        CellData& cellData4 = problem().variables().cellData(eIdxGlobal4);
        Scalar temp1 = globalPos * this->gravity_;
        Scalar temp2 = globalPosNeighbor * this->gravity_;
        Scalar temp3 = globalPos3 * this->gravity_;
        Scalar temp4 = globalPos4 * this->gravity_;

        for (int phaseIdx = 0; phaseIdx < NumPhases; ++phaseIdx)
        {
            potential[phaseIdx] = (cellDataI.pressure(phaseIdx)-temp1*rhoMean[phaseIdx]) * T[0]
                             + (cellDataJ.pressure(phaseIdx)-temp2*rhoMean[phaseIdx]) * T[1]
                             + (cellData3.pressure(phaseIdx)-temp3*rhoMean[phaseIdx]) * T[2]
                             + (cellData4.pressure(phaseIdx)-temp4*rhoMean[phaseIdx]) * T[3];
        }

        // regard more interaction regions, if there are more
        if(interactionRegions != 1)
        {
            for(int banana = 1; banana < interactionRegions; banana ++)
            {
                // get data for second interaction region
                problem().variables().getMpfaData3D(intersection, additionalT,
                                globalPosAdditional3, eIdxGlobalAdditional3,
                                globalPosAdditional4, eIdxGlobalAdditional4 ,
                                banana); // offset for second interaction region

                Scalar gravityContributionAdditonal
                    = temp1 * additionalT[0] + temp2 * additionalT[1]
                        + globalPosAdditional3*this->gravity_ * additionalT[2]
                        + globalPosAdditional4*this->gravity_ * additionalT[3];
                CellData& cellDataA3 = problem().variables().cellData(eIdxGlobalAdditional3);
                CellData& cellDataA4 = problem().variables().cellData(eIdxGlobalAdditional4);

                for (int phaseIdx = 0; phaseIdx < NumPhases; ++phaseIdx)
                {
                    potential[phaseIdx] += cellDataI.pressure(phaseIdx) * additionalT[0]
                                     + cellDataJ.pressure(phaseIdx) * additionalT[1]
                                     + cellDataA3.pressure(phaseIdx) * additionalT[2]
                                     + cellDataA4.pressure(phaseIdx) * additionalT[3];
                    potential[phaseIdx] -= gravityContributionAdditonal * rhoMean[phaseIdx];
                }
            }
        }

        // initialize convenience shortcuts
        PhaseVector lambda(0.), dV(0.), gV(0.);

        //do the upwinding of the mobility depending on the phase potentials
        std::vector<const CellData*> upwindCellData(NumPhases);

        for (int phaseIdx = 0; phaseIdx < NumPhases; ++phaseIdx)
        {
            int eqIdx = phaseIdx + 1;
            if (potential[phaseIdx] > 0.)
                upwindCellData[phaseIdx] = &cellDataI;
            else if (potential[phaseIdx] < 0.)
                upwindCellData[phaseIdx] = &cellDataJ;
            else
            {
                if(cellDataI.isUpwindCell(intersection.indexInInside(), eqIdx))
                    upwindCellData[phaseIdx] = &cellDataI;
                else if(cellDataJ.isUpwindCell(intersection.indexInOutside(), eqIdx))
                    upwindCellData[phaseIdx] = &cellDataJ;
                //else
                //  upwinding is not done!
            }

            //perform upwinding if desired
            if(!upwindCellData[phaseIdx])
            {
                Scalar averagedMassFraction[NumComponents];
                for (int compIdx = 0; compIdx < NumComponents; ++compIdx)
                    averagedMassFraction[compIdx]
                         = harmonicMean(cellDataI.massFraction(phaseIdx, compIdx), cellDataJ.massFraction(phaseIdx, compIdx));
                Scalar averageDensity = harmonicMean(cellDataI.density(phaseIdx), cellDataJ.density(phaseIdx));

                for (int compIdx = 0; compIdx < NumComponents; ++compIdx)
                {
                    dV[phaseIdx] += dv_dC[compIdx] * averagedMassFraction[compIdx];
                    gV[phaseIdx] += graddv_dC[compIdx] * averagedMassFraction[compIdx];
                }
                dV[phaseIdx] *= averageDensity;
                gV[phaseIdx] *= averageDensity;
                lambda[phaseIdx] = harmonicMean(cellDataI.mobility(phaseIdx), cellDataJ.mobility(phaseIdx));
            }
            else //perform upwinding
            {
                for (int compIdx = 0; compIdx < NumComponents; ++compIdx)
                {
                     dV[phaseIdx] += dv_dC[compIdx] * upwindCellData[phaseIdx]->massFraction(phaseIdx, compIdx);
                     gV[phaseIdx] += graddv_dC[compIdx] * upwindCellData[phaseIdx]->massFraction(phaseIdx, compIdx);
                }
                lambda[phaseIdx] = upwindCellData[phaseIdx]->mobility(phaseIdx);
                dV[phaseIdx] *= upwindCellData[phaseIdx]->density(phaseIdx);
                gV[phaseIdx] *= upwindCellData[phaseIdx]->density(phaseIdx);
            }
        }

    /* compute matrix entry: advective fluxes */
    /* extend T with other matrix entries and assemble to A_    */
    this->A_[eIdxGlobalI][eIdxGlobalI] += (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * T[0];
    this->A_[eIdxGlobalI][eIdxGlobalJ] += (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * T[1];
    this->A_[eIdxGlobalI][eIdxGlobal3] += (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * T[2];
    this->A_[eIdxGlobalI][eIdxGlobal4] += (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * T[3];


    // add gravity to RHS vector
    Scalar gravityContribution = temp1 * T[0] + temp2 * T[1] + temp3 * T[2] + temp4 * T[3];
    this->f_[eIdxGlobalI] += (rhoMean[wPhaseIdx] * lambda[wPhaseIdx] * dV[wPhaseIdx]
                              + rhoMean[nPhaseIdx] * lambda[nPhaseIdx] * dV[nPhaseIdx]) * gravityContribution;

    // weithing accounts for the fraction of the subcontrol volume
    Scalar weightingFactor = volume / perimeter;    // transforms flux through area A -> V * A/perimeter
    if(enableVolumeIntegral_) // switch off volume integral for mpfa case
    {
        // correct for area integral
        this->A_[eIdxGlobalI][eIdxGlobalI] -=
                weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * T[0];
        this->A_[eIdxGlobalI][eIdxGlobalJ] -=
                weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * T[1];
        this->A_[eIdxGlobalI][eIdxGlobal3] -=
                weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * T[2];
        this->A_[eIdxGlobalI][eIdxGlobal4] -=
                weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * T[3];

        // add gravity to RHS vector
        this->f_[eIdxGlobalI] -= weightingFactor * gravityContribution *
                (rhoMean[wPhaseIdx] * lambda[wPhaseIdx] * gV[wPhaseIdx] + rhoMean[nPhaseIdx] * lambda[nPhaseIdx] * gV[nPhaseIdx]);
    }

    // capillary pressure flux
    Scalar pcGradient = cellDataI.capillaryPressure() * T[0]
                       + cellDataJ.capillaryPressure() * T[1]
                       + cellData3.capillaryPressure() * T[2]
                       + cellData4.capillaryPressure() * T[3];

    if (this->pressureType == pw)
        pcGradient *= + lambda[nPhaseIdx] * dV[nPhaseIdx]
                        - enableVolumeIntegral_ * weightingFactor * lambda[nPhaseIdx] * gV[nPhaseIdx];
    else if (this->pressureType == pn)
        pcGradient *= - lambda[wPhaseIdx] * dV[wPhaseIdx]
                        + enableVolumeIntegral_ * weightingFactor * lambda[wPhaseIdx] * gV[wPhaseIdx];

    this->f_[eIdxGlobalI] += pcGradient;

    // regard more interaction regions, if there are more
    if(interactionRegions != 1)
    {
        for(int banana = 1; banana < interactionRegions; banana ++)
        {
            // get data for second interaction region
            problem().variables().getMpfaData3D(intersection, additionalT,
                            globalPosAdditional3, eIdxGlobalAdditional3,
                            globalPosAdditional4, eIdxGlobalAdditional4 ,
                            banana); // offset for second interaction region

            Scalar gravityContributionAdditonal
                = temp1 * additionalT[0] + temp2 * additionalT[1]
                    + globalPosAdditional3*this->gravity_ * additionalT[2]
                    + globalPosAdditional4*this->gravity_ * additionalT[3];
            CellData& cellDataA3 = problem().variables().cellData(eIdxGlobalAdditional3);
            CellData& cellDataA4 = problem().variables().cellData(eIdxGlobalAdditional4);

            /* compute matrix entry: advective fluxes */
            /* extend T with other matrix entries and assemble to A_    */
            this->A_[eIdxGlobalI][eIdxGlobalI] +=
                    (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * additionalT[0];
            this->A_[eIdxGlobalI][eIdxGlobalJ] +=
                    (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * additionalT[1];
            this->A_[eIdxGlobalI][eIdxGlobalAdditional3] +=
                    (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * additionalT[2];
            this->A_[eIdxGlobalI][eIdxGlobalAdditional4] +=
                    (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * additionalT[3];


            // add gravity to RHS vector
            this->f_[eIdxGlobalI] += (rhoMean[wPhaseIdx] * lambda[wPhaseIdx] * dV[wPhaseIdx]
                                     + rhoMean[nPhaseIdx] * lambda[nPhaseIdx] * dV[nPhaseIdx]) * gravityContributionAdditonal;

            // weithing accounts for the fraction of the subcontrol volume
            if(enableVolumeIntegral_) // switch off volume integral for mpfa case
            {
                // correct for area integral
                this->A_[eIdxGlobalI][eIdxGlobalI] -=
                        weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * additionalT[0];
                this->A_[eIdxGlobalI][eIdxGlobalJ] -=
                        weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * additionalT[1];
                this->A_[eIdxGlobalI][eIdxGlobalAdditional3] -=
                        weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * additionalT[2];
                this->A_[eIdxGlobalI][eIdxGlobalAdditional4] -=
                        weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * additionalT[3];

                // add gravity to RHS vector
                this->f_[eIdxGlobalI] -= weightingFactor * gravityContribution *
                        (rhoMean[wPhaseIdx] * lambda[wPhaseIdx] * gV[wPhaseIdx] + rhoMean[nPhaseIdx] * lambda[nPhaseIdx] * gV[nPhaseIdx]);
            }

            // capillary pressure flux
            Scalar addPcGradient = cellDataI.capillaryPressure() * additionalT[0]
                                 + cellDataJ.capillaryPressure() * additionalT[1]
                                 + cellDataA3.capillaryPressure() * additionalT[2]
                                 + cellDataA4.capillaryPressure() * additionalT[3];

            if (this->pressureType == pw)
                addPcGradient *= + lambda[nPhaseIdx] * dV[nPhaseIdx]
                               - enableVolumeIntegral_ * weightingFactor * lambda[nPhaseIdx] * gV[nPhaseIdx];
            else if (this->pressureType == pn)
                addPcGradient *= - lambda[wPhaseIdx] * dV[wPhaseIdx]
                               + enableVolumeIntegral_ * weightingFactor * lambda[wPhaseIdx] * gV[wPhaseIdx];

            this->f_[eIdxGlobalI] += addPcGradient;
        }
    }
}


/*!
 * \brief Compute single-phase flux through an irregular interface using a \a mpfa method
 *
 * A mpfa l-method is applied to calculate fluxes near hanging nodes for
 * multiphysics models, using:
 * \f[
      - \lambda_{\alpha}
        \left( \sum_k \tau_{2k} p^t_{\alpha,k} + \varrho_{\alpha} \sum_k \tau_{2k} \mathbf{g}^T \mathbf{x}_{k} \right)
                \sum_{\kappa} X^{\kappa}_{\alpha}
    \f]
 *
 * Depending on the applied interaction volumes (e.g. on the boundary only
 * a single interaction volume might be applied), the stencil also spans
 * over neighboring cells. The contribution in other cells than I or J make
 * it necessary that the matrix and rhs entries are filled up within this function.
 * \param isIt Iterator of the intersection between cell I and J
 * \param cellDataI Data of cell I
 */
template<class TypeTag>
void FV3dPressure2P2CAdaptive<TypeTag>::get1pMpfaFlux(const IntersectionIterator& isIt,
                                                    const CellData& cellDataI)
{
    const auto& intersection = *isIt;

    // acess Cell I
    auto elementI = intersection.inside();
    int eIdxGlobalI = problem().variables().index(elementI);

    // get global coordinate of cell center
    const GlobalPosition& globalPos = elementI.geometry().center();

    // access neighbor
    auto neighbor = intersection.outside();
    int eIdxGlobalJ = problem().variables().index(neighbor);
    CellData& cellDataJ = problem().variables().cellData(eIdxGlobalJ);

    // gemotry info of neighbor
    const GlobalPosition& globalPosNeighbor = neighbor.geometry().center();

    // due to "safety cell" around subdomain, both cells I and J
    // have single-phase conditions, although one is in 2p domain.
    using std::min;
    int phaseIdx = min(cellDataI.subdomain(), cellDataJ.subdomain());

    // get average density for gravity flux
    Scalar rhoMean = 0.5 * (cellDataI.density(phaseIdx) + cellDataJ.density(phaseIdx));

    // 1p => no pC => only 1 pressure, potential
    Scalar potential = 0.;

        // Prepare MPFA
        /** get geometrical Info, transmissibility matrix */
        GlobalPosition globalPos3(0.);
        int eIdxGlobal3=-1;
        GlobalPosition globalPos4(0.);
        int eIdxGlobal4=-1;
        TransmissivityMatrix T(0.);

            // prepare second interaction region
            GlobalPosition globalPosAdditional3(0.);
            int eIdxGlobalAdditional3=-1;
            GlobalPosition globalPosAdditional4(0.);
            int eIdxGlobalAdditional4=-1;

            TransmissivityMatrix additionalT(0.);

        int interactionRegions
            = problem().variables().getMpfaData3D(intersection, T,
                                    globalPos3, eIdxGlobal3, globalPos4, eIdxGlobal4  );
        if (interactionRegions == 0)
            interactionRegions = problem().pressureModel().computeTransmissibilities(isIt,T,
                    globalPos3, eIdxGlobal3,  globalPos4, eIdxGlobal4  );
        if(!interactionRegions)
            Dune::dgrave << "something went wrong getting mpfa data on cell " << eIdxGlobalI << std::endl;

        // shortcurts mpfa case
        CellData& cellData3 = problem().variables().cellData(eIdxGlobal3);
        CellData& cellData4 = problem().variables().cellData(eIdxGlobal4);
        Scalar temp1 = globalPos * this->gravity_;
        Scalar temp2 = globalPosNeighbor * this->gravity_;
        Scalar temp3 = globalPos3 * this->gravity_;
        Scalar temp4 = globalPos4 * this->gravity_;

        potential = (cellDataI.pressure(phaseIdx)-temp1*rhoMean) * T[0]
                     + (cellDataJ.pressure(phaseIdx)-temp2*rhoMean) * T[1]
                     + (cellData3.pressure(phaseIdx)-temp3*rhoMean) * T[2]
                     + (cellData4.pressure(phaseIdx)-temp4*rhoMean) * T[3];

        // regard more interaction regions, if there are more
        if(interactionRegions != 1)
        {
            for(int banana = 1; banana < interactionRegions; banana ++)
            {
                // get data for second interaction region
                problem().variables().getMpfaData3D(intersection, additionalT,
                                globalPosAdditional3, eIdxGlobalAdditional3,
                                globalPosAdditional4, eIdxGlobalAdditional4 ,
                                banana); // offset for second interaction region

                Scalar gravityContributionAdditonal
                    = temp1 * additionalT[0] + temp2 * additionalT[1]
                        + globalPosAdditional3*this->gravity_ * additionalT[2]
                        + globalPosAdditional4*this->gravity_ * additionalT[3];
                CellData& cellDataA3 = problem().variables().cellData(eIdxGlobalAdditional3);
                CellData& cellDataA4 = problem().variables().cellData(eIdxGlobalAdditional4);

                potential += cellDataI.pressure(phaseIdx) * additionalT[0]
                             + cellDataJ.pressure(phaseIdx) * additionalT[1]
                             + cellDataA3.pressure(phaseIdx) * additionalT[2]
                             + cellDataA4.pressure(phaseIdx) * additionalT[3];
                potential -= gravityContributionAdditonal * rhoMean;
            }
        }

    // initialize convenience shortcuts
    Scalar lambda(0.);

    //do the upwinding of the mobility depending on the phase potentials
    if (potential >= 0.)
        lambda = cellDataI.mobility(phaseIdx);
    else
        lambda = cellDataJ.mobility(phaseIdx);

    /** compute matrix entry: advective fluxes */
    /* extend T with other matrix entries and assemble to A_    */
    this->A_[eIdxGlobalI][eIdxGlobalI] += lambda * T[0];
    this->A_[eIdxGlobalI][eIdxGlobalJ] += lambda * T[1];
    this->A_[eIdxGlobalI][eIdxGlobal3] += lambda * T[2];
    this->A_[eIdxGlobalI][eIdxGlobal4] += lambda * T[3];

    // add gravity to RHS vector
    Scalar gravityContribution = temp1 * T[0] + temp2 * T[1] + temp3 * T[2] + temp4 * T[3];
    this->f_[eIdxGlobalI] += lambda * rhoMean * gravityContribution;

    // regard more interaction regions, if there are more
    if(interactionRegions != 1)
    {
        for(int banana = 1; banana < interactionRegions; banana ++)
        {
            // get data for second interaction region
            problem().variables().getMpfaData3D(intersection, additionalT,
                            globalPosAdditional3, eIdxGlobalAdditional3,
                            globalPosAdditional4, eIdxGlobalAdditional4 ,
                            banana); // offset for second interaction region

            Scalar gravityContributionAdditonal
                = temp1 * additionalT[0] + temp2 * additionalT[1]
                    + globalPosAdditional3*this->gravity_ * additionalT[2]
                    + globalPosAdditional4*this->gravity_ * additionalT[3];

            /** compute matrix entry: advective fluxes */
            /* extend T with other matrix entries and assemble to A_    */
            this->A_[eIdxGlobalI][eIdxGlobalI] += lambda * additionalT[0];
            this->A_[eIdxGlobalI][eIdxGlobalJ] += lambda * additionalT[1];
            this->A_[eIdxGlobalI][eIdxGlobalAdditional3] += lambda * additionalT[2];
            this->A_[eIdxGlobalI][eIdxGlobalAdditional4] += lambda * additionalT[3];

            // add gravity to RHS vector
            this->f_[eIdxGlobalI] += lambda * rhoMean* gravityContributionAdditonal;
        }
    }
}

//! constitutive functions are updated once if new concentrations are calculated
//! and stored in the variables container
/*!
 * In contrast to the standard sequential 2p2c model, this method also holds routines
 * to adapt the subdomain. The subdomain indicates weather we are in 1p domain (value = 1)
 * or in the two phase subdomain (value = 2).
 * After the grid adaption, update of the subdomain is explicitly avoided.
 */
template<class TypeTag>
void FV3dPressure2P2CAdaptive<TypeTag>::updateMaterialLaws(bool fromPostTimestep)
{
//#define noMultiphysics
#ifdef noMultiphysics
    FVPressure2P2C<TypeTag>::updateMaterialLaws();
    return;
#endif

    if(!fromPostTimestep) // this happens after grid update
    {
        Scalar maxError = 0.;
        // iterate through leaf grid an evaluate c0 at cell center
        for (const auto& element : elements(problem().gridView()))
        {
            int eIdxGlobal = problem().variables().index(element);
            CellData& cellData = problem().variables().cellData(eIdxGlobal);

            if(cellData.fluidStateType() == 0) // i.e. it is complex
                problem().pressureModel().updateMaterialLawsInElement(element, fromPostTimestep);
            else
                problem().pressureModel().update1pMaterialLawsInElement(element, cellData, fromPostTimestep);
            using std::max;
            maxError = max(maxError, fabs(cellData.volumeError()));
        }
        this->maxError_ = maxError/problem().timeManager().timeStepSize();
    }
    else
    {
        // resize multiphysics vectors
        FVPressure2P2CMultiPhysics<TypeTag>::nextSubdomain.resize(problem().gridView().size(0));

        FVPressure2P2CMultiPhysics<TypeTag>::updateMaterialLaws();
    }
}

/*!
 * \brief Computes the transmissibility coefficients for the MPFA-l method in 3D
 *
 * For faces with a hanging node in 3D, there are four sub-faces. The first subface
 * contains a unique interaction volume, with is directly calculated by this method.
 * For the remainder of the sub-faces, the interaction volumes are build and calculated
 * by the common mpfa-l--implementation of the 2p models. The latter is established via the
 * protected method FV3dPressure2P2CAdaptive::transmissibilityAdapter_().
 * The calculated Transmissivity Matrices are (along with some geometric information)
 * stored for later use in Variableclass2p2cadaptive .
 *
 * \param isIt Iterator to the current intersection
 * \param T Transmissitivity matrix of the first unique interaction volume
 * \param globalPos4 Position of the 3rd cell (with local Idx 4) of the unique interaction volume
 * \param eIdxGlobal4 Index of the 3rd cell (with local Idx 4) of the unique interaction volume
 * \param globalPos6 Position of the 4th cell (with local Idx 6) of the unique interaction volume
 * \param eIdxGlobal6 Index of the 4th cell (with local Idx 6) of the unique interaction volume
 */
template <class TypeTag>
int FV3dPressure2P2CAdaptive<TypeTag>::computeTransmissibilities(const IntersectionIterator& isIt,
        TransmissivityMatrix& T,
        GlobalPosition& globalPos4,
        int& eIdxGlobal4,
        GlobalPosition& globalPos6,
        int& eIdxGlobal6)
{
    const auto& intersection = *isIt;

    // get geometry information of cellI = cell1, cellJ = cell2
    auto element = intersection.inside();
    auto neighbor = intersection.outside();
    GlobalPosition globalPos1 = element.geometry().center();
    GlobalPosition globalPos2 = neighbor.geometry().center();
    DimMatrix K1(problem().spatialParams().intrinsicPermeability(element));
    DimMatrix K2(problem().spatialParams().intrinsicPermeability(neighbor));

    // determine ID of intersection seen from larger cell
    int intersectionID = 0;
    if(intersection.inside().level() < intersection.outside().level())
        intersectionID = problem().grid().localIdSet().subId(element,
                intersection.indexInInside(), 1);
    else
        DUNE_THROW(Dune::NotImplemented, " ABORT, transmiss calculated from wrong side!!");

    std::vector<int> localIrregularCells = irregularCellMap_[intersectionID];

    /** 1) get geometrical information of interaction triangle   */
    // geometry and Data of face IJ in nomenclature of mpfa
    GlobalPosition globalPosFace12 = intersection.geometry().center();  // = 'x'1
    GlobalPosition outerNormaln12 = intersection.centerUnitOuterNormal();

    /**** a) search for intersection 24 and 26 ************/
    auto face24=isIt; // as long as face24 = isIt, it is still not found!
    auto face26=isIt; // as long as face26 = isIt, it is still not found!

    const auto isEndIt = problem().gridView().iend(neighbor);
    for (auto isIt2 = problem().gridView().ibegin(neighbor); isIt2 != isEndIt; ++isIt2)
    {
        const auto& intersection2 = *isIt2;

        // continue if no neighbor or arrived at intersection
        if(!(intersection2.neighbor()) || intersection2.outside() == element)
            continue;

        int currentNeighbor = problem().variables().index(intersection2.outside());

        // have we found face24?
        if (find(localIrregularCells.begin(), localIrregularCells.end(),
                currentNeighbor) != localIrregularCells.end() && face24==isIt)
            face24 = isIt2;
        else if (find(localIrregularCells.begin(), localIrregularCells.end(),
                currentNeighbor) != localIrregularCells.end() && face26==isIt)
        {
            // we now found both intersections, but we have to investigate orientation:
            // Calculate the vector product of the normals in current orientation
            GlobalPosition vectorProduct = crossProduct(face24->centerUnitOuterNormal(), intersection2.centerUnitOuterNormal());
            // the vector product of the two isNormals for the desired configuration points
            // in the direction of outerNormaln12 => ScalarProduct > 0.
            if((vectorProduct * outerNormaln12) > 0.)
                face26 = isIt2;
            else
            {
                // revert orientation: switch face24 and face26
                face26 = face24;
                face24 = isIt2;
            }
        }
    }

    // get information of cell4
    globalPos4 = face24->outside().geometry().center();
    eIdxGlobal4 = problem().variables().index(face24->outside());
    GlobalPosition outerNormaln24 = face24->centerUnitOuterNormal();
    // get absolute permeability of neighbor cell 3
    DimMatrix K4(problem().spatialParams().intrinsicPermeability(face24->outside()));

    // get information of cell6
    globalPos6 = face26->outside().geometry().center();
    eIdxGlobal6 = problem().variables().index(face26->outside());
    GlobalPosition outerNormaln26 = face26->centerUnitOuterNormal();
    // get absolute permeability of neighbor cell 3
    DimMatrix K6(problem().spatialParams().intrinsicPermeability(face26->outside()));

    /**** b) Get Points on the edges (in plane of isIt), 'x'4 and 'x'5 ***********/
    int localFace12 = intersection.indexInOutside(); // isIt is face12
    int localFace24 = face24->indexInInside();
    int localFace26 = face26->indexInInside();

    const auto refElement = referenceElement(neighbor);

    //find 'x'5 = edgeCoord1226
    int edge1226;
    // search through edges of face 12
    for(int nectarine=0; nectarine < refElement.size(localFace12, 1, dim-1); nectarine++)
    {
        // get local Idx of edge on face 12
        int localEdgeOn12 = refElement.subEntity(localFace12, 1, nectarine, dim-1);
        // search through edges of face 26
        for(int plum = 0; plum < refElement.size(localFace26, 1,dim-1); plum++)
        {
//            int localEdge26 = refElement.subEntity(localFace26, 1, plum, dim-1);
            if(refElement.subEntity(localFace12, 1, nectarine, dim-1)
                    == refElement.subEntity(localFace26, 1, plum, dim-1))
            {
                edge1226 = localEdgeOn12;
                break;
            }
        }
    }
    GlobalPosition edgeCoord1226 =  // 'x'5
            neighbor.geometry().global(refElement.position(edge1226, dim-1));

    //find 'x'4 = edgeCoord1224
    int edge1224;
    // search through edges of face 12
    for(int nectarine=0; nectarine < refElement.size(localFace12, 1, dim-1); nectarine++)
    {
        // get local Idx of edge on face 12
        int localEdgeOn12 = refElement.subEntity(localFace12, 1, nectarine, dim-1);
        // search through edges of face 24
        for(int plum = 0; plum < refElement.size(localFace24, 1, dim-1); plum++)
        {
            int localEdge24 = refElement.subEntity(localFace24, 1, plum, dim-1);
            if(localEdgeOn12 == localEdge24)
            {
                edge1224 = localEdgeOn12;
                break;
            }
        }
    }
    GlobalPosition edgeCoord1224 =  // 'x'4
            neighbor.geometry().global(refElement.position(edge1224, dim-1));

    //find 'x'6 = edgeCoord2426
    int edge2426;
    // search through edges of face 24
    for(int nectarine=0; nectarine < refElement.size(localFace24, 1, dim-1); nectarine++)
    {
        // get local Idx of edge on face 24
        int localEdgeOn24 = refElement.subEntity(localFace24, 1, nectarine, dim-1);
        // search through edges of face 26
        for(int plum = 0; plum < refElement.size(localFace26, 1, dim-1); plum++)
        {
            int localEdge26 = refElement.subEntity(localFace26, 1, plum, dim-1);
            if(localEdgeOn24 == localEdge26)
            {
                edge2426 = localEdgeOn24;
                break;
            }
        }
    }
    GlobalPosition edgeCoord2426 =   // 'x'6
            neighbor.geometry().global(refElement.position(edge2426, dim-1));

    /** 2) Calculate omega, chi for matrices  **/
    // center of face in global coordinates, i.e., the midpoint of face 'isIt24'
    GlobalPosition globalPosFace24 = face24->geometry().center();
    GlobalPosition globalPosFace26 = face26->geometry().center();

    // get face volumes
    //Area12 = | 'x'1->'x'5  x  'x'1->'x'4 |
    Scalar subFaceArea12 = crossProduct(edgeCoord1226-globalPosFace12, edgeCoord1224-globalPosFace12).two_norm();

    //Area24 = | 'x'2->'x'4  x  'x'2->'x'6 |
    Scalar subFaceArea24 = crossProduct(edgeCoord1224-globalPosFace24, edgeCoord2426-globalPosFace24).two_norm();

    //Area26 = | 'x'3->'x'5  x  'x'3->'x'6 |
    Scalar subFaceArea26 = crossProduct(edgeCoord1226-globalPosFace26, edgeCoord2426-globalPosFace26).two_norm();

    // compute normal vectors for case 2 (idx1, idx2, idx4, idx6)
    GlobalPosition nu11C2 = crossProduct(edgeCoord1226-globalPos1, globalPosFace12-globalPos1);
    GlobalPosition nu12C2 = crossProduct(edgeCoord1224-globalPos1, edgeCoord1226-globalPos1);
    GlobalPosition nu13C2 = crossProduct(globalPosFace12-globalPos1, edgeCoord1224-globalPos1);

    GlobalPosition nu21C2 = crossProduct(globalPosFace26-globalPos2, globalPosFace24-globalPos2);
    GlobalPosition nu22C2 = crossProduct(globalPosFace12-globalPos2, globalPosFace26-globalPos2);
    GlobalPosition nu23C2 = crossProduct(globalPosFace24-globalPos2, globalPosFace12-globalPos2);

    GlobalPosition nu41C2 = crossProduct(edgeCoord2426-globalPos4, edgeCoord1224-globalPos4);
    GlobalPosition nu42C2 = crossProduct(globalPosFace24-globalPos4, edgeCoord2426-globalPos4);
    GlobalPosition nu43C2 = crossProduct(edgeCoord1224-globalPos4, globalPosFace24-globalPos4);

    GlobalPosition nu61C2 = crossProduct(globalPosFace26-globalPos6, edgeCoord1226-globalPos6);
    GlobalPosition nu62C2 = crossProduct(edgeCoord2426-globalPos6, globalPosFace26-globalPos6);
    GlobalPosition nu63C2 = crossProduct(edgeCoord1226-globalPos6, edgeCoord2426-globalPos6);

    // compute T, i.e., the volume of the subtetrahedron
    Scalar T1C2 = (globalPosFace12-globalPos1) * crossProduct(edgeCoord1224-globalPos1, edgeCoord1226-globalPos1);
    Scalar T2C2 = (globalPosFace12-globalPos2) * crossProduct(globalPosFace26-globalPos2, globalPosFace24-globalPos2);
    Scalar T4C2 = (globalPosFace24-globalPos4) * crossProduct(edgeCoord2426-globalPos4, edgeCoord1224-globalPos4);
    Scalar T6C2 = (globalPosFace26-globalPos6) * crossProduct(edgeCoord1226-globalPos6, edgeCoord2426-globalPos6);

    // compute components needed for flux calculation, denoted as 'omega' and 'r'
    GlobalPosition K1nu11C2(0);
    K1.mv(nu11C2, K1nu11C2);
    GlobalPosition K1nu12C2(0);
    K1.mv(nu12C2, K1nu12C2);
    GlobalPosition K1nu13C2(0);
    K1.mv(nu13C2, K1nu13C2);

    GlobalPosition K2nu21C2(0);
    K2.mv(nu21C2, K2nu21C2);
    GlobalPosition K2nu22C2(0);
    K2.mv(nu22C2, K2nu22C2);
    GlobalPosition K2nu23C2(0);
    K2.mv(nu23C2, K2nu23C2);

    GlobalPosition K4nu41C2(0);
    K4.mv(nu41C2, K4nu41C2);
    GlobalPosition K4nu42C2(0);
    K4.mv(nu42C2, K4nu42C2);
    GlobalPosition K4nu43C2(0);
    K4.mv(nu43C2, K4nu43C2);

    GlobalPosition K6nu61C2(0);
    K6.mv(nu61C2, K6nu61C2);
    GlobalPosition K6nu62C2(0);
    K6.mv(nu62C2, K6nu62C2);
    GlobalPosition K6nu63C2(0);
    K6.mv(nu63C2, K6nu63C2);

    // upwinding of lambda is not included here to "recycle" it-> only possible if lambda does
    // not vary harshly. If it would, refinement indicator would not allow different grid levels here.
    Scalar omega111C2 = /*lambda12C2 */ (outerNormaln12 * K1nu11C2) * subFaceArea12/T1C2;
    Scalar omega112C2 = /*lambda12C2 */ (outerNormaln12 * K1nu12C2) * subFaceArea12/T1C2;
    Scalar omega113C2 = /*lambda12C2 */ (outerNormaln12 * K1nu13C2) * subFaceArea12/T1C2;

    Scalar omega121C2 = /*lambda21C2 */ (outerNormaln12 * K2nu21C2) * subFaceArea12/T2C2;
    Scalar omega122C2 = /*lambda21C2 */ (outerNormaln12 * K2nu22C2) * subFaceArea12/T2C2;
    Scalar omega123C2 = /*lambda21C2 */ (outerNormaln12 * K2nu23C2) * subFaceArea12/T2C2;

    Scalar omega221C2 = /*lambda24C2 */ (outerNormaln24 * K2nu21C2) * subFaceArea24/T2C2;
    Scalar omega222C2 = /*lambda24C2 */ (outerNormaln24 * K2nu22C2) * subFaceArea24/T2C2;
    Scalar omega223C2 = /*lambda24C2 */ (outerNormaln24 * K2nu23C2) * subFaceArea24/T2C2;

    Scalar omega241C2 = /*lambda42C2 */ (outerNormaln24 * K4nu41C2) * subFaceArea24/T4C2;
    Scalar omega242C2 = /*lambda42C2 */ (outerNormaln24 * K4nu42C2) * subFaceArea24/T4C2;
    Scalar omega243C2 = /*lambda42C2 */ (outerNormaln24 * K4nu43C2) * subFaceArea24/T4C2;

    Scalar omega321C2 = /*lambda26C2 */ (outerNormaln26 * K2nu21C2) * subFaceArea26/T2C2;
    Scalar omega322C2 = /*lambda26C2 */ (outerNormaln26 * K2nu22C2) * subFaceArea26/T2C2;
    Scalar omega323C2 = /*lambda26C2 */ (outerNormaln26 * K2nu23C2) * subFaceArea26/T2C2;

    Scalar omega361C2 = /*lambda62C2 */ (outerNormaln26 * K6nu61C2) * subFaceArea26/T6C2;
    Scalar omega362C2 = /*lambda62C2 */ (outerNormaln26 * K6nu62C2) * subFaceArea26/T6C2;
    Scalar omega363C2 = /*lambda62C2 */ (outerNormaln26 * K6nu63C2) * subFaceArea26/T6C2;

    Scalar r211C2 = (nu21C2 * (edgeCoord1224-globalPos2))/T2C2;
    Scalar r212C2 = (nu21C2 * (edgeCoord1226-globalPos2))/T2C2;
    Scalar r213C2 = (nu21C2 * (edgeCoord2426-globalPos2))/T2C2;

    Scalar r221C2 = (nu22C2 * (edgeCoord1224-globalPos2))/T2C2;
    Scalar r222C2 = (nu22C2 * (edgeCoord1226-globalPos2))/T2C2;
    Scalar r223C2 = (nu22C2 * (edgeCoord2426-globalPos2))/T2C2;

    Scalar r231C2 = (nu23C2 * (edgeCoord1224-globalPos2))/T2C2;
    Scalar r232C2 = (nu23C2 * (edgeCoord1226-globalPos2))/T2C2;
    Scalar r233C2 = (nu23C2 * (edgeCoord2426-globalPos2))/T2C2;


    /** 3) Calculate A, B, C, D and solve for T **/
    // compute transmissibility matrix T = CA^{-1}B+D
    DimMatrix C(0), A(0);
    Dune::FieldMatrix<Scalar,dim,dim+1> D(0), B(0);

    // evaluate matrix C, D, A, B
    C[0][0] = -omega121C2;
    C[0][1] = -omega122C2;
    C[0][2] = -omega123C2;
    C[1][0] = -omega221C2;
    C[1][1] = -omega222C2;
    C[1][2] = -omega223C2;
    C[2][0] = -omega321C2;
    C[2][1] = -omega322C2;
    C[2][2] = -omega323C2;

    D[0][1] = omega121C2 + omega122C2 + omega123C2;
    D[1][1] = omega221C2 + omega222C2 + omega223C2;
    D[2][1] = omega321C2 + omega322C2 + omega323C2;

    A[0][0] = omega121C2 - omega112C2 - omega111C2*r211C2 - omega113C2*r212C2;
    A[0][1] = omega122C2 - omega111C2*r221C2 - omega113C2*r222C2;
    A[0][2] = omega123C2 - omega111C2*r231C2 - omega113C2*r232C2;
    A[1][0] = omega221C2 - omega242C2*r211C2 - omega243C2*r213C2;
    A[1][1] = omega222C2 - omega241C2 - omega242C2*r221C2 - omega243C2*r223C2;
    A[1][2] = omega223C2 - omega242C2*r231C2 - omega243C2*r233C2;
    A[2][0] = omega321C2 - omega361C2*r213C2 - omega362C2*r212C2;
    A[2][1] = omega322C2 - omega361C2*r223C2 - omega362C2*r222C2;
    A[2][2] = omega323C2 - omega363C2 - omega361C2*r233C2 - omega362C2*r232C2;

    B[0][0] = -omega111C2 - omega112C2 - omega113C2;
    B[0][1] = omega121C2 + omega122C2 + omega123C2 + omega111C2*(1.0 - r211C2 - r221C2 -r231C2)
            + omega113C2*(1.0 - r212C2 - r222C2 - r232C2);
    B[1][1] = omega221C2 + omega222C2 + omega223C2 + omega242C2*(1.0 - r211C2 - r221C2 - r231C2)
            + omega243C2*(1.0 - r213C2 - r223C2 - r233C2);
    B[1][2] = -omega241C2 - omega242C2 - omega243C2;
    B[2][1] = omega321C2 + omega322C2 + omega323C2 + omega361C2*(1.0 - r213C2 - r223C2 - r233C2)
            + omega362C2*(1.0 - r212C2 -r222C2 -r232C2);
    B[2][3] = -omega361C2 - omega362C2 -omega363C2;

//    printmatrix(std::cout, A, "matrix A ", "row", 11, 3);
//    printmatrix(std::cout, B, "matrix B ", "row", 11, 3);
//    printmatrix(std::cout, C, "matrix C ", "row", 11, 3);
//    printmatrix(std::cout, D, "matrix D ", "row", 11, 3);

    // compute T
    A.invert();
    D += B.leftmultiply(C.rightmultiply(A));
    T = D[0];

    /** 4) Store and/or calculate additional half edges *******/
    if(maxInteractionVolumes ==1)
    {
        T *= intersection.geometry().volume()/subFaceArea12 ;

        // set your map entry
        problem().variables().storeMpfaData3D(intersection, T, globalPos4, eIdxGlobal4, globalPos6, eIdxGlobal6);
        return 1; // indicates that only 1 interaction region was regarded
    }
    else
    {
        int countInteractionRegions = 1;
        //initialize additional transmissitivity matrix
        TransmissivityMatrix additionalT(0.);

        auto outerCorner = intersection.inside().template subEntity<dim>(0); //initialize with rubbish
        // prepare additonal pointer to cells
        auto additional2 = intersection.inside(); //initialize with something wrong!
        auto additional3 = intersection.inside();
        int caseL = -2;

        /**** 2nd interaction region: get corner of interest ************/
        // search through corners of large cell with isIt
        int localIdxLarge = searchCommonVertex_(intersection, outerCorner);

        //in the parallel case, skip all border entities
        #if HAVE_MPI
        if (problem().gridView().comm().size() > 1)
            if(outerCorner.partitionType() != Dune::InteriorEntity)
                caseL = -1; // abort this specific interaction volume
        #endif

        // get Interaction Volume object
        int vIdxGlobal = problem().variables().index(outerCorner);
        InteractionVolume& interactionVolume
                        = interactionVolumesContainer_->interactionVolume(vIdxGlobal);

        // abort if we are on boundary
        if(interactionVolume.isBoundaryInteractionVolume())
            caseL = -1;

        // use Markus mpfa-l implementation
        int subVolumeFaceIdx = -1;
        bool properFluxDirection = true;
        if(caseL == -2)
        {
            //get Hanging Node type: no hanging node would mean -1
            int hangingNodeType = interactionVolume.getHangingNodeType();

            if(hangingNodeType == InteractionVolume::noHangingNode)
                subVolumeFaceIdx = interactionVolumesContainer_->getMpfaCase8cells(isIt, localIdxLarge, interactionVolume, properFluxDirection);
            else if(hangingNodeType == InteractionVolume::sixSmallCells)
                subVolumeFaceIdx = interactionVolumesContainer_->getMpfaCase6cells(isIt, interactionVolume, properFluxDirection);
            else
                Dune::dgrave << "HangingType " << hangingNodeType << " not implemented " << std::endl;

            caseL = this->transmissibilityAdapter_(isIt, interactionVolume, subVolumeFaceIdx,
                                    properFluxDirection, additional2, additional3, additionalT);
        }

        //store everything as second interaction region
        if(caseL != -1) //check if we regard 2 interaction regions
        {
            problem().variables().storeMpfaData3D(intersection, additionalT,
                    additional2.geometry().center(), problem().variables().index(additional2),
                    additional3.geometry().center(), problem().variables().index(additional3),
                    1); // offset for second interaction region
            countInteractionRegions++;
        }

        /***** 3rd and 4th interaction region *******/
        if(maxInteractionVolumes>2)
        {
            // loop through remaining 2 points
            std::vector<int> diagonal;
            for(int verticeSmall = 0; verticeSmall < intersection.outside().subEntities(dim); ++verticeSmall)
            {
                auto vSmall = intersection.outside().template subEntity<dim>(verticeSmall);

                //in the parallel case, skip all border entities
                #if HAVE_MPI
                if (problem().gridView().comm().size() > 1)
                    if(vSmall.partitionType() != Dune::InteriorEntity)
                        continue;
                #endif

                // get postion as seen from element
                GlobalPosition vertexOnElement
                    = refElement.position(verticeSmall, dim);

                for (int indexOnFace = 0; indexOnFace < 4; indexOnFace++)
                {
                    // get position as seen from interface
                    GlobalPosition vertexOnInterface
                        = intersection.geometryInOutside().corner(indexOnFace);

                    if(vSmall != outerCorner
                            && ((vertexOnInterface - vertexOnElement).two_norm()<1e-5))
                    {
                        int vIdxGlobal2 = problem().variables().index(vSmall);
                        // acess interactionVolume
                        InteractionVolume& interactionVolume2
                            = interactionVolumesContainer_->interactionVolume(vIdxGlobal2);
                        if(interactionVolume2.isBoundaryInteractionVolume())
                            continue;

                        int hangingNodeType = interactionVolume2.getHangingNodeType();
                        // reset flux direction indicator
                        properFluxDirection = true;

                        if(hangingNodeType != InteractionVolume::fourSmallCellsFace)
                        {
                            diagonal.push_back(problem().variables().index(vSmall));
                            // a) take interaction volume and determine fIdx
                            if(hangingNodeType == InteractionVolume::noHangingNode)
                            {
                                //TODO determine current localIdxLarge!!!!
                                Dune::dgrave << " args, noHangingNode on additional interaction region";
    //                            subVolumeFaceIdx = getMpfaCase8cells_(isIt, localIdxLarge, interactionVolume2, properFluxDirection);
                            }
                            else if(hangingNodeType == InteractionVolume::sixSmallCells)
                                subVolumeFaceIdx = interactionVolumesContainer_->getMpfaCase6cells(isIt,
                                                                        interactionVolume2, properFluxDirection);
                            else
                                subVolumeFaceIdx = interactionVolumesContainer_->getMpfaCase2or4cells(isIt,
                                                                        interactionVolume2, properFluxDirection);

                            // b) calculate T, eIdxGlobal3+4
                            caseL = this->transmissibilityAdapter_(isIt, interactionVolume2, subVolumeFaceIdx,
                                                        properFluxDirection, additional2, additional3, additionalT);

                            // c) store it
                            //store everything as third/4th interaction region
                            if(caseL != -1) //check if we regard this interaction region
                            {
                                problem().variables().storeMpfaData3D(intersection, additionalT,
                                        additional2.geometry().center(), problem().variables().index(additional2),
                                        additional3.geometry().center(), problem().variables().index(additional3),
                                        countInteractionRegions); // offset for this interaction region
                                countInteractionRegions++;
                            }
                        }
                        //else we would have type 1, which is handled explicitly in the beginning of this method
                    }
                }
            }
        }

        // set your map entry
        problem().variables().storeMpfaData3D(intersection, T, globalPos4, eIdxGlobal4, globalPos6, eIdxGlobal6);

        // determine weights
        Scalar weight = intersection.geometry().volume()/subFaceArea12; // =4 for 1 interaction region
        weight /= static_cast<Scalar>(countInteractionRegions);

        // perform weighting of transmissibilies
        problem().variables().performTransmissitivityWeighting(intersection, weight);

        return countInteractionRegions;
    }

    return 0;
}

/*!
 * \brief Adapter to use the general implementation of the mpfa-l for the compositional models
 *
 * Depending on the subVolumeFaceIdx, the appropriate method in
 * FvMpfaL2dTransmissibilityCalculator (potentially specifying certain cases)
 * gets called and the transmissibility and geometric information of the applied additional
 * cells of the interaction regions are passed back.
 *
* \param isIt Iterator to the current intersection
* \param interactionVolume The current interaction Volume object of interest
* \param subVolumeFaceIdx The local index of the intersection of interest in the interaction volume
* \param properFluxDirection True if the intersection normal coincides
*           with the local indexing in the interaction volume
* \param[out] additional2 The 3rd cell's element in the interaction volume
* \param[out] additional3 The 4th cell's element in the interaction volume
* \param[out] additionalT Transmissitivity matrix calculated
*/
template<class TypeTag>
int FV3dPressure2P2CAdaptive<TypeTag>::transmissibilityAdapter_(const IntersectionIterator& isIt,
                                InteractionVolume& interactionVolume,
                                const int& subVolumeFaceIdx,
                                bool properFluxDirection,
                                Element& additional2,
                                Element& additional3,
                                TransmissivityMatrix& additionalT)
{
    // abort if we are on boundary
    if(interactionVolume.isBoundaryInteractionVolume())
        return -1;

    //get Hanging Node type: no hanging node would mean -1
    int hangingNodeType = interactionVolume.getHangingNodeType();

    // prepare necessary variables etc
    Dune::FieldMatrix<Scalar, dim, 2 * dim - dim + 1> T(0);
    int caseL = -1;
    Dune::FieldVector<Scalar, dim> unity(1.);
    std::vector<Dune::FieldVector<Scalar, dim> > lambda
        = {unity, unity, unity, unity, unity, unity, unity, unity};
    Dune::FieldVector<bool, 4> useCases(false);

    /*********** calculate transmissibility for that case ******/
    switch(subVolumeFaceIdx)
    {
        case 0:{
            caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                    lambda, 0, 1, 2, 3, 4, 5);

            additionalT = T[0];

            // determine cells for T-matrix entry 2 and 3
            switch (caseL)
            {
            case 1:{
                additional2 = interactionVolume.getSubVolumeElement(2);
                additional3 = interactionVolume.getSubVolumeElement(4);
            break; }
            case 2:{
                additional2 = interactionVolume.getSubVolumeElement(3);
                additional3 = interactionVolume.getSubVolumeElement(5);
            break; }
            case 3:{
                additional2 = interactionVolume.getSubVolumeElement(3);
                additional3 = interactionVolume.getSubVolumeElement(4);
            break; }
            case 4:{
                additional2 = interactionVolume.getSubVolumeElement(2);
                additional3 = interactionVolume.getSubVolumeElement(5);
            break; }
            }
            break;
        }

        case 1:{
            if (hangingNodeType == InteractionVolume::twoSmallCells
                    || hangingNodeType == InteractionVolume::fourSmallCellsDiag )
            {
                useCases[0] = true;
                useCases[1] = false;
                useCases[2] = false;
                if(hangingNodeType != InteractionVolume::twoSmallCells)
                    useCases[3] = true; //differs from 2p Implementation
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                    lambda, 1, 3, 0, 2, 5, 7, useCases);
            }
            else
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                    lambda, 1, 3, 0, 2, 5, 7);

            additionalT = T[0];

            // determine cells for T-matrix entry 2 and 3
            switch (caseL)
            {
                case 1:
                {
                    additional2 = interactionVolume.getSubVolumeElement(0);
                    additional3 = interactionVolume.getSubVolumeElement(5);
                    break;}
                case 2:
                {
                    additional2 = interactionVolume.getSubVolumeElement(2);
                    additional3 = interactionVolume.getSubVolumeElement(7);
                    break;}
                case 3:
                {
                    additional2 = interactionVolume.getSubVolumeElement(2);
                    additional3 = interactionVolume.getSubVolumeElement(5);
                    break;}
                case 4:
                {
                    additional2 = interactionVolume.getSubVolumeElement(0);
                    additional3 = interactionVolume.getSubVolumeElement(7);
                    break;}
            }
            break;
        }

        case 2:
        {
            assert (hangingNodeType != InteractionVolume::twoSmallCells
                        && hangingNodeType != InteractionVolume::fourSmallCellsDiag);

            caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                    lambda, 3, 2, 1, 0, 7, 6);

            additionalT = T[0];

            // determine cells for T-matrix entry 2 and 3
            switch (caseL)
            {
                case 1:
                {
                    additional2 = interactionVolume.getSubVolumeElement(1);
                    additional3 = interactionVolume.getSubVolumeElement(7);
                    break;}
                case 2:
                {
                    additional2 = interactionVolume.getSubVolumeElement(0);
                    additional3 = interactionVolume.getSubVolumeElement(6);
                    break;}
                case 3:
                {
                    additional2 = interactionVolume.getSubVolumeElement(0);
                    additional3 = interactionVolume.getSubVolumeElement(7);
                    break;}
                case 4:
                {
                    additional2 = interactionVolume.getSubVolumeElement(1);
                    additional3 = interactionVolume.getSubVolumeElement(6);
                    break;}
            }
            break;
        }

        case 3:{
            if (hangingNodeType == InteractionVolume::twoSmallCells
                    || hangingNodeType == InteractionVolume::fourSmallCellsDiag)
            {
                useCases[0] = false;
                useCases[1] = true;
                if(hangingNodeType != InteractionVolume::twoSmallCells)
                    useCases[2] = true; //TODO: warning, hacked from markus!!
                useCases[3] = false;
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T,
                                    interactionVolume, lambda, 2, 0, 3, 1, 6, 4, useCases);
            }
            else
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                    lambda, 2, 0, 3, 1, 6, 4);

            additionalT = T[0];

            // determine cells for T-matrix entry 2 and 3
            switch (caseL)
            {
                case 1:
                {
                    additional2 = interactionVolume.getSubVolumeElement(3);
                    additional3 = interactionVolume.getSubVolumeElement(6);
                    break;}
                case 2:
                {
                    additional2 = interactionVolume.getSubVolumeElement(1);
                    additional3 = interactionVolume.getSubVolumeElement(4);
                    break;}
                case 3:
                {
                    additional2 = interactionVolume.getSubVolumeElement(1);
                    additional3 = interactionVolume.getSubVolumeElement(6);
                    break;}
                case 4:
                {
                    additional2 = interactionVolume.getSubVolumeElement(3);
                    additional3 = interactionVolume.getSubVolumeElement(4);
                    break;}
            }
            break;
        }

        case 4:{
            if (hangingNodeType == InteractionVolume::fourSmallCellsEdge) // should never happen, cause TPFA should be applied
            {
                assert(problem().variables().index(interactionVolume.getSubVolumeElement(4))
                        != problem().variables().index(interactionVolume.getSubVolumeElement(5))); // else there would not be a subVolFaceIdx 4
                Dune::dgrave << " SubVolumeFace4 in hanging node type 3 should be modelled by"
                                << " Tpfa!!" <<std::endl; // TODO: or use 1,5,3,4 and case 4
                return -1;
            }
            else
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 5, 4, 7, 6, 1, 0);

            additionalT = T[0];

            // determine cells for T-matrix entry 2 and 3
            switch (caseL)
            {
            case 1:{
                additional2 = interactionVolume.getSubVolumeElement(7);
                additional3 = interactionVolume.getSubVolumeElement(1);
            break; }
            case 2:{
                additional2 = interactionVolume.getSubVolumeElement(6);
                additional3 = interactionVolume.getSubVolumeElement(0);
            break; }
            case 3:{
                additional2 = interactionVolume.getSubVolumeElement(6);
                additional3 = interactionVolume.getSubVolumeElement(1);
            break; }
            case 4:{
                additional2 = interactionVolume.getSubVolumeElement(7);
                additional3 = interactionVolume.getSubVolumeElement(0);
            break; }
            }
            break;
        }

        case 5:{
            if (hangingNodeType == InteractionVolume::fourSmallCellsEdge)
            {
                assert (problem().variables().index(interactionVolume.getSubVolumeElement(4))
                        != problem().variables().index(interactionVolume.getSubVolumeElement(6))); // else there would not be a subVolFaceIdx 5
                Dune::dgrave << " SubVolumeFace5 in hanging node type 3 should be modelled by"
                                    << " Tpfa!!" <<std::endl; // TODO: or use 1,5,7,0 and case 3
                return -1;
            }
            else if (hangingNodeType == InteractionVolume::sixSmallCells)
            {
                useCases[0] = false;
                useCases[1] = true;
                useCases[2] = true;
                useCases[3] = false;

                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 7, 5, 6, 4, 3, 1, useCases);
            }
            else if (hangingNodeType == InteractionVolume::fourSmallCellsDiag)
            {
                useCases[0] = true;
                useCases[1] = false;
                useCases[2] = false;
                useCases[3] = true;
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda,  7, 5, 6, 4, 3, 1, useCases);
            }
            else
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 7, 5, 6, 4, 3, 1);

            additionalT = T[0];

            // determine cells for T-matrix entry 2 and 3
            switch (caseL)
            {
                case 1:
                {
                    additional2 = interactionVolume.getSubVolumeElement(6);
                    additional3 = interactionVolume.getSubVolumeElement(3);
                    break;}
                case 2:
                {
                    additional2 = interactionVolume.getSubVolumeElement(4);
                    additional3 = interactionVolume.getSubVolumeElement(1);
                    break;}
                case 3:
                {
                    additional2 = interactionVolume.getSubVolumeElement(4);
                    additional3 = interactionVolume.getSubVolumeElement(3);
                    break;}
                case 4:
                {
                    additional2 = interactionVolume.getSubVolumeElement(6);
                    additional3 = interactionVolume.getSubVolumeElement(1);
                    break;}
            }
            break;
        }

        case 6:{
            if (hangingNodeType == InteractionVolume::fourSmallCellsEdge)
            {
                assert (problem().variables().index(interactionVolume.getSubVolumeElement(4))
                        != problem().variables().index(interactionVolume.getSubVolumeElement(5))); // else there would not be a subVolFaceIdx 4
                Dune::dgrave << " SubVolumeFace6 in hanging node type 3 should be modelled by"
                        << " Tpfa!!" <<std::endl; // TODO: or use 2,6,0,7 and case 4
                return -1;
            }
            else
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 6, 7, 4, 5, 2, 3);

            additionalT = T[0];

            // determine cells for T-matrix entry 2 and 3
            switch (caseL)
            {
                case 1:
                {
                    additional2 = interactionVolume.getSubVolumeElement(4);
                    additional3 = interactionVolume.getSubVolumeElement(2);
                    break;}
                case 2:
                {
                    additional2 = interactionVolume.getSubVolumeElement(5);
                    additional3 = interactionVolume.getSubVolumeElement(3);
                    break;}
                case 3:
                {
                    additional2 = interactionVolume.getSubVolumeElement(5);
                    additional3 = interactionVolume.getSubVolumeElement(2);
                    break;}
                case 4:
                {
                    additional2 = interactionVolume.getSubVolumeElement(4);
                    additional3 = interactionVolume.getSubVolumeElement(3);
                    break;}
            }
            break;
        }

        case 7:
        {
            if (hangingNodeType == InteractionVolume::fourSmallCellsEdge)
            {
                assert (problem().variables().index(interactionVolume.getSubVolumeElement(4))
                        != problem().variables().index(interactionVolume.getSubVolumeElement(6))); // else there would not be a subVolFaceIdx 5
                Dune::dgrave << " SubVolumeFace5 in hanging node type 3 should be modelled by"
                                    << " Tpfa!!" <<std::endl; // TODO: or use 4,0,6,1 and case 4
                return -1;
            }
            else if (hangingNodeType == InteractionVolume::sixSmallCells)
            {
                useCases[0] = true;
                useCases[1] = false;
                useCases[2] = false;
                useCases[3] = true;

                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                    lambda, 4, 6, 5, 7, 0, 2, useCases);
            }
            else if (hangingNodeType == InteractionVolume::fourSmallCellsDiag)
            {
                useCases[0] = false;
                useCases[1] = true;
                useCases[2] = true;
                useCases[3] = false;
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                    lambda, 4, 6, 5, 7, 0, 2, useCases);
            }
            else
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                    lambda, 4, 6, 5, 7, 0, 2);

            additionalT = T[0];

            // determine cells for T-matrix entry 2 and 3
            switch (caseL)
            {
                case 1:
                {
                    additional2 = interactionVolume.getSubVolumeElement(5);
                    additional3 = interactionVolume.getSubVolumeElement(0);
                    break;}
                case 2:
                {
                    additional2 = interactionVolume.getSubVolumeElement(7);
                    additional3 = interactionVolume.getSubVolumeElement(2);
                    break;}
                case 3:
                {
                    additional2 = interactionVolume.getSubVolumeElement(7);
                    additional3 = interactionVolume.getSubVolumeElement(0);
                    break;}
                case 4:
                {
                    additional2 = interactionVolume.getSubVolumeElement(5);
                    additional3 = interactionVolume.getSubVolumeElement(2);
                    break;}
            }
            break;
        }

        case 8:
        {
            if(hangingNodeType == InteractionVolume::noHangingNode
               || hangingNodeType == InteractionVolume::sixSmallCells)
            {
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 4, 0, 6, 2, 5, 1);
            }
            else if (hangingNodeType == InteractionVolume::twoSmallCells
                            || hangingNodeType == InteractionVolume::fourSmallCellsFace)
            {
                caseL = mpfal3DTransmissibilityCalculator_.transmissibilityCaseTwo(T, interactionVolume,
                                        lambda, 4, 0, 2, 1);
            }
            else if (hangingNodeType == InteractionVolume::fourSmallCellsDiag
                            || (hangingNodeType == InteractionVolume::fourSmallCellsEdge
                                && (problem().variables().index(interactionVolume.getSubVolumeElement(4))
                                    != problem().variables().index(interactionVolume.getSubVolumeElement(6)))) )
            {
                useCases[0] = false;
                useCases[1] = true;
                useCases[2] = false;
                useCases[3] = true;

                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 4, 0, 6, 2, 5, 1, useCases);
            }
            else if (hangingNodeType == InteractionVolume::fourSmallCellsEdge)
            {
                useCases[0] = false;
                useCases[1] = true;
                useCases[2] = true;
                useCases[3] = false;

                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 4, 0, 6, 2, 5, 1, useCases);
            }

            additionalT = T[0];

            // determine cells for T-matrix entry 2 and 3
            switch (caseL)
            {
                case 1:
                {
                    additional2 = interactionVolume.getSubVolumeElement(6);
                    additional3 = interactionVolume.getSubVolumeElement(5);
                    break;}
                case 2:
                {
                    additional2 = interactionVolume.getSubVolumeElement(2);
                    additional3 = interactionVolume.getSubVolumeElement(1);
                    break;}
                case 3:
                {
                    additional2 = interactionVolume.getSubVolumeElement(2);
                    additional3 = interactionVolume.getSubVolumeElement(5);
                    break;}
                case 4:
                {
                    additional2 = interactionVolume.getSubVolumeElement(6);
                    additional3 = interactionVolume.getSubVolumeElement(1);
                    break;}
            }
            break;
        }

        case 9:
        {
            if(hangingNodeType == InteractionVolume::noHangingNode
               || hangingNodeType == InteractionVolume::sixSmallCells)
            {
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 1, 5, 3, 7, 0, 4);
            }
            else if (hangingNodeType == InteractionVolume::twoSmallCells || hangingNodeType == InteractionVolume::fourSmallCellsFace)
            {
                caseL = mpfal3DTransmissibilityCalculator_.transmissibilityCaseOne(T, interactionVolume,
                                        lambda, 1, 5, 3, 0);
            }
            else if (hangingNodeType == InteractionVolume::fourSmallCellsDiag
                    || (hangingNodeType == InteractionVolume::fourSmallCellsEdge
                        &&(problem().variables().index(interactionVolume.getSubVolumeElement(4))
                            != problem().variables().index(interactionVolume.getSubVolumeElement(6))) ))
            {
                useCases[0] = true;
                useCases[1] = false;
                useCases[2] = true;
                useCases[3] = false;

                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 1, 5, 3, 7, 0, 4, useCases);
            }
            else if (hangingNodeType == InteractionVolume::fourSmallCellsEdge
                    &&(problem().variables().index(interactionVolume.getSubVolumeElement(4))
                            != problem().variables().index(interactionVolume.getSubVolumeElement(5))) )
            {
                useCases[0] = true;
                useCases[1] = false;
                useCases[2] = false;
                useCases[3] = true;

                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 1, 5, 3, 7, 0, 4, useCases);
            }
            else
                Dune::dgrave << " Missing case for subVolFaceIdx 9!!" <<std::endl;

            additionalT = T[0];

            // determine cells for T-matrix entry 2 and 3
            switch (caseL)
            {
                case 1:
                {
                    additional2 = interactionVolume.getSubVolumeElement(3);
                    additional3 = interactionVolume.getSubVolumeElement(0);
                    break;}
                case 2:
                {
                    additional2 = interactionVolume.getSubVolumeElement(7);
                    additional3 = interactionVolume.getSubVolumeElement(4);
                    break;}
                case 3:
                {
                    additional2 = interactionVolume.getSubVolumeElement(7);
                    additional3 = interactionVolume.getSubVolumeElement(0);
                    break;}
                case 4:
                {
                    additional2 = interactionVolume.getSubVolumeElement(3);
                    additional3 = interactionVolume.getSubVolumeElement(4);
                    break;}
            }
            break;
        }

        case 10:
        {
            if (hangingNodeType == InteractionVolume::fourSmallCellsFace)
                caseL = mpfal3DTransmissibilityCalculator_.transmissibilityCaseTwo(T, interactionVolume,
                                        lambda, 7, 3, 1, 2);
            else if ((hangingNodeType == InteractionVolume::fourSmallCellsEdge
                      &&(problem().variables().index(interactionVolume.getSubVolumeElement(4))
                            != problem().variables().index(interactionVolume.getSubVolumeElement(6))) )
                      || hangingNodeType == InteractionVolume::sixSmallCells)
            {
                useCases[0] = false;
                useCases[1] = true;
                useCases[2] = false;
                useCases[3] = true;

                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 7, 3, 5, 1, 6, 2, useCases);
            }
            else if (hangingNodeType == InteractionVolume::fourSmallCellsEdge)
            {
                useCases[0] = false;
                useCases[1] = true;
                useCases[2] = true;
                useCases[3] = false;

                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 7, 3, 5, 1, 6, 2, useCases);
            }
            else if (hangingNodeType == InteractionVolume::fourSmallCellsDiag)
            {
                useCases[0] = true;
                useCases[1] = false;
                useCases[2] = true;
                useCases[3] = false;

                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 7, 3, 5, 1, 6, 2, useCases);
            }
            else
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 7, 3, 5, 1, 6, 2);

            additionalT = T[0];

            // determine cells for T-matrix entry 2 and 3
            switch (caseL)
            {
                case 1:
                {
                    additional2 = interactionVolume.getSubVolumeElement(5);
                    additional3 = interactionVolume.getSubVolumeElement(6);
                    break;}
                case 2:
                {
                    additional2 = interactionVolume.getSubVolumeElement(1);
                    additional3 = interactionVolume.getSubVolumeElement(2);
                    break;}
                case 3:
                {
                    additional2 = interactionVolume.getSubVolumeElement(1);
                    additional3 = interactionVolume.getSubVolumeElement(6);
                    break;}
                case 4:
                {
                    additional2 = interactionVolume.getSubVolumeElement(5);
                    additional3 = interactionVolume.getSubVolumeElement(2);
                    break;}
            }
            break;
        }

        case 11:
        {
            if(!interactionVolume.isHangingNodeVolume())
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 2, 6, 0, 4, 3, 7);
            else if (hangingNodeType == InteractionVolume::fourSmallCellsFace)
                caseL = mpfal3DTransmissibilityCalculator_.transmissibilityCaseOne(T, interactionVolume,
                                        lambda, 2, 6, 0, 3);
            else if ((hangingNodeType == InteractionVolume::fourSmallCellsEdge
                            && (problem().variables().index(interactionVolume.getSubVolumeElement(4))
                                != problem().variables().index(interactionVolume.getSubVolumeElement(6))))
                    || hangingNodeType == InteractionVolume::sixSmallCells)
            {
                useCases[0] = true;
                useCases[1] = false;
                useCases[2] = true;
                useCases[3] = false;

                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 2, 6, 0, 4, 3, 7, useCases);
            }
            else if (hangingNodeType == InteractionVolume::fourSmallCellsEdge
                        && (problem().variables().index(interactionVolume.getSubVolumeElement(4))
                            != problem().variables().index(interactionVolume.getSubVolumeElement(5))) )
            {
                useCases[0] = true;
                useCases[1] = false;
                useCases[2] = false;
                useCases[3] = true;

                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 2, 6, 0, 4, 3, 7, useCases);
            }
            else if (hangingNodeType == InteractionVolume::fourSmallCellsDiag)
            {
                useCases[0] = false;
                useCases[1] = true;
                useCases[2] = false;
                useCases[3] = true;

                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 2, 6, 0, 4, 3, 7, useCases);
            }

            additionalT = T[0];

            // determine cells for T-matrix entry 2 and 3
            switch (caseL)
            {
                case 1:
                {
                    additional2 = interactionVolume.getSubVolumeElement(0);
                    additional3 = interactionVolume.getSubVolumeElement(3);
                    break;}
                case 2:
                {
                    additional2 = interactionVolume.getSubVolumeElement(4);
                    additional3 = interactionVolume.getSubVolumeElement(7);
                    break;}
                case 3:
                {
                    additional2 = interactionVolume.getSubVolumeElement(4);
                    additional3 = interactionVolume.getSubVolumeElement(3);
                    break;}
                case 4:
                {
                    additional2 = interactionVolume.getSubVolumeElement(0);
                    additional3 = interactionVolume.getSubVolumeElement(7);
                    break;}
            }
            break;}
    }

    //if necessary, correct for reversely calculated flux direction
    if(!properFluxDirection)
    {
        // mpfal3DTransmissibilityCalculator_.setTransChoiceThreshold(transChoiceThreshold_);
        // a) reverse direction
        additionalT *= -1;
        // b) swap cell I and J
        using std::swap;
        swap(additionalT[0], additionalT[1]);
    }

    return caseL;
}


}//end namespace Dumux
#endif
