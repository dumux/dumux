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
#ifndef DUMUX_FV3DPRESSURE2P2C_ADAPTIVE_HH
#define DUMUX_FV3DPRESSURE2P2C_ADAPTIVE_HH

// dune environent:
#include <dune/common/version.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

// dumux environment
#include <dumux/decoupled/2p2c/fvpressure2p2cmultiphysics.hh>
#include <dumux/common/math.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/decoupled/2p2c/2p2cadaptiveproperties.hh>

// include pressure model from Markus
#include <dumux/decoupled/common/fv/mpfa/fvmpfaproperties.hh>
#include <dumux/decoupled/2p2c/fvmpfal3d2p2cinteractionvolumecontaineradaptive.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/lmethod/fvmpfal3dtransmissibilitycalculator.hh>

/**
 * @file
 * @brief  Finite Volume Diffusion Model
 * @author Benjamin Faigle, Bernd Flemisch, Jochen Fritz, Markus Wolff
 */

namespace Dumux
{
namespace Properties
{
SET_TYPE_PROP(DecoupledTwoPTwoCAdaptive, MPFAInteractionVolume, Dumux::FvMpfaL3dInteractionVolumeAdaptive<TypeTag>);
SET_TYPE_PROP(DecoupledTwoPTwoCAdaptive, MPFAInteractionVolumeContainer, Dumux::FvMpfaL3d2P2CInteractionVolumeContainerAdaptive<TypeTag>);
}

//! The finite volume model for the solution of the compositional pressure equation
/*! \ingroup Adaptive2p2c
 *  Provides a Finite Volume implementation for the pressure equation of a compressible
 *  system with two components. An IMPES-like method is used for the sequential
 *  solution of the problem.  Diffusion is neglected, capillarity can be regarded.
 *  Isothermal conditions and local thermodynamic
 *  equilibrium are assumed.  Gravity is included.
 *  \f[
         c_{total}\frac{\partial p}{\partial t} + \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}} \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{\alpha} \bf{v}_{\alpha}\right)
          = \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}} q^{\kappa},
 *  \f]
 *  where \f$\bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 *  \f$ c_{total} \f$ represents the total compressibility, for constant porosity this yields \f$ - \frac{\partial V_{total}}{\partial p_{\alpha}} \f$,
 *  \f$p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability, \f$ \lambda_{\alpha} \f$ the phase mobility,
 *  \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and \f$ C^{\kappa} \f$ the total Component concentration.
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
    typedef typename GET_PROP_TYPE(TypeTag, PressureModel) Implementation;
    typedef FVPressure2P2CMultiPhysics<TypeTag> ParentType;
    typedef FVPressure<TypeTag> BaseType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;
    typedef typename SpatialParams::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld,
        NumPhases = GET_PROP_VALUE(TypeTag, NumPhases), NumComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureN,
        pglobal = Indices::pressureGlobal,
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

    // typedefs to abbreviate several dune classes...
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::ReferenceElements<Scalar, dim> ReferenceElementContainer;
    typedef Dune::ReferenceElement<Scalar, dim> ReferenceElement;

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::template Codim<dim>::EntityPointer VertexPointer;
    typedef typename GridView::Intersection Intersection;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    // convenience shortcuts for Vectors/Matrices
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar,dim+1> TransmissivityMatrix;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    typedef Dune::FieldVector<Scalar, GET_PROP_VALUE(TypeTag, NumPhases)> PhaseVector;
    typedef Dune::FieldVector<Scalar, GET_PROP_VALUE(TypeTag, NumComponents)> ComponentVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    // the typenames used for the stiffness matrix and solution vector
    typedef typename GET_PROP_TYPE(TypeTag, PressureCoefficientMatrix) Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, PressureRHSVector) RHSVector;

    // Dumux MPFA typedefs
    typedef typename GET_PROP_TYPE(TypeTag, MPFAInteractionVolumeContainer) InteractionVolumeContainer;
    typedef typename InteractionVolumeContainer::InteractionVolume InteractionVolume;

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
        // update of interaction Volumes if grid was changed
        if (true or problem().gridAdapt().wasAdapted())
        {
            if(enableMPFA && maxInteractionVolumes>1)
            {
                if(!interactionVolumesContainer_)
                    interactionVolumesContainer_ =
                            new InteractionVolumeContainer(problem());

                interactionVolumesContainer_->update();
            }
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
    };

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

        ElementIterator eItEnd = problem().gridView().template end<0>();
        for (ElementIterator eIt = problem().gridView().template begin<0>(); eIt != eItEnd; ++eIt)
        {
            // get the global index of the cell
            int globalIdxI = problem().variables().index(*eIt);

            // assemble interior element contributions
            if (eIt->partitionType() == Dune::InteriorEntity)
            {
                this->pressure()[globalIdxI]
                      = problem().variables().cellData(globalIdxI).pressure(this->pressureType);
            }
        }
#if HAVE_MPI
    // communicate updated values
    typedef typename GET_PROP(TypeTag, SolutionTypes) SolutionTypes;
    typedef typename SolutionTypes::ElementMapper ElementMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PressureSolutionVector) PressureSolution;
    typedef VectorExchange<ElementMapper, PressureSolution> DataHandle;

        DataHandle dataHandle(problem().variables().elementMapper(), this->pressure());
        problem().gridView().template communicate<DataHandle>(dataHandle,
                                                            Dune::InteriorBorder_All_Interface,
                                                            Dune::ForwardCommunication);
#endif
    }

    //! Constructs a FVPressure2P2C object
    /**
     * \param problem a problem class object
     */
    FV3dPressure2P2CAdaptive(Problem& problem) : FVPressure2P2CMultiPhysics<TypeTag>(problem),
            enableMPFA(false),
            interactionVolumesContainer_(0), mpfal3DTransmissibilityCalculator_(problem)
    {
        enableVolumeIntegral_ = this->enableVolumeIntegral;
        enableMPFA= GET_PARAM_FROM_GROUP(TypeTag,bool, GridAdapt, EnableMultiPointFluxApproximation);

        maxInteractionVolumes = GET_PARAM_FROM_GROUP(TypeTag,int, GridAdapt, MaxInteractionVolumes);
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    {   return *static_cast<Implementation *>(this);}

    //! \copydoc Dumux::IMPETProblem::asImp_()
    const Implementation &asImp_() const
    {   return *static_cast<const Implementation *>(this);}

    int searchCommonVertex_(const Intersection& is, VertexPointer& vertexPointer)
    {
        /******* get corner of interest ************/
        // search through corners of large cell with isIt
        int localIdxLarge = 0;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
        for(localIdxLarge = 0; localIdxLarge < is.inside()->subEntities(dim); ++localIdxLarge)
#else
        for(localIdxLarge = 0; localIdxLarge<is.inside()->template count<dim>(); ++localIdxLarge)
#endif
        {
            const VertexPointer vPtrLarge = is.inside()->template subEntity<dim>(localIdxLarge);

            // search through corners of small cell with isIt
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            for(int verticeSmall = 0; verticeSmall<is.outside()->subEntities(dim); ++verticeSmall)
#else
            for(int verticeSmall = 0; verticeSmall<is.outside()->template count<dim>(); ++verticeSmall)
#endif
            {
                const VertexPointer vPtrSmall = is.outside()->template subEntity<dim>(verticeSmall);

                if(problem().variables().index(*vPtrSmall) == problem().variables().index(*vPtrLarge) )
                {
                    vertexPointer = vPtrSmall;
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
                                    ElementPointer& additional2,
                                    ElementPointer& additional3,
                                    TransmissivityMatrix& additionalT);

    std::map<int, std::vector<int> > irregularCellMap_; //!< Container to store all cell's Indice with a hanging node
    bool enableVolumeIntegral_; //!< Calculates the volume integral (on by default)
    bool enableMPFA; //!< Enables mpfa on hanging nodes (on by default)
    int maxInteractionVolumes; //!< Maximum number of interaction volumes considered (4 by default)

    //! A pointer to the adaptive interaction volumes container
    InteractionVolumeContainer* interactionVolumesContainer_;
    //! The common implementation to calculate the Transmissibility with the mpfa-L-method
    Dumux::FvMpfaL3dTransmissibilityCalculator<TypeTag> mpfal3DTransmissibilityCalculator_;
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
    ElementIterator eItEnd = problem().gridView().template end<0> ();
    for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem().variables().index(*eIt);
        CellData& cellDataI = problem().variables().cellData(globalIdxI);

        // initialize row size
        int rowSize = 1;

        // set perimeter to zero
        cellDataI.perimeter() = 0;

        // prepare storage for all found 3rd cells
        std::vector<int> foundAdditionals;

        int numberOfIntersections = 0;
        // run through all intersections with neighbors
        IntersectionIterator isItEnd = problem().gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem().gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            cellDataI.perimeter()
                    += isIt->geometry().volume();
            numberOfIntersections++;
            if (isIt->neighbor())
            {
                rowSize++;

                // special treatment for hanging nodes in the mpfa case
                if (enableMPFA && (eIt->level() < isIt->outside()->level()))
                {
                    // each cell might have 4 hanging nodes with 4 irregular neighbors each
                    // get global ID of Interface from larger cell
                    int intersectionID = problem().grid().localIdSet().subId(*eIt,
                            isIt->indexInInside(), 1);
                    //index outside
                    int globalIdxJ = problem().variables().index(*isIt->outside());

                    // add Entry of current neighbor cell to the IS seen from large cell
                    irregularCellMap_[intersectionID].push_back(globalIdxJ);
                }
            }
        }
        cellDataI.fluxData().resize(numberOfIntersections);
        this->A_.incrementrowsize(globalIdxI, rowSize);
    } //end first loop (that already reserved enough space for the MPFA connections on hanging nodes)

    // second loop to determine which matrix entries will be occupied
    if(enableMPFA && maxInteractionVolumes>1)
    {
        //prepare map for additional cell-connection through mpfa
        std::multimap<int, int> addionalRelations;
        typedef std::pair<int,int> IntPair;
        std::pair<std::multimap<int,int>::iterator,std::multimap<int,int>::iterator> range;
        std::multimap<int,int>::iterator rangeIt;

        // second loop for further sub-faces
        ElementIterator eItEnd = problem().gridView().template end<0> ();
        for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
        {
            // cell index
            int globalIdxI = problem().variables().index(*eIt);
            // run through all intersections with neighbors
            IntersectionIterator isItEnd = problem().gridView().iend(*eIt);
            for (IntersectionIterator isIt = problem().gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
            {
                if (isIt->neighbor())
                {
                    //index outside
                    int globalIdxJ = problem().variables().index(*isIt->outside());

                    // if mpfa is used, more entries might be needed if all interactionRegions are regarded
                    if (isIt->outside()->level() > eIt->level()) //look from larger cell
                    {
                        VertexPointer outerCornerPtr(isIt->inside()->template subEntity<dim>(0)); //initialize with rubbish
                        // prepare additional pointer to cells
                        ElementPointer additional2(isIt->inside()); //initialize with something wrong!
                        ElementPointer additional3(isIt->inside());

                        // Prepare MPFA
                        /** get geometric Info, transmissibility matrix */
                        GlobalPosition globalPos3(0.);
                        int globalIdx3=-1;
                        GlobalPosition globalPos4(0.);
                        int globalIdx4=-1;
                        TransmissivityMatrix T(0.);

                        int interactionRegions
                            = problem().variables().getMpfaData3D(*isIt, T,
                                                    globalPos3, globalIdx3, globalPos4, globalIdx4  );
                        if (interactionRegions == 0)
                            interactionRegions = problem().pressureModel().computeTransmissibilities(isIt,T,
                                    globalPos3, globalIdx3,  globalPos4, globalIdx4  );
                        if(!interactionRegions)
                            Dune::dgrave << "something went wrong getting mpfa data on cell " << globalIdxI << std::endl;
                        if (interactionRegions == 1) // no second subface
                            continue;

                        //loop over all found interaction regions
                        for (int cocumber=1; cocumber<interactionRegions; cocumber++ )
                        {
                            problem().variables().getMpfaData3D(*isIt, T,
                                   globalPos3, globalIdx3, globalPos4, globalIdx4, cocumber);
                            // indices
                            int additionalIdx2 = globalIdx3;
                            int additionalIdx3 = globalIdx4;

                            bool addIndex = true;

                            // check if additionals are "normal" neighbors of I
                            bool additional2isNeighbor(false), additional3isNeighbor(false);
                            // run through all intersections with neighbors if eIt
                            for (IntersectionIterator isItcheck = problem().gridView().ibegin(*eIt);
                                    isItcheck != problem().gridView().iend(*eIt); ++isItcheck)
                            {
                                if (isItcheck->neighbor())
                                {
                                    if(additionalIdx2==problem().variables().index(*isItcheck->outside()))
                                        additional2isNeighbor = true;
                                    if(additionalIdx3 == problem().variables().index(*isItcheck->outside()))
                                        additional3isNeighbor = true;
                                }

                            }
                            // if not "normal" neighbor, increase row size
                            if(!additional2isNeighbor)
                            {
                                // check if relation not already added
                                IntPair intPair(globalIdxI,additionalIdx2);
                                if(globalIdxI > additionalIdx2)
                                    std::swap(intPair.first, intPair.second);
                                range = addionalRelations.equal_range(intPair.first);
                                for (rangeIt=range.first; range.first!=range.second
                                                            and rangeIt!=range.second; ++rangeIt)
                                    if((*rangeIt).second == intPair.second)
                                        addIndex = false;
                                if(addIndex)
                                {
                                    this->A_.incrementrowsize(globalIdxI);
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
                                IntPair intPair(globalIdxI,additionalIdx3);
                                if(globalIdxI > additionalIdx3)
                                    std::swap(intPair.first, intPair.second);
                                range = addionalRelations.equal_range(intPair.first);
                                for (rangeIt=range.first; range.first!=range.second
                                                            and rangeIt!=range.second; ++rangeIt)
                                    if((*rangeIt).second == intPair.second)
                                        addIndex = false;
                                if(addIndex)
                                {
                                    this->A_.incrementrowsize(globalIdxI);
                                    // add space for additional itsself
                                    this->A_.incrementrowsize(additionalIdx3);
                                    // mark relation as added
                                    addionalRelations.insert(intPair);
                                }
                            }

                            //reset bools for J
                            additional2isNeighbor = additional3isNeighbor = false;
                            // run through all intersections with neighbors of J
                            for (IntersectionIterator isItcheck = problem().gridView().ibegin(*isIt->outside());
                                    isItcheck != problem().gridView().iend(*isIt->outside()); ++isItcheck)
                            {
                                if (isItcheck->neighbor())
                                {
                                    if(additionalIdx2 == problem().variables().index(*isItcheck->outside()))
                                        additional2isNeighbor = true;
                                    if(additionalIdx3 == problem().variables().index(*isItcheck->outside()))
                                        additional3isNeighbor = true;
                                }
                            }

                            // if not "normal" neighbor, increase row size
                            if(!additional2isNeighbor)
                            {
                                addIndex = true;
                                // check if relation not already added
                                IntPair intPair(globalIdxJ,additionalIdx2);
                                if(globalIdxJ > additionalIdx2)
                                    std::swap(intPair.first, intPair.second);
                                range = addionalRelations.equal_range(intPair.first);
                                for (rangeIt=range.first; range.first!=range.second
                                                            and rangeIt!=range.second; ++rangeIt)
                                    if((*rangeIt).second == intPair.second)
                                        addIndex = false;
                                if(addIndex)
                                {
                                    this->A_.incrementrowsize(globalIdxJ);
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
                                IntPair intPair(globalIdxJ,additionalIdx3);
                                if(globalIdxJ > additionalIdx3)
                                    std::swap(intPair.first, intPair.second);
                                range = addionalRelations.equal_range(intPair.first);
                                for (rangeIt=range.first; range.first!=range.second
                                                            and rangeIt!=range.second; ++rangeIt)
                                    if((*rangeIt).second == intPair.second)
                                        addIndex = false;
                                if(addIndex)
                                {
                                    this->A_.incrementrowsize(globalIdxJ);
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
    return;
}

//!Determine position of matrix entries
/* Method adds TPFA and MPFA matrix entries
 */
template<class TypeTag>
void FV3dPressure2P2CAdaptive<TypeTag>::initializeMatrixIndices()
{
    // determine position of matrix entries
    ElementIterator eItEnd = problem().gridView().template end<0> ();
    for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = problem().variables().index(*eIt);

        // add diagonal index
        this->A_.addindex(globalIdxI, globalIdxI);

        // run through all intersections with neighbors
        IntersectionIterator isItEnd = problem().gridView().template iend(*eIt);
        for (IntersectionIterator isIt = problem().gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer outside = isIt->outside();
                int globalIdxJ = problem().variables().index(*outside);

                // add off diagonal index
                this->A_.addindex(globalIdxI, globalIdxJ);

                // special treatment for hanging nodes in the mpfa case
                if (enableMPFA && (eIt->level() < isIt->outside()->level()))
                {
                    // prepare stuff to enter transmissibility calculation
                    GlobalPosition globalPos3(0.);
                    int globalIdx3=-1;
                    GlobalPosition globalPos4(0.);
                    int globalIdx4=-1;
                    TransmissivityMatrix T(0.);
                    TransmissivityMatrix additionalT(0.);

                    int interactionRegions
                        = problem().variables().getMpfaData3D(*isIt, T, globalPos3, globalIdx3, globalPos4, globalIdx4  );
                    if (interactionRegions == 0)
                        interactionRegions = problem().pressureModel().computeTransmissibilities(isIt,T,
                                globalPos3, globalIdx3,  globalPos4, globalIdx4  );

                    for (int cocumber=1; cocumber<interactionRegions; cocumber++ )
                    {
                        problem().variables().getMpfaData3D(*isIt, T,
                               globalPos3, globalIdx3, globalPos4, globalIdx4, cocumber);

                        // add off diagonal index in both directions!!
                        this->A_.addindex(globalIdxI, globalIdx3);
                        this->A_.addindex(globalIdx3, globalIdxI);
                        this->A_.addindex(globalIdxI, globalIdx4);
                        this->A_.addindex(globalIdx4, globalIdxI);
                        this->A_.addindex(globalIdxJ, globalIdx3);
                        this->A_.addindex(globalIdx3, globalIdxJ);
                        this->A_.addindex(globalIdxJ, globalIdx4);
                        this->A_.addindex(globalIdx4, globalIdxJ);
                    }
                }
            }
    }
    return;
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

    ElementIterator eItEnd = problem().gridView().template end<0>();
    for (ElementIterator eIt = problem().gridView().template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // get the global index of the cell
        int globalIdxI = problem().variables().index(*eIt);

        // assemble interior element contributions
        if (eIt->partitionType() == Dune::InteriorEntity)
        {
                // get the cell data
            CellData& cellDataI = problem().variables().cellData(globalIdxI);

            Dune::FieldVector<Scalar, 2> entries(0.);

            /*****  source term ***********/
#ifndef noMultiphysics
            if(cellDataI.subdomain() != 2)
                problem().pressureModel().get1pSource(entries,*eIt, cellDataI);
            else
#endif
                problem().pressureModel().getSource(entries,*eIt, cellDataI, first);

            this->f_[globalIdxI] += entries[rhs];

            /*****  flux term ***********/
            // iterate over all faces of the cell
            IntersectionIterator isItEnd = problem().gridView().template iend(*eIt);
            for (IntersectionIterator isIt = problem().gridView().template ibegin(*eIt); isIt != isItEnd; ++isIt)
            {
                /************* handle interior face *****************/
                if (isIt->neighbor())
                {

                    ElementPointer elementNeighbor = isIt->outside();
                    int globalIdxJ = problem().variables().index(*elementNeighbor);
                    //check for hanging nodes
                    //take a hanging node never from the element with smaller level!
                    bool haveSameLevel = (eIt->level() == elementNeighbor->level());
                    // calculate only from one side, but add matrix entries for both sides
                    // the last condition is needed to properly assemble in the presence 
                    // of ghost elements
                    if (GET_PROP_VALUE(TypeTag, VisitFacesOnlyOnce) 
                        && (globalIdxI > globalIdxJ) && haveSameLevel
                        && elementNeighbor->partitionType() == Dune::InteriorEntity)
                        continue;

                    entries = 0;
                    //check for hanging nodes
                    if(!haveSameLevel && enableMPFA)
                    {
                        if (cellDataI.subdomain() != 2
                                or problem().variables().cellData(globalIdxJ).subdomain() != 2) // cell in the 1p domain
                            asImp_().get1pMpfaFlux(isIt, cellDataI);
                        else
                            asImp_().getMpfaFlux(isIt, cellDataI);
                    }
                    else
                    {
                        CellData cellDataJ = problem().variables().cellData(globalIdxJ);
                        if (cellDataI.subdomain() != 2
                                or problem().variables().cellData(globalIdxJ).subdomain() != 2) // cell in the 1p domain
                            asImp_().get1pFlux(entries, *isIt, cellDataI);
                        else
                            asImp_().getFlux(entries, *isIt, cellDataI, first);


                        //set right hand side
                        this->f_[globalIdxI] -= entries[rhs];

                        // set diagonal entry
                        this->A_[globalIdxI][globalIdxI] += entries[matrix];

                            // set off-diagonal entry
                        this->A_[globalIdxI][globalIdxJ] -= entries[matrix];

                        // The second condition is needed to not spoil the ghost element entries
                        if (GET_PROP_VALUE(TypeTag, VisitFacesOnlyOnce) 
                            && elementNeighbor->partitionType() == Dune::InteriorEntity)
                        {
                            this->f_[globalIdxJ] += entries[rhs];
                            this->A_[globalIdxJ][globalIdxJ] += entries[matrix];
                            this->A_[globalIdxJ][globalIdxI] -= entries[matrix];
                        }
                    }
                } // end neighbor

                /************* boundary face ************************/
                else
                {
                    entries = 0;
                    if (cellDataI.subdomain() != 2) //the current cell in the 1p domain
                        asImp_().get1pFluxOnBoundary(entries, *isIt, cellDataI);
                    else
                        asImp_().getFluxOnBoundary(entries, *isIt, cellDataI, first);

                    //set right hand side
                    this->f_[globalIdxI] += entries[rhs];
                    // set diagonal entry
                    this->A_[globalIdxI][globalIdxI] += entries[matrix];
                }
            } //end interfaces loop

            /*****  storage term ***********/
            if (cellDataI.subdomain() != 2) //the current cell in the 1p domain
                asImp_().get1pStorage(entries, *eIt, cellDataI);
            else
                asImp_().getStorage(entries, *eIt, cellDataI, first);

            this->f_[globalIdxI] += entries[rhs];
            // set diagonal entry
            this->A_[globalIdxI][globalIdxI] += entries[matrix];
        }
        // assemble overlap and ghost element contributions
        else
        {
            this->A_[globalIdxI] = 0.0;
            this->A_[globalIdxI][globalIdxI] = 1.0;
            this->f_[globalIdxI] = this->pressure()[globalIdxI];
        }
    } // end grid traversal

    return;
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
    // acess Cell I
    ElementPointer elementPointerI = isIt->inside();
    int globalIdxI = problem().variables().index(*elementPointerI);

    // get global coordinate of cell center
    const GlobalPosition& globalPos = elementPointerI->geometry().center();

    // cell volume & perimeter, assume linear map here
    Scalar volume = elementPointerI->geometry().volume();
    Scalar perimeter = cellDataI.perimeter();

    // get absolute permeability
    DimMatrix permeabilityI(problem().spatialParams().intrinsicPermeability(*elementPointerI));

    // access neighbor
    ElementPointer neighborPointer = isIt->outside();
    int globalIdxJ = problem().variables().index(*neighborPointer);
    CellData& cellDataJ = problem().variables().cellData(globalIdxJ);

    // gemotry info of neighbor
    const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

    // distance vector between barycenters
    GlobalPosition distVec = globalPosNeighbor - globalPos;

    // compute distance between cell centers
    Scalar dist = distVec.two_norm();

    GlobalPosition unitDistVec(distVec);
    unitDistVec /= dist;

    DimMatrix permeabilityJ
        = problem().spatialParams().intrinsicPermeability(*neighborPointer);

    // compute vectorized permeabilities
    DimMatrix meanPermeability(0);
    Dumux::harmonicMeanMatrix(meanPermeability, permeabilityI, permeabilityJ);

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
        asImp_().volumeDerivatives(globalPosNeighbor, *neighborPointer);

    ComponentVector dv_dC(0.), graddv_dC(0.);
    for (int compIdx = 0; compIdx < NumComponents; ++compIdx)
    {
        dv_dC[compIdx]= (cellDataJ.dv(compIdx) // dV/dm1= dv/dC^1
                + cellDataI.dv(compIdx)) * 0.5;
        graddv_dC[compIdx] = (cellDataJ.dv(compIdx)
                                - cellDataI.dv(compIdx)) / dist;
    }
//                    potentialW = problem().variables().potentialWetting(globalIdxI, isIndex);
//                    potentialNW = problem().variables().potentialNonwetting(globalIdxI, isIndex);
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
        int globalIdx3=-1;
        GlobalPosition globalPos4(0.);
        int globalIdx4=-1;
        TransmissivityMatrix T(0.);

            // prepare second interaction region
            GlobalPosition globalPosAdditional3(0.);
            int globalIdxAdditional3=-1;
            GlobalPosition globalPosAdditional4(0.);
            int globalIdxAdditional4=-1;

            TransmissivityMatrix additionalT(0.);

        int interactionRegions
            = problem().variables().getMpfaData3D(*isIt, T,
                                    globalPos3, globalIdx3, globalPos4, globalIdx4  );
        if (interactionRegions == 0)
            interactionRegions = problem().pressureModel().computeTransmissibilities(isIt,T,
                    globalPos3, globalIdx3,  globalPos4, globalIdx4  );
        if(!interactionRegions)
            Dune::dgrave << "something went wrong getting mpfa data on cell " << globalIdxI << std::endl;

        // shortcurts mpfa case
        CellData& cellData3 = problem().variables().cellData(globalIdx3);
        CellData& cellData4 = problem().variables().cellData(globalIdx4);
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
                problem().variables().getMpfaData3D(*isIt, additionalT,
                                globalPosAdditional3, globalIdxAdditional3,
                                globalPosAdditional4, globalIdxAdditional4 ,
                                banana); // offset for second interaction region

                Scalar gravityContributionAdditonal
                    = temp1 * additionalT[0] + temp2 * additionalT[1]
                        + globalPosAdditional3*this->gravity_ * additionalT[2]
                        + globalPosAdditional4*this->gravity_ * additionalT[3];
                CellData& cellDataA3 = problem().variables().cellData(globalIdxAdditional3);
                CellData& cellDataA4 = problem().variables().cellData(globalIdxAdditional4);

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
                if(cellDataI.isUpwindCell(isIt->indexInInside(), eqIdx))
                    upwindCellData[phaseIdx] = &cellDataI;
                else if(cellDataJ.isUpwindCell(isIt->indexInOutside(), eqIdx))
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
    this->A_[globalIdxI][globalIdxI] += (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * T[0];
    this->A_[globalIdxI][globalIdxJ] += (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * T[1];
    this->A_[globalIdxI][globalIdx3] += (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * T[2];
    this->A_[globalIdxI][globalIdx4] += (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * T[3];


    // add gravity to RHS vector
    Scalar gravityContribution = temp1 * T[0] + temp2 * T[1] + temp3 * T[2] + temp4 * T[3];
    this->f_[globalIdxI] += (rhoMean[wPhaseIdx] * lambda[wPhaseIdx] * dV[wPhaseIdx]
                              + rhoMean[nPhaseIdx] * lambda[nPhaseIdx] * dV[nPhaseIdx]) * gravityContribution;

    // weithing accounts for the fraction of the subcontrol volume
    Scalar weightingFactor = volume / perimeter;    // transforms flux through area A -> V * A/perimeter
    if(enableVolumeIntegral_) // switch off volume integral for mpfa case
    {
        // correct for area integral
        this->A_[globalIdxI][globalIdxI] -=
                weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * T[0];
        this->A_[globalIdxI][globalIdxJ] -=
                weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * T[1];
        this->A_[globalIdxI][globalIdx3] -=
                weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * T[2];
        this->A_[globalIdxI][globalIdx4] -=
                weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * T[3];

        // add gravity to RHS vector
        this->f_[globalIdxI] -= weightingFactor * gravityContribution *
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

    this->f_[globalIdxI] += pcGradient;

    // regard more interaction regions, if there are more
    if(interactionRegions != 1)
    {
        for(int banana = 1; banana < interactionRegions; banana ++)
        {
            // get data for second interaction region
            problem().variables().getMpfaData3D(*isIt, additionalT,
                            globalPosAdditional3, globalIdxAdditional3,
                            globalPosAdditional4, globalIdxAdditional4 ,
                            banana); // offset for second interaction region

            Scalar gravityContributionAdditonal
                = temp1 * additionalT[0] + temp2 * additionalT[1]
                    + globalPosAdditional3*this->gravity_ * additionalT[2]
                    + globalPosAdditional4*this->gravity_ * additionalT[3];
            CellData& cellDataA3 = problem().variables().cellData(globalIdxAdditional3);
            CellData& cellDataA4 = problem().variables().cellData(globalIdxAdditional4);

            /* compute matrix entry: advective fluxes */
            /* extend T with other matrix entries and assemble to A_    */
            this->A_[globalIdxI][globalIdxI] +=
                    (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * additionalT[0];
            this->A_[globalIdxI][globalIdxJ] +=
                    (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * additionalT[1];
            this->A_[globalIdxI][globalIdxAdditional3] +=
                    (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * additionalT[2];
            this->A_[globalIdxI][globalIdxAdditional4] +=
                    (lambda[wPhaseIdx] * dV[wPhaseIdx] + lambda[nPhaseIdx] * dV[nPhaseIdx]) * additionalT[3];


            // add gravity to RHS vector
            this->f_[globalIdxI] += (rhoMean[wPhaseIdx] * lambda[wPhaseIdx] * dV[wPhaseIdx]
                                     + rhoMean[nPhaseIdx] * lambda[nPhaseIdx] * dV[nPhaseIdx]) * gravityContributionAdditonal;

            // weithing accounts for the fraction of the subcontrol volume
            Scalar weightingFactor = volume / perimeter;    // transforms flux through area A -> V * A/perimeter
            if(enableVolumeIntegral_) // switch off volume integral for mpfa case
            {
                // correct for area integral
                this->A_[globalIdxI][globalIdxI] -=
                        weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * additionalT[0];
                this->A_[globalIdxI][globalIdxJ] -=
                        weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * additionalT[1];
                this->A_[globalIdxI][globalIdxAdditional3] -=
                        weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * additionalT[2];
                this->A_[globalIdxI][globalIdxAdditional4] -=
                        weightingFactor * (lambda[wPhaseIdx] * gV[wPhaseIdx] + lambda[nPhaseIdx] * gV[nPhaseIdx]) * additionalT[3];

                // add gravity to RHS vector
                this->f_[globalIdxI] -= weightingFactor * gravityContribution *
                        (rhoMean[wPhaseIdx] * lambda[wPhaseIdx] * gV[wPhaseIdx] + rhoMean[nPhaseIdx] * lambda[nPhaseIdx] * gV[nPhaseIdx]);
            }

            // capillary pressure flux
            Scalar pcGradient = cellDataI.capillaryPressure() * additionalT[0]
                               + cellDataJ.capillaryPressure() * additionalT[1]
                               + cellDataA3.capillaryPressure() * additionalT[2]
                               + cellDataA4.capillaryPressure() * additionalT[3];

            if (this->pressureType == pw)
                pcGradient *= + lambda[nPhaseIdx] * dV[nPhaseIdx]
                               - enableVolumeIntegral_ * weightingFactor * lambda[nPhaseIdx] * gV[nPhaseIdx];
            else if (this->pressureType == pn)
                pcGradient *= - lambda[wPhaseIdx] * dV[wPhaseIdx]
                               + enableVolumeIntegral_ * weightingFactor * lambda[wPhaseIdx] * gV[wPhaseIdx];

            this->f_[globalIdxI] += pcGradient;
        }
    }
}

//! Compute single-phase flux through an irregular interface using a \a mpfa method
/*! A mpfa l-method is applied to calculate fluxes near hanging nodes for
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
    // acess Cell I
    ElementPointer elementPointerI = isIt->inside();
    int globalIdxI = problem().variables().index(*elementPointerI);

    // get global coordinate of cell center
    const GlobalPosition& globalPos = elementPointerI->geometry().center();

    // access neighbor
    ElementPointer neighborPointer = isIt->outside();
    int globalIdxJ = problem().variables().index(*neighborPointer);
    CellData& cellDataJ = problem().variables().cellData(globalIdxJ);

    // gemotry info of neighbor
    const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().center();

    // due to "safety cell" around subdomain, both cells I and J
    // have single-phase conditions, although one is in 2p domain.
    int phaseIdx = std::min(cellDataI.subdomain(), cellDataJ.subdomain());

    // get average density for gravity flux
    Scalar rhoMean = 0.5 * (cellDataI.density(phaseIdx) + cellDataJ.density(phaseIdx));

    // 1p => no pC => only 1 pressure, potential
    Scalar potential = 0.;

        // Prepare MPFA
        /** get geometrical Info, transmissibility matrix */
        GlobalPosition globalPos3(0.);
        int globalIdx3=-1;
        GlobalPosition globalPos4(0.);
        int globalIdx4=-1;
        TransmissivityMatrix T(0.);

            // prepare second interaction region
            GlobalPosition globalPosAdditional3(0.);
            int globalIdxAdditional3=-1;
            GlobalPosition globalPosAdditional4(0.);
            int globalIdxAdditional4=-1;

            TransmissivityMatrix additionalT(0.);

        int interactionRegions
            = problem().variables().getMpfaData3D(*isIt, T,
                                    globalPos3, globalIdx3, globalPos4, globalIdx4  );
        if (interactionRegions == 0)
            interactionRegions = problem().pressureModel().computeTransmissibilities(isIt,T,
                    globalPos3, globalIdx3,  globalPos4, globalIdx4  );
        if(!interactionRegions)
            Dune::dgrave << "something went wrong getting mpfa data on cell " << globalIdxI << std::endl;

        // shortcurts mpfa case
        CellData& cellData3 = problem().variables().cellData(globalIdx3);
        CellData& cellData4 = problem().variables().cellData(globalIdx4);
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
                problem().variables().getMpfaData3D(*isIt, additionalT,
                                globalPosAdditional3, globalIdxAdditional3,
                                globalPosAdditional4, globalIdxAdditional4 ,
                                banana); // offset for second interaction region

                Scalar gravityContributionAdditonal
                    = temp1 * additionalT[0] + temp2 * additionalT[1]
                        + globalPosAdditional3*this->gravity_ * additionalT[2]
                        + globalPosAdditional4*this->gravity_ * additionalT[3];
                CellData& cellDataA3 = problem().variables().cellData(globalIdxAdditional3);
                CellData& cellDataA4 = problem().variables().cellData(globalIdxAdditional4);

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
    this->A_[globalIdxI][globalIdxI] += lambda * T[0];
    this->A_[globalIdxI][globalIdxJ] += lambda * T[1];
    this->A_[globalIdxI][globalIdx3] += lambda * T[2];
    this->A_[globalIdxI][globalIdx4] += lambda * T[3];

    // add gravity to RHS vector
    Scalar gravityContribution = temp1 * T[0] + temp2 * T[1] + temp3 * T[2] + temp4 * T[3];
    this->f_[globalIdxI] += lambda * rhoMean * gravityContribution;

    // regard more interaction regions, if there are more
    if(interactionRegions != 1)
    {
        for(int banana = 1; banana < interactionRegions; banana ++)
        {
            // get data for second interaction region
            problem().variables().getMpfaData3D(*isIt, additionalT,
                            globalPosAdditional3, globalIdxAdditional3,
                            globalPosAdditional4, globalIdxAdditional4 ,
                            banana); // offset for second interaction region

            Scalar gravityContributionAdditonal
                = temp1 * additionalT[0] + temp2 * additionalT[1]
                    + globalPosAdditional3*this->gravity_ * additionalT[2]
                    + globalPosAdditional4*this->gravity_ * additionalT[3];

            /** compute matrix entry: advective fluxes */
            /* extend T with other matrix entries and assemble to A_    */
            this->A_[globalIdxI][globalIdxI] += lambda * additionalT[0];
            this->A_[globalIdxI][globalIdxJ] += lambda * additionalT[1];
            this->A_[globalIdxI][globalIdxAdditional3] += lambda * additionalT[2];
            this->A_[globalIdxI][globalIdxAdditional4] += lambda * additionalT[3];

            // add gravity to RHS vector
            this->f_[globalIdxI] += lambda * rhoMean* gravityContributionAdditonal;
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
        ElementIterator eItEnd = problem().gridView().template end<0> ();
        for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
        {
            int globalIdx = problem().variables().index(*eIt);
            CellData& cellData = problem().variables().cellData(globalIdx);

            if(cellData.fluidStateType() == 0) // i.e. it is complex
                problem().pressureModel().updateMaterialLawsInElement(*eIt, fromPostTimestep);
            else
                problem().pressureModel().update1pMaterialLawsInElement(*eIt, cellData, fromPostTimestep);

            maxError = std::max(maxError, fabs(cellData.volumeError()));
        }
        this->maxError_ = maxError/problem().timeManager().timeStepSize();
        return;
    }
    else
    {
        // resize multiphysics vectors
        FVPressure2P2CMultiPhysics<TypeTag>::nextSubdomain.resize(problem().gridView().size(0));

        FVPressure2P2CMultiPhysics<TypeTag>::updateMaterialLaws();
    }
}

//! Computes the transmissibility coefficients for the MPFA-l method in 3D
/*! For faces with a hanging node in 3D, there are four sub-faces. The first subface
 * contains a unique interaction volume, with is directly calculated by this method.
 * For the remainder of the sub-faces, the interaction volumes are build and calculated
 * by the common mpfa-l--implementation of the 2p models. The latter is established via the
 * protected method FV3dPressure2P2CAdaptive::transmissibilityAdapter_().
 * The calculated Transmissivity Matrices are (along with some geometric information)
 * stored for later use in Dumux::Variableclass2p2cadaptive .
 *
* \param isIt Iterator to the current intersection
* \param T Transmissitivity matrix of the first unique interaction volume
* \param globalPos4 Position of the 3rd cell (with local Idx 4) of the unique interaction volume
* \param globalIdx4 Index of the 3rd cell (with local Idx 4) of the unique interaction volume
* \param globalPos6 Position of the 4th cell (with local Idx 6) of the unique interaction volume
* \param globalIdx6 Index of the 4th cell (with local Idx 6) of the unique interaction volume
*/
template <class TypeTag>
int FV3dPressure2P2CAdaptive<TypeTag>::computeTransmissibilities(const IntersectionIterator& isIt,
        TransmissivityMatrix& T,
        GlobalPosition& globalPos4,
        int& globalIdx4,
        GlobalPosition& globalPos6,
        int& globalIdx6)
{
    // get geometry information of cellI = cell1, cellJ = cell2
    ElementPointer eIt = isIt->inside();
//    int globalIdxI =  problem().variables().index(*eIt);
    ElementPointer neighborPointer = isIt->outside();
    GlobalPosition globalPos1 = eIt->geometry().center();
    GlobalPosition globalPos2 = neighborPointer->geometry().center();
    DimMatrix K1(problem().spatialParams().intrinsicPermeability(*eIt));
    DimMatrix K2(problem().spatialParams().intrinsicPermeability(*neighborPointer));
//    int globalIdxJ =  problem().variables().index(*isIt->outside());

    // determine ID of intersection seen from larger cell
    int intersectionID = 0;
    if(isIt->inside()->level() < isIt->outside()->level())
        intersectionID = problem().grid().localIdSet().subId(*eIt,
                isIt->indexInInside(), 1);
    else
        DUNE_THROW(Dune::NotImplemented, " ABORT, transmiss calculated from wrong side!!");

    std::vector<int> localIrregularCells = irregularCellMap_[intersectionID];

    /** 1) get geometrical information of interaction triangle   */
    // geometry and Data of face IJ in nomenclature of mpfa
    GlobalPosition globalPosFace12 = isIt->geometry().center();  // = 'x'1
    GlobalPosition outerNormaln12 = isIt->centerUnitOuterNormal();

    /**** a) search for intersection 24 and 26 ************/
    IntersectionIterator face24=isIt; // as long as face24 = isIt, it is still not found!
    IntersectionIterator face26=isIt; // as long as face26 = isIt, it is still not found!

    IntersectionIterator nextIsEnd = problem().gridView().iend(*neighborPointer);
    for (IntersectionIterator isIt2 = problem().gridView().ibegin(*neighborPointer); isIt2 != nextIsEnd; ++isIt2)
    {
        // continue if no neighbor or arrived at intersection
        if(!(isIt2->neighbor()) or isIt2->outside() == eIt)
            continue;

        int currentNeighbor = problem().variables().index(*isIt2->outside());

        // have we found face24?
        if (find(localIrregularCells.begin(), localIrregularCells.end(),
                currentNeighbor) != localIrregularCells.end() && face24==isIt)
            face24 = isIt2;
        else if (find(localIrregularCells.begin(), localIrregularCells.end(),
                currentNeighbor) != localIrregularCells.end() && face26==isIt)
        {
            // we now found both intersections, but we have to investigate orientation:
            // Calculate the vector product of the normals in current orientation
            GlobalPosition vectorProduct = crossProduct(face24->centerUnitOuterNormal(), isIt2->centerUnitOuterNormal());
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
    globalPos4 = face24->outside()->geometry().center();
    globalIdx4 = problem().variables().index(*(face24->outside()));
    GlobalPosition outerNormaln24 = face24->centerUnitOuterNormal();
    // get absolute permeability of neighbor cell 3
    DimMatrix K4(problem().spatialParams().intrinsicPermeability(*(face24->outside())));

    // get information of cell6
    globalPos6 = face26->outside()->geometry().center();
    globalIdx6 = problem().variables().index(*(face26->outside()));
    GlobalPosition outerNormaln26 = face26->centerUnitOuterNormal();
    // get absolute permeability of neighbor cell 3
    DimMatrix K6(problem().spatialParams().intrinsicPermeability(*(face26->outside())));

    /**** b) Get Points on the edges (in plane of isIt), 'x'4 and 'x'5 ***********/
    int localFace12 = isIt->indexInOutside(); // isIt is face12
    int localFace24 = face24->indexInInside();
    int localFace26 = face26->indexInInside();

    const ReferenceElement& referenceElement = ReferenceElementContainer::general(neighborPointer->geometry().type());
    //find 'x'5 = edgeCoord1226
    int edge1226;
    // search through edges of face 12
    for(int nectarine=0; nectarine < referenceElement.size(localFace12, 1, dim-1); nectarine++)
    {
        // get local Idx of edge on face 12
        int localEdgeOn12 = referenceElement.subEntity(localFace12, 1, nectarine, dim-1);
        // search through edges of face 26
        for(int plum = 0; plum < referenceElement.size(localFace26, 1,dim-1); plum++)
        {
//            int localEdge26 = referenceElement.subEntity(localFace26, 1, plum, dim-1);
            if(referenceElement.subEntity(localFace12, 1, nectarine, dim-1)
                    == referenceElement.subEntity(localFace26, 1, plum, dim-1))
            {
                edge1226 = localEdgeOn12;
                break;
            }
        }
    }
    GlobalPosition edgeCoord1226 =  // 'x'5
            neighborPointer->geometry().global(referenceElement.position(edge1226, dim-1));

    //find 'x'4 = edgeCoord1224
    int edge1224;
    // search through edges of face 12
    for(int nectarine=0; nectarine < referenceElement.size(localFace12, 1, dim-1); nectarine++)
    {
        // get local Idx of edge on face 12
        int localEdgeOn12 = referenceElement.subEntity(localFace12, 1, nectarine, dim-1);
        // search through edges of face 24
        for(int plum = 0; plum < referenceElement.size(localFace24, 1, dim-1); plum++)
        {
            int localEdge24 = referenceElement.subEntity(localFace24, 1, plum, dim-1);
            if(localEdgeOn12 == localEdge24)
            {
                edge1224 = localEdgeOn12;
                break;
            }
        }
    }
    GlobalPosition edgeCoord1224 =  // 'x'4
            neighborPointer->geometry().global(referenceElement.position(edge1224, dim-1));

    //find 'x'6 = edgeCoord2426
    int edge2426;
    // search through edges of face 24
    for(int nectarine=0; nectarine < referenceElement.size(localFace24, 1, dim-1); nectarine++)
    {
        // get local Idx of edge on face 24
        int localEdgeOn24 = referenceElement.subEntity(localFace24, 1, nectarine, dim-1);
        // search through edges of face 26
        for(int plum = 0; plum < referenceElement.size(localFace26, 1, dim-1); plum++)
        {
            int localEdge26 = referenceElement.subEntity(localFace26, 1, plum, dim-1);
            if(localEdgeOn24 == localEdge26)
            {
                edge2426 = localEdgeOn24;
                break;
            }
        }
    }
    GlobalPosition edgeCoord2426 =   // 'x'6
            neighborPointer->geometry().global(referenceElement.position(edge2426, dim-1));

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
    	T *= isIt->geometry().volume()/subFaceArea12 ;

        // set your map entry
        problem().variables().storeMpfaData3D(*isIt, T, globalPos4, globalIdx4, globalPos6, globalIdx6);
        return 1; // indicates that only 1 interaction region was regarded
    }
    else
    {
        int countInteractionRegions = 1;
        //initialize additional transmissitivity matrix
        TransmissivityMatrix additionalT(0.);

        VertexPointer outerCornerPtr(isIt->inside()->template subEntity<dim>(0)); //initialize with rubbish
        // prepare additonal pointer to cells
        ElementPointer additional2(isIt->inside()); //initialize with something wrong!
        ElementPointer additional3(isIt->inside());
        int caseL = -2;

        /**** 2nd interaction region: get corner of interest ************/
        // search through corners of large cell with isIt
        int localIdxLarge = searchCommonVertex_(*isIt, outerCornerPtr);

        //in the parallel case, skip all border entities
        #if HAVE_MPI
        if (problem().gridView().comm().size() > 1)
            if(outerCornerPtr->partitionType() != Dune::InteriorEntity)
                caseL = -1; // abort this specific interaction volume
        #endif

        // get Interaction Volume object
        int vIdxGlobal = problem().variables().index(*outerCornerPtr);
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
            problem().variables().storeMpfaData3D(*isIt, additionalT,
                    additional2->geometry().center(), problem().variables().index(*additional2),
                    additional3->geometry().center(), problem().variables().index(*additional3),
                    1); // offset for second interaction region
            countInteractionRegions++;
        }

        /***** 3rd and 4th interaction region *******/
        if(maxInteractionVolumes>2)
        {
            // loop through remaining 2 points
            std::vector<int> diagonal;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            for(int verticeSmall = 0; verticeSmall < isIt->outside()->subEntities(dim); ++verticeSmall)
#else
            for(int verticeSmall = 0; verticeSmall<isIt->outside()->template count<dim>(); ++verticeSmall)
#endif
            {
                const VertexPointer vPtrSmall = isIt->outside()->template subEntity<dim>(verticeSmall);

                //in the parallel case, skip all border entities
                #if HAVE_MPI
                if (problem().gridView().comm().size() > 1)
                    if(vPtrSmall->partitionType() != Dune::InteriorEntity)
                        continue;
                #endif

                // get postion as seen from element
                GlobalPosition vertexOnElement
                    = referenceElement.position(verticeSmall, dim);

                for (int indexOnFace = 0; indexOnFace < 4; indexOnFace++)
                {
                    // get position as seen from interface
                    GlobalPosition vertexOnInterface
                        = isIt->geometryInOutside().corner(indexOnFace);

                    if(vPtrSmall != outerCornerPtr
                            && ((vertexOnInterface - vertexOnElement).two_norm()<1e-5))
                    {
                        int vIdxGlobal = problem().variables().index(*vPtrSmall);
                        // acess interactionVolume
                        InteractionVolume& interactionVolume
                            = interactionVolumesContainer_->interactionVolume(vIdxGlobal);
                        if(interactionVolume.isBoundaryInteractionVolume())
                            continue;

                        int hangingNodeType = interactionVolume.getHangingNodeType();
                        // reset flux direction indicator
                        properFluxDirection = true;

                        if(hangingNodeType != InteractionVolume::fourSmallCellsFace)
                        {
                            diagonal.push_back(problem().variables().index(*vPtrSmall));
                            // a) take interaction volume and determine faceIdx
                            if(hangingNodeType == InteractionVolume::noHangingNode)
                            {
                                //TODO determine current localIdxLarge!!!!
                                Dune::dgrave << " args, noHangingNode on additional interaction region";
    //                            subVolumeFaceIdx = getMpfaCase8cells_(isIt, localIdxLarge, interactionVolume, properFluxDirection);
                            }
                            else if(hangingNodeType == InteractionVolume::sixSmallCells)
                                subVolumeFaceIdx = interactionVolumesContainer_->getMpfaCase6cells(isIt,
                                                                        interactionVolume, properFluxDirection);
                            else
                                subVolumeFaceIdx = interactionVolumesContainer_->getMpfaCase2or4cells(isIt,
                                                                        interactionVolume, properFluxDirection);

                            // b) calculate T, globalIdx3+4
                            caseL = this->transmissibilityAdapter_(isIt, interactionVolume, subVolumeFaceIdx,
                                                        properFluxDirection, additional2, additional3, additionalT);

                            // c) store it
                            //store everything as third/4th interaction region
                            if(caseL != -1) //check if we regard this interaction region
                            {
                                problem().variables().storeMpfaData3D(*isIt, additionalT,
                                        additional2->geometry().center(), problem().variables().index(*additional2),
                                        additional3->geometry().center(), problem().variables().index(*additional3),
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
        problem().variables().storeMpfaData3D(*isIt, T, globalPos4, globalIdx4, globalPos6, globalIdx6);

        // determine weights
        Scalar weight = isIt->geometry().volume()/subFaceArea12; // =4 for 1 interaction region
        weight /= static_cast<Scalar>(countInteractionRegions);

        // perform weighting of transmissibilies
        problem().variables().performTransmissitivityWeighting(*isIt, weight);

        return countInteractionRegions;
    }

    return 0;
}

//! Adapter to use the general implementation of the mpfa-l for the compositional models
/*! Depending on the subVolumeFaceIdx, the appropriate method in
 * Dumux::FvMpfaL2dTransmissibilityCalculator (potentially specifying certain cases)
 * gets called and the transmissibility and geometric information of the applied additional
 * cells of the interaction regions are passed back.
 *
* \param isIt Iterator to the current intersection
* \param interactionVolume The current interaction Volume object of interest
* \param subVolumeFaceIdx The local index of the intersection of interest in the interaction volume
* \param properFluxDirection True if the intersection normal coincides
*           with the local indexing in the interaction volume
* \param[out] additional2 Pointer to the 3rd cell's element in the interaction volume
* \param[out] additional3 Pointer to the 4th cell's element in the interaction volume
* \param[out] additionalT Transmissitivity matrix calculated
*/
template<class TypeTag>
int FV3dPressure2P2CAdaptive<TypeTag>::transmissibilityAdapter_(const IntersectionIterator& isIt,
                                InteractionVolume& interactionVolume,
                                const int& subVolumeFaceIdx,
                                bool properFluxDirection,
                                ElementPointer& additional2,
                                ElementPointer& additional3,
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

//    static Scalar transChoiceThreshold_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, MPFA, TransmissibilityCriterionThreshold);
//    if(!properFluxDirection)
//        mpfal3DTransmissibilityCalculator_.setTransChoiceThreshold(-transChoiceThreshold_);

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
                assert(problem().variables().index(*(interactionVolume.getSubVolumeElement(4)))
                        != problem().variables().index(*(interactionVolume.getSubVolumeElement(5)))); // else there would not be a subVolFaceIdx 4
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
                assert (problem().variables().index(*(interactionVolume.getSubVolumeElement(4)))
                        != problem().variables().index(*(interactionVolume.getSubVolumeElement(6)))); // else there would not be a subVolFaceIdx 5
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
                assert (problem().variables().index(*(interactionVolume.getSubVolumeElement(4)))
                        != problem().variables().index(*(interactionVolume.getSubVolumeElement(5)))); // else there would not be a subVolFaceIdx 4
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
                assert (problem().variables().index(*(interactionVolume.getSubVolumeElement(4)))
                        != problem().variables().index(*(interactionVolume.getSubVolumeElement(6)))); // else there would not be a subVolFaceIdx 5
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
                                or hangingNodeType == InteractionVolume::sixSmallCells)
                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 4, 0, 6, 2, 5, 1);
            else if (hangingNodeType == InteractionVolume::twoSmallCells
                            || hangingNodeType == InteractionVolume::fourSmallCellsFace)
            {
                caseL = mpfal3DTransmissibilityCalculator_.transmissibilityCaseTwo(T, interactionVolume,
                                        lambda, 4, 0, 2, 1);
            }
            else if (hangingNodeType == InteractionVolume::fourSmallCellsDiag
                            || (hangingNodeType == InteractionVolume::fourSmallCellsEdge
                                && (problem().variables().index(*(interactionVolume.getSubVolumeElement(4)))
                                    != problem().variables().index(*(interactionVolume.getSubVolumeElement(6))))) )
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
                                or hangingNodeType == InteractionVolume::sixSmallCells)
            caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 1, 5, 3, 7, 0, 4);
            else if (hangingNodeType == InteractionVolume::twoSmallCells || hangingNodeType == InteractionVolume::fourSmallCellsFace)
            {
                caseL = mpfal3DTransmissibilityCalculator_.transmissibilityCaseOne(T, interactionVolume,
                                        lambda, 1, 5, 3, 0);
            }
            else if (hangingNodeType == InteractionVolume::fourSmallCellsDiag
                    || (hangingNodeType == InteractionVolume::fourSmallCellsEdge
                        &&(problem().variables().index(*(interactionVolume.getSubVolumeElement(4)))
                            != problem().variables().index(*(interactionVolume.getSubVolumeElement(6)))) ))
            {
                useCases[0] = true;
                useCases[1] = false;
                useCases[2] = true;
                useCases[3] = false;

                caseL = mpfal3DTransmissibilityCalculator_.transmissibility(T, interactionVolume,
                                        lambda, 1, 5, 3, 7, 0, 4, useCases);
            }
            else if (hangingNodeType == InteractionVolume::fourSmallCellsEdge
                    &&(problem().variables().index(*(interactionVolume.getSubVolumeElement(4)))
                            != problem().variables().index(*(interactionVolume.getSubVolumeElement(5)))) )
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
                      &&(problem().variables().index(*(interactionVolume.getSubVolumeElement(4)))
                            != problem().variables().index(*(interactionVolume.getSubVolumeElement(6)))) )
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
                            && (problem().variables().index(*(interactionVolume.getSubVolumeElement(4)))
                                != problem().variables().index(*(interactionVolume.getSubVolumeElement(6)))))
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
                        && (problem().variables().index(*(interactionVolume.getSubVolumeElement(4)))
                            != problem().variables().index(*(interactionVolume.getSubVolumeElement(5)))) )
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
        std::swap(additionalT[0], additionalT[1]);
    }

    return caseL;
}


}//end namespace Dumux
#endif
