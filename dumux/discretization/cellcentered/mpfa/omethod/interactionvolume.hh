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
 * \brief Base classes for interaction volume of mpfa methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_INTERACTIONVOLUME_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_INTERACTIONVOLUME_HH

#include <dumux/common/math.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/interactionvolumebase.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

#include "interactionvolumeseed.hh"
#include "localsubcontrolentities.hh"

namespace Dumux
{
//! Specialization of the interaction volume traits class for the mpfa-o method
template<class TypeTag>
class CCMpfaOInteractionVolumeTraits : public CCMpfaInteractionVolumeTraitsBase<TypeTag>
{
    using BaseTraits = CCMpfaInteractionVolumeTraitsBase<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    using BoundaryInteractionVolume = CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethod>;

    using PositionVector = std::vector<GlobalPosition>;
    using Matrix = Dune::DynamicMatrix<Scalar>;
    using Vector = typename Matrix::row_type;

    using LocalScvType = CCMpfaOLocalScv<TypeTag>;
    using LocalScvfType = CCMpfaOLocalScvf<TypeTag>;
    using typename BaseTraits::LocalIndexSet;
    using typename BaseTraits::GlobalIndexSet;
    using Seed = CCMpfaOInteractionVolumeSeed<GlobalIndexSet, LocalIndexSet, dim, dimWorld>;
};

//! Forward declaration of the mpfa-o interaction volume
template<class TypeTag, class Traits, class Implementation>
class CCMpfaOInteractionVolume;

/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volumes of the mpfa-o method.
 *        We introduce one more level of inheritance here because the o-method with
 *        full pressure support uses the mpfa-o interaction volume but uses a different
 *        traits class.
 */
template<class TypeTag>
class CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethod>
          : public CCMpfaOInteractionVolume<TypeTag,
                                            CCMpfaOInteractionVolumeTraits<TypeTag>,
                                            CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethod>>
{
    using ThisType = CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethod>;
    using TraitsType = CCMpfaOInteractionVolumeTraits<TypeTag>;
    using ParentType = CCMpfaOInteractionVolume<TypeTag, TraitsType, ThisType>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using IVSeed = typename TraitsType::Seed;
public:
    // state the traits class type
    using Traits = TraitsType;
};

template<class TypeTag, class Traits, class Implementation>
class CCMpfaOInteractionVolume : public CCMpfaInteractionVolumeBase<TypeTag, Traits>
{
    // The actual implementation has to be friend
    friend Implementation;

    using ParentType = CCMpfaInteractionVolumeBase<TypeTag, Traits>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using InteriorBoundaryData = typename GET_PROP_TYPE(TypeTag, InteriorBoundaryData);

    static constexpr bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);
    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using Vector = typename Traits::Vector;
    using Matrix = typename Traits::Matrix;
    using Tensor = typename Traits::Tensor;

    using LocalScvType = typename Traits::LocalScvType;
    using LocalScvfType = typename Traits::LocalScvfType;

    using LocalIndexType = typename Traits::LocalIndexType;
    using LocalIndexSet = typename Traits::LocalIndexSet;
    using GlobalIndexSet = typename Traits::GlobalIndexSet;
    using PositionVector = typename Traits::PositionVector;
    using DataHandle = typename Traits::DataHandle;
    using Seed = typename Traits::Seed;

public:

    using typename ParentType::GlobalLocalFaceDataPair;
    using typename ParentType::LocalFaceData;

    //! Sets up the local scope for a given seed
    //! This function has to be called before using the IV!
    void bind(const Seed& seed,
              const Problem& problem,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              DataHandle& dataHandle)
    {
        problemPtr_ = &problem;
        fvGeometryPtr_ = &fvGeometry;
        elemVolVarsPtr_ = &elemVolVars;
        onDomainOrInteriorBoundary_ = seed.onDomainOrInteriorBoundary();

        // set up the local scope of this interaction volume
        setLocalScope_(seed, dataHandle);

        // initialize the vector containing the neumann fluxes
        asImp_().assembleNeumannFluxVector();
    }

    //! solves for the transmissibilities subject to a given tensor
    template<typename GetTensorFunction>
    void solveLocalSystem(const GetTensorFunction& getTensor, DataHandle& dataHandle)
    {
        // if only dirichlet faces are present, assemble T_ directly
        if (numFluxFaces_ == 0)
            assemblePureDirichletSystem_(getTensor, dataHandle.T());
        else
        {
            // assemble
            assembleLocalMatrices_(getTensor);

            // solve
            A_.invert();

            // set up T-matrix
            auto& CA = dataHandle.CA();
            auto& T = dataHandle.T();
            CA = C_.rightmultiply(A_);
            T = multiplyMatrices(CA, B_);
            T += D_;
        }

        // set vol vars stencil & positions pointer in handle
        dataHandle.setVolVarsStencilPointer(volVarsStencil());
        dataHandle.setVolVarsPositionsPointer(volVarsPositions());

        // store A-1B only when gradient reconstruction is necessary
        if (requireABMatrix_())
            dataHandle.AB() = B_.leftmultiply(A_);

        // on surface grids, additionally prepare the outside transmissibilities
        if (dim < dimWorld)
            computeOutsideTransmissibilities_(dataHandle);
    }

    //! Gets the transmissibilities for a sub-control volume face within the interaction volume.
    //! specialization for dim == dimWorld
    template<int d = dim, int dw = dimWorld>
    typename std::enable_if< (d == dw), const Vector& >::type
    getTransmissibilities(const SubControlVolumeFace& scvf, const LocalFaceData& localFaceData, const DataHandle& dataHandle) const
    { return dataHandle.T()[localFaceData.localScvfIndex]; }

    //! Gets the transmissibilities for a sub-control volume face within the interaction volume.
    //! specialization for dim < dimWorld.
    template<int d = dim, int dw = dimWorld>
    typename std::enable_if< (d < dw), const Vector& >::type
    getTransmissibilities(const SubControlVolumeFace& scvf, const LocalFaceData& localFaceData, const DataHandle& dataHandle) const
    {
        // If we come from the inside, simply return tij
        if (!localFaceData.isOutside)
            return dataHandle.T()[localFaceData.localScvfIndex];
        else
            return dataHandle.outsideTij()[this->findIndexInVector(outsideScvfIndices_, scvf.index())];
    }

    //! Returns the vector of coefficients with which the vector of (interior) neumann boundary conditions
    //! have to be multiplied in order to transform them on the scvf this cache belongs to
    const Vector& getNeumannFluxTransformationCoefficients(const SubControlVolumeFace& scvf, const LocalFaceData& localFaceData, const DataHandle& dataHandle) const
    { return dataHandle.CA()[localFaceData.localScvfIndex]; }

    //! Returns the boundary neumann flux for a given face
    Scalar getNeumannFlux(const SubControlVolumeFace& scvf, const LocalFaceData& localFaceData, const DataHandle& dataHandle, unsigned int eqIdx) const
    {
        if (!onDomainOrInteriorBoundary() || useTpfaBoundary || fluxScvfIndexSet().size() == 0 )
            return 0.0;

        // Do the scalar product CAinv_*neumannFluxes[eqIdx]
        assert(dataHandle.CA()[localFaceData.localScvfIndex].size() == neumannFluxes_.size() &&
               "Number of neumann flux entries does not match with size of coefficent vector!");

        Scalar flux(0.0);
        for (unsigned int i = 0; i < neumannFluxes_.size(); ++i)
            flux += dataHandle.CA()[localFaceData.localScvfIndex][i]*neumannFluxes_[i][eqIdx];
        return flux;
    }

    // Per default we do not add additional terms for interior Neumann boundaries
    // Overload this function to realize interior membranes etc...
    template<typename GetTensorFunction>
    Scalar interiorNeumannTerm(const GetTensorFunction& getTensor,
                               const Element& element,
                               const LocalScvfType& localScvf,
                               const InteriorBoundaryData& data) const
    { return 0.0; }

    //! This function sets up the vector which stores the neumann boundary conditions in the IV
    void assembleNeumannFluxVector()
    {
        // initialize the neumann fluxes vector to zero
        neumannFluxes_.resize(fluxFaceIndexSet_.size(), PrimaryVariables(0.0));

        if (!onDomainOrInteriorBoundary() || useTpfaBoundary)
            return;

        LocalIndexType fluxFaceIdx = 0;
        for (auto localFluxFaceIdx : fluxFaceIndexSet_)
        {
            const auto& curLocalScvf = localScvf(localFluxFaceIdx);
            const auto faceType = curLocalScvf.faceType();

            if (faceType == MpfaFaceTypes::neumann || faceType == MpfaFaceTypes::interiorNeumann)
            {
                const auto& element = localElement(curLocalScvf.insideLocalScvIndex());
                const auto& globalScvf = fvGeometry().scvf(curLocalScvf.insideGlobalScvfIndex());
                auto neumannFlux = problem().neumann(element, this->fvGeometry(), this->elemVolVars(), globalScvf);
                neumannFlux *= globalScvf.area();
                neumannFlux *= elemVolVars()[globalScvf.insideScvIdx()].extrusionFactor();

                // The flux is assumed to be prescribed in the form of -D*gradU
                neumannFluxes_[fluxFaceIdx] = neumannFlux;
            }

            fluxFaceIdx++;
        }
    }

    //! obtain the local data object for a given global scvf
    const LocalFaceData& getLocalFaceData(const SubControlVolumeFace& scvf) const
    { return globalLocalScvfPairedData_[this->findIndexInVector(globalScvfIndices_, scvf.index())].second; }

    //! returns whether or not a domain or interior boundary is in the iv
    bool onDomainOrInteriorBoundary() const
    { return onDomainOrInteriorBoundary_; }

    //! returns the vector of data pairs storing the local face data for each global scvf
    const std::vector<GlobalLocalFaceDataPair>& globalLocalScvfPairedData() const
    { return globalLocalScvfPairedData_; }

    //! returns the volume variables stencil of this interaction volume
    const GlobalIndexSet& volVarsStencil() const
    { return volVarsStencil_; }

    //! returns positions (cells and dirichlet faces) the volvars in this interaction volume
    const PositionVector& volVarsPositions() const
    { return volVarsPositions_; }

    //! returns the grid element corresponding to a given local (in the iv) scv idx
    const Element& localElement(const LocalIndexType localScvIdx) const
    { return localElements_[localScvIdx]; }

    //! returns the local scvf entity corresponding to a given local (in the iv) scvf idx
    const LocalScvfType& localScvf(const LocalIndexType localScvfIdx) const
    { return localScvfs_[localScvfIdx]; }

    //! returns the local scv entity corresponding to a given local (in the iv) scv idx
    const LocalScvType& localScv(const LocalIndexType localScvIdx) const
    { return localScvs_[localScvIdx]; }

    //! returns the (local) index set of scvfs that have associated intermediate unknowns
    const LocalIndexSet& fluxScvfIndexSet() const
    { return fluxFaceIndexSet_; }

    //! returns the (local) index set of scvfs that are Dirichlet boundary scvfs
    const LocalIndexSet& dirichletScvfIndexSet() const
    { return dirichletFaceIndexSet_; }

    //! returns the (local) index set of scvfs that are interior Dirichlet boundary scvfs
    const LocalIndexSet& interiorDirichletScvfIndexSet() const
    { return interiorDirichletFaceIndexSet_; }

    //! returns the (local) index set of all scvfs that are on interior boundaries
    const LocalIndexSet& interiorBoundaryScvfIndexSet() const
    { return interiorBoundaryFaceIndexSet_; }

    //! returns the container storing the data on interior boundaries
    const std::vector<InteriorBoundaryData>& interiorBoundaryData() const
    { return interiorBoundaryData_; }

    //! returns a reference to the problem to be solved
    const Problem& problem() const
    { return *problemPtr_; }

    //! returns a reference to the fvGeometry object
    const FVElementGeometry& fvGeometry() const
    { return *fvGeometryPtr_; }

    //! returns a reference to the element volume variables
    const ElementVolumeVariables& elemVolVars() const
    { return *elemVolVarsPtr_; }

private:
    //! returns a boolean whether or not the AB matrix has to be passed to the handles
    bool requireABMatrix_() const
    {
        static const bool requireAB = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddVelocity) || dim < dimWorld;
        return requireAB;
    }

    //! clears all the containers
    void reset_()
    {
        localElements_.clear();
        localScvs_.clear();
        localScvfs_.clear();
        globalLocalScvfPairedData_.clear();
        globalScvfIndices_.clear();
        outsideScvfIndices_.clear();
        dirichletFaceIndexSet_.clear();
        interiorDirichletFaceIndexSet_.clear();
        interiorBoundaryFaceIndexSet_.clear();
        fluxFaceIndexSet_.clear();
    }

    //! sets up the local scope for a given seed
    //! the local scv and scvf entities are constructed etc...
    void setLocalScope_(const Seed& seed, DataHandle& dataHandle)
    {
        //! clear previous data
        reset_();

        //! set sizes that can be determined already
        numFaces_ = seed.scvfSeeds().size();
        const auto numLocalScvs = seed.scvSeeds().size();
        const auto numGlobalScvfs = seed.globalScvfIndices().size();
        const auto maxNumVolVars = numLocalScvs + numFaces_;

        //! reserve memory for local entities
        localElements_.reserve(numLocalScvs);
        localScvs_.reserve(numLocalScvs);
        localScvfs_.reserve(numFaces_);
        globalLocalScvfPairedData_.reserve(numGlobalScvfs);
        globalScvfIndices_.reserve(numGlobalScvfs);

        // store outside scvf idx set on surface grids
        if (dim < dimWorld)
            outsideScvfIndices_.reserve(numGlobalScvfs);

        //! prepare interior boundary datadata handle
        interiorBoundaryData_.reserve(numFaces_);

        // prepare the stencil and the positions of the volvars
        volVarsPositions_.reserve(maxNumVolVars);
        volVarsStencil_ = seed.globalScvIndices();
        volVarsStencil_.reserve(maxNumVolVars);

        //! prepare temporary index sets
        dirichletFaceIndexSet_.reserve(numFaces_);
        interiorDirichletFaceIndexSet_.reserve(numFaces_);
        interiorBoundaryFaceIndexSet_.reserve(numFaces_);
        fluxFaceIndexSet_.reserve(numFaces_);

        // set up quantities related to sub-control volumes
        for (const auto& scvSeed : seed.scvSeeds())
        {
            auto element = problem().model().globalFvGeometry().element(scvSeed.globalIndex());
            localScvs_.emplace_back(problem(), element, fvGeometry(), scvSeed);
            localElements_.emplace_back(std::move(element));
            volVarsPositions_.push_back(localScvs_.back().center());
        }

        // set up quantitites related to sub-control volume faces
        LocalIndexType localFaceIdx = 0;
        for (const auto& scvfSeed : seed.scvfSeeds())
        {
            const auto faceType = scvfSeed.faceType();

            // we have to use the "inside" scv face here
            const auto& scvf = fvGeometry().scvf(scvfSeed.insideGlobalScvfIndex());
            localScvfs_.emplace_back(problem(), localElements_[scvfSeed.insideLocalScvIndex()], scvfSeed, scvf);

            // create global/local face data for this face
            // we simultaneously store the corresponding global scvf indices (allows global to local mapping later)
            globalLocalScvfPairedData_.emplace_back(&scvf, LocalFaceData(localFaceIdx, scvfSeed.insideLocalScvIndex(), false));
            globalScvfIndices_.push_back(scvf.index());

            // set data depending on the face type
            // obtain the local scvf entity just inserted
            const auto& curLocalScvf = localScvfs_.back();

            // interior faces are flux faces
            // also, set up outside global/local data for all the neighbors
            if (faceType == MpfaFaceTypes::interior)
            {
                for (unsigned int i = 0; i < curLocalScvf.outsideGlobalScvfIndices().size(); ++i)
                {
                    const auto& outsideScvf = fvGeometry().scvf(curLocalScvf.outsideGlobalScvfIndex(i));
                    globalLocalScvfPairedData_.emplace_back(&outsideScvf, LocalFaceData(localFaceIdx, scvfSeed.outsideLocalScvIndex(i), true));
                    globalScvfIndices_.push_back(outsideScvf.index());
                    if (dim < dimWorld)
                        outsideScvfIndices_.push_back(outsideScvf.index());
                }
                fluxFaceIndexSet_.push_back(localFaceIdx++);
            }
            // dirichlet faces are in the "stencil"
            else if (faceType == MpfaFaceTypes::dirichlet)
            {
                volVarsStencil_.push_back(curLocalScvf.outsideGlobalScvIndex());
                volVarsPositions_.push_back(curLocalScvf.ip());
                dirichletFaceIndexSet_.push_back(localFaceIdx++);
            }
            // neumann faces have an unknown associated with them
            else if (faceType == MpfaFaceTypes::neumann)
            {
                fluxFaceIndexSet_.push_back(localFaceIdx++);
            }
            // interior neumann faces additionally produce interior boundary data
            else if (faceType == MpfaFaceTypes::interiorNeumann)
            {
                interiorBoundaryData_.push_back(InteriorBoundaryData(problem(),
                                                                     curLocalScvf.insideGlobalScvIndex(),
                                                                     curLocalScvf.insideGlobalScvfIndex(),
                                                                     fluxFaceIndexSet_.size(),
                                                                     faceType));
                fluxFaceIndexSet_.push_back(localFaceIdx);
                interiorBoundaryFaceIndexSet_.push_back(localFaceIdx++);
            }
            // as well as interior Dirichlet boundaries
            else if (faceType == MpfaFaceTypes::interiorDirichlet)
            {
                interiorBoundaryData_.push_back(InteriorBoundaryData(problem(),
                                                                     curLocalScvf.insideGlobalScvIndex(),
                                                                     curLocalScvf.insideGlobalScvfIndex(),
                                                                     interiorDirichletFaceIndexSet_.size(),
                                                                     faceType));
                interiorDirichletFaceIndexSet_.push_back(localFaceIdx);
                interiorBoundaryFaceIndexSet_.push_back(localFaceIdx++);
            }
            else
                DUNE_THROW(Dune::InvalidStateException, "Face type can not be handled by the mpfa o-method.");
        }

        // set local sizes
        numFluxFaces_ = fluxFaceIndexSet_.size();
        numPotentials_ = volVarsStencil_.size() + interiorDirichletFaceIndexSet_.size();

        // resize the matrices in the data handle
        dataHandle.resizeT(numFaces_, numPotentials_);
        dataHandle.resizeCA(numFaces_, numFluxFaces_);
        if (requireABMatrix_())
            dataHandle.resizeAB(numFluxFaces_, numPotentials_);

        // resize the local matrices
        A_.resize(numFluxFaces_, numFluxFaces_);
        B_.resize(numFluxFaces_, numPotentials_);
        C_.resize(numFaces_, numFluxFaces_);
        D_.resize(numFaces_, numPotentials_);

        // on surface grids, resize the vector containing the "outside" transmissibilities
        if (dim < dimWorld)
            dataHandle.resizeOutsideTij(outsideScvfIndices_.size(), numPotentials_);
    }

    //! Assembles the local matrices that define the local system of equations and flux expressions
    template<typename GetTensorFunction>
    void assembleLocalMatrices_(const GetTensorFunction& getTensor)
    {
        // Xi factor for the coupling on interior neumann boundary facets
        static const auto xi = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Mpfa, Xi);

        // reset matrices
        A_ = 0.0;
        B_ = 0.0;
        C_ = 0.0;
        D_ = 0.0;

        // get no of Dirichlet faces
        const auto numLocalScvs = localScvs_.size();
        const auto numDirichletScvfs = dirichletScvfIndexSet().size();

        // reserve space for the omegas
        wijk_.resize(localScvfs_.size());

        // loop over the local faces
        for (unsigned int rowIdx = 0; rowIdx < numFaces_; ++rowIdx)
        {
            const auto& curLocalScvf = localScvf(rowIdx);
            const auto faceType = curLocalScvf.faceType();
            const bool curIsOnBoundary = curLocalScvf.boundary();
            const bool hasUnknown = faceType != MpfaFaceTypes::dirichlet && faceType != MpfaFaceTypes::interiorDirichlet;
            const auto curIdxInFluxFaces = hasUnknown ? this->findIndexInVector(fluxScvfIndexSet(), rowIdx) : -1;

            // get diffusion tensor in "positive" sub volume
            const auto posLocalScvIdx = curLocalScvf.insideLocalScvIndex();
            const auto& posLocalScv = localScv(posLocalScvIdx);
            const auto& posGlobalScv = fvGeometry().scv(posLocalScv.globalIndex());
            const auto& posVolVars = elemVolVars()[posGlobalScv];
            const auto& element = localElement(posLocalScvIdx);
            const auto tensor = getTensor(problem(), element, posVolVars, fvGeometry(), posGlobalScv);

            // the omega factors of the "positive" sub volume
            auto posWijk = calculateOmegas_(posLocalScv, curLocalScvf.unitOuterNormal(), curLocalScvf.area(), tensor);
            posWijk *= posVolVars.extrusionFactor();

            // Check the local directions of the positive sub volume
            for (int localDir = 0; localDir < dim; localDir++)
            {
                const auto otherLocalScvfIdx = posLocalScv.localScvfIndex(localDir);
                const auto& otherLocalScvf = localScvf(otherLocalScvfIdx);
                const auto otherFaceType = otherLocalScvf.faceType();

                // First, add the entries associated with unknown face pressures
                if (otherFaceType != MpfaFaceTypes::dirichlet && otherFaceType != MpfaFaceTypes::interiorDirichlet)
                {
                    // we need the index of the current local scvf in the flux face indices
                    auto otherIdxInFluxFaces = this->findIndexInVector(fluxScvfIndexSet(), otherLocalScvfIdx);

                    // this creates an entry in matrix C
                    C_[rowIdx][otherIdxInFluxFaces] -= posWijk[localDir];

                    // proceed depending on if the current face has an unknown
                    if (hasUnknown)
                    {
                        if (faceType == MpfaFaceTypes::interiorNeumann)
                        {
                            // on interior neumann faces, apply xi factor. However, if this is a boundary face at the same time, don't!
                            A_[curIdxInFluxFaces][otherIdxInFluxFaces] -= curIsOnBoundary ? posWijk[localDir] : xi*posWijk[localDir];

                            // maybe add terms stemming from the interior boundary (in case of facet coupling or membrane modelling)
                            if (otherIdxInFluxFaces == curIdxInFluxFaces)
                            {
                                const auto& data = interiorBoundaryData_[this->findIndexInVector(interiorBoundaryScvfIndexSet(), otherLocalScvfIdx)];
                                A_[curIdxInFluxFaces][otherIdxInFluxFaces] += asImp_().interiorNeumannTerm(getTensor, element, otherLocalScvf, data);
                            }
                        }
                        // this means we are on an interior face
                        else
                            A_[curIdxInFluxFaces][otherIdxInFluxFaces] -= posWijk[localDir];
                    }
                }
                else if (otherFaceType == MpfaFaceTypes::dirichlet)
                {
                    // the current face is a Dirichlet face and creates entries in D & eventually B
                    auto otherIdxInDiriFaces = this->findIndexInVector(dirichletScvfIndexSet(), otherLocalScvfIdx);

                    D_[rowIdx][numLocalScvs + otherIdxInDiriFaces] -= posWijk[localDir];
                    if (hasUnknown)
                    {
                        if (faceType == MpfaFaceTypes::interiorNeumann && !curIsOnBoundary)
                            B_[curIdxInFluxFaces][numLocalScvs + otherIdxInDiriFaces] += xi*posWijk[localDir];
                        else
                            B_[curIdxInFluxFaces][numLocalScvs + otherIdxInDiriFaces] += posWijk[localDir];
                    }
                }
                else if (otherFaceType == MpfaFaceTypes::interiorDirichlet)
                {
                    // the current face is an interior Dirichlet face and creates entries in D & eventually B
                    auto otherIdxInInteriorDiriFaces = this->findIndexInVector(interiorDirichletScvfIndexSet(), otherLocalScvfIdx);

                    D_[rowIdx][numLocalScvs + numDirichletScvfs + otherIdxInInteriorDiriFaces] -= posWijk[localDir];
                    if (hasUnknown)
                    {
                        if (faceType == MpfaFaceTypes::interiorNeumann && !curIsOnBoundary)
                            B_[curIdxInFluxFaces][numLocalScvs + numDirichletScvfs + otherIdxInInteriorDiriFaces] += xi*posWijk[localDir];
                        else
                            B_[curIdxInFluxFaces][numLocalScvs + numDirichletScvfs + otherIdxInInteriorDiriFaces] += posWijk[localDir];
                    }
                }

                // add entries related to pressures at the scv centers (dofs)
                D_[rowIdx][posLocalScvIdx] += posWijk[localDir];

                if (hasUnknown)
                {
                    if (faceType == MpfaFaceTypes::interiorNeumann && !curIsOnBoundary)
                        B_[curIdxInFluxFaces][posLocalScvIdx] -= xi*posWijk[localDir];
                    else
                        B_[curIdxInFluxFaces][posLocalScvIdx] -= posWijk[localDir];
                }
            }

            // store the omegas
            wijk_[rowIdx].emplace_back(std::move(posWijk));

            // If we are on an interior or interior neumann face, add values from negative sub volume
            if (!curIsOnBoundary && hasUnknown)
            {
                // loop over all the outside neighbors of this face and add entries
                unsigned int indexInOutsideData = 0;
                for (auto negLocalScvIdx : curLocalScvf.outsideLocalScvIndices())
                {
                    const auto& negLocalScv = localScv(negLocalScvIdx);
                    const auto& negGlobalScv = fvGeometry().scv(negLocalScv.globalIndex());
                    const auto& negVolVars = elemVolVars()[negGlobalScv];
                    const auto& negElement = localElement(negLocalScvIdx);
                    const auto negTensor = getTensor(problem(), negElement, negVolVars, fvGeometry(), negGlobalScv);

                    // the omega factors of the "negative" sub volume
                    DimVector negWijk;

                    // if dim < dimWorld, use outside normal vector
                    if (dim < dimWorld)
                    {
                        // outside scvf
                        const auto& outsideScvf = fvGeometry().scvf(curLocalScvf.outsideGlobalScvfIndex(indexInOutsideData));
                        auto negNormal = outsideScvf.unitOuterNormal();
                        negNormal *= -1.0;
                        negWijk = calculateOmegas_(negLocalScv, negNormal, curLocalScvf.area(), negTensor);
                    }
                    else
                        negWijk = calculateOmegas_(negLocalScv, curLocalScvf.unitOuterNormal(), curLocalScvf.area(), negTensor);

                    // scale by extrusion factpr
                    negWijk *= negVolVars.extrusionFactor();

                    // Check local directions of negative sub volume
                    for (int localDir = 0; localDir < dim; localDir++)
                    {
                        const auto otherLocalScvfIdx = negLocalScv.localScvfIndex(localDir);
                        const auto& otherLocalScvf = localScvf(otherLocalScvfIdx);
                        const auto otherFaceType = otherLocalScvf.faceType();

                        if (otherFaceType != MpfaFaceTypes::dirichlet && otherFaceType != MpfaFaceTypes::interiorDirichlet)
                        {
                            if (faceType == MpfaFaceTypes::interiorNeumann)
                                A_[curIdxInFluxFaces][this->findIndexInVector(fluxScvfIndexSet(), otherLocalScvfIdx)] -= (1-xi)*negWijk[localDir];
                            else
                                A_[curIdxInFluxFaces][this->findIndexInVector(fluxScvfIndexSet(), otherLocalScvfIdx)] += negWijk[localDir];
                        }
                        else if (otherFaceType == MpfaFaceTypes::interiorDirichlet)
                        {
                            const auto otherIdxInInteriorDiriFaces = this->findIndexInVector(interiorDirichletScvfIndexSet(), otherLocalScvfIdx);
                            if (faceType == MpfaFaceTypes::interiorNeumann)
                                B_[curIdxInFluxFaces][numLocalScvs + numDirichletScvfs + otherIdxInInteriorDiriFaces] += (1-xi)*negWijk[localDir];
                            else
                                B_[curIdxInFluxFaces][numLocalScvs + numDirichletScvfs + otherIdxInInteriorDiriFaces] -= negWijk[localDir];
                        }
                        else // Dirichlet face
                        {
                            if (faceType == MpfaFaceTypes::interiorNeumann)
                                B_[curIdxInFluxFaces][numLocalScvs + this->findIndexInVector(dirichletScvfIndexSet(), otherLocalScvfIdx)] += (1-xi)*negWijk[localDir];
                            else
                                B_[curIdxInFluxFaces][numLocalScvs + this->findIndexInVector(dirichletScvfIndexSet(), otherLocalScvfIdx)] -= negWijk[localDir];
                        }

                        // add entries to matrix B
                        if (faceType == MpfaFaceTypes::interiorNeumann)
                            B_[curIdxInFluxFaces][negLocalScvIdx] -= (1-xi)*negWijk[localDir];
                        else
                            B_[curIdxInFluxFaces][negLocalScvIdx] += negWijk[localDir];
                    }

                    // store the omegas
                    wijk_[rowIdx].emplace_back(std::move(negWijk));

                    // increment counter in outside data
                    indexInOutsideData++;
                }
            }
        }
    }

    //! for interaction volumes that have only dirichlet scvfs,
    //! the transmissibility matrix can be assembled directly
    template<typename GetTensorFunction>
    void assemblePureDirichletSystem_(const GetTensorFunction& getTensor, Matrix& T)
    {
        // reset the transmissibility matrix beforehand
        T = 0.0;

        const auto numLocalScvs = localScvs_.size();
        const auto numInteriorDirichletFaces = interiorDirichletScvfIndexSet().size();

        // Loop over all the faces, in this case these are all dirichlet boundaries
        for (unsigned int rowIdx = 0; rowIdx < numFaces_; ++rowIdx)
        {
            const auto& curLocalScvf = localScvf(rowIdx);

            // get diffusion tensor in "positive" sub volume
            const auto posLocalScvIdx = curLocalScvf.insideLocalScvIndex();
            const auto& posLocalScv = localScv(posLocalScvIdx);
            const auto& posGlobalScv = fvGeometry().scv(posLocalScv.globalIndex());
            const auto& posVolVars = elemVolVars()[posGlobalScv];
            const auto element = localElement(posLocalScvIdx);
            const auto tensor = getTensor(problem(), element, posVolVars, fvGeometry(), posGlobalScv);

            // the omega factors of the "positive" sub volume
            auto posWijk = calculateOmegas_(posLocalScv, curLocalScvf.unitOuterNormal(), curLocalScvf.area(), tensor);
            posWijk *= posVolVars.extrusionFactor();

            for (int localDir = 0; localDir < dim; localDir++)
            {
                // When interior boundaries are disabled, all faces will be of dirichlet type
                if (!enableInteriorBoundaries)
                {
                    const auto otherLocalScvfIdx = posLocalScv.localScvfIndex(localDir);
                    const auto otherIdxInDiriFaces = this->findIndexInVector(dirichletScvfIndexSet(), otherLocalScvfIdx);
                    T[rowIdx][numLocalScvs + otherIdxInDiriFaces] -= posWijk[localDir];
                    T[rowIdx][posLocalScvIdx] += posWijk[localDir];
                }
                else
                {
                    const auto otherLocalScvfIdx = posLocalScv.localScvfIndex(localDir);
                    const auto& curLocalScvf = localScvf(otherLocalScvfIdx);
                    const auto curFaceType = curLocalScvf.faceType();

                    const auto otherIdxInDiriFaces = curFaceType == MpfaFaceTypes::dirichlet ?
                                                   this->findIndexInVector(dirichletScvfIndexSet(), otherLocalScvfIdx) :
                                                   numInteriorDirichletFaces + this->findIndexInVector(interiorDirichletScvfIndexSet(), otherLocalScvfIdx);

                    T[rowIdx][numLocalScvs + otherIdxInDiriFaces] -= posWijk[localDir];
                    T[rowIdx][posLocalScvIdx] += posWijk[localDir];
                }
            }
        }
    }

    //! computes the transmissibilities associated with "outside" faces on surface grids
    void computeOutsideTransmissibilities_(DataHandle& dataHandle) const
    {
        if (!(dim < dimWorld))
            DUNE_THROW(Dune::InvalidStateException, "transmissibilities for outside scvfs can only be computed for surface grids!");

        // get the local scv and iterate over local coordinates
        const auto numLocalScvs = localScvs_.size();
        const auto numDirichletScvfs = dirichletScvfIndexSet().size();

        for (const auto& globalLocalData : globalLocalScvfPairedData())
        {
            const auto& localFaceData = globalLocalData.second;

            // continue only for "outside" faces
            if (!localFaceData.isOutside)
                continue;

            const auto& posLocalScv = localScv(localFaceData.localScvIndex);
            const auto& posLocalScvf = localScvf(localFaceData.localScvfIndex);

            const auto idxInOutside = this->findIndexInVector(posLocalScvf.outsideLocalScvIndices(), localFaceData.localScvIndex);
            const auto& wijk = wijk_[localFaceData.localScvfIndex][idxInOutside+1];

            // store the calculated transmissibilities in the data handle
            auto& tij = dataHandle.outsideTij()[this->findIndexInVector(outsideScvfIndices_, globalLocalData.first->index())];
            tij = 0.0;

            for (int localDir = 0; localDir < dim; localDir++)
            {
                const auto localScvfIdx = posLocalScv.localScvfIndex(localDir);
                const auto faceType = localScvf(localScvfIdx).faceType();
                if (faceType != MpfaFaceTypes::dirichlet && faceType != MpfaFaceTypes::interiorDirichlet)
                {
                    const auto fluxFaceIndex = this->findIndexInVector(fluxScvfIndexSet(), localScvfIdx);
                    auto tmp = dataHandle.AB()[fluxFaceIndex];
                    tmp *= wijk[localDir];

                    tij -= tmp;
                }
                else if (faceType == MpfaFaceTypes::dirichlet)
                {
                    const auto idxInDiriFaces = this->findIndexInVector(dirichletScvfIndexSet(), localScvfIdx);
                    tij[numLocalScvs + idxInDiriFaces] -= wijk[localDir];
                }
                else if (faceType == MpfaFaceTypes::interiorDirichlet)
                {
                    const auto idxInInteriorDiriFaces = this->findIndexInVector(interiorDirichletScvfIndexSet(), localScvfIdx);
                    tij[numLocalScvs + numDirichletScvfs + idxInInteriorDiriFaces] -= wijk[localDir];
                }

                // add entry from the scv unknown
                tij[localFaceData.localScvIndex] += wijk[localDir];
            }
        }
    }

    // calculates n_i^T*K_j*nu_k
    DimVector calculateOmegas_(const LocalScvType& localScv,
                               const GlobalPosition normal,
                               const Scalar area,
                               const Tensor& T) const
    {
        // make sure we have positive definite diffsion tensors
        assert(this->tensorIsPositiveDefinite(T) && "only positive definite tensors can be handled by mpfa methods");

        DimVector wijk;
        GlobalPosition tmp;
        for (int dir = 0; dir < dim; ++dir)
        {
            T.mv(localScv.innerNormal(dir), tmp);
            wijk[dir] = tmp*normal;
        }
        wijk *= area;
        wijk /= localScv.detX();

        return wijk;
    }

    // calculates n_i^T*K_j*nu_k
    DimVector calculateOmegas_(const LocalScvType& localScv,
                               const GlobalPosition normal,
                               const Scalar area,
                               const Scalar t) const
    {
        // make sure we have positive diffusion coefficients
        assert(t > 0.0 && "non-positive diffusion coefficients cannot be handled by mpfa methods");

        DimVector wijk;
        GlobalPosition tmp(normal);
        tmp *= t;

        for (int dir = 0; dir < dim; ++dir)
            wijk[dir] = tmp*localScv.innerNormal(dir);
        wijk *= area;
        wijk /= localScv.detX();

        return wijk;
    }

    Implementation& asImp_()
    { return static_cast<Implementation&> (*this); }

    const Implementation& asImp_() const
    { return static_cast<const Implementation&> (*this); }

    const Problem* problemPtr_;
    const FVElementGeometry* fvGeometryPtr_;
    const ElementVolumeVariables* elemVolVarsPtr_;

    // tells us if the IV touches any kind of boundary
    bool onDomainOrInteriorBoundary_;

    // Variables defining the local scope
    std::vector<Element> localElements_;
    std::vector<LocalScvType> localScvs_;
    std::vector<LocalScvfType> localScvfs_;
    std::vector<GlobalLocalFaceDataPair> globalLocalScvfPairedData_;
    GlobalIndexSet globalScvfIndices_;
    GlobalIndexSet outsideScvfIndices_;

    LocalIndexSet fluxFaceIndexSet_;
    LocalIndexSet dirichletFaceIndexSet_;
    LocalIndexSet interiorDirichletFaceIndexSet_;
    LocalIndexSet interiorBoundaryFaceIndexSet_;

    // the stencil and the vol var positions
    GlobalIndexSet volVarsStencil_;
    PositionVector volVarsPositions_;

    // container with data on interior boundaries
    std::vector<InteriorBoundaryData> interiorBoundaryData_;

    // sizes involved in the local matrices
    unsigned int numFaces_;
    unsigned int numFluxFaces_;
    unsigned int numPotentials_;

    // The omega factors and the matrix A⁻1*B are stored
    // in order to recover the transmissibilities of outside faces on network grids
    std::vector< std::vector< DimVector > > wijk_;
    Matrix CAinv_;

    // Matrices involved in transmissibility calculations
    Matrix A_;
    Matrix B_;
    Matrix C_;
    Matrix D_;

    // stores all the neumann fluxes appearing in this interaction volume
    std::vector<PrimaryVariables> neumannFluxes_;
};

} // end namespace

#endif
