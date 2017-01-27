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
template<class TypeTag, class Traits> class CCMpfaOInteractionVolume;

/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volumes of the mpfa-o method.
 *        We introduce one more level of inheritance here because the o-method with
 *        full pressure support uses the mpfa-o interaction volume but uses a different
 *        traits class.
 */
template<class TypeTag>
class CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethod> : public CCMpfaOInteractionVolume<TypeTag, CCMpfaOInteractionVolumeTraits<TypeTag>>
{
    using Traits = CCMpfaOInteractionVolumeTraits<TypeTag>;
    using ParentType = CCMpfaOInteractionVolume<TypeTag, Traits>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using IVSeed = typename Traits::Seed;
public:
    CCMpfaInteractionVolumeImplementation(const IVSeed& seed,
                                          const Problem& problem,
                                          const FVElementGeometry& fvGeometry,
                                          const ElementVolumeVariables& elemVolVars)
    : ParentType(seed, problem, fvGeometry, elemVolVars)
    {}

};

template<class TypeTag, class Traits>
class CCMpfaOInteractionVolume : public CCMpfaInteractionVolumeBase<TypeTag, Traits>
{
    // The interaction volume implementation has to be friend,
    // because some methods use the mpfa-o interaction volume as base
    friend typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using Implementation = typename GET_PROP_TYPE(TypeTag, InteractionVolume);
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
    static constexpr bool facetCoupling = GET_PROP_VALUE(TypeTag, MpfaFacetCoupling);
    static constexpr bool enableInteriorBoundaries = GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using DynamicVector = typename Traits::Vector;
    using DynamicMatrix = typename Traits::Matrix;
    using Tensor = typename Traits::Tensor;

    using LocalScvType = typename Traits::LocalScvType;
    using LocalScvfType = typename Traits::LocalScvfType;

public:
    using typename ParentType::LocalIndexType;
    using typename ParentType::LocalIndexSet;
    using typename ParentType::LocalFaceData;
    using typename ParentType::GlobalIndexSet;
    using typename ParentType::PositionVector;
    using typename ParentType::Seed;

    CCMpfaOInteractionVolume(const Seed& seed,
                             const Problem& problem,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars)
    : problemPtr_(&problem),
      fvGeometryPtr_(&fvGeometry),
      elemVolVarsPtr_(&elemVolVars),
      onDomainOrInteriorBoundary_(seed.onDomainOrInteriorBoundary()),
      globalScvfIndices_(seed.globalScvfIndices())
    {
        // create local sub control entities from the seed
        createLocalEntities_(seed);

        // reserve memory
        auto numLocalScvs = localScvs_.size();
        auto numLocalScvfs = localScvfs_.size();
        auto maxNumVolVars = numLocalScvs + numLocalScvfs;
        volVarsStencil_ = seed.globalScvIndices(); // boundary vol vars are placed at the end
        volVarsStencil_.reserve(maxNumVolVars);
        volVarsPositions_.reserve(maxNumVolVars);
        dirichletFaceIndexSet_.reserve(numLocalScvfs);
        interiorDirichletFaceIndexSet_.reserve(numLocalScvfs);
        interiorBoundaryFaceIndexSet_.reserve(numLocalScvfs);
        interiorBoundaryData_.reserve(numLocalScvfs);
        fluxFaceIndexSet_.reserve(numLocalScvfs);

        // the positions where the vol vars are defined at (required for the gravitational acceleration)
        for (auto&& localScv : localScvs_)
            volVarsPositions_.push_back(localScv.center());

        // maybe add dirichlet vol var indices and set up local index sets of flux and dirichlet faces
        LocalIndexType localScvfIdx = 0;
        for (auto&& localScvf : localScvfs_)
        {
            auto faceType = localScvf.faceType();

            // eventually add vol var index and corresponding position
            if (faceType == MpfaFaceTypes::dirichlet)
            {
                volVarsStencil_.push_back(localScvf.outsideGlobalScvIndex());
                volVarsPositions_.push_back(localScvf.ip());
                dirichletFaceIndexSet_.push_back(localScvfIdx++);
            }
            else if (faceType == MpfaFaceTypes::interior || faceType == MpfaFaceTypes::neumann)
            {
                fluxFaceIndexSet_.push_back(localScvfIdx++);
            }
            else if (faceType == MpfaFaceTypes::interiorNeumann)
            {
                interiorBoundaryData_.push_back(InteriorBoundaryData(problem,
                                                                     localScvf.insideGlobalScvIndex(),
                                                                     localScvf.insideGlobalScvfIndex(),
                                                                     fluxFaceIndexSet_.size(),
                                                                     faceType));
                fluxFaceIndexSet_.push_back(localScvfIdx);
                interiorBoundaryFaceIndexSet_.push_back(localScvfIdx++);
            }
            else if(faceType == MpfaFaceTypes::interiorDirichlet)
            {
                interiorBoundaryData_.push_back(InteriorBoundaryData(problem,
                                                                     localScvf.insideGlobalScvIndex(),
                                                                     localScvf.insideGlobalScvfIndex(),
                                                                     interiorDirichletFaceIndexSet_.size(),
                                                                     faceType));
                interiorDirichletFaceIndexSet_.push_back(localScvfIdx);
                interiorBoundaryFaceIndexSet_.push_back(localScvfIdx++);
            }
            else
                DUNE_THROW(Dune::InvalidStateException, "Unknown mpfa face type.");
        }

        // initialize the vector containing the neumann fluxes
        assembleNeumannFluxVector_();
    }

    template<typename GetTensorFunction>
    void solveLocalSystem(const GetTensorFunction& getTensor)
    {
        const auto numFluxFaces = fluxScvfIndexSet_().size();

        // if only dirichlet faces are present, assemble T_ directly
        if (numFluxFaces == 0)
            return assemblePureDirichletSystem_(getTensor);

        const auto numFaces = localScvfs_.size();
        const auto numPotentials = volVarsStencil().size() + interiorDirichletScvfIndexSet_().size();

        // the local matrices
        DynamicMatrix A(numFluxFaces, numFluxFaces, 0.0);
        DynamicMatrix B(numFluxFaces, numPotentials, 0.0);
        DynamicMatrix C(numFaces, numFluxFaces, 0.0);
        DynamicMatrix D(numFaces, numPotentials, 0.0);

        assembleLocalMatrices_(getTensor, A, B, C, D);

        // solve local system and store matrices
        DynamicMatrix copy(B);
        A.invert();
        AinvB_ = B.leftmultiply(A);
        CAinv_ = C.rightmultiply(A);
        T_ = multiplyMatrices(CAinv_, copy);
        T_ += D;
    }

    //! gets local data on a global scvf within the interaction volume
    //! specialization for dim = dimWorld
    template<int d = dim, int dw = dimWorld>
    typename std::enable_if< (d == dw), LocalFaceData >::type
    getLocalFaceData(const SubControlVolumeFace& scvf) const
    {
        auto scvfGlobalIdx = scvf.index();

        // if interior boundaries are disabled, do simplified search
        if (!enableInteriorBoundaries)
        {
            LocalIndexType localIdx = 0;
            for (auto&& localScvf : localScvfs_)
            {
                if (localScvf.insideGlobalScvfIndex() == scvfGlobalIdx)
                    return LocalFaceData(localIdx, localScvf.insideLocalScvIndex(), false);

                if (!localScvf.boundary() && localScvf.outsideGlobalScvfIndex() == scvfGlobalIdx)
                    return LocalFaceData(localIdx, localScvf.outsideLocalScvIndex(), true);

                localIdx++;
            }
        }
        else
        {
            // first check the inside data of the faces
            LocalIndexType localIdx = 0;
            for (auto&& localScvf : localScvfs_)
            {
                if (localScvf.insideGlobalScvfIndex() == scvfGlobalIdx)
                    return LocalFaceData(localIdx, localScvf.insideLocalScvIndex(), false);
                localIdx++;
            }

            // then check the outside data of the faces
            localIdx = 0;
            for (auto&& localScvf : localScvfs_)
            {
                if (!localScvf.boundary() && localScvf.outsideGlobalScvfIndex() == scvfGlobalIdx)
                    return LocalFaceData(localIdx, localScvf.outsideLocalScvIndex(), true);
                localIdx++;
            }
        }

        DUNE_THROW(Dune::InvalidStateException, "Could not find the local scv face in the interaction volume for the given scvf with index: " << scvf.index());
    }

    //! gets local data on a global scvf within the interaction volume
    //! specialization for dim < dimWorld
    template<int d = dim, int dw = dimWorld>
    typename std::enable_if< (d < dw), LocalFaceData >::type
    getLocalFaceData(const SubControlVolumeFace& scvf) const
    {
        const auto scvfGlobalIdx = scvf.index();

        // if interior boundaries are disabled, do simplified search
        if (!enableInteriorBoundaries)
        {
            LocalIndexType localFaceIdx = 0;
            for (auto&& localScvf : localScvfs_)
            {
                if (localScvf.insideGlobalScvfIndex() == scvfGlobalIdx)
                    return LocalFaceData(localFaceIdx, localScvf.insideLocalScvIndex(), false);

                if (!localScvf.boundary() && MpfaHelper::contains(localScvf.outsideGlobalScvfIndices(), scvfGlobalIdx))
                {
                    if (localScvf.outsideLocalScvIndices().size() == 1)
                        return LocalFaceData(localFaceIdx, localScvf.outsideLocalScvIndex(), true);

                    // Now we have to find wich local scv is the one the global scvf is embedded in
                    const auto scvGlobalIdx = scvf.insideScvIdx();
                    for (auto localScvIdx : localScvf.outsideLocalScvIndices())
                        if (localScv_(localScvIdx).globalIndex() == scvGlobalIdx)
                            return LocalFaceData(localFaceIdx, localScvIdx, true);

                    DUNE_THROW(Dune::InvalidStateException, "Could not find the local scv in the interaction volume for the given scvf with index: " << scvfGlobalIdx);
                }

                localFaceIdx++;
            }
        }
        else
        {
            // first check the inside data of the faces
            LocalIndexType localFaceIdx = 0;
            for (auto&& localScvf : localScvfs_)
            {
                if (localScvf.insideGlobalScvfIndex() == scvfGlobalIdx)
                    return LocalFaceData(localFaceIdx, localScvf.insideLocalScvIndex(), false);
                localFaceIdx++;
            }

            // then check the outside data
            localFaceIdx = 0;
            for (auto&& localScvf : localScvfs_)
            {
                if (!localScvf.boundary() && MpfaHelper::contains(localScvf.outsideGlobalScvfIndices(), scvfGlobalIdx))
                {
                    if (localScvf.outsideLocalScvIndices().size() == 1)
                        return LocalFaceData(localFaceIdx, localScvf.outsideLocalScvIndex(), true);

                    // Now we have to find wich local scv is the one the global scvf is embedded in
                    const auto scvGlobalIdx = scvf.insideScvIdx();
                    for (auto localScvIdx : localScvf.outsideLocalScvIndices())
                        if (localScv_(localScvIdx).globalIndex() == scvGlobalIdx)
                            return LocalFaceData(localFaceIdx, localScvIdx, true);

                    DUNE_THROW(Dune::InvalidStateException, "Could not find the local scv in the interaction volume for the given scvf with index: " << scvfGlobalIdx);
                }

                localFaceIdx++;
            }
        }

        DUNE_THROW(Dune::InvalidStateException, "Could not find the local scv face in the interaction volume for the given scvf with index: " << scvfGlobalIdx);
    }

    //! Gets the transmissibilities for a sub-control volume face within the interaction volume.
    //! specialization for dim == dimWorld
    template<int d = dim, int dw = dimWorld>
    typename std::enable_if< (d == dw), DynamicVector >::type
    getTransmissibilities(const LocalFaceData& localFaceData) const
    {
        if (localFaceData.isOutside)
        {
            auto tij = T_[localFaceData.localScvfIndex];
            tij *= -1.0;
            return tij;
        }
        else
            return T_[localFaceData.localScvfIndex];
    }

    //! Gets the transmissibilities for a sub-control volume face within the interaction volume.
    //! specialization for dim < dimWorld.
    template<int d = dim, int dw = dimWorld>
    typename std::enable_if< (d < dw), DynamicVector >::type
    getTransmissibilities(const LocalFaceData& localFaceData) const
    {
        // If we come from the inside, simply return tij
        if (!localFaceData.isOutside)
            return T_[localFaceData.localScvfIndex];

        // compute the outside transmissibilities
        DynamicVector tij(volVarsStencil().size() + interiorDirichletScvfIndexSet_().size(), 0.0);

        // get the local scv and iterate over local coordinates
        const auto numLocalScvs = localScvs_.size();
        const auto numDirichletScvfs = dirichletScvfIndexSet_().size();
        const auto& localScv = localScv_(localFaceData.localScvIndex);
        const auto& localScvf = localScvf_(localFaceData.localScvfIndex);

        const auto idxInOutside = this->findIndexInVector(localScvf.outsideLocalScvIndices(), localFaceData.localScvIndex);
        const auto& wijk = wijk_[localFaceData.localScvfIndex][idxInOutside+1];
        for (int localDir = 0; localDir < dim; localDir++)
        {
            const auto localScvfIdx = localScv.localScvfIndex(localDir);
            const auto& localScvf = localScvf_(localScvfIdx);

            const auto faceType = localScvf.faceType();
            if (faceType != MpfaFaceTypes::dirichlet && faceType != MpfaFaceTypes::interiorDirichlet)
            {
                const auto fluxFaceIndex = this->findIndexInVector(fluxScvfIndexSet_(), localScvfIdx);
                auto tmp = AinvB_[fluxFaceIndex];
                tmp *= wijk[localDir];

                tij += tmp;
            }
            else if (faceType == MpfaFaceTypes::dirichlet)
            {
                const auto idxInDiriFaces = this->findIndexInVector(dirichletScvfIndexSet_(), localScvfIdx);
                tij[numLocalScvs + idxInDiriFaces] += wijk[localDir];
            }
            else if (faceType == MpfaFaceTypes::interiorDirichlet)
            {
                const auto idxInInteriorDiriFaces = this->findIndexInVector(interiorDirichletScvfIndexSet_(), localScvfIdx);
                tij[numLocalScvs + numDirichletScvfs + idxInInteriorDiriFaces] += wijk[localDir];
            }

            // add entry from the scv unknown
            tij[localFaceData.localScvIndex] -= wijk[localDir];
        }

        return tij;
    }

    //! Returns the vector of coefficients with which the vector of neumann boundary conditions
    //! has to be multiplied in order to transform them on the scvf this cache belongs to
    DynamicVector getNeumannFluxTransformationCoefficients(const LocalFaceData& localFaceData) const
    {
        auto cij = CAinv_[localFaceData.localScvfIndex];
        if (localFaceData.isOutside)
            cij *= -1.0;
        return cij;
    }

    Scalar getNeumannFlux(const LocalFaceData& localFaceData, unsigned int eqIdx) const
    {
        if (!onDomainOrInteriorBoundary() || useTpfaBoundary || fluxScvfIndexSet_().size() == 0 )
            return 0.0;

        // Do the scalar product CAinv_*neumannFluxes[eqIdx]
        assert(CAinv_[localFaceData.localScvfIndex].size() == neumannFluxes_.size() &&
               "Number of columns of matrix does not correspond to number entries in vector!");

        Scalar flux(0.0);
        for (unsigned int i = 0; i < neumannFluxes_.size(); ++i)
            flux += CAinv_[localFaceData.localScvfIndex][i] * neumannFluxes_[i][eqIdx];

        // flip sign if we are coming from the outside
        if (localFaceData.isOutside)
            return -1.0*flux;

        return flux;
    }

    bool onDomainOrInteriorBoundary() const
    { return onDomainOrInteriorBoundary_; }

    const GlobalIndexSet& volVarsStencil() const
    { return volVarsStencil_; }

    const PositionVector& volVarsPositions() const
    { return volVarsPositions_; }

    const GlobalIndexSet& globalScvfs() const
    { return globalScvfIndices_; }

    const std::vector<InteriorBoundaryData>& interiorBoundaryData() const
    { return interiorBoundaryData_; }

private:

    const LocalScvfType& localScvf_(const LocalIndexType localScvfIdx) const
    { return localScvfs_[localScvfIdx]; }

    const LocalScvType& localScv_(const LocalIndexType localScvIdx) const
    { return localScvs_[localScvIdx]; }

    const LocalIndexSet& fluxScvfIndexSet_() const
    { return fluxFaceIndexSet_; }

    const LocalIndexSet& dirichletScvfIndexSet_() const
    { return dirichletFaceIndexSet_; }

    const LocalIndexSet& interiorDirichletScvfIndexSet_() const
    { return interiorDirichletFaceIndexSet_; }

    const LocalIndexSet& interiorBoundaryScvfIndexSet_() const
    { return interiorBoundaryFaceIndexSet_; }

    const Element& localElement_(const LocalIndexType localScvIdx) const
    { return localElements_[localScvIdx]; }

    void createLocalEntities_(const Seed& seed)
    {
        const auto numScvs = seed.scvSeeds().size();
        const auto numScvfs = seed.scvfSeeds().size();

        localElements_.reserve(numScvs);
        localScvs_.reserve(numScvs);
        localScvfs_.reserve(numScvfs);

        for (auto&& scvSeed : seed.scvSeeds())
        {
            auto element = problem_().model().globalFvGeometry().element(scvSeed.globalIndex());
            localScvs_.emplace_back(LocalScvType(problem_(), element, fvGeometry_(), scvSeed));
            localElements_.emplace_back(std::move(element));
        }

        for (auto&& scvfSeed : seed.scvfSeeds())
        {
            // we have to use the "inside" scv face here
            auto&& scvf = fvGeometry_().scvf(scvfSeed.insideGlobalScvfIndex());
            localScvfs_.emplace_back(LocalScvfType(scvfSeed, scvf));
        }
    }

    void assembleNeumannFluxVector_()
    {
        // initialize the neumann fluxes vector to zero
        neumannFluxes_.resize(fluxFaceIndexSet_.size(), PrimaryVariables(0.0));

        if (!onDomainOrInteriorBoundary() || useTpfaBoundary)
            return;

        LocalIndexType fluxFaceIdx = 0;
        for (auto localFluxFaceIdx : fluxFaceIndexSet_)
        {
            const auto& localScvf = localScvf_(localFluxFaceIdx);
            const auto faceType = localScvf.faceType();

            if (faceType == MpfaFaceTypes::neumann || (faceType == MpfaFaceTypes::interiorNeumann && !facetCoupling))
            {
                const auto& element = localElement_(localScvf.insideLocalScvIndex());
                const auto& globalScvf = fvGeometry_().scvf(localScvf.insideGlobalScvfIndex());
                auto neumannFlux = problem_().neumann(element, this->fvGeometry_(), this->elemVolVars_(), globalScvf);
                neumannFlux *= globalScvf.area();
                neumannFlux *= elemVolVars_()[globalScvf.insideScvIdx()].extrusionFactor();

                // The flux is assumed to be prescribed in the form of -D*gradU
                neumannFluxes_[fluxFaceIdx] = neumannFlux;
            }

            fluxFaceIdx++;
        }
    }

    template<typename GetTensorFunction>
    void assembleLocalMatrices_(const GetTensorFunction& getTensor,
                                DynamicMatrix& A,
                                DynamicMatrix& B,
                                DynamicMatrix& C,
                                DynamicMatrix& D)
    {
        static const auto xi = GET_PARAM_FROM_GROUP(TypeTag, Scalar, Mpfa, Xi);
        const auto numLocalScvs = localScvs_.size();
        const auto numDirichletScvfs = dirichletScvfIndexSet_().size();

        // reserve space for the omegas
        wijk_.resize(localScvfs_.size());

        // loop over the local faces
        unsigned int rowIdx = 0;
        for (auto&& localScvf : localScvfs_)
        {
            const auto faceType = localScvf.faceType();
            const bool hasUnknown = faceType != MpfaFaceTypes::dirichlet && faceType != MpfaFaceTypes::interiorDirichlet;
            const LocalIndexType idxInFluxFaces = hasUnknown ? this->findIndexInVector(fluxScvfIndexSet_(), rowIdx) : -1;

            // get diffusion tensor in "positive" sub volume
            const auto posLocalScvIdx = localScvf.insideLocalScvIndex();
            const auto& posLocalScv = localScv_(posLocalScvIdx);
            const auto& posGlobalScv = fvGeometry_().scv(posLocalScv.globalIndex());
            const auto& posVolVars = elemVolVars_()[posGlobalScv];
            const auto& element = localElement_(posLocalScvIdx);
            const auto tensor = getTensor(element, posVolVars, posGlobalScv);

            // the omega factors of the "positive" sub volume
            auto posWijk = calculateOmegas_(posLocalScv, localScvf.unitOuterNormal(), localScvf.area(), tensor);
            posWijk *= posVolVars.extrusionFactor();

            // Check the local directions of the positive sub volume
            for (int localDir = 0; localDir < dim; localDir++)
            {
                const auto curLocalScvfIdx = posLocalScv.localScvfIndex(localDir);
                const auto& curLocalScvf = localScvf_(curLocalScvfIdx);
                const auto curFaceType = curLocalScvf.faceType();

                // First, add the entries associated with unknown face pressures
                if (curFaceType != MpfaFaceTypes::dirichlet && curFaceType != MpfaFaceTypes::interiorDirichlet)
                {
                    // we need the index of the current local scvf in the flux face indices
                    auto curIdxInFluxFaces = this->findIndexInVector(fluxScvfIndexSet_(), curLocalScvfIdx);

                    // this creates an entry in matrix C
                    C[rowIdx][curIdxInFluxFaces] += posWijk[localDir];

                    // proceed depending on if the current face has an unknown
                    if (hasUnknown)
                    {
                        if (faceType == MpfaFaceTypes::interiorNeumann)
                        {
                            // on interior neumann faces, apply xi factor
                            // However, if this is a boundary face at the same time, don't!
                            A[idxInFluxFaces][curIdxInFluxFaces] += localScvf.globalScvf().boundary() ?
                                                                    posWijk[localDir] : xi*posWijk[localDir];

                            // add values from other domain in case of facet coupling
                            if (facetCoupling && curIdxInFluxFaces == idxInFluxFaces)
                            {
                                // get interior boundary data
                                const auto& data = interiorBoundaryData_[this->findIndexInVector(interiorBoundaryScvfIndexSet_(), curLocalScvfIdx)];

                                // get the volvars on the actual interior boundary face
                                const auto facetVolVars = data.facetVolVars(fvGeometry_());

                                // calculate "leakage factor"
                                const auto n = curLocalScvf.unitOuterNormal();
                                const auto v = [&] ()
                                        {
                                            auto res = n;
                                            res *= -0.5*facetVolVars.extrusionFactor();
                                            res -= curLocalScvf.ip();
                                            res += curLocalScvf.globalScvf().facetCorner();
                                            res /= res.two_norm2();
                                            return res;
                                        } ();

                                // substract from matrix entry
                                A[idxInFluxFaces][curIdxInFluxFaces] -= MpfaHelper::nT_M_v(n, posVolVars.permeability(), v);
                            }
                        }
                        // this means we are on an interior face
                        else
                            A[idxInFluxFaces][curIdxInFluxFaces] += posWijk[localDir];
                    }
                }
                else if (curFaceType == MpfaFaceTypes::dirichlet)
                {
                    // the current face is a Dirichlet face and creates entries in D & eventually B
                    auto curIdxInDiriFaces = this->findIndexInVector(dirichletScvfIndexSet_(), curLocalScvfIdx);

                    D[rowIdx][numLocalScvs + curIdxInDiriFaces] += posWijk[localDir];
                    if (hasUnknown)
                        B[idxInFluxFaces][numLocalScvs + curIdxInDiriFaces] -= posWijk[localDir];
                }
                else if (curFaceType == MpfaFaceTypes::interiorDirichlet)
                {
                    // the current face is an interior Dirichlet face and creates entries in D & eventually B
                    auto curIdxInInteriorDiriFaces = this->findIndexInVector(interiorDirichletScvfIndexSet_(), curLocalScvfIdx);

                    D[rowIdx][numLocalScvs + numDirichletScvfs + curIdxInInteriorDiriFaces] += posWijk[localDir];
                    if (hasUnknown)
                        B[idxInFluxFaces][numLocalScvs + numDirichletScvfs + curIdxInInteriorDiriFaces] -= posWijk[localDir];
                }

                // add entries related to pressures at the scv centers (dofs)
                D[rowIdx][posLocalScvIdx] -= posWijk[localDir];

                if (hasUnknown)
                {
                    if (faceType == MpfaFaceTypes::interiorNeumann && !localScvf.globalScvf().boundary())
                        B[idxInFluxFaces][posLocalScvIdx] += xi*posWijk[localDir];
                    else
                        B[idxInFluxFaces][posLocalScvIdx] += posWijk[localDir];
                }
            }

            // store the omegas
            wijk_[rowIdx].emplace_back(std::move(posWijk));

            // If we are on an interior or interior neumann face, add values from negative sub volume
            if (!localScvf.globalScvf().boundary() &&
                (faceType == MpfaFaceTypes::interior || faceType == MpfaFaceTypes::interiorNeumann))
            {
                // loop over all the outside neighbors of this face and add entries
                unsigned int indexInOutsideData = 0;
                for (auto negLocalScvIdx : localScvf.outsideLocalScvIndices())
                {
                    const auto& negLocalScv = localScv_(negLocalScvIdx);
                    const auto& negGlobalScv = fvGeometry_().scv(negLocalScv.globalIndex());
                    const auto& negVolVars = elemVolVars_()[negGlobalScv];
                    const auto& negElement = localElement_(negLocalScvIdx);
                    const auto negTensor = getTensor(negElement, negVolVars, negGlobalScv);

                    // the omega factors of the "negative" sub volume
                    DimVector negWijk;

                    // if dim < dimWorld, use outside normal vector
                    if (dim < dimWorld)
                    {
                        // outside scvf
                        const auto& outsideScvf = fvGeometry_().scvf(localScvf.outsideGlobalScvfIndex(indexInOutsideData));
                        auto negNormal = outsideScvf.unitOuterNormal();
                        negNormal *= -1.0;
                        negWijk = calculateOmegas_(negLocalScv, negNormal, localScvf.area(), negTensor);
                    }
                    else
                        negWijk = calculateOmegas_(negLocalScv, localScvf.unitOuterNormal(), localScvf.area(), negTensor);

                    // scale by extrusion factpr
                    negWijk *= negVolVars.extrusionFactor();

                    // Check local directions of negative sub volume
                    for (int localDir = 0; localDir < dim; localDir++)
                    {
                        const auto curLocalScvfIdx = negLocalScv.localScvfIndex(localDir);
                        const auto& curLocalScvf = localScvf_(curLocalScvfIdx);
                        const auto curFaceType = curLocalScvf.faceType();

                        if (curFaceType != MpfaFaceTypes::dirichlet && curFaceType != MpfaFaceTypes::interiorDirichlet)
                        {
                            if (faceType == MpfaFaceTypes::interiorNeumann)
                                A[idxInFluxFaces][this->findIndexInVector(fluxScvfIndexSet_(), curLocalScvfIdx)] -= (1-xi)*negWijk[localDir];
                            else
                                A[idxInFluxFaces][this->findIndexInVector(fluxScvfIndexSet_(), curLocalScvfIdx)] -= negWijk[localDir];
                        }
                        else if (curFaceType == MpfaFaceTypes::interiorDirichlet)
                        {
                            const auto idxInInteriorDiriFaces = this->findIndexInVector(interiorDirichletScvfIndexSet_(), curLocalScvfIdx);
                            B[idxInFluxFaces][numLocalScvs + numDirichletScvfs + idxInInteriorDiriFaces] += negWijk[localDir];
                        }
                        else // Dirichlet face
                            B[idxInFluxFaces][numLocalScvs + this->findIndexInVector(dirichletScvfIndexSet_(), curLocalScvfIdx)] += negWijk[localDir];

                        // add entries to matrix B
                        if (faceType == MpfaFaceTypes::interiorNeumann)
                            B[idxInFluxFaces][negLocalScvIdx] -= (1-xi)*negWijk[localDir];
                        else
                            B[idxInFluxFaces][negLocalScvIdx] -= negWijk[localDir];
                    }

                    // store the omegas (negative, because of normal vector sign switch)
                    negWijk *= -1.0;
                    wijk_[rowIdx].emplace_back(std::move(negWijk));

                    // increment counter in outside data
                    indexInOutsideData++;
                }
            }
            // go to the next face
            rowIdx++;
        }
    }

    template<typename GetTensorFunction>
    void assemblePureDirichletSystem_(const GetTensorFunction& getTensor)
    {
        const auto numLocalScvs = localScvs_.size();
        const auto numFaces = localScvfs_.size();
        const auto numInteriorDirichletFaces = interiorDirichletScvfIndexSet_().size();
        const auto numPotentials = volVarsStencil().size() + numInteriorDirichletFaces;

        // resize matrices, only T_ will have entries
        T_.resize(numFaces, numPotentials, 0.0);
        AinvB_.resize(0, 0);
        CAinv_.resize(0, 0);

        // resize the omegas
        wijk_.resize(numFaces);

        // Loop over all the faces, in this case these are all dirichlet boundaries
        LocalIndexType rowIdx = 0;
        for (auto&& localScvf : localScvfs_)
        {
            // get diffusion tensor in "positive" sub volume
            const auto posLocalScvIdx = localScvf.insideLocalScvIndex();
            const auto& posLocalScv = localScv_(posLocalScvIdx);
            const auto& posGlobalScv = fvGeometry_().scv(posLocalScv.globalIndex());
            const auto& posVolVars = elemVolVars_()[posGlobalScv];
            const auto element = localElement_(posLocalScvIdx);
            const auto tensor = getTensor(element, posVolVars, posGlobalScv);

            // the omega factors of the "positive" sub volume
            auto posWijk = calculateOmegas_(posLocalScv, localScvf.unitOuterNormal(), localScvf.area(), tensor);
            posWijk *= posVolVars.extrusionFactor();

            for (int localDir = 0; localDir < dim; localDir++)
            {
                // When interior boundaries are disabled, all faces will be of dirichlet type
                if (!enableInteriorBoundaries)
                {
                    const auto curLocalScvfIdx = posLocalScv.localScvfIndex(localDir);
                    const auto curIdxInDiriFaces = this->findIndexInVector(dirichletScvfIndexSet_(), curLocalScvfIdx);
                    T_[rowIdx][numLocalScvs + curIdxInDiriFaces] += posWijk[localDir];
                    T_[rowIdx][posLocalScvIdx] -= posWijk[localDir];
                }
                else
                {
                    const auto curLocalScvfIdx = posLocalScv.localScvfIndex(localDir);
                    const auto& curLocalScvf = localScvf_(curLocalScvfIdx);
                    const auto curFaceType = curLocalScvf.faceType();

                    const auto curIdxInDiriFaces = curFaceType == MpfaFaceTypes::dirichlet ?
                                                   this->findIndexInVector(dirichletScvfIndexSet_(), curLocalScvfIdx) :
                                                   numInteriorDirichletFaces + this->findIndexInVector(interiorDirichletScvfIndexSet_(), curLocalScvfIdx);

                    T_[rowIdx][numLocalScvs + curIdxInDiriFaces] += posWijk[localDir];
                    T_[rowIdx][posLocalScvIdx] -= posWijk[localDir];
                }
            }

            // store the omegas
            wijk_[rowIdx].emplace_back(std::move(posWijk));

            // go to the next face
            rowIdx++;
        }
    }

    // TODO: how to do the assertion of positive coefficients for tensors?
    DimVector calculateOmegas_(const LocalScvType& localScv,
                               const GlobalPosition normal,
                               const Scalar area,
                               const Tensor& T) const
    {
        DimVector wijk;
        GlobalPosition tmp;
        for (int dir = 0; dir < dim; ++dir)
        {
            T.mv(localScv.innerNormal(dir), tmp);
            wijk[dir] = tmp*normal;
        }
        wijk *= area;
        wijk /= localScv.detX();
        wijk *= -1.0;

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
        wijk *= -1.0;

        return wijk;
    }

    const Problem& problem_() const
    { return *problemPtr_; }

    const FVElementGeometry& fvGeometry_() const
    { return *fvGeometryPtr_; }

    const ElementVolumeVariables& elemVolVars_() const
    { return *elemVolVarsPtr_; }

    const Problem* problemPtr_;
    const FVElementGeometry* fvGeometryPtr_;
    const ElementVolumeVariables* elemVolVarsPtr_;

    bool onDomainOrInteriorBoundary_;

    std::vector<Element> localElements_;
    std::vector<LocalScvType> localScvs_;
    std::vector<LocalScvfType> localScvfs_;

    GlobalIndexSet globalScvfIndices_;
    GlobalIndexSet volVarsStencil_;
    PositionVector volVarsPositions_;

    LocalIndexSet fluxFaceIndexSet_;
    LocalIndexSet dirichletFaceIndexSet_;
    LocalIndexSet interiorDirichletFaceIndexSet_;
    LocalIndexSet interiorBoundaryFaceIndexSet_;
    std::vector<InteriorBoundaryData> interiorBoundaryData_;

    std::vector< std::vector< DimVector > > wijk_;

    DynamicMatrix T_;
    DynamicMatrix AinvB_;
    DynamicMatrix CAinv_;

    std::vector<PrimaryVariables> neumannFluxes_;
};

} // end namespace

#endif
