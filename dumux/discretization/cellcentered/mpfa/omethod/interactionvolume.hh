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
#include <dumux/discretization/cellcentered/mpfa/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>

#include "localsubcontrolentities.hh"
#include "interactionvolumeindexset.hh"

namespace Dumux
{
//! Specialization of the interaction volume traits class for the mpfa-o method
template<class TypeTag>
class CCMpfaOInteractionVolumeTraits : public CCMpfaInteractionVolumeTraitsBase<TypeTag>
{
    using BaseTraits = CCMpfaInteractionVolumeTraitsBase<TypeTag>;
    using NodalIndexSet = typename CCMpfaDualGridIndexSet<TypeTag>::NodalIndexSet;

public:
    using SecondaryInteractionVolume = CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethod>;

    using typename BaseTraits::DynamicLocalIndexContainer;
    using typename BaseTraits::DynamicGlobalIndexContainer;
    using IndexSet = CCMpfaOInteractionVolumeIndexSet<NodalIndexSet, DynamicGlobalIndexContainer, DynamicLocalIndexContainer>;

    // In the o-scheme, matrix & vector types are dynamic types
    using typename BaseTraits::DynamicVector;
    using typename BaseTraits::DynamicMatrix;
    using Vector = DynamicVector;
    using Matrix = DynamicMatrix;

    using LocalScvType = CCMpfaOInteractionVolumeLocalScv<TypeTag, IndexSet>;
    using LocalScvfType = CCMpfaOInteractionVolumeLocalScvf<TypeTag>;
};

//! Forward declaration of the mpfa-o interaction volume
template<class TypeTag, class Traits>
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
                                            CCMpfaOInteractionVolumeTraits<TypeTag>>
{
    using TraitsType = CCMpfaOInteractionVolumeTraits<TypeTag>;
public:
    // state the traits class type
    using Traits = TraitsType;
};

template<class TypeTag, class Traits>
class CCMpfaOInteractionVolume : public CCMpfaInteractionVolumeBase<TypeTag, Traits>
{
    using ParentType = CCMpfaInteractionVolumeBase<TypeTag, Traits>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using DualGridNodalIndexSet = typename CCMpfaDualGridIndexSet<TypeTag>::NodalIndexSet;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using Vector = typename Traits::DynamicVector;
    using Matrix = typename Traits::DynamicMatrix;
    using Tensor = typename Traits::Tensor;

    using LocalScvType = typename Traits::LocalScvType;
    using LocalScvfType = typename Traits::LocalScvfType;

    using IndexSet = typename Traits::IndexSet;
    using LocalIndexContainer = typename Traits::DynamicLocalIndexContainer;
    using LocalIndexType = typename LocalIndexContainer::value_type;
    using GlobalIndexContainer = typename Traits::DynamicGlobalIndexContainer;
    using DataHandle = typename Traits::DataHandle;
public:

    using typename ParentType::LocalFaceData;
    using typename ParentType::DirichletDataContainer;
    using typename ParentType::LocalFaceDataContainer;

    //! Sets up the local scope for a given iv index set!
    void setUpLocalScope(const IndexSet& indexSet,
                         const Problem& problem,
                         const FVElementGeometry& fvGeometry)
    {
        //! store a pointer to the index set
        indexSetPtr_ = &indexSet;

        //! clear previous data
        clear_();

        //! number of interaction-volume-local faces
        numFaces_ = indexSet.numFaces();

        //! number of interaction-volume-local (and node-local) scvs
        const auto& scvIndices = indexSet.globalScvIndices();
        const auto numLocalScvs = indexSet.numScvs();

        //! number of global scvfs appearing in this interaction volume
        const auto numGlobalScvfs = indexSet.nodalIndexSet().numScvfs();

        //! reserve memory for local entities
        elements_.reserve(numLocalScvs);
        scvs_.reserve(numLocalScvs);
        scvfs_.reserve(numFaces_);
        dirichletData_.reserve(numFaces_);
        localFaceData_.reserve(numGlobalScvfs);

        // set up quantities related to sub-control volumes
        for (LocalIndexType scvIdxLocal = 0; scvIdxLocal < numLocalScvs; scvIdxLocal++)
        {
            const auto scvIdxGlobal = scvIndices[scvIdxLocal];
            scvs_.emplace_back(fvGeometry, fvGeometry.scv(scvIdxGlobal), scvIdxLocal, indexSet);
            elements_.emplace_back(fvGeometry.fvGridGeometry().element(scvIdxGlobal));
        }

        // keep track of the number of unknowns etc
        numUnknowns_ = 0;
        numOutsideFaces_ = 0;
        numPotentials_ = numLocalScvs;

        // set up quantitites related to sub-control volume faces
        for (LocalIndexType faceIdxLocal = 0; faceIdxLocal < numFaces_; ++faceIdxLocal)
        {
            const auto scvfIdxGlobal = indexSet.scvfIdxGlobal(faceIdxLocal);
            const auto& neighborScvIndicesLocal = indexSet.neighboringLocalScvIndices(faceIdxLocal);
            const auto insideLocalScvIdx = neighborScvIndicesLocal[0];

            // we have to use the "inside" scv face here
            const auto& scvf = fvGeometry.scvf(scvfIdxGlobal);

            // create local face data object for this face
            localFaceData_.emplace_back(faceIdxLocal, insideLocalScvIdx, scvf.index());

            // create iv-local scvf object
            if (scvf.boundary())
            {
                const auto insideElement = elements_[insideLocalScvIdx];
                const auto bcTypes = problem.boundaryTypes(insideElement, scvf);

                if (bcTypes.hasOnlyDirichlet())
                {
                    scvfs_.emplace_back(scvf, neighborScvIndicesLocal, /*isDirichlet*/true, numPotentials_++);
                    dirichletData_.emplace_back(scvf.outsideScvIdx(), scvf.ipGlobal());
                }
                else
                    scvfs_.emplace_back(scvf, neighborScvIndicesLocal, /*isDirichlet*/false, numUnknowns_++);
            }
            else
            {
                scvfs_.emplace_back(scvf, neighborScvIndicesLocal, /*isDirichlet*/false, numUnknowns_++);

                // add local face data object for the outside faces
                for (unsigned int i = 1; i < neighborScvIndicesLocal.size(); ++i)
                {
                    const auto outsideLocalScvIdx = neighborScvIndicesLocal[i];

                    // loop over scvfs in outside scv until we find the one coinciding with current scvf
                    for (int coord = 0; coord < dim; ++coord)
                    {
                        if (indexSet.scvfIdxLocal(outsideLocalScvIdx, coord) == faceIdxLocal)
                        {
                            const auto globalScvfIdx = indexSet.nodalIndexSet().scvfIdxGlobal(outsideLocalScvIdx, coord);
                            const auto& flipScvf = fvGeometry.scvf(globalScvfIdx);
                            localFaceData_.emplace_back(faceIdxLocal,         //! iv-local scvf idx
                                                        outsideLocalScvIdx,   //! iv-local scv index
                                                        numOutsideFaces_++,    //! iv-local index in outside faces
                                                        i-1,                  //! scvf-local index in outside faces
                                                        flipScvf.index());   //! global scvf index
                        }
                    }
                }
            }
        }

        // resize the local matrices
        A_.resize(numUnknowns_, numUnknowns_);
        B_.resize(numUnknowns_, numPotentials_);
        C_.resize(numFaces_, numUnknowns_);
        D_.resize(numFaces_, numPotentials_);
    }

    //! sets the sizes of the corresponding matrices in the data handle
    void prepareDataHandle(DataHandle& dataHandle)
    {
      // resize the transmissibility matrix in the data handle
      dataHandle.resizeT(numFaces_, numPotentials_);

      // resize possible additional containers in the data handle
      if (requireABMatrix_()) dataHandle.resizeAB(numUnknowns_, numPotentials_);
      if (dim < dimWorld) dataHandle.resizeOutsideTij(numOutsideFaces_, numPotentials_);
    }

    //! solves for the transmissibilities subject to a given tensor
    template<typename GetTensorFunction>
    void solveLocalSystem(const GetTensorFunction& getTensor,
                          const Problem& problem,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          DataHandle& dataHandle)
    {
        // if only dirichlet faces are present, assemble T_ directly
        if (numUnknowns_ == 0)
            assemblePureDirichletSystem_(getTensor, problem, fvGeometry, elemVolVars, dataHandle.T());
        else
        {
            // assemble
            assembleLocalMatrices_(getTensor, problem, fvGeometry, elemVolVars);

            // solve
            A_.invert();

            // T = C*A^-1*B + D
            dataHandle.T() = multiplyMatrices(C_.rightmultiply(A_), B_);
            dataHandle.T() += D_;

            // store A-1B only when gradient reconstruction is necessary
            if (requireABMatrix_())
                dataHandle.AB() = B_.leftmultiply(A_);
        }

        // // set vol vars stencil & positions pointer in handle
        dataHandle.setVolVarsStencilPointer(indexSet().globalScvIndices());
        dataHandle.setDirichletDataPointer(dirichletData_);

        // on surface grids, additionally prepare the outside transmissibilities
        if (dim < dimWorld)
            computeOutsideTransmissibilities_(dataHandle);
    }

    //! obtain the local data object for a given global scvf
    const LocalFaceData& getLocalFaceData(const SubControlVolumeFace& scvf) const
    {
        //! find corresponding entry in the local face data container
        const auto scvfIdxGlobal = scvf.index();
        auto it = std::find_if(localFaceData_.begin(),
                               localFaceData_.end(),
                               [scvfIdxGlobal] (const LocalFaceData& d) { return d.globalScvfIndex() == scvfIdxGlobal; });
        assert(it != localFaceData_.end() && "Could not find the local face data corresponding to the given scvf");
        return localFaceData_[std::distance(localFaceData_.begin(), it)];
    }

    //! returns the grid element corresponding to a given iv-local scv idx
    const Element& element(const LocalIndexType ivLocalScvIdx) const
    { return elements_[ivLocalScvIdx]; }

    //! returns the local scvf entity corresponding to a given iv-local scvf idx
    const LocalScvfType& localScvf(const LocalIndexType ivLocalScvfIdx) const
    { return scvfs_[ivLocalScvfIdx]; }

    //! returns the local scv entity corresponding to a given iv-local scv idx
    const LocalScvType& localScv(const LocalIndexType ivLocalScvfIdx) const
    { return scvs_[ivLocalScvfIdx]; }

    //! returns a reference to the container with the data on Dirichlet boundaries
    const DirichletDataContainer& dirichletData() const
    { return dirichletData_; }

    //! returns a reference to the container with the local face data
    const std::vector<LocalFaceData>& localFaceData() const
    { return localFaceData_; }

    //! returns a reference to the index set of this iv
    const IndexSet& indexSet() const
    { return *indexSetPtr_; }

    //! returns the number of interaction volumes living around a vertex
    //! the mpfa-o scheme always constructs one iv per vertex
    static std::size_t numInteractionVolumesAtVertex(const DualGridNodalIndexSet& nodalIndexSet)
    { return 1; }

    //! adds the iv index sets living around a vertex to a given container
    //! and stores the the corresponding index in a map for each scvf
    template<class IvIndexSetContainer, class ScvfIndexMap>
    static void addInteractionVolumeIndexSets(IvIndexSetContainer& ivIndexSetContainer,
                                              ScvfIndexMap& scvfIndexMap,
                                              const DualGridNodalIndexSet& nodalIndexSet)
    {
        // the global index of the iv index set that is about to be created
        const auto curGlobalIndex = ivIndexSetContainer.size();

        // make the one index set for this node
        ivIndexSetContainer.emplace_back(nodalIndexSet);

        // store the index mapping
        for (const auto scvfIdx : nodalIndexSet.globalScvfIndices())
            scvfIndexMap[scvfIdx] = curGlobalIndex;
    }

private:
    //! returns a boolean whether or not the AB matrix has to be passed to the handles
    bool requireABMatrix_() const
    {
        static const bool requireAB = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Vtk.AddVelocity") || dim < dimWorld;
        return requireAB;
    }

    //! clears all the containers
    void clear_()
    {
        elements_.clear();
        scvs_.clear();
        scvfs_.clear();
        localFaceData_.clear();
        dirichletData_.clear();
    }

    //! Assembles the local matrices that define the local system of equations and flux expressions
    template<typename GetTensorFunction>
    void assembleLocalMatrices_(const GetTensorFunction& getTensor,
                                const Problem& problem,
                                const FVElementGeometry& fvGeometry,
                                const ElementVolumeVariables& elemVolVars)
    {
        // reset matrices
        A_ = 0.0;
        B_ = 0.0;
        C_ = 0.0;
        D_ = 0.0;

        // reserve space for the omegas
        wijk_.resize(scvfs_.size());

        // loop over the local faces
        for (unsigned int faceIdx = 0; faceIdx < numFaces_; ++faceIdx)
        {
            const auto& curLocalScvf = localScvf(faceIdx);
            const auto& curGlobalScvf = fvGeometry.scvf(curLocalScvf.globalScvfIndex());
            const auto curIsDirichlet = curLocalScvf.isDirichlet();
            const auto curLocalDofIdx = curLocalScvf.localDofIndex();

            // get diffusion tensor in "positive" sub volume
            const auto& neighborScvIndices = curLocalScvf.neighboringLocalScvIndices();
            const auto posLocalScvIdx = neighborScvIndices[0];
            const auto& posLocalScv = localScv(posLocalScvIdx);
            const auto& posGlobalScv = fvGeometry.scv(posLocalScv.globalScvIndex());
            const auto& posVolVars = elemVolVars[posGlobalScv];
            const auto& posElement = element(posLocalScvIdx);
            const auto tensor = getTensor(problem, posElement, posVolVars, fvGeometry, posGlobalScv);

            // the omega factors of the "positive" sub volume
            auto posWijk = calculateOmegas_(posLocalScv, curGlobalScvf.unitOuterNormal(), curGlobalScvf.area(), tensor);
            posWijk *= posVolVars.extrusionFactor();

            // go over the coordinate directions in the positive sub volume
            for (unsigned int localDir = 0; localDir < dim; localDir++)
            {
                const auto otherLocalScvfIdx = posLocalScv.scvfIdxLocal(localDir);
                const auto& otherLocalScvf = localScvf(otherLocalScvfIdx);
                const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();

                // if we are not on a Dirichlet face, add entries associated with unknown face pressures
                // i.e. in matrix C and maybe A (if current face is not a Dirichlet face)
                if (!otherLocalScvf.isDirichlet())
                {
                    C_[faceIdx][otherLocalDofIdx] -= posWijk[localDir];
                    if (!curIsDirichlet)
                        A_[curLocalDofIdx][otherLocalDofIdx] -= posWijk[localDir];
                }
                // the current face is a Dirichlet face and creates entries in D & maybe B
                else
                {
                    D_[faceIdx][otherLocalDofIdx] -= posWijk[localDir];
                    if (!curIsDirichlet)
                        B_[curLocalDofIdx][otherLocalDofIdx] += posWijk[localDir];
                }

                // add entries related to pressures at the scv centers (dofs)
                const auto posScvLocalDofIdx = posLocalScv.localDofIndex();
                D_[faceIdx][posScvLocalDofIdx] += posWijk[localDir];

                if (!curIsDirichlet)
                    B_[curLocalDofIdx][posScvLocalDofIdx] -= posWijk[localDir];
            }

            // store the omegas
            wijk_[faceIdx].emplace_back(std::move(posWijk));

            // If we are on an interior face, add values from negative sub volume
            if (!curGlobalScvf.boundary())
            {
                // loop over all the outside neighbors of this face and add entries
                for (unsigned int idxInOutside = 0; idxInOutside < curGlobalScvf.numOutsideScvs(); ++idxInOutside)
                {
                    const auto negLocalScvIdx = neighborScvIndices[idxInOutside+1];
                    const auto& negLocalScv = localScv(negLocalScvIdx);
                    const auto& negGlobalScv = fvGeometry.scv(negLocalScv.globalScvIndex());
                    const auto& negVolVars = elemVolVars[negGlobalScv];
                    const auto& negElement = element(negLocalScvIdx);
                    const auto negTensor = getTensor(problem, negElement, negVolVars, fvGeometry, negGlobalScv);

                    // the omega factors of the "negative" sub volume
                    DimVector negWijk;

                    // if dim < dimWorld, use outside normal vector
                    if (dim < dimWorld)
                    {
                        const auto& flipScvf = fvGeometry.flipScvf(curGlobalScvf.index(), idxInOutside);
                        auto negNormal = flipScvf.unitOuterNormal();
                        negNormal *= -1.0;
                        negWijk = calculateOmegas_(negLocalScv, negNormal, curGlobalScvf.area(), negTensor);
                    }
                    else
                        negWijk = calculateOmegas_(negLocalScv, curGlobalScvf.unitOuterNormal(), curGlobalScvf.area(), negTensor);

                    // scale by extrusion factpr
                    negWijk *= negVolVars.extrusionFactor();

                    // go over the coordinate directions in the positive sub volume
                    for (int localDir = 0; localDir < dim; localDir++)
                    {
                        const auto otherLocalScvfIdx = negLocalScv.scvfIdxLocal(localDir);
                        const auto& otherLocalScvf = localScvf(otherLocalScvfIdx);
                        const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();

                        if (!otherLocalScvf.isDirichlet())
                            A_[curLocalDofIdx][otherLocalDofIdx] += negWijk[localDir];
                        else
                            B_[curLocalDofIdx][otherLocalDofIdx] -= negWijk[localDir];

                        // add entries to matrix B
                        B_[curLocalDofIdx][negLocalScv.localDofIndex()] += negWijk[localDir];
                    }

                    // store the omegas
                    wijk_[faceIdx].emplace_back(std::move(negWijk));
                }
            }
        }
    }

    //! for interaction volumes that have only dirichlet scvfs,
    //! the transmissibility matrix can be assembled directly
    template<typename GetTensorFunction>
    void assemblePureDirichletSystem_(const GetTensorFunction& getTensor,
                                      const Problem& problem,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVolumeVariables& elemVolVars,
                                      Matrix& T)
    {
        // reset the transmissibility matrix beforehand
        T = 0.0;

        // Loop over all the faces, in this case these are all dirichlet boundaries
        for (unsigned int faceIdx = 0; faceIdx < numFaces_; ++faceIdx)
        {
            const auto& curLocalScvf = localScvf(faceIdx);
            const auto& curGlobalScvf = fvGeometry.scvf(curLocalScvf.globalScvfIndex());

            // get diffusion tensor in "positive" sub volume
            const auto& neighborScvIndices = curLocalScvf.neighboringLocalScvIndices();
            const auto posLocalScvIdx = neighborScvIndices[0];
            const auto& posLocalScv = localScv(posLocalScvIdx);
            const auto& posGlobalScv = fvGeometry.scv(posLocalScv.globalScvIndex());
            const auto& posVolVars = elemVolVars[posGlobalScv];
            const auto& posElement = element(posLocalScvIdx);
            const auto tensor = getTensor(problem, posElement, posVolVars, fvGeometry, posGlobalScv);

            // the omega factors of the "positive" sub volume
            auto posWijk = calculateOmegas_(posLocalScv, curGlobalScvf.unitOuterNormal(), curGlobalScvf.area(), tensor);
            posWijk *= posVolVars.extrusionFactor();

            const auto posScvLocalDofIdx = posLocalScv.localDofIndex();
            for (int localDir = 0; localDir < dim; localDir++)
            {
                const auto otherLocalScvfIdx = posLocalScv.scvfIdxLocal(localDir);
                const auto& otherLocalScvf = localScvf(otherLocalScvfIdx);
                const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();
                T[faceIdx][otherLocalDofIdx] -= posWijk[localDir];
                T[faceIdx][posScvLocalDofIdx] += posWijk[localDir];
            }
        }
    }

    //! computes the transmissibilities associated with "outside" faces on surface grids
    void computeOutsideTransmissibilities_(DataHandle& dataHandle) const
    {
        assert(dim < dimWorld && "only for dim < dimWorld the outside transmissiblity container has the right size");

        for (const auto& localFaceData : localFaceData_)
        {
            //! continue only for "outside" faces
            if (!localFaceData.isOutside()) continue;

            const auto localScvIdx = localFaceData.ivLocalInsideScvIndex();
            const auto localScvfIdx = localFaceData.ivLocalScvfIndex();
            const auto& posLocalScv = localScv(localScvIdx);
            const auto& wijk = wijk_[localScvfIdx][localFaceData.scvfLocalOutsideScvfIndex() + 1];

            //! store the calculated transmissibilities in the data handle
            auto& tij = dataHandle.outsideTij()[localFaceData.ivLocalOutsideScvfIndex()];

            //! reset transmissibility vector
            tij = 0.0;

            //! add contributions from all local directions
            for (int localDir = 0; localDir < dim; localDir++)
            {
                //! the scvf corresponding to this local direction in the scv
                const auto& curLocalScvf = localScvf(posLocalScv.scvfIdxLocal(localDir));

                //! on interior faces the coefficients of the AB matrix come into play
                if (!curLocalScvf.isDirichlet())
                {
                    auto tmp = dataHandle.AB()[curLocalScvf.localDofIndex()];
                    tmp *= wijk[localDir];
                    tij -= tmp;
                }
                else
                    tij[curLocalScvf.localDofIndex()] -= wijk[localDir];

                // add entry from the scv unknown
                tij[localScvIdx] += wijk[localDir];
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

    const IndexSet* indexSetPtr_;

    // Variables defining the local scope
    std::vector<Element> elements_;
    std::vector<LocalScvType> scvs_;
    std::vector<LocalScvfType> scvfs_;
    std::vector<LocalFaceData> localFaceData_;
    DirichletDataContainer dirichletData_;

    // sizes involved in the local matrices
    unsigned int numFaces_;
    unsigned int numUnknowns_;
    unsigned int numPotentials_;
    unsigned int numOutsideFaces_;

    // The omega factors and the matrix A⁻1*B are stored
    // in order to recover the transmissibilities of outside faces on network grids
    std::vector< std::vector< DimVector > > wijk_;
    Matrix CAinv_;

    // Matrices involved in transmissibility calculations
    Matrix A_;
    Matrix B_;
    Matrix C_;
    Matrix D_;
};

} // end namespace

#endif
