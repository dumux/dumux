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
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_FPS_INTERACTIONVOLUME_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_FPS_INTERACTIONVOLUME_HH

#include <dumux/common/math.hh>

#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>

#include <dumux/discretization/cellcentered/mpfa/interactionvolumebase.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

#include "localsubcontrolentities.hh"

namespace Dumux
{
//! Specialization of the interaction volume traits class
template<class TypeTag>
class CCMpfaOFpsInteractionVolumeTraits : public CCMpfaOInteractionVolumeTraits<TypeTag>
{
public:
    // Interior boundaries can not yet be handled by the currend o-method fps implementation
    // In that case we use the o-interactionvolume, otherwise we use its own interaction volumes at the boundary
    using BoundaryInteractionVolume = typename std::conditional<GET_PROP_VALUE(TypeTag, EnableInteriorBoundaries),
                                                                CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethod>,
                                                                CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethodFps>>::type;

    // The local sub-control volume type differs from the standard mpfa-o method
    using LocalScvType = CCMpfaOFpsLocalScv<TypeTag>;
};
/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volumes of the mpfa-o method with full pressure support.
 */
template<class TypeTag>
class CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethodFps>
          : public CCMpfaOInteractionVolume<TypeTag,
                                            CCMpfaOFpsInteractionVolumeTraits<TypeTag>,
                                            CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethodFps>>
{
    using Traits = CCMpfaOFpsInteractionVolumeTraits<TypeTag>;
    using ThisType = CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethodFps>;
    using ParentType = CCMpfaOInteractionVolume<TypeTag, Traits, ThisType>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using LocalScvType = typename Traits::LocalScvType;
    using LocalScvfType = typename Traits::LocalScvfType;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using CoordScalar = typename GridView::ctype;
    using FeCache = Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1>;
    using FeLocalBasis = typename FeCache::FiniteElementType::Traits::LocalBasisType;
    using ShapeJacobian = typename FeLocalBasis::Traits::JacobianType;
    using ShapeValue = typename Dune::FieldVector<Scalar, 1>;
    using JacobianInverseTransposed = typename LocalScvType::Geometry::JacobianInverseTransposed;

    using LocalPosition = typename LocalScvType::Geometry::LocalCoordinate;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using DynamicVector = typename Traits::Vector;
    using DynamicMatrix = typename Traits::Matrix;
    using Tensor = typename Traits::Tensor;

    enum
    {
        jacRows = JacobianInverseTransposed::rows,
        jacCols = JacobianInverseTransposed::cols
    };

    // stores the matrices required to calculate Tij
    struct LocalMatrixContainer
    {
        DynamicMatrix AL; //! coefficients of unknown face pressures (LHS)
        DynamicMatrix AR; //! coefficients of unknown face pressures (RHS)
        DynamicMatrix BL; //! coefficients of cell/Dirichlet pressures (LHS)
        DynamicMatrix BR; //! coefficients of cell/Dirichlet pressures (RHS)

        // the matrices for the expression of the fluxes
        DynamicMatrix AF; //! coefficients for the unknown face pressures
        DynamicMatrix BF; //! coefficients for the cell/Dirichlet pressures
    };

public:
    using typename ParentType::LocalIndexType;
    using typename ParentType::LocalFaceData;
    using typename ParentType::Seed;

    CCMpfaInteractionVolumeImplementation(const Seed& seed,
                                          const Problem& problem,
                                          const FVElementGeometry& fvGeometry,
                                          const ElementVolumeVariables& elemVolVars)
    : ParentType(seed, problem, fvGeometry, elemVolVars),
      c_(GET_PARAM_FROM_GROUP(TypeTag, Scalar, Mpfa, C)),
      p_(GET_PARAM_FROM_GROUP(TypeTag, Scalar, Mpfa, P)),
      divEqIdx_(this->fluxScvfIndexSet_().size())
    {
        if (dim == 3)
            DUNE_THROW(Dune::NotImplemented, "Fps scheme in 3d");

        //! add entry to the vector of neumann fluxes
        addAuxiliaryCellNeumannFlux_();
    }

    template<typename GetTensorFunction>
    void solveLocalSystem(const GetTensorFunction& getTensor)
    {
        const std::size_t numFluxFaces = this->fluxScvfIndexSet_().size();
        const std::size_t numUnknowns = numFluxFaces + 1;
        const std::size_t numFaces = this->localScvfs_.size();
        const std::size_t numPotentials = this->volVarsStencil().size();

        // instantiate and resize the local matrices
        LocalMatrixContainer mc;
        mc.AL.resize(numUnknowns, numUnknowns, 0.0);
        mc.AR.resize(numUnknowns, numUnknowns, 0.0);
        mc.BL.resize(numUnknowns, numPotentials, 0.0);
        mc.BR.resize(numUnknowns, numPotentials, 0.0);
        mc.AF.resize(numFaces, numUnknowns, 0.0);
        mc.BF.resize(numFaces, numPotentials, 0.0);

        // assemble the local eq system and matrices
        assembleLocalMatrices_(getTensor, mc);

        // solve local system and store transmissibility matrix
        mc.AL -= mc.AR;
        mc.BR -= mc.BL;
        mc.AL.invert();
        this->CAinv_ = Dumux::multiplyMatrices(mc.AF, mc.AL);
        this->T_ = Dumux::multiplyMatrices(mc.AF, Dumux::multiplyMatrices(mc.AL, mc.BR));
        this->T_ += mc.BF;
    }

private:

    void addAuxiliaryCellNeumannFlux_()
    {
        if (!this->onDomainOrInteriorBoundary() || GET_PROP_VALUE(TypeTag, UseTpfaBoundary))
            return;

        // add one entry to the vector of neumann fluxes
        auto& neumannFluxes = this->neumannFluxes_;
        neumannFluxes.resize(this->fluxScvfIndexSet_().size() + 1);
        neumannFluxes[neumannFluxes.size()-1] = std::accumulate(neumannFluxes.begin(), neumannFluxes.end()-1, PrimaryVariables(0.0));
    }

    template<typename GetTensorFunction>
    void assembleLocalMatrices_(const GetTensorFunction& getTensor, LocalMatrixContainer& mc)
    {
        // loop over the local faces
        LocalIndexType localFaceIdx = 0;
        for (const auto& localScvf : this->localScvfs_)
        {
            if (localScvf.boundary())
                assemblePositiveScv(getTensor, localScvf, localFaceIdx, mc, true);
            else
            {
                assemblePositiveScv(getTensor, localScvf, localFaceIdx, mc);
                assembleNegativeScv(getTensor, localScvf, localFaceIdx, mc);
            }

            // go to the next face
            localFaceIdx++;
        }
    }

    template<typename GetTensorFunction>
    void assemblePositiveScv(const GetTensorFunction& getTensor,
                             const LocalScvfType& localScvf,
                             LocalIndexType localScvfIdx,
                             LocalMatrixContainer& mc,
                             bool boundary = false)
    {
        // get diffusion tensor in "positive" sub volume
        auto localScvIdx = localScvf.insideLocalScvIndex();
        auto&& localScv = this->localScv_(localScvIdx);
        auto&& globalScv = this->fvGeometry_().scv(localScv.globalIndex());
        auto&& element = this->localElement_(localScvIdx);
        auto D = makeTensor_(getTensor(this->problem_(), element, this->elemVolVars_()[globalScv], this->fvGeometry_(), globalScv));
        // the local finite element basis
        const auto& localBasis = feCache_.get(localScv.geometry().type()).localBasis();

        // the normal and local integration point
        // On the ref element, normal vector points in the direction of local coordinate
        auto normalDir = localScv.getScvfIdxInScv(localScvfIdx);
        auto ipLocal = localScv.geometry().local(localScvf.ip());

        // find normal coordinate direction and integration point for divergence condition
        LocalIndexType divEqNormalDir = 1 - normalDir;
        LocalPosition divEqIpLocal(0.0);
        divEqIpLocal[divEqNormalDir] = divEqNormalDir == 1 ? c_ : 1.0 - (1.0-c_)*p_;
        divEqIpLocal[normalDir] = divEqNormalDir == 1 ? c_ + (1.0-c_)*p_ : c_;

        // does the face has an unknown associated with it?
        bool isFluxFace = localScvf.faceType() != MpfaFaceTypes::dirichlet;

        // assemble coefficients for the face fluxes
        addFaceFluxCoefficients_(localScv, localBasis, D, localScvfIdx, ipLocal, normalDir, mc, isFluxFace);
        // assemble matrix entries for the condition of zero divergence
        addDivEquationCoefficients_(localScv, localBasis, D, divEqIpLocal, divEqNormalDir, mc);

        // on dirichlet boundary faces, add coefficients for the boundary fluxes
        if (boundary && !isFluxFace)
        {
            LocalPosition bcIpLocal(0.0);
            bcIpLocal[normalDir] = 1.0;
            bcIpLocal[divEqNormalDir] = normalDir == 1 ? 1.0 - (1.0-c_)*p_ : c_ + (1.0-c_)*p_;
            addDivEquationCoefficients_(localScv, localBasis, D, bcIpLocal, normalDir, mc, boundary);
        }
    }

    template<typename GetTensorFunction>
    void assembleNegativeScv(const GetTensorFunction& getTensor,
                             const LocalScvfType& localScvf,
                             LocalIndexType localScvfIdx,
                             LocalMatrixContainer& mc)
    {
        // get diffusion tensor in "negative" sub volume
        for (auto localScvIdx : localScvf.outsideLocalScvIndices())
        {
            auto&& localScv = this->localScv_(localScvIdx);
            auto&& globalScv = this->fvGeometry_().scv(localScv.globalIndex());
            auto&& element = this->localElement_(localScvIdx);;
            auto D = makeTensor_(getTensor(this->problem_(), element, this->elemVolVars_()[globalScv], this->fvGeometry_(), globalScv));

            // the local finite element bases of the scvs
            const auto& localBasis = feCache_.get(localScv.geometry().type()).localBasis();

            // the normal and local integration point
            // On the ref element, normal vector points in the direction of local coordinate
            auto normalDir = localScv.getScvfIdxInScv(localScvfIdx);
            auto ipLocal = localScv.geometry().local(localScvf.ip());

            // find normals and integration points in the two scvs for condition of zero divergence
            LocalIndexType divEqNormalDir = 1 - normalDir;
            LocalPosition divEqIpLocal(0.0);
            divEqIpLocal[divEqNormalDir] = divEqNormalDir == 1 ? c_ : 1.0 - (1.0-c_)*p_;
            divEqIpLocal[normalDir] = divEqNormalDir == 1 ? c_ + (1.0-c_)*p_ : c_;

            // does the face has an unknown associated with it?
            bool isFluxFace = localScvf.faceType() != MpfaFaceTypes::dirichlet;

            // assemble coefficients for the face fluxes
            addFaceFluxCoefficients_(localScv, localBasis, D, localScvfIdx, ipLocal, normalDir, mc, isFluxFace, true);

            // assemble matrix entries for the condition of zero divergence
            addDivEquationCoefficients_(localScv, localBasis, D, divEqIpLocal, divEqNormalDir, mc);
        }
    }

    void addFaceFluxCoefficients_(const LocalScvType& localScv,
                                  const FeLocalBasis& localBasis,
                                  const Tensor& D,
                                  LocalIndexType rowIdx,
                                  const LocalPosition& ipLocal,
                                  LocalIndexType normalDir,
                                  LocalMatrixContainer& mc,
                                  bool isFluxEq,
                                  bool isRHS = false)
    {
        // In case we're on a flux continuity face, get local index
        LocalIndexType eqSystemIdx = isFluxEq ? this->findIndexInVector(this->fluxScvfIndexSet_(), rowIdx) : -1;

        // Fluxes stemming from the RHS have to have the opposite sign
        Scalar factor = isRHS ? 1.0 : -1.0;
        DynamicMatrix& A = isRHS ? mc.AR : mc.AL;
        DynamicMatrix& B = isRHS ? mc.BR : mc.BL;

        // evaluate shape functions gradients at the ip
        std::vector<ShapeJacobian> shapeJacobian;
        localBasis.evaluateJacobian(ipLocal, shapeJacobian);

        // get the inverse jacobian and its transposed
        const auto jacInvT = localScv.geometry().jacobianInverseTransposed(ipLocal);
        const auto jacInv = Dumux::getTransposed(jacInvT);

        // transform the diffusion tensor into local coordinates
        auto DJinvT = jacInvT.leftmultiplyany(D);
        auto localD = DJinvT.leftmultiplyany(jacInv);
        localD *= localScv.geometry().integrationElement(ipLocal);

        // add matrix entries for the pressure in the cell center
        auto cellPressureIdx = this->findIndexInVector(this->volVarsStencil(), localScv.globalIndex());
        Scalar bi0 = factor*(localD[normalDir]*shapeJacobian[0][0]);
        if (!isRHS) mc.BF[rowIdx][cellPressureIdx] += bi0;
        if (isFluxEq) B[eqSystemIdx][cellPressureIdx] += bi0;

        // Add entries from the local scv faces
        for (int localDir = 0; localDir < dim; localDir++)
        {
            auto localScvfIdx = localScv.localScvfIndex(localDir);
            if (this->localScvf_(localScvfIdx).faceType() != MpfaFaceTypes::dirichlet)
            {
                Scalar aij = factor*(localD[normalDir]*shapeJacobian[localDir+1][0]);
                auto colIdx = this->findIndexInVector(this->fluxScvfIndexSet_(), localScvfIdx);
                if (!isRHS) mc.AF[rowIdx][colIdx] += aij;
                if (isFluxEq) A[eqSystemIdx][colIdx] += aij;
            }
            else
            {
                Scalar bij = factor*(localD[normalDir]*shapeJacobian[localDir+1][0]);
                auto colIdx = this->localScvs_.size() + this->findIndexInVector(this->dirichletScvfIndexSet_(), localScvfIdx);
                if (!isRHS) mc.BF[rowIdx][colIdx] += bij;
                if (isFluxEq) B[eqSystemIdx][colIdx] += bij;
            }
        }

        // add entry from the vertex pressure
        Scalar ain = factor*(localD[normalDir]*shapeJacobian[3][0]);
        if (!isRHS) mc.AF[rowIdx][divEqIdx_] += ain;
        if (isFluxEq) A[eqSystemIdx][divEqIdx_] += ain;
    }

    void addDivEquationCoefficients_(const LocalScvType& localScv,
                                     const FeLocalBasis& localBasis,
                                     const Tensor& D,
                                     const LocalPosition& ipLocal,
                                     LocalIndexType normalDir,
                                     LocalMatrixContainer& mc,
                                     bool isBoundary = false)
    {
        // fluxes on the auxiliary volume have to be scaled
        static const Scalar auxArea = 1.0 - c_;
        Scalar factor = isBoundary ? -1.0*auxArea : auxArea;

        // evaluate shape functions gradients at the ip
        std::vector<ShapeJacobian> shapeJacobian;
        localBasis.evaluateJacobian(ipLocal, shapeJacobian);

        // get the inverse jacobian and its transposed
        const auto jacInvT = localScv.geometry().jacobianInverseTransposed(ipLocal);
        const auto jacInv = Dumux::getTransposed(jacInvT);

        // transform the diffusion tensor into local coordinates
        auto DJinvT = jacInvT.leftmultiplyany(D);
        auto localD = DJinvT.leftmultiplyany(jacInv);
        localD *= localScv.geometry().integrationElement(ipLocal);

        // add matrix entries for the pressure in the cell center
        auto cellPressureIdx = this->findIndexInVector(this->volVarsStencil(), localScv.globalIndex());
        mc.BL[divEqIdx_][cellPressureIdx] += factor*(localD[normalDir]*shapeJacobian[0][0]);

        // Add entries from the local scv faces
        for (int localDir = 0; localDir < dim; localDir++)
        {
            auto localScvfIdx = localScv.localScvfIndex(localDir);

            if (this->localScvf_(localScvfIdx).faceType() != MpfaFaceTypes::dirichlet)
            {
                auto colIdx = this->findIndexInVector(this->fluxScvfIndexSet_(), localScvfIdx);
                mc.AL[divEqIdx_][colIdx] += factor*(localD[normalDir]*shapeJacobian[localDir+1][0]);
            }
            else
            {
                auto colIdx = this->localScvs_.size() + this->findIndexInVector(this->dirichletScvfIndexSet_(), localScvfIdx);
                mc.BL[divEqIdx_][colIdx] += factor*(localD[normalDir]*shapeJacobian[localDir+1][0]);
            }
        }

        // add entry from the vertex pressure
        mc.AL[divEqIdx_][divEqIdx_] += factor*(localD[normalDir]*shapeJacobian[3][0]);
    }

    // TODO: how to do the assertion of positive coefficients for tensors?
    Tensor makeTensor_(Tensor&& tensor) const
    { return std::move(tensor); }

    // turns a scalar into a tensor
    Tensor makeTensor_(Scalar&& t) const
    {
        // make sure we have positive diffusion coefficients
        assert(t > 0.0 && "non-positive diffusion coefficients cannot be handled by mpfa methods");

        Tensor T(0.0);
        for (int i = 0; i < dimWorld; ++i)
            T[i][i] = t;

        return T;
    }

    const FeCache feCache_;

    const Scalar c_;
    const Scalar p_;
    const LocalIndexType divEqIdx_;
};

} // end namespace

#endif
