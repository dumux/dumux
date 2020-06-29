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
 * \ingroup FacetCoupling
 * \copydoc Dumux::MpfaOFacetCouplingInteractionVolumeAssembler
 */
#ifndef DUMUX_MULTIDOMAIN_FACET_CC_MPFA_O_LOCAL_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_FACET_CC_MPFA_O_LOCAL_ASSEMBLER_HH

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/localassemblerbase.hh>
#include <dumux/discretization/cellcentered/mpfa/localassemblerhelper.hh>
#include <dumux/discretization/cellcentered/mpfa/computetransmissibility.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Specialization of the interaction volume-local
 *        assembler class for the schemes using an mpfa-o type
 *        assembly in the context of facet coupling.
 *
 * \tparam P The problem type
 * \tparam EG The element finite volume geometry
 * \tparam EV The element volume variables type
 */
template< class P, class EG, class EV >
class MpfaOFacetCouplingInteractionVolumeAssembler
: public InteractionVolumeAssemblerBase< P, EG, EV >
{
    using ParentType = InteractionVolumeAssemblerBase< P, EG, EV >;
    using Helper = InteractionVolumeAssemblerHelper;
    using Extrusion = Extrusion_t<typename EG::GridGeometry>;

    template< class IV >
    using Scalar = typename IV::Traits::MatVecTraits::FaceVector::value_type;

public:
    //! Pull up constructor of the base class
    using ParentType::ParentType;

    /*!
     * \brief Assembles the matrices involved in the flux expressions
     *        and the local system of equations within an interaction volume.
     *
     * \param handle The data handle in which the matrices are stored
     * \param iv The interaction volume
     * \param getT Lambda to evaluate the scv-wise tensors
     * \param wijZeroThresh Threshold below which transmissibilities are taken to be zero.
     *                      On the basis of this threshold, trivial (0 = 0) rows in the
     *                      A matrix are identified and modified accordingly in order to
     *                      avoid ending up with singular matrices. This can occur when the
     *                      tensor is zero in some cells.
     */
    template< class DataHandle, class IV, class TensorFunc >
    void assembleMatrices(DataHandle& handle, IV& iv, const TensorFunc& getT, Scalar<IV> wijZeroThresh = 0.0)
    {
        assembleLocalMatrices_(handle.A(), handle.AB(), handle.CA(), handle.T(), handle.omegas(), iv, getT, wijZeroThresh);

        // maybe solve the local system
        if (iv.numUnknowns() > 0)
            Helper::solveLocalSystem(this->fvGeometry(), handle, iv);
    }

    /*!
     * \brief Assembles the vector of primary (cell) unknowns and (maybe)
     *        Dirichlet boundary conditions within an interaction volume.
     *
     * \param handle The data handle in which the vector is stored
     * \param iv The interaction volume
     * \param getU Lambda to obtain the desired cell/Dirichlet value from vol vars
     */
    template< class DataHandle, class IV, class GetU >
    void assembleU(DataHandle& handle, const IV& iv, const GetU& getU)
    {
        auto& u = handle.uj();
        Helper::resizeVector(u, iv.numKnowns());

        // put the cell unknowns first ...
        typename IV::Traits::IndexSet::LocalIndexType i = 0;
        for (; i < iv.numScvs(); i++)
            u[i] = getU( this->elemVolVars()[iv.localScv(i).gridScvIndex()] );

        // ... then put facet unknowns ...
        for (const auto& data : iv.interiorBoundaryData())
        {
            const auto& scvf = this->fvGeometry().scvf(data.scvfIndex());
            const auto element = this->fvGeometry().gridGeometry().element(scvf.insideScvIdx());
            u[i++] = getU( this->problem().couplingManager().getLowDimVolVars(element, scvf) );
        }

        // ... then put the Dirichlet unknowns
        for (const auto& data : iv.dirichletData())
            u[i++] = getU( this->elemVolVars()[data.volVarIndex()] );
    }

    /*!
     * \brief Assembles the gravitational flux contributions on the scvfs within an
     *        interaction volume.
     *
     * \param handle The data handle in which the vector is stored
     * \param iv The interaction volume
     * \param getRho Lambda to obtain the density from volume variables
     */
    template< class DataHandle, class IV, class GetRho >
    void assembleGravity(DataHandle& handle, const IV& iv, const GetRho& getRho)
    {
        using GridView = typename IV::Traits::GridView;
        static constexpr int dim = GridView::dimension;
        static constexpr int dimWorld = GridView::dimensionworld;
        static constexpr bool isSurfaceGrid = dim < dimWorld;

        // resize the gravity vectors
        auto& g = handle.g();
        auto& deltaG = handle.deltaG();
        auto& outsideG = handle.gOutside();
        Helper::resizeVector(g, iv.numFaces());
        Helper::resizeVector(deltaG, iv.numUnknowns());
        if (isSurfaceGrid)
            Helper::resizeVector(outsideG, iv.numFaces());

        //! For each face, we...
        //! - arithmetically average the phase densities
        //! - compute the term \f$ \alpha := \mathbf{A} \rho \ \mathbf{n}^T \mathbf{K} \mathbf{g} \f$ in each neighboring cell
        //! - compute \f$ \alpha^* = \sum{\alpha_{outside, i}} - \alpha_{inside} \f$
        using Scalar = typename IV::Traits::MatVecTraits::TMatrix::value_type;
        using LocalIndexType = typename IV::Traits::IndexSet::LocalIndexType;

        // xi factor for coupling conditions
        static const Scalar xi = getParamFromGroup<Scalar>(this->problem().paramGroup(), "FacetCoupling.Xi", 1.0);

        for (LocalIndexType faceIdx = 0; faceIdx < iv.numFaces(); ++faceIdx)
        {
            // gravitational acceleration on this face
            const auto& curLocalScvf = iv.localScvf(faceIdx);
            const auto& curGlobalScvf = this->fvGeometry().scvf(curLocalScvf.gridScvfIndex());
            const auto& gravity = this->problem().spatialParams().gravity(curGlobalScvf.ipGlobal());
            const auto curIsInteriorBoundary = curLocalScvf.isOnInteriorBoundary();
            const Scalar curXiFactor = curIsInteriorBoundary ? (curGlobalScvf.boundary() ? 1.0 : xi) : 1.0;

            // get permeability tensor in "positive" sub volume
            const auto& neighborScvIndices = curLocalScvf.neighboringLocalScvIndices();
            const auto& posGlobalScv = this->fvGeometry().scv(iv.localScv(neighborScvIndices[0]).gridScvIndex());
            const auto& posVolVars = this->elemVolVars()[posGlobalScv];
            const auto alpha_inside = posVolVars.extrusionFactor()*vtmv(curGlobalScvf.unitOuterNormal(),
                                                                        posVolVars.permeability(),
                                                                        gravity);

            const auto numOutsideFaces = !curGlobalScvf.boundary() ? curGlobalScvf.numOutsideScvs() : 0;
            using OutsideAlphaStorage = std::conditional_t< isSurfaceGrid,
                                                            std::vector<Scalar>,
                                                            Dune::ReservedVector<Scalar, 1> >;
            OutsideAlphaStorage alpha_outside; alpha_outside.resize(numOutsideFaces);
            std::fill(alpha_outside.begin(), alpha_outside.end(), 0.0);
            Scalar rho;

            if (isSurfaceGrid && !curIsInteriorBoundary)
                Helper::resizeVector(outsideG[faceIdx], numOutsideFaces);

            if (!curLocalScvf.isDirichlet())
            {
                const auto localDofIdx = curLocalScvf.localDofIndex();
                const Scalar curOneMinusXi = curIsInteriorBoundary ?  -(1.0 - xi) : 1.0;

                rho = getRho(posVolVars);
                deltaG[localDofIdx] = 0.0;

                if (!curGlobalScvf.boundary())
                {
                    for (unsigned int idxInOutside = 0; idxInOutside < curGlobalScvf.numOutsideScvs(); ++idxInOutside)
                    {
                        // obtain outside tensor
                        const auto negLocalScvIdx = neighborScvIndices[idxInOutside+1];
                        const auto& negGlobalScv = this->fvGeometry().scv(iv.localScv(negLocalScvIdx).gridScvIndex());
                        const auto& negVolVars = this->elemVolVars()[negGlobalScv];
                        const auto& flipScvf = !isSurfaceGrid ? curGlobalScvf
                                                              : this->fvGeometry().flipScvf(curGlobalScvf.index(), idxInOutside);

                        alpha_outside[idxInOutside] = negVolVars.extrusionFactor()*vtmv(flipScvf.unitOuterNormal(),
                                                                                        negVolVars.permeability(),
                                                                                        gravity);
                        if (isSurfaceGrid)
                            alpha_outside[idxInOutside] *= -1.0;

                        if (!curIsInteriorBoundary)
                            rho += getRho(negVolVars);

                        deltaG[localDofIdx] += curOneMinusXi*alpha_outside[idxInOutside];
                    }
                }

                if (curIsInteriorBoundary)
                {
                    const auto posElement = this->fvGeometry().gridGeometry().element(posGlobalScv.elementIndex());
                    const auto& facetVolVars = this->problem().couplingManager().getLowDimVolVars(posElement, curGlobalScvf);
                    rho += getRho(facetVolVars);
                    rho /= 2.0;
                    const auto alphaFacet = posVolVars.extrusionFactor()*vtmv(curGlobalScvf.unitOuterNormal(),
                                                                              facetVolVars.permeability(),
                                                                              gravity);
                    deltaG[localDofIdx] += alphaFacet;
                }
                else
                    rho /= numOutsideFaces + 1;

                deltaG[localDofIdx] -= curXiFactor*alpha_inside;
                deltaG[localDofIdx] *= rho*Extrusion::area(curGlobalScvf);
            }
            // use average density between facet and cell density on interior Dirichlet boundaries
            else if (curIsInteriorBoundary)
            {
                const auto posElement = this->fvGeometry().gridGeometry().element(posGlobalScv.elementIndex());
                const auto& facetVolVars = this->problem().couplingManager().getLowDimVolVars(posElement, curGlobalScvf);
                rho = 0.5*(getRho(facetVolVars) + getRho(posVolVars));
            }
            // use density resulting from Dirichlet BCs
            else
                rho = getRho(this->elemVolVars()[curGlobalScvf.outsideScvIdx()]);

            // add "inside" & "outside" alphas to gravity containers
            g[faceIdx] = alpha_inside*rho*Extrusion::area(curGlobalScvf);

            if (isSurfaceGrid && !curIsInteriorBoundary)
            {
                unsigned int i = 0;
                for (const auto& alpha : alpha_outside)
                    outsideG[faceIdx][i++] = alpha*rho*Extrusion::area(curGlobalScvf);
            }
        }

        // add iv-wide contributions to gravity vectors
        handle.CA().umv(deltaG, g);
        if (isSurfaceGrid)
        {
            using FaceVector = typename IV::Traits::MatVecTraits::FaceVector;
            FaceVector AG;
            Helper::resizeVector(AG, iv.numUnknowns());
            handle.A().mv(deltaG, AG);

            // compute gravitational accelerations
            for (const auto& localFaceData : iv.localFaceData())
            {
                // continue only for "outside" faces
                if (!localFaceData.isOutsideFace())
                    continue;

                const auto localScvIdx = localFaceData.ivLocalInsideScvIndex();
                const auto localScvfIdx = localFaceData.ivLocalScvfIndex();
                const auto idxInOutside = localFaceData.scvfLocalOutsideScvfIndex();
                const auto& posLocalScv = iv.localScv(localScvIdx);
                const auto& wijk = handle.omegas()[localScvfIdx][idxInOutside+1];

                // add contributions from all local directions
                for (LocalIndexType localDir = 0; localDir < dim; localDir++)
                {
                    // the scvf corresponding to this local direction in the scv
                    const auto& curLocalScvf = iv.localScvf(posLocalScv.localScvfIndex(localDir));
                    if (!curLocalScvf.isDirichlet())
                        outsideG[localScvfIdx][idxInOutside] -= wijk[localDir]*AG[curLocalScvf.localDofIndex()];
                }
            }
        }
    }

private:
    /*!
     * \copydoc Dumux::MpfaOInteractionVolumeAssembler::assembleLocalMatrices_
     */
    template< class IV, class TensorFunc >
    void assembleLocalMatrices_(typename IV::Traits::MatVecTraits::AMatrix& A,
                                typename IV::Traits::MatVecTraits::BMatrix& B,
                                typename IV::Traits::MatVecTraits::CMatrix& C,
                                typename IV::Traits::MatVecTraits::DMatrix& D,
                                typename IV::Traits::MatVecTraits::OmegaStorage& wijk,
                                IV& iv, const TensorFunc& getT,
                                Scalar<IV> wijZeroThresh)
    {
        using Scalar = typename IV::Traits::MatVecTraits::TMatrix::value_type;
        using LocalIndexType = typename IV::Traits::IndexSet::LocalIndexType;
        static constexpr int dim = IV::Traits::GridView::dimension;
        static constexpr int dimWorld = IV::Traits::GridView::dimensionworld;

        // xi factor for coupling conditions
        static const Scalar xi = getParamFromGroup<Scalar>(this->problem().paramGroup(), "FacetCoupling.Xi", 1.0);

        // On surface grids only xi = 1.0 can be used, as the coupling condition
        // for xi != 1.0 does not generalize for surface grids where the normal
        // vectors of the inside/outside elements have different orientations.
        if (dim < dimWorld)
            if (Dune::FloatCmp::ne(xi, 1.0, 1e-6))
                DUNE_THROW(Dune::InvalidStateException, "Xi != 1.0 cannot be used on surface grids");

        std::vector< std::pair<LocalIndexType, LocalIndexType> > faceMarkers;
        Helper::resizeVector(wijk, iv.numFaces());

        // if only Dirichlet faces are present in the iv,
        // the matrices A, B & C are undefined and D = T
        if (iv.numUnknowns() == 0)
        {
            // resize & reset D matrix
            Helper::resizeMatrix(D, iv.numFaces(), iv.numKnowns()); D = 0.0;

            // Loop over all the faces, in this case these are all dirichlet boundaries
            for (LocalIndexType faceIdx = 0; faceIdx < iv.numFaces(); ++faceIdx)
            {
                const auto& curLocalScvf = iv.localScvf(faceIdx);
                const auto& curGlobalScvf = this->fvGeometry().scvf(curLocalScvf.gridScvfIndex());
                const auto& neighborScvIndices = curLocalScvf.neighboringLocalScvIndices();

                // get tensor in "positive" sub volume
                const auto& posLocalScv = iv.localScv(neighborScvIndices[0]);
                const auto& posGlobalScv = this->fvGeometry().scv(posLocalScv.gridScvIndex());
                const auto& posVolVars = this->elemVolVars()[posGlobalScv];
                const auto tensor = getT(posVolVars);

                // the omega factors of the "positive" sub volume
                Helper::resizeVector(wijk[faceIdx], /*no outside scvs present*/1);
                wijk[faceIdx][0] = computeMpfaTransmissibility<EG>(posLocalScv, curGlobalScvf, tensor, posVolVars.extrusionFactor());

                const auto posScvLocalDofIdx = posLocalScv.localDofIndex();
                for (LocalIndexType localDir = 0; localDir < dim; localDir++)
                {
                    const auto& otherLocalScvf = iv.localScvf( posLocalScv.localScvfIndex(localDir) );
                    const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();
                    D[faceIdx][otherLocalDofIdx] -= wijk[faceIdx][0][localDir];
                    D[faceIdx][posScvLocalDofIdx] += wijk[faceIdx][0][localDir];
                }
            }
        }
        else
        {
            // resize & reset matrices
            Helper::resizeMatrix(A, iv.numUnknowns(), iv.numUnknowns()); A = 0.0;
            Helper::resizeMatrix(B, iv.numUnknowns(), iv.numKnowns());   B = 0.0;
            Helper::resizeMatrix(C, iv.numFaces(), iv.numUnknowns());    C = 0.0;
            Helper::resizeMatrix(D, iv.numFaces(), iv.numKnowns());      D = 0.0;

            for (LocalIndexType faceIdx = 0; faceIdx < iv.numFaces(); ++faceIdx)
            {
                const auto& curLocalScvf = iv.localScvf(faceIdx);
                const auto& curGlobalScvf = this->fvGeometry().scvf(curLocalScvf.gridScvfIndex());
                const auto curIsDirichlet = curLocalScvf.isDirichlet();
                const auto curLocalDofIdx = curLocalScvf.localDofIndex();
                const auto curIsInteriorBoundary = curLocalScvf.isOnInteriorBoundary();
                const Scalar curXiFactor = curIsInteriorBoundary ? (curGlobalScvf.boundary() ? 1.0 : xi) : 1.0;

                // get tensor in "positive" sub volume
                const auto& neighborScvIndices = curLocalScvf.neighboringLocalScvIndices();
                const auto& posLocalScv = iv.localScv(neighborScvIndices[0]);
                const auto& posGlobalScv = this->fvGeometry().scv(posLocalScv.gridScvIndex());
                const auto& posVolVars = this->elemVolVars()[posGlobalScv];
                const auto& posElement = iv.element(neighborScvIndices[0]);
                const auto tensor = getT(posVolVars);

                // the omega factors of the "positive" sub volume
                Helper::resizeVector(wijk[faceIdx], neighborScvIndices.size());
                wijk[faceIdx][0] = computeMpfaTransmissibility<EG>(posLocalScv, curGlobalScvf, tensor, posVolVars.extrusionFactor());

                using std::abs;
                bool isZeroWij = false;

                // go over the coordinate directions in the positive sub volume
                for (unsigned int localDir = 0; localDir < dim; localDir++)
                {
                    const auto& otherLocalScvf = iv.localScvf( posLocalScv.localScvfIndex(localDir) );
                    const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();

                    if (otherLocalDofIdx == curLocalDofIdx)
                    {
                        if (abs(wijk[faceIdx][0][localDir]) <= wijZeroThresh)
                        {
                            if (!curIsDirichlet)
                            {
                                isZeroWij = true;
                                faceMarkers.emplace_back( std::make_pair(curLocalDofIdx, faceIdx) );
                            }
                        }
                    }

                    // if we are not on a Dirichlet face, add entries associated with unknown face pressures
                    // i.e. in matrix C and maybe A (if current face is not a Dirichlet face)
                    if (!otherLocalScvf.isDirichlet())
                    {
                        C[faceIdx][otherLocalDofIdx] -= wijk[faceIdx][0][localDir];
                        if (!curIsDirichlet)
                            A[curLocalDofIdx][otherLocalDofIdx] -= curXiFactor*wijk[faceIdx][0][localDir];
                    }
                    // the current face is a Dirichlet face and creates entries in D & maybe B
                    else
                    {
                        D[faceIdx][otherLocalDofIdx] -= wijk[faceIdx][0][localDir];
                        if (!curIsDirichlet)
                            B[curLocalDofIdx][otherLocalDofIdx] += curXiFactor*wijk[faceIdx][0][localDir];
                    }

                    // add entries related to pressures at the scv centers (dofs)
                    const auto posScvLocalDofIdx = posLocalScv.localDofIndex();
                    D[faceIdx][posScvLocalDofIdx] += wijk[faceIdx][0][localDir];

                    if (!curIsDirichlet)
                        B[curLocalDofIdx][posScvLocalDofIdx] -= curXiFactor*wijk[faceIdx][0][localDir];
                }

                // handle the entries related to the coupled facet element
                if (curIsInteriorBoundary && !curIsDirichlet)
                {
                    const auto facetVolVars = this->problem().couplingManager().getLowDimVolVars(posElement, curGlobalScvf);
                    const auto facetTensor = getT(facetVolVars);

                    // On surface grids we use the square root of the extrusion factor as approximation of the aperture
                    using std::sqrt;
                    const auto wFacet = 2.0*Extrusion::area(curGlobalScvf)*posVolVars.extrusionFactor()
                                           *vtmv(curGlobalScvf.unitOuterNormal(), facetTensor, curGlobalScvf.unitOuterNormal())
                                           / (dim < dimWorld ? sqrt(facetVolVars.extrusionFactor()) : facetVolVars.extrusionFactor());

                    A[curLocalDofIdx][curLocalDofIdx] -= wFacet;
                    B[curLocalDofIdx][curLocalScvf.coupledFacetLocalDofIndex()] -= wFacet;

                    // check for zero transmissibilities (skip if inside has been zero already)
                    if (!isZeroWij && abs(wFacet) <= wijZeroThresh)
                    {
                        isZeroWij = true;
                        faceMarkers.emplace_back( std::make_pair(curLocalDofIdx, faceIdx) );
                    }
                }

                // If we are on an interior face (which isn't coupled via Dirichlet Bs), add values from negative sub volume
                if (!curGlobalScvf.boundary() && !curIsDirichlet)
                {
                    // the minus ensures the right sign of the coupling condition in the matrices
                    const Scalar curOneMinusXi = curIsInteriorBoundary ?  -(1.0 - xi) : 1.0;

                    // loop over all the outside neighbors of this face and add entries
                    for (unsigned int idxInOutside = 0; idxInOutside < curGlobalScvf.numOutsideScvs(); ++idxInOutside)
                    {
                        const auto idxOnScvf = idxInOutside+1;
                        const auto& negLocalScv = iv.localScv( neighborScvIndices[idxOnScvf] );
                        const auto& negGlobalScv = this->fvGeometry().scv(negLocalScv.gridScvIndex());
                        const auto& negVolVars = this->elemVolVars()[negGlobalScv];
                        const auto negTensor = getT(negVolVars);

                        // On surface grids, use outside face for "negative" transmissibility calculation
                        const auto& scvf = dim < dimWorld ? this->fvGeometry().flipScvf(curGlobalScvf.index(), idxInOutside)
                                                          : curGlobalScvf;
                        wijk[faceIdx][idxOnScvf] = computeMpfaTransmissibility<EG>(negLocalScv, scvf, negTensor, negVolVars.extrusionFactor());

                        // flip sign on surface grids (since we used the "outside" normal)
                        if (dim < dimWorld)
                            wijk[faceIdx][idxOnScvf] *= -1.0;

                        // go over the coordinate directions in the negative sub volume
                        for (int localDir = 0; localDir < dim; localDir++)
                        {
                            const auto otherLocalScvfIdx = negLocalScv.localScvfIndex(localDir);
                            const auto& otherLocalScvf = iv.localScvf(otherLocalScvfIdx);
                            const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();

                            // check for zero transmissibilities (skip if inside has been zero already)
                            if (otherLocalDofIdx == curLocalDofIdx && !isZeroWij)
                                if (abs(wijk[faceIdx][idxOnScvf][localDir]) <= wijZeroThresh)
                                    faceMarkers.emplace_back( std::make_pair(curLocalDofIdx, faceIdx) );

                            if (!otherLocalScvf.isDirichlet())
                                A[curLocalDofIdx][otherLocalDofIdx] += curOneMinusXi*wijk[faceIdx][idxOnScvf][localDir];
                            else
                                B[curLocalDofIdx][otherLocalDofIdx] -= curOneMinusXi*wijk[faceIdx][idxOnScvf][localDir];

                            // add entries to matrix B
                            B[curLocalDofIdx][negLocalScv.localDofIndex()] += curOneMinusXi*wijk[faceIdx][idxOnScvf][localDir];
                        }
                    }
                }
            }
        }
    }
};

} // end namespace

#endif
