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
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \ingroup CCMpfaDiscretization
 * \copydoc Dumux::MpfaOFacetCouplingInteractionVolumeAssembler
 */
#ifndef DUMUX_MULTIDOMAIN_FACET_CC_MPFA_O_LOCAL_ASSEMBLER_HH
#define DUMUX_MULTIDOMAIN_FACET_CC_MPFA_O_LOCAL_ASSEMBLER_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/math.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/localassemblerbase.hh>
#include <dumux/discretization/cellcentered/mpfa/localassemblerhelper.hh>
#include <dumux/discretization/cellcentered/mpfa/computetransmissibility.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \ingroup CCMpfaDiscretization
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
     */
    template< class DataHandle, class IV, class TensorFunc >
    void assembleMatrices(DataHandle& handle, IV& iv, const TensorFunc& getT)
    {
        assembleLocalMatrices_(handle.A(), handle.AB(), handle.CA(), handle.T(), handle.omegas(), iv, getT);

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
            const auto element = this->fvGeometry().fvGridGeometry().element(scvf.insideScvIdx());
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
        DUNE_THROW(Dune::NotImplemented, "mpfa facet coupling with gravity");
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
                                IV& iv, const TensorFunc& getT)
    {
        using Scalar = typename IV::Traits::MatVecTraits::TMatrix::value_type;
        using LocalIndexType = typename IV::Traits::IndexSet::LocalIndexType;
        static constexpr int dim = IV::Traits::GridView::dimension;
        static constexpr int dimWorld = IV::Traits::GridView::dimensionworld;

        // xi factor for coupling conditions
        static const Scalar xi = getParamFromGroup<Scalar>(this->problem().paramGroup(), "FacetCoupling.Xi", 1.0);

        // resize omegas
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
                const auto& posElement = iv.element(neighborScvIndices[0]);
                const auto tensor = getT(this->problem(), posElement, posVolVars, this->fvGeometry(), posGlobalScv);

                // the omega factors of the "positive" sub volume
                Helper::resizeVector(wijk[faceIdx], /*no outside scvs present*/1);
                wijk[faceIdx][0] = computeMpfaTransmissibility(posLocalScv, curGlobalScvf, tensor, posVolVars.extrusionFactor());

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
                const auto tensor = getT(this->problem(), posElement, posVolVars, this->fvGeometry(), posGlobalScv);

                // the omega factors of the "positive" sub volume
                Helper::resizeVector(wijk[faceIdx], neighborScvIndices.size());
                wijk[faceIdx][0] = computeMpfaTransmissibility(posLocalScv, curGlobalScvf, tensor, posVolVars.extrusionFactor());

                // go over the coordinate directions in the positive sub volume
                for (unsigned int localDir = 0; localDir < dim; localDir++)
                {
                    const auto& otherLocalScvf = iv.localScvf( posLocalScv.localScvfIndex(localDir) );
                    const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();

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
                    const auto facetTensor = this->problem().couplingManager().getLowDimTensor(posElement, curGlobalScvf, getT);

                    // On surface grids we use the square root of the extrusion factor as approximation of the aperture
                    using std::sqrt;
                    const auto wFacet = 2.0*curGlobalScvf.area()*posVolVars.extrusionFactor()
                                           *vtmv(curGlobalScvf.unitOuterNormal(), facetTensor, curGlobalScvf.unitOuterNormal())
                                           / (dim < dimWorld ? sqrt(facetVolVars.extrusionFactor()) : facetVolVars.extrusionFactor());

                    A[curLocalDofIdx][curLocalDofIdx] -= wFacet;
                    B[curLocalDofIdx][curLocalScvf.coupledFacetLocalDofIndex()] -= wFacet;
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
                        const auto& negElement = iv.element( neighborScvIndices[idxOnScvf] );
                        const auto negTensor = getT(this->problem(), negElement, negVolVars, this->fvGeometry(), negGlobalScv);

                        // On surface grids, use outside face for "negative" transmissibility calculation
                        const auto& scvf = dim < dimWorld ? this->fvGeometry().flipScvf(curGlobalScvf.index(), idxInOutside)
                                                          : curGlobalScvf;
                        wijk[faceIdx][idxOnScvf] = computeMpfaTransmissibility(negLocalScv, scvf, negTensor, negVolVars.extrusionFactor());

                        // flip sign on surface grids (since we used the "outside" normal)
                        if (dim < dimWorld)
                            wijk[faceIdx][idxOnScvf] *= -1.0;

                        // go over the coordinate directions in the negative sub volume
                        for (int localDir = 0; localDir < dim; localDir++)
                        {
                            const auto otherLocalScvfIdx = negLocalScv.localScvfIndex(localDir);
                            const auto& otherLocalScvf = iv.localScvf(otherLocalScvfIdx);
                            const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();

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
