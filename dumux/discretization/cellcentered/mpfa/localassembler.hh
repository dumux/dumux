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
 * \ingroup CCMpfaDiscretization
 * \brief Defines the general interface of classes used for the assembly
 *        of the local systems of equations involved in the transmissibility
 *        computaion in mpfa schemes.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_LOCAL_ASSEMBLER_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_LOCAL_ASSEMBLER_HH

#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

namespace Dumux
{
//! Forward declaration of the implementation
template< class P, class EG, class EV, MpfaMethods M > class InteractionVolumeAssemblerImpl;

//! Alias to select the right implementation.
template< class P, class EG, class EV, MpfaMethods M >
using InteractionVolumeAssembler = InteractionVolumeAssemblerImpl< P, EG, EV, M >;

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Defines the general interface of the local assembler
 *        classes for the assembly of the interaction volume-local
 *        transmissibility matrix. Specializations have to be provided
 *        for the available interaction volume implementations. these
 *        should derive from this base clases.
 *
 * \tparam P The problem type
 * \tparam EG The element finite volume geometry
 * \tparam EV The element volume variables type
 */
template< class P, class EG, class EV >
class InteractionVolumeAssemblerBase
{
    using Problem = P;
    using FVElementGeometry = EG;
    using ElementVolumeVariables = EV;

    typedef std::true_type yes;
    typedef std::false_type no;

    //! Determines whether or not a matrix has a resize() function
    template<typename M>
    struct matrix_has_resize_method
    {
    private:
        // resize function is called with two indices for matrices
        template<typename U> static auto test(int) -> decltype(std::declval<U>().resize(0, 0), yes());
        template<typename> static no test(...);
    public:
        static constexpr bool value = std::is_same<decltype(test<M>(0)), yes>::value;
    };

    //! determines whether or not a vector has a resize() function
    template<typename V>
    struct vector_has_resize_method
    {
    private:
        // resize function is called with one index for vectors
        template<typename U> static auto test(int) -> decltype(std::declval<U>().resize(0), yes());
        template<typename> static no test(...);
    public:
        static constexpr bool value = std::is_same<decltype(test<V>(0)), yes>::value;
    };

 public:
    /*!
     * \brief The constructor.
     *        Sets pointers to the objects required for a subsequent call to assemble().
     *
     * \param problem The problem to be solved (boundary/initial conditions etc.)
     * \param fvGeometry The local view on the finite volume grid geometry
     * \param elemVolVars The local view on the primary/secondary variables
     */
    InteractionVolumeAssemblerBase(const Problem& problem,
                                   const FVElementGeometry& fvGeometry,
                                   const ElementVolumeVariables& elemVolVars)
    {
        problemPtr_ = &problem;
        fvGeometryPtr_ = &fvGeometry;
        elemVolVarsPtr_ = &elemVolVars;
    }

    // return functions to the local views & problem
    const Problem& problem() const { return *problemPtr_; }
    const FVElementGeometry& fvGeometry() const { return *fvGeometryPtr_; }
    const ElementVolumeVariables& elemVolVars() const { return *elemVolVarsPtr_; }

    /*!
     * \brief Assembles the matrices involved in the flux
     *        expressions and the local system of equations
     *        within an mpfa interaction volume.
     *
     * \tparam IV The interaction volume type implementation
     * \tparam TensorFunc Lambda to obtain the tensor w.r.t.
     *                    which the local system is to be solved
     *
     * \param handle The data handle in which the matrices are stored
     * \param iv The interaction volume
     * \param getT Lambda to evaluate the scv-wise tensors
     */
    template< class DataHandle, class IV, class TensorFunc >
    void assembleMatrices(DataHandle& handle, IV& iv, const TensorFunc& getT)
    {
        DUNE_THROW(Dune::NotImplemented, "Implementation does not provide a assembleMatrices() function");
    }

    /*!
     * \brief Assembles the vector of primary (cell) unknowns and (maybe)
     *        Dirichlet boundary conditions within an interaction volume.
     *
     * \tparam IV The interaction volume type implementation
     * \tparam GetU Lambda to obtain the cell unknowns from grid indices
     *
     * \param handle The data handle in which the vector is stored
     * \param iv The mpfa-o interaction volume
     * \param getU Lambda to obtain the desired cell/Dirichlet value from vol vars
     */
    template< class DataHandle, class IV, class GetU >
    void assembleU(DataHandle& handle, const IV& iv, const GetU& getU)
    {
        DUNE_THROW(Dune::NotImplemented, "Implementation does not provide a assemble() function for the cell/Dirichlet unknowns");
    }

    /*!
     * \brief Assembles the gravitational flux contributions on the scvfs within an
     *        interaction volume.
     *
     * \param handle The data handle in which the vector is stored
     * \param iv The mpfa-o interaction volume
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
        resizeVector_(g, iv.numFaces());
        resizeVector_(deltaG, iv.numUnknowns());
        if (isSurfaceGrid)
            resizeVector_(outsideG, iv.numFaces());

        // we require the CA matrix to have the correct size already
        assert(CA.rows() == iv.numFaces() && CA.cols() == iv.numUnknowns());

        //! For each face, we...
        //! - arithmetically average the phase densities
        //! - compute the term \f$ \alpha := \mathbf{A} \rho \ \mathbf{n}^T \mathbf{K} \mathbf{g} \f$ in each neighboring cell
        //! - compute \f$ \alpha^* = \sum{\alpha_{outside, i}} - \alpha_{inside} \f$
        using Scalar = typename IV::Traits::MatVecTraits::TMatrix::value_type;
        using LocalIndexType = typename IV::Traits::IndexSet::LocalIndexType;

        for (LocalIndexType faceIdx = 0; faceIdx < iv.numFaces(); ++faceIdx)
        {
            // gravitational acceleration on this face
            const auto& curLocalScvf = iv.localScvf(faceIdx);
            const auto& curGlobalScvf = fvGeometry().scvf(curLocalScvf.gridScvfIndex());
            const auto& gravity = problem().gravityAtPos(curGlobalScvf.ipGlobal());

            // get permeability tensor in "positive" sub volume
            const auto& neighborScvIndices = curLocalScvf.neighboringLocalScvIndices();
            const auto& posGlobalScv = fvGeometry().scv(iv.localScv(neighborScvIndices[0]).gridScvIndex());
            const auto& posVolVars = elemVolVars()[posGlobalScv];
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

            if (isSurfaceGrid)
                resizeVector_(outsideG[faceIdx], numOutsideFaces);

            if (!curLocalScvf.isDirichlet())
            {
                const auto localDofIdx = curLocalScvf.localDofIndex();

                rho = getRho(posVolVars);
                deltaG[localDofIdx] = 0.0;

                if (!curGlobalScvf.boundary())
                {
                    for (unsigned int idxInOutside = 0; idxInOutside < curGlobalScvf.numOutsideScvs(); ++idxInOutside)
                    {
                        // obtain outside tensor
                        const auto negLocalScvIdx = neighborScvIndices[idxInOutside+1];
                        const auto& negGlobalScv = fvGeometry().scv(iv.localScv(negLocalScvIdx).gridScvIndex());
                        const auto& negVolVars = elemVolVars()[negGlobalScv];
                        const auto& flipScvf = !isSurfaceGrid ? curGlobalScvf
                                                              : fvGeometry().flipScvf(curGlobalScvf.index(), idxInOutside);

                        alpha_outside[idxInOutside] = negVolVars.extrusionFactor()*vtmv(flipScvf.unitOuterNormal(),
                                                                                        negVolVars.permeability(),
                                                                                        gravity);
                        if (isSurfaceGrid)
                            alpha_outside[idxInOutside] *= -1.0;

                        rho += getRho(negVolVars);
                        deltaG[localDofIdx] += alpha_outside[idxInOutside];
                    }
                }

                rho /= numOutsideFaces + 1;
                deltaG[localDofIdx] -= alpha_inside;
                deltaG[localDofIdx] *= rho*curGlobalScvf.area();
            }
            // use density resulting from Dirichlet BCs
            else
                rho = getRho(elemVolVars()[curGlobalScvf.outsideScvIdx()]);

            // add "inside" & "outside" alphas to gravity containers
            g[faceIdx] = alpha_inside*rho*curGlobalScvf.area();

            if (isSurfaceGrid)
            {
                unsigned int i = 0;
                for (const auto& alpha : alpha_outside)
                    outsideG[faceIdx][i++] = alpha*rho*curGlobalScvf.area();
            }
        }

        // add iv-wide contributions to gravity vectors
        handle.CA().umv(deltaG, g);
        if (isSurfaceGrid)
        {
            using FaceVector = typename IV::Traits::MatVecTraits::FaceVector;
            FaceVector AG;
            resizeVector_(AG, iv.numUnknowns());
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
                const auto& wijk = iv.omegas()[localScvfIdx][idxInOutside+1];

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

protected:
    //! resizes a matrix to the given sizes (specialization for dynamic matrix type)
    template< class Matrix,
              class size_type,
              std::enable_if_t<matrix_has_resize_method<Matrix>::value, int> = 0 >
    void resizeMatrix_(Matrix& M, size_type rows, size_type cols)
    {
        M.resize(rows, cols);
    }

    //! resizes a matrix to the given sizes (specialization for static matrix type - do nothing)
    template< class Matrix,
              class size_type,
              std::enable_if_t<!matrix_has_resize_method<Matrix>::value, int> = 0 >
    void resizeMatrix_(Matrix& M, size_type rows, size_type cols)
    {}

    //! resizes a vector to the given size (specialization for dynamic matrix type)
    template< class Vector,
              class size_type,
              std::enable_if_t<vector_has_resize_method<Vector>::value, int> = 0 >
    void resizeVector_(Vector& v, size_type size)
    {
        v.resize(size);
    }

    //! resizes a vector to the given size (specialization for static vector type - do nothing)
    template< class Vector,
              class size_type,
              std::enable_if_t<!vector_has_resize_method<Vector>::value, int> = 0 >
    void resizeVector_(Vector& v, size_type rows)
    {}

private:
    // pointers to the data required for assembly
    const Problem* problemPtr_;
    const FVElementGeometry* fvGeometryPtr_;
    const ElementVolumeVariables* elemVolVarsPtr_;
};

} // end namespace Dumux

//! include all specializations for different mpfa schemes
#include <dumux/discretization/cellcentered/mpfa/omethod/localassembler.hh>

#endif
