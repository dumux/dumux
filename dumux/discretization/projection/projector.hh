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
 * \ingroup Discretization
 * \brief Contains functionality for L2-projections from one function
 *        space into another, which can live both on the same or different
 *        grids of potentially different dimensionality.
 * \note In the case of the function space bases living on the same grid
 *       an optimized implementation could be implemented. Currently, we
 *       require a glue object containing the intersections between two
 *       grids to be provided to the free factory functions.
 */
#ifndef DUMUX_DISCRETIZATION_PROJECTOR_HH
#define DUMUX_DISCRETIZATION_PROJECTOR_HH

#include <algorithm>
#include <string>
#include <utility>
#include <type_traits>

#include <dune/common/timer.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/promotiontraits.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/jacobianpattern.hh>
#include <dumux/discretization/functionspacebasis.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Does an L2-projection from one discrete function space
 *        into another. The convenience functions makeProjectorPair
 *        or makeProjector can be used to create such a projection.
 */
template<class ScalarType>
class Projector
{
    using MatrixBlockType = Dune::FieldMatrix<ScalarType, 1, 1>;

public:
    //! Export the scalar type
    using Scalar = ScalarType;
    //! Export the type of the projection matrices
    using Matrix = Dune::BCRSMatrix< MatrixBlockType >;

    //! Parameters that can be passed to project()
    struct Params
    {
        std::size_t maxIterations{100};
        Scalar residualReduction{1e-13};
    };

    //! delete default constructor
    Projector() = delete;

    /*!
     * \brief Constructor. Receives the mass and projection
     *        matrix that define the linear system describing
     *        the L2-projection from a function space into another.
     */
    Projector(Matrix&& massMatrix, Matrix&& projectionMatrix)
    : massMat_(std::move(massMatrix))
    , projMat_(std::move(projectionMatrix))
    {}

    /*!
     * \brief Project a solution u into up
     * \param u The solution living on the domain space
     * \return The projection of u into the target space
     */
    template< class BlockType, std::enable_if_t<std::is_convertible<BlockType, ScalarType>::value, int> = 0 >
    Dune::BlockVector<BlockType> project(const Dune::BlockVector<BlockType>& u, const Params& params = Params{}) const
    {
        // be picky about size of u
        if ( u.size() != projMat_.M())
            DUNE_THROW(Dune::InvalidStateException, "Vector size mismatch" );

        Dune::BlockVector<BlockType> up(massMat_.N());

        auto rhs = up;
        projMat_.mv(u, rhs);

        SSORCGBackend solver;
        solver.setResidualReduction(params.residualReduction);
        solver.setMaxIter(params.maxIterations);
        solver.solve(massMat_, up, rhs);

        return up;
    }

    /*!
     * \brief Project a solution u into up
     * \param u The solution living on the domain space
     * \return The projection of u into the target space
     */
    template< class BlockType, std::enable_if_t<!std::is_convertible<BlockType, ScalarType>::value, int> = 0 >
    Dune::BlockVector<BlockType> project(const Dune::BlockVector<BlockType>& u, const Params& params = Params{}) const
    {
        Dune::BlockVector<BlockType> result(massMat_.N());

        for (int pvIdx = 0; pvIdx < BlockType::size(); ++pvIdx)
        {
            Dune::BlockVector<Dune::FieldVector<Scalar, 1>> tmp(u.size());
            std::transform(u.begin(), u.end(), tmp.begin(), [pvIdx] (const auto& v) { return v[pvIdx]; });

            const auto p = project(tmp, params);
            for (std::size_t i = 0; i < p.size(); ++i)
                result[i][pvIdx] = p[i];
        }

        return result;
    }

    /*!
     * \brief Returns the default parameters.
     */
    static Params defaultParams()
    { return {}; }

private:
    Matrix massMat_;
    Matrix projMat_;
};

/*!
 * \brief Traits class stating the type of projector between to bases
 */
template<class FEBasisDomain, class FEBasisTarget>
class ProjectorTraits
{
    using FiniteElementDomain = typename FEBasisDomain::LocalView::Tree::FiniteElement;
    using FiniteElementTarget = typename FEBasisTarget::LocalView::Tree::FiniteElement;
    using ScalarDomain = typename FiniteElementDomain::Traits::LocalBasisType::Traits::RangeFieldType;
    using ScalarTarget = typename FiniteElementTarget::Traits::LocalBasisType::Traits::RangeFieldType;
    using Scalar = typename Dune::PromotionTraits<ScalarDomain, ScalarTarget>::PromotedType;
public:
    using Projector = Dumux::Projector< Scalar >;
};


// Projector construction details
namespace Impl {

/*!
 * \brief Creates a projector class between two function space bases
 * \tparam doBidirectional If false, the backward projection matrix is not assembled
 */
template<bool doBidirectional, class FEBasisDomain, class FEBasisTarget, class GlueType>
std::pair< typename ProjectorTraits<FEBasisDomain, FEBasisTarget>::Projector,
           typename ProjectorTraits<FEBasisTarget, FEBasisDomain>::Projector >
makeProjectorPair(const FEBasisDomain& feBasisDomain,
                  const FEBasisTarget& feBasisTarget,
                  const GlueType& glue)
{
    // we assume that target dim <= domain dimension
    static constexpr int domainDim = FEBasisDomain::GridView::dimension;
    static constexpr int targetDim = FEBasisTarget::GridView::dimension;
    static_assert(targetDim <= domainDim, "This expects target dim < domain dim, please swap arguments");

    using ForwardProjector = typename ProjectorTraits<FEBasisDomain, FEBasisTarget>::Projector;
    using BackwardProjector = typename ProjectorTraits<FEBasisTarget, FEBasisDomain>::Projector;

    using ForwardProjectionMatrix = typename ForwardProjector::Matrix;
    using BackwardProjectionMatrix = typename BackwardProjector::Matrix;

    auto domainLocalView = feBasisDomain.localView();
    auto targetLocalView = feBasisTarget.localView();

    // determine mass matrix patterns (standard FE scheme pattern)
    Dune::MatrixIndexSet backwardPatternM, forwardPatternM;
    forwardPatternM = getFEJacobianPattern(feBasisTarget);
    if (doBidirectional) backwardPatternM = getFEJacobianPattern(feBasisDomain);

    // determine projection matrix patterns
    Dune::MatrixIndexSet backwardPatternP, forwardPatternP;
    forwardPatternP.resize(feBasisTarget.size(), feBasisDomain.size());
    if (doBidirectional) backwardPatternP.resize(feBasisDomain.size(), feBasisTarget.size());

    using std::max;
    unsigned int maxBasisOrder = 0;
    for (const auto& is : intersections(glue))
    {
        // "outside" element is of target domain
        // since target dim <= domain dim there is maximum one!
        targetLocalView.bind( is.outside(0) );
        const auto& targetLocalBasis = targetLocalView.tree().finiteElement().localBasis();

        for (unsigned int nIdx = 0; nIdx < is.neighbor(/*targetSide*/1); ++nIdx)
        {
            domainLocalView.bind( is.inside(nIdx) );
            const auto& domainLocalBasis = domainLocalView.tree().finiteElement().localBasis();

            // keep track of maximum basis order (used in integration)
            maxBasisOrder = max(maxBasisOrder, max(domainLocalBasis.order(), targetLocalBasis.order()));

            for (unsigned int i = 0; i < domainLocalBasis.size(); ++i)
                for (unsigned int j = 0; j < targetLocalBasis.size(); ++j)
                {
                    forwardPatternP.add(targetLocalView.index(j), domainLocalView.index(i));
                    if (doBidirectional) backwardPatternP.add(domainLocalView.index(i), targetLocalView.index(j));
                }
        }
    }

    // assemble matrices
    ForwardProjectionMatrix forwardM, forwardP;
    forwardPatternM.exportIdx(forwardM); forwardM = 0.0;
    forwardPatternP.exportIdx(forwardP); forwardP = 0.0;

    BackwardProjectionMatrix backwardM, backwardP;
    if (doBidirectional)
    {
        backwardPatternM.exportIdx(backwardM); backwardM = 0.0;
        backwardPatternP.exportIdx(backwardP); backwardP = 0.0;
    }

    for (const auto& is : intersections(glue))
    {
        const auto& targetElement = is.outside(0);
        const auto& targetElementGeometry = targetElement.geometry();

        targetLocalView.bind( targetElement );
        const auto& targetLocalBasis = targetLocalView.tree().finiteElement().localBasis();

        // perform integration over intersection geometry
        using IsGeometry = typename std::decay_t<decltype(is.geometry())>;
        using ctype = typename IsGeometry::ctype;

        const auto& isGeometry = is.geometry();
        const int intOrder = maxBasisOrder + 1;
        const auto& quad = Dune::QuadratureRules<ctype, IsGeometry::mydimension>::rule(isGeometry.type(), intOrder);
        for (auto&& qp : quad)
        {
            const auto weight = qp.weight();
            const auto ie = isGeometry.integrationElement(qp.position());
            const auto globalPos = isGeometry.global(qp.position());

            std::vector< Dune::FieldVector<ctype, 1> > targetShapeVals;
            targetLocalBasis.evaluateFunction(targetElementGeometry.local(globalPos), targetShapeVals);

            // mass matrix entries target domain
            for (unsigned int i = 0; i < targetLocalBasis.size(); ++i)
            {
                const auto dofIdxI = targetLocalView.index(i);
                forwardM[dofIdxI][dofIdxI] += ie*weight*targetShapeVals[i]*targetShapeVals[i];

                for (unsigned int j = i+1; j < targetLocalBasis.size(); ++j)
                {
                    const auto dofIdxJ = targetLocalView.index(j);
                    const auto value = ie*weight*targetShapeVals[i]*targetShapeVals[j];
                    forwardM[dofIdxI][dofIdxJ] += value;
                    forwardM[dofIdxJ][dofIdxI] += value;
                }
            }

            // If targetDim < domainDim, there can be several "neighbors" if
            // targetElement is aligned with a facet of domainElement. In
            // this case make sure the basis functions are not added
            // multiple times! (division by numNeighbors)
            const auto numNeighbors = is.neighbor(/*targetSide*/1);
            for (unsigned int nIdx = 0; nIdx < numNeighbors; ++nIdx)
            {
                const auto& domainElement = is.inside(nIdx);
                domainLocalView.bind( domainElement );
                const auto& domainLocalBasis = domainLocalView.tree().finiteElement().localBasis();

                std::vector< Dune::FieldVector<ctype, 1> > domainShapeVals;
                domainLocalBasis.evaluateFunction(domainElement.geometry().local(globalPos), domainShapeVals);

                // add entries in matrices
                for (unsigned int i = 0; i < domainLocalBasis.size(); ++i)
                {
                    const auto dofIdxDomain = domainLocalView.index(i);
                    const auto domainShapeVal = domainShapeVals[i];
                    if (doBidirectional)
                    {
                        backwardM[dofIdxDomain][dofIdxDomain] += ie*weight*domainShapeVal*domainShapeVal;

                        for (unsigned int j = i+1; j < domainLocalBasis.size(); ++j)
                        {
                            const auto dofIdxDomainJ = domainLocalView.index(j);
                            const auto value = ie*weight*domainShapeVal*domainShapeVals[j];
                            backwardM[dofIdxDomain][dofIdxDomainJ] += value;
                            backwardM[dofIdxDomainJ][dofIdxDomain] += value;
                        }
                    }

                    for (unsigned int j = 0; j < targetLocalBasis.size(); ++j)
                    {
                        const auto dofIdxTarget = targetLocalView.index(j);
                        const auto entry = ie*weight*domainShapeVal*targetShapeVals[j];

                        forwardP[dofIdxTarget][dofIdxDomain] += entry/numNeighbors;
                        if (doBidirectional)
                            backwardP[dofIdxDomain][dofIdxTarget] += entry;
                    }
                }
            }
        }
    }

    // On those dofs that to not take part in any intersection we will have zeroes
    // in the mass matrices. The right hand side will be zero because the projection
    // matrix has no entries for those dofs as there was no intersection. Thus, we set
    // 1.0 on the diagonal for those dofs, such that after projection, the projected
    // solution is 0 on them.
    for (std::size_t dofIdx = 0; dofIdx < forwardM.N(); ++dofIdx)
        if (forwardM[dofIdx][dofIdx] == 0.0)
            forwardM[dofIdx][dofIdx] = 1.0;

    if (doBidirectional)
    {
        for (std::size_t dofIdx = 0; dofIdx < backwardM.N(); ++dofIdx)
            if (backwardM[dofIdx][dofIdx] == 0.0)
                backwardM[dofIdx][dofIdx] = 1.0;
    }

    ForwardProjector fw(std::move(forwardM), std::move(forwardP));
    BackwardProjector bw(std::move(backwardM), std::move(backwardP));

    return std::make_pair(std::move(fw), std::move(bw));
}

} // end namespace Implementation


/*!
 * \ingroup Discretization
 * \brief Creates a pair of projectors between the space with
 *        basis feBasisDomain to the space with basis feBasisTarget.
 * \param feBasisDomain The domain finite element space basis
 * \param feBasisTarget The target finite element space basis
 * \param glue The glue object containing the intersections between the grids.
 * \return A pair of projectors where the first is the forward
 *         projector from the space with basis feBasisDomain to
 *         the space with basis feBasisTarget and the second
 *         does the backward projection.
 */
template< class FEBasisDomain, class FEBasisTarget, class GlueType >
std::pair< typename ProjectorTraits<FEBasisDomain, FEBasisTarget>::Projector,
           typename ProjectorTraits<FEBasisTarget, FEBasisDomain>::Projector >
makeProjectorPair(const FEBasisDomain& feBasisDomain,
                  const FEBasisTarget& feBasisTarget,
                  GlueType glue)
{
    // we assume that target dim <= domain dimension
    static constexpr int domainDim = FEBasisDomain::GridView::dimension;
    static constexpr int targetDim = FEBasisTarget::GridView::dimension;
    static_assert(targetDim <= domainDim, "makeProjectorPair() expects targetDim < domainDim, please swap arguments");

    return Impl::makeProjectorPair<true>(feBasisDomain, feBasisTarget, glue);
}

/*!
 * \ingroup Discretization
 * \brief Creates a forward projector from the space feBasisDomain
 *        to the space with basis feBasisTarget.
 * \param feBasisDomain The domain finite element space basis
 * \param feBasisTarget The target finite element space basis
 * \param glue The glue object containing the intersections between the grids.
 * \return The forward projector from the space with basis feBasisDomain
 *         to the space with basis feBasisTarget.
 */
template< class FEBasisDomain, class FEBasisTarget, class GlueType >
typename ProjectorTraits<FEBasisDomain, FEBasisTarget>::Projector
makeProjector(const FEBasisDomain& feBasisDomain,
              const FEBasisTarget& feBasisTarget,
              GlueType glue)
{
    // we assume that target dim <= domain dimension
    static constexpr int domainDim = FEBasisDomain::GridView::dimension;
    static constexpr int targetDim = FEBasisTarget::GridView::dimension;
    static_assert(targetDim <= domainDim, "makeProjectorPair() expects targetDim < domainDim, please swap arguments");

    return Impl::makeProjectorPair<false>(feBasisDomain, feBasisTarget, glue).first;
}

} // end namespace Dumux

#endif // DUMUX_PROJECTOR_HH
