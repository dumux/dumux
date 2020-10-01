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
    , numDofsTarget_(massMat_.N())
    {
        if (massMat_.N() != projMat_.N())
            DUNE_THROW(Dune::InvalidStateException, "Matrix row size mismatch: " << massMat_.N() << " vs " << projMat_.N());
    }

    /*!
     * \brief Constructor for projection into a target space that occupies
     *        a larger geometric region than the domain space. In this case,
     *        the mass matrix can be chosen such that is is solved for only
     *        those dofs which will be populated with values from the projection.
     *        This requires an additional index map that maps the entries of the
     *        projected solution into the solution vector for the target space.
     *        Furthermore, the number of degrees of freedom must be specified to
     *        set up the coefficient vector with correct size for the target space.
     */
    Projector(Matrix&& massMatrix,
              Matrix&& projectionMatrix,
              std::vector<std::size_t>&& indexMap,
              std::size_t numDofsTarget)
    : massMat_(std::move(massMatrix))
    , projMat_(std::move(projectionMatrix))
    , indexMapTarget_(std::move(indexMap))
    , numDofsTarget_(numDofsTarget)
    {
        if (indexMapTarget_.size() != massMat_.N())
            DUNE_THROW(Dune::InvalidStateException, "Target index map size mismatch: " << indexMapTarget_.size() << " vs " << massMat_.N());

        if (massMat_.N() != projMat_.N())
            DUNE_THROW(Dune::InvalidStateException, "Matrix row size mismatch: " << massMat_.N() << " vs " << projMat_.N());

        if (*std::max_element(indexMapTarget_.begin(), indexMapTarget_.end()) > numDofsTarget_)
            DUNE_THROW(Dune::InvalidStateException, "Index map exceeds provided number of dofs in target domain!");
    }

    /*!
     * \brief Project a solution u into up
     * \param u The solution living on the domain space
     * \param params Optional parameters for mass matrix solve
     * \return The projection of u into the target space
     */
    template< class BlockType, std::enable_if_t<std::is_convertible<BlockType, ScalarType>::value, int> = 0 >
    Dune::BlockVector<BlockType> project(const Dune::BlockVector<BlockType>& u, const Params& params = Params{}) const
    {
        // be picky about size of u
        if ( u.size() != projMat_.M())
            DUNE_THROW(Dune::InvalidStateException, "Vector size mismatch");

        Dune::BlockVector<BlockType> up(massMat_.N());

        auto rhs = up;
        projMat_.mv(u, rhs);

        SSORCGBackend solver;
        solver.setResidualReduction(params.residualReduction);
        solver.setMaxIter(params.maxIterations);
        solver.solve(massMat_, up, rhs);

        // if target space occupies a larger region, fill missing entries with zero
        if (!indexMapTarget_.empty())
        {
            Dune::BlockVector<BlockType> result(numDofsTarget_);

            result = 0.0;
            for (std::size_t i = 0; i < indexMapTarget_.size(); ++i)
                result[indexMapTarget_[i]] = up[i];

            return result;
        }

        return up;
    }

    /*!
     * \brief Project a solution u into up
     * \param u The solution living on the domain space
     * \param params Optional parameters for mass matrix solve
     * \return The projection of u into the target space
     */
    template< class BlockType, std::enable_if_t<!std::is_convertible<BlockType, ScalarType>::value, int> = 0 >
    Dune::BlockVector<BlockType> project(const Dune::BlockVector<BlockType>& u, const Params& params = Params{}) const
    {
        Dune::BlockVector<BlockType> result(numDofsTarget_);

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

    std::vector<std::size_t> indexMapTarget_;
    std::size_t numDofsTarget_;
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
namespace Detail {

/*!
 * \brief Reduces a mass matrix and projection matrix such that they are composed
 *        of only those dofs that actually take part in the projection. Simultaneously,
 *        a container with the index map into the complete target space is filled so that
 *        the entries after projection can be assigned to the corresponding dof in the
 *        overall target space.
 */
template<class Matrix>
void setupReducedMatrices(const Matrix& massMatrix, const Matrix& projMatrix, const std::vector<bool>& dofIsVoid,
                          Matrix& reducedM, Matrix& reducedP, std::vector<std::size_t>& expansionMap)
{
    const std::size_t numNonVoidDofs = std::count_if(dofIsVoid.begin(), dofIsVoid.end(), [] (bool v) { return !v; });

    // reduce matrices to only dofs that take part and create index map
    std::vector<std::size_t> reductionMap(massMatrix.N());
    expansionMap.resize(numNonVoidDofs);

    std::size_t idxInReducedSpace = 0;
    for (std::size_t dofIdx = 0; dofIdx < dofIsVoid.size(); ++dofIdx)
        if (!dofIsVoid[dofIdx])
        {
            reductionMap[dofIdx] = idxInReducedSpace;
            expansionMap[idxInReducedSpace] = dofIdx;
            idxInReducedSpace++;
        }

    // create reduced M/P matrix
    Dune::MatrixIndexSet patternMReduced, patternPReduced;
    patternMReduced.resize(numNonVoidDofs, numNonVoidDofs);
    patternPReduced.resize(numNonVoidDofs, projMatrix.M());
    for (auto rowIt = massMatrix.begin(); rowIt != massMatrix.end(); ++rowIt)
        if (!dofIsVoid[rowIt.index()])
        {
            const auto reducedRowIdx = reductionMap[rowIt.index()];
            for (auto colIt = (*rowIt).begin(); colIt != (*rowIt).end(); ++colIt)
                if (!dofIsVoid[colIt.index()])
                    patternMReduced.add(reducedRowIdx, reductionMap[colIt.index()]);
        }

    for (auto rowIt = projMatrix.begin(); rowIt != projMatrix.end(); ++rowIt)
        if (!dofIsVoid[rowIt.index()])
        {
            const auto reducedRowIdx = reductionMap[rowIt.index()];
            for (auto colIt = (*rowIt).begin(); colIt != (*rowIt).end(); ++colIt)
                patternPReduced.add(reducedRowIdx, colIt.index());
        }

    patternMReduced.exportIdx(reducedM);
    patternPReduced.exportIdx(reducedP);

    // fill matrix entries
    for (auto rowIt = massMatrix.begin(); rowIt != massMatrix.end(); ++rowIt)
        if (!dofIsVoid[rowIt.index()])
        {
            const auto reducedRowIdx = reductionMap[rowIt.index()];
            for (auto colIt = (*rowIt).begin(); colIt != (*rowIt).end(); ++colIt)
                if (!dofIsVoid[colIt.index()])
                    reducedM[reducedRowIdx][reductionMap[colIt.index()]] = *colIt;
        }

    for (auto rowIt = projMatrix.begin(); rowIt != projMatrix.end(); ++rowIt)
        if (!dofIsVoid[rowIt.index()])
        {
            const auto reducedRowIdx = reductionMap[rowIt.index()];
            for (auto colIt = (*rowIt).begin(); colIt != (*rowIt).end(); ++colIt)
                reducedP[reducedRowIdx][colIt.index()] = *colIt;
        }
}

/*!
 * \brief Creates the matrices underlying l2-projections
 * \tparam doBidirectional If false, the backward projection matrix is not assembled
 * \param feBasisDomain The basis to the domain finite element space
 * \param feBasisTarget The basis to the target finite element space
 * \param glue The glue object containing the intersections between the two grids
 * \param treatDiagonalZeroes If true, zero entries on the diagonal of the matrices
 *        that appear if the two domains occupy different geometric regions (and some
 *        dofs to not take part in the projection as a result) are substituted by ones.
 *        This substitution will lead to those dofs being mapped to zeroes in the target space.
 * \returns An std::pair of projection matrices, where the first entry stores the
 *          matrices of the forward projection and the second entry stores those
 *          of the backward projection. The entries of the returned pair are itself
 *          std::pairs which store the mass matrix in the first and the projection
 *          matrix in the second entry.
 */
template<bool doBidirectional, class FEBasisDomain, class FEBasisTarget, class GlueType>
auto createProjectionMatrices(const FEBasisDomain& feBasisDomain,
                              const FEBasisTarget& feBasisTarget,
                              const GlueType& glue,
                              bool treatDiagonalZeroes = true)
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
        // since target dim <= domain dim there is maximum one!
        targetLocalView.bind( is.targetEntity(0) );
        const auto& targetLocalBasis = targetLocalView.tree().finiteElement().localBasis();

        for (unsigned int nIdx = 0; nIdx < is.numDomainNeighbors(); ++nIdx)
        {
            domainLocalView.bind( is.domainEntity(nIdx) );
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
        const auto& targetElement = is.targetEntity(0);
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
            const auto numNeighbors = is.numDomainNeighbors();
            for (unsigned int nIdx = 0; nIdx < numNeighbors; ++nIdx)
            {
                const auto& domainElement = is.domainEntity(nIdx);
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

    // maybe treat zeroes on the diagonal
    if (treatDiagonalZeroes)
    {
        for (std::size_t dofIdxTarget = 0; dofIdxTarget < forwardM.N(); ++dofIdxTarget)
            if (forwardM[dofIdxTarget][dofIdxTarget] == 0.0)
                forwardM[dofIdxTarget][dofIdxTarget] = 1.0;

        if (doBidirectional)
        {
            for (std::size_t dofIdxDomain = 0; dofIdxDomain < backwardM.N(); ++dofIdxDomain)
                if (backwardM[dofIdxDomain][dofIdxDomain] == 0.0)
                    backwardM[dofIdxDomain][dofIdxDomain] = 1.0;
        }
    }

    return std::make_pair( std::make_pair(std::move(forwardM), std::move(forwardP)),
                           std::make_pair(std::move(backwardM), std::move(backwardP)) );
}

/*!
 * \brief Creates a projector class between two function space bases
 * \tparam doBidirectional If false, the backward projection matrix is not assembled
 * \returns an std::pair with the forward and backward projector
 */
template<bool doBidirectional, class FEBasisDomain, class FEBasisTarget, class GlueType>
auto makeProjectorPair(const FEBasisDomain& feBasisDomain,
                       const FEBasisTarget& feBasisTarget,
                       const GlueType& glue)
{
    using ForwardProjector = typename ProjectorTraits<FEBasisDomain, FEBasisTarget>::Projector;
    using BackwardProjector = typename ProjectorTraits<FEBasisTarget, FEBasisDomain>::Projector;

    using ForwardProjectionMatrix = typename ForwardProjector::Matrix;
    using BackwardProjectionMatrix = typename BackwardProjector::Matrix;

    auto projectionMatrices = createProjectionMatrices<doBidirectional>(feBasisDomain, feBasisTarget, glue, false);
    auto& forwardMatrices = projectionMatrices.first;
    auto& backwardMatrices = projectionMatrices.second;

    auto& forwardM = forwardMatrices.first;
    auto& forwardP = forwardMatrices.second;

    auto& backwardM = backwardMatrices.first;
    auto& backwardP = backwardMatrices.second;

    // determine the dofs that do not take part in intersections
    std::vector<bool> isVoidTarget(forwardM.N(), false);
    for (std::size_t dofIdxTarget = 0; dofIdxTarget < forwardM.N(); ++dofIdxTarget)
        if (forwardM[dofIdxTarget][dofIdxTarget] == 0.0)
            isVoidTarget[dofIdxTarget] = true;

    std::vector<bool> isVoidDomain;
    if (doBidirectional)
    {
        isVoidDomain.resize(backwardM.N(), false);
        for (std::size_t dofIdxDomain = 0; dofIdxDomain < backwardM.N(); ++dofIdxDomain)
            if (backwardM[dofIdxDomain][dofIdxDomain] == 0.0)
                isVoidDomain[dofIdxDomain] = true;
    }

    const bool hasVoidTarget = std::any_of(isVoidTarget.begin(), isVoidTarget.end(), [] (bool v) { return v; });
    const bool hasVoidDomain = std::any_of(isVoidDomain.begin(), isVoidDomain.end(), [] (bool v) { return v; });
    if (!hasVoidDomain && !hasVoidTarget)
    {
        return std::make_pair(ForwardProjector(std::move(forwardM), std::move(forwardP)),
                              BackwardProjector(std::move(backwardM), std::move(backwardP)));
    }
    else if (!hasVoidDomain && hasVoidTarget)
    {
        std::vector<std::size_t> expansionMapTarget;
        ForwardProjectionMatrix forwardMReduced, forwardPReduced;
        setupReducedMatrices(forwardM, forwardP, isVoidTarget,
                             forwardMReduced, forwardPReduced, expansionMapTarget);

        return std::make_pair( ForwardProjector(std::move(forwardMReduced),
                                                std::move(forwardPReduced),
                                                std::move(expansionMapTarget),
                                                forwardM.N()),
                               BackwardProjector(std::move(backwardM), std::move(backwardP)) );
    }
    else if (hasVoidDomain && !hasVoidTarget)
    {
        if (doBidirectional)
        {
            std::vector<std::size_t> expansionMapDomain;
            BackwardProjectionMatrix backwardMReduced, backwardPReduced;
            setupReducedMatrices(backwardM, backwardP, isVoidDomain,
                                 backwardMReduced, backwardPReduced, expansionMapDomain);

            return std::make_pair( ForwardProjector(std::move(forwardM), std::move(forwardP)),
                                   BackwardProjector(std::move(backwardMReduced),
                                                     std::move(backwardPReduced),
                                                     std::move(expansionMapDomain),
                                                     backwardM.N()) );
        }
        else
            return std::make_pair( ForwardProjector(std::move(forwardM), std::move(forwardP)),
                                   BackwardProjector(std::move(backwardM), std::move(backwardP)) );
    }
    else
    {
        std::vector<std::size_t> expansionMapTarget;
        ForwardProjectionMatrix forwardMReduced, forwardPReduced;
        setupReducedMatrices(forwardM, forwardP, isVoidTarget,
                             forwardMReduced, forwardPReduced, expansionMapTarget);

        if (doBidirectional)
        {
            std::vector<std::size_t> expansionMapDomain;
            BackwardProjectionMatrix backwardMReduced, backwardPReduced;
            setupReducedMatrices(backwardM, backwardP, isVoidDomain,
                                 backwardMReduced, backwardPReduced, expansionMapDomain);

            return std::make_pair( ForwardProjector(std::move(forwardMReduced),
                                                    std::move(forwardPReduced),
                                                    std::move(expansionMapTarget),
                                                    forwardM.N()),
                                   BackwardProjector(std::move(backwardMReduced),
                                                     std::move(backwardPReduced),
                                                     std::move(expansionMapDomain),
                                                     backwardM.N()) );
        }
        else
            return std::make_pair( ForwardProjector(std::move(forwardMReduced),
                                                    std::move(forwardPReduced),
                                                    std::move(expansionMapTarget),
                                                    forwardM.N()),
                                   BackwardProjector(std::move(backwardM), std::move(backwardP)) );
    }
}

} // end namespace Detail


/*!
 * \ingroup Discretization
 * \brief Creates a pair of projectors between the space with
 *        basis feBasisDomain to the space with basis feBasisTarget.
 * \param feBasisDomain The domain finite element space basis
 * \param feBasisTarget The target finite element space basis
 * \param glue The glue object containing the intersections between the grids.
 * \return An std::pair of projectors where the first is the forward
 *         projector from the space with basis feBasisDomain to
 *         the space with basis feBasisTarget and the second
 *         does the backward projection.
 */
template< class FEBasisDomain, class FEBasisTarget, class GlueType >
auto makeProjectorPair(const FEBasisDomain& feBasisDomain,
                       const FEBasisTarget& feBasisTarget,
                       GlueType glue)
{
    // we assume that target dim <= domain dimension
    static constexpr int domainDim = FEBasisDomain::GridView::dimension;
    static constexpr int targetDim = FEBasisTarget::GridView::dimension;
    static_assert(targetDim <= domainDim, "makeProjectorPair() expects targetDim < domainDim, please swap arguments");

    return Detail::makeProjectorPair<true>(feBasisDomain, feBasisTarget, glue);
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
auto makeProjector(const FEBasisDomain& feBasisDomain,
                   const FEBasisTarget& feBasisTarget,
                   GlueType glue)
{
    // we assume that target dim <= domain dimension
    static constexpr int domainDim = FEBasisDomain::GridView::dimension;
    static constexpr int targetDim = FEBasisTarget::GridView::dimension;
    static_assert(targetDim <= domainDim, "makeProjectorPair() expects targetDim < domainDim, please swap arguments");

    return Detail::makeProjectorPair<false>(feBasisDomain, feBasisTarget, glue).first;
}

/*!
 * \brief Creates the matrices underlying l2-projections
 * \param feBasisDomain The basis to the domain finite element space
 * \param feBasisTarget The basis to the target finite element space
 * \param glue The glue object containing the intersections between the two grids
 * \returns An std::pair of projection matrices, where the first entry stores the
 *          matrices of the forward projection and the second entry stores those
 *          of the backward projection. The entries of the returned pair are itself
 *          std::pairs which store the mass matrix in the first and the projection
 *          matrix in the second entry.
 */
template< class FEBasisDomain, class FEBasisTarget, class GlueType >
auto makeProjectionMatricesPair(const FEBasisDomain& feBasisDomain,
                                const FEBasisTarget& feBasisTarget,
                                GlueType glue)
{
    // we assume that target dim <= domain dimension
    static constexpr int domainDim = FEBasisDomain::GridView::dimension;
    static constexpr int targetDim = FEBasisTarget::GridView::dimension;
    static_assert(targetDim <= domainDim, "makeProjectionMatrixPair() expects targetDim < domainDim, please swap arguments");

    return Detail::createProjectionMatrices<true>(feBasisDomain, feBasisTarget, glue);
}

/*!
 * \brief Creates the matrices underlying l2-projections
 * \param feBasisDomain The basis to the domain finite element space
 * \param feBasisTarget The basis to the target finite element space
 * \param glue The glue object containing the intersections between the two grids
 * \returns An std::pair of matrices, which store the mass matrix in the first and the
 *          projection matrix in the second entry.
 */
template< class FEBasisDomain, class FEBasisTarget, class GlueType >
auto makeProjectionMatrices(const FEBasisDomain& feBasisDomain,
                            const FEBasisTarget& feBasisTarget,
                            GlueType glue)
{
    // we assume that target dim <= domain dimension
    static constexpr int domainDim = FEBasisDomain::GridView::dimension;
    static constexpr int targetDim = FEBasisTarget::GridView::dimension;
    static_assert(targetDim <= domainDim, "makeProjectionMatrixPair() expects targetDim < domainDim, please swap arguments");

    return Detail::createProjectionMatrices<false>(feBasisDomain, feBasisTarget, glue).first;
}

} // end namespace Dumux

#endif
