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
 * \brief L2-projections of analytic functions into a given function space
 */
#ifndef DUMUX_DISCRETIZATION_L2_PROJECTION_HH
#define DUMUX_DISCRETIZATION_L2_PROJECTION_HH
#if HAVE_DUNE_FUNCTIONS

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/functions/gridviewfunction.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/jacobianpattern.hh>

namespace Dumux {


template <class FEBasis>
class L2Projection
{
    static constexpr int dim = FEBasis::GridView::dimension;
    using FiniteElement = typename FEBasis::LocalView::Tree::FiniteElement;
    using Scalar = typename FiniteElement::Traits::LocalBasisType::Traits::RangeFieldType;
    using ShapeValue = typename FiniteElement::Traits::LocalBasisType::Traits::RangeType;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, 1, 1>>;
public:
    using CoefficientVector = Dune::BlockVector<Dune::FieldVector<Scalar, 1>>;

    //! Parameters that can be passed to project()
    struct Params
    {
        std::size_t maxIterations{100};
        Scalar residualReduction{1e-13};
    };

    L2Projection(const FEBasis& feBasis)
    : feBasis_(feBasis)
    , massMatrix_(createMassMatrix_(feBasis))
    {}

    template <class Function>
    CoefficientVector project(Function&& function, const Params& params = Params{})
    {
        CoefficientVector projection, rhs;
        projection.resize(feBasis_.size());
        rhs.resize(feBasis_.size());
        rhs = 0.0;

        // assemble right hand size
        auto localFunc = localFunction(Dune::Functions::makeGridViewFunction(function, feBasis_.gridView()));
        auto localView = feBasis_.localView();
        for (const auto& element : elements(feBasis_.gridView()))
        {
            localView.bind(element);
            localFunc.bind(element);

            const auto& localFiniteElement = localView.tree().finiteElement();
            const int order = dim*localFiniteElement.localBasis().order();
            const auto& quad = Dune::QuadratureRules<Scalar, dim>::rule(element.type(), order);
            const auto geometry = element.geometry();

            for (auto&& qp : quad)
            {
                const auto weight = qp.weight();
                const auto ie = geometry.integrationElement(qp.position());
                const auto globalPos = geometry.global(qp.position());

                std::vector<ShapeValue> shapeValues;
                localFiniteElement.localBasis().evaluateFunction(geometry.local(globalPos), shapeValues);
                const auto functionValue = localFunc(qp.position());

                for (int i = 0; i < localFiniteElement.localBasis().size(); ++i)
                {
                    const auto globalI = localView.index(i);
                    rhs[globalI] += ie*weight*shapeValues[i]*functionValue;
                }
            }
        }

        // solve projection
        SSORCGBackend solver;
        solver.setResidualReduction(params.residualReduction);
        solver.setMaxIter(params.maxIterations);
        solver.solve(massMatrix_, projection, rhs);

        return projection;
    }

private:
    Matrix createMassMatrix_(const FEBasis& feBasis) const
    {
        Matrix massMatrix;

        auto pattern = getFEJacobianPattern(feBasis);
        pattern.exportIdx(massMatrix);

        auto localView = feBasis.localView();
        for (const auto& element : elements(feBasis.gridView()))
        {
            localView.bind(element);

            const auto& localFiniteElement = localView.tree().finiteElement();
            const int order = 2*(dim*localFiniteElement.localBasis().order()-1);
            const auto& quad = Dune::QuadratureRules<Scalar, dim>::rule(element.type(), order);
            const auto geometry = element.geometry();

            for (auto&& qp : quad)
            {
                const auto weight = qp.weight();
                const auto ie = geometry.integrationElement(qp.position());
                const auto globalPos = geometry.global(qp.position());

                std::vector<ShapeValue> shapeValues;
                localFiniteElement.localBasis().evaluateFunction(geometry.local(globalPos), shapeValues);

                for (int i = 0; i < localFiniteElement.localBasis().size(); ++i)
                {

                    const auto globalI = localView.index(i);
                    massMatrix[globalI][globalI] += ie*weight*shapeValues[i]*shapeValues[i];

                    for (int j = i+1; j < localFiniteElement.localBasis().size(); ++j)
                    {
                        const auto globalJ = localView.index(j);
                        const auto value = ie*weight*shapeValues[i]*shapeValues[j];
                        massMatrix[globalI][globalJ] += value;
                        massMatrix[globalJ][globalI] += value;
                    }
                }
            }
        }

        return massMatrix;
    }

    const FEBasis& feBasis_;
    const Matrix massMatrix_;
};

} // end namespace Dumux

#endif // HAVE_DUNE_FUNCTIONS
#endif
