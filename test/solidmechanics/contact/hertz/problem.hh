// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Two elastic half disks pressed together by rigid frictionless platens
 *
 * The upper half disk (centre (0, R)) and the lower half disk (centre (0, -R)) touch at the
 * origin. Their flat edges y = R and y = -R are moved towards each other by the approach
 * \f$ \Delta \f$ each and are free to slide; the axis x = 0 is a symmetry line. The material is
 * neo-Hookean, \f$ \mathbf{P} = \mu(\mathbf{F} - \mathbf{F}^{-T}) + \lambda \ln J\, \mathbf{F}^{-T} \f$,
 * with the Lame parameters of the Young's modulus and Poisson ratio.
 */
#ifndef DUMUX_TEST_CONTACT_HERTZ_PROBLEM_HH
#define DUMUX_TEST_CONTACT_HERTZ_PROBLEM_HH

#include <cmath>
#include <vector>

#include <dune/common/fmatrix.hh>

#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

template<class TypeTag>
class HertzProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using GridIndex = typename IndexTraits<GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<Dumux::DirichletConstraintInfo<numEq>, PrimaryVariables, GridIndex>;
    static constexpr int dim = GridView::dimension;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    HertzProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , radius_(getParam<Scalar>("Problem.Radius"))
    , youngsModulus_(getParam<Scalar>("Problem.YoungsModulus"))
    , poissonRatio_(getParam<Scalar>("Problem.PoissonRatio"))
    {
        const Scalar eps = 1e-8*radius_;
        std::vector<bool> visited(gridGeometry->numDofs(), false);
        auto elemDisc = localView(*gridGeometry);
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            elemDisc.bind(element);
            for (const auto& localDof : localDofs(elemDisc))
            {
                const GridIndex dofIdx = localDof.dofIndex();
                if (visited[dofIdx])
                    continue;
                visited[dofIdx] = true;
                const auto pos = ipData(elemDisc, localDof).global();
                const int side = pos[1] > radius_ - eps ? 1 : (pos[1] < -radius_ + eps ? -1 : 0);
                const bool onAxis = std::abs(pos[0]) < eps;
                if (side != 0 || onAxis)
                    constrainedDofs_.push_back({dofIdx, side, onAxis});
            }
        }
        setApproach(0.0);
    }

    //! Move both platens towards the origin by delta
    void setApproach(Scalar delta)
    {
        constraints_.clear();
        for (const auto& dof : constrainedDofs_)
        {
            typename DirichletConstraintData::ConstraintInfo info;
            PrimaryVariables values(0.0);
            if (dof.onAxis)
                info.set(0);
            if (dof.side != 0)
            {
                info.set(1);
                values[1] = -dof.side*delta;
            }
            constraints_.push_back(DirichletConstraintData{std::move(info), values, dof.index});
        }
    }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    Tensor firstPiolaKirchhoffStressTensor(const Tensor& F) const
    {
        const Scalar mu = youngsModulus_/(2.0*(1.0 + poissonRatio_));
        const Scalar lambda = youngsModulus_*poissonRatio_/((1.0 + poissonRatio_)*(1.0 - 2.0*poissonRatio_));
        auto FinvT = F;
        FinvT.invert();
        FinvT = getTransposed(FinvT);
        auto P = F;
        P -= FinvT;
        P *= mu;
        P.axpy(lambda*std::log(F.determinant()), FinvT);
        return P;
    }

    Scalar radius() const
    { return radius_; }

    Scalar youngsModulus() const
    { return youngsModulus_; }

    Scalar poissonRatio() const
    { return poissonRatio_; }

private:
    struct ConstrainedDof
    {
        GridIndex index;
        int side;
        bool onAxis;
    };

    Scalar radius_, youngsModulus_, poissonRatio_;
    std::vector<ConstrainedDof> constrainedDofs_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
