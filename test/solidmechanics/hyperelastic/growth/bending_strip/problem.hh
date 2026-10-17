// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief A strip bent into a circular arc by growth that varies linearly through its thickness
 *
 * The strip \f$ [0, L] \times [0, h] \f$ (in 3D \f$ [0, L] \times [0, b] \times [0, h] \f$) grows
 * along its axis \f$ X \f$ with the stretch \f$ 1 + \kappa T \f$, where \f$ T \f$ is the distance
 * from the mid-surface in the thickness direction (the last coordinate),
 * \f[ \mathbf{F}_\mathrm{g} = \mathbf{I} + \kappa T\, \mathbf{e}_X \otimes \mathbf{e}_X. \f]
 * This growth is compatible: the deformation
 * \f[ x_X = (\kappa^{-1} + T) \sin(\kappa X), \quad
 *     x_T = h/2 + (\kappa^{-1} + T) \cos(\kappa X) - \kappa^{-1}, \f]
 * which bends the strip into an arc of radius \f$ \kappa^{-1} \f$ about the mid-surface, has
 * \f$ \mathbf{F} = \mathbf{R}\mathbf{F}_\mathrm{g} \f$ with a rotation \f$ \mathbf{R} \f$, so it is
 * stress free for any objective constitutive law. It satisfies the clamped end at \f$ X = 0 \f$
 * and the traction-free remaining boundary. The material is compressible neo-Hookean.
 */
#ifndef DUMUX_TEST_HYPERELASTIC_GROWTH_BENDING_STRIP_PROBLEM_HH
#define DUMUX_TEST_HYPERELASTIC_GROWTH_BENDING_STRIP_PROBLEM_HH

#include <cmath>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/solidmechanics/hyperelastic/growth.hh>

namespace Dumux {

template<class TypeTag>
class BendingStripProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using DirichletConstraintData = Dumux::DirichletConstraintData<Dumux::DirichletConstraintInfo<numEq>, PrimaryVariables,
                                                                   typename IndexTraits<GridView>::GridIndex>;

    static constexpr int dim = GridView::dimension;
    static constexpr int thicknessIdx = dim - 1;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;
    using Vector = Dune::FieldVector<Scalar, dim>;

public:
    BendingStripProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        mu_ = getParam<Scalar>("Problem.ShearModulus");
        lambda_ = getParam<Scalar>("Problem.LameLambda");
        thickness_ = gridGeometry->bBoxMax()[thicknessIdx] - gridGeometry->bBoxMin()[thicknessIdx];

        auto elemDisc = localView(*gridGeometry);
        std::vector<bool> isConstrained(gridGeometry->numDofs(), false);
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            elemDisc.bind(element);
            for (const auto& boundaryFace : boundaryFaces(elemDisc))
                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                {
                    const auto globalPos = ipData(elemDisc, localDof).global();
                    if (globalPos[0] < gridGeometry->bBoxMin()[0] + eps_ && !isConstrained[localDof.dofIndex()])
                    {
                        isConstrained[localDof.dofIndex()] = true;
                        typename DirichletConstraintData::ConstraintInfo info;
                        info.setAll();
                        constraints_.push_back(DirichletConstraintData{std::move(info), PrimaryVariables(0.0), localDof.dofIndex()});
                    }
                }
        }
    }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    void setCurvature(const Scalar curvature)
    { curvature_ = curvature; }

    template<class IpData>
    Tensor firstPiolaKirchhoffStressTensor(const Tensor& F, const ElementDiscretization&, const IpData& ipData) const
    {
        return Growth::firstPiolaKirchhoffStressTensor(F, growthTensor(ipData.global()),
                                                       [&](const Tensor& Fe) { return neoHookeanStress_(Fe); });
    }

    Tensor growthTensor(const GlobalPosition& pos) const
    {
        Vector axis(0.0);
        axis[0] = 1.0;
        return Growth::transverselyIsotropicTensor(axis, 1.0 + curvature_*thicknessCoordinate_(pos), Scalar(1.0));
    }

    Vector exactDisplacement(const GlobalPosition& pos) const
    {
        using std::sin; using std::cos;
        Vector u(0.0);
        if (curvature_ == 0.0)
            return u;

        const auto radius = 1.0/curvature_ + thicknessCoordinate_(pos);
        const auto angle = curvature_*pos[0];
        u[0] = radius*sin(angle) - pos[0];
        const auto midSurface = this->gridDiscretization().bBoxMin()[thicknessIdx] + 0.5*thickness_;
        u[thicknessIdx] = midSurface + radius*cos(angle) - 1.0/curvature_ - pos[thicknessIdx];
        return u;
    }

    Tensor exactDisplacementGradient(const GlobalPosition& pos) const
    {
        using std::sin; using std::cos;
        Tensor gradU(0.0);
        const auto stretch = 1.0 + curvature_*thicknessCoordinate_(pos);
        const auto angle = curvature_*pos[0];
        gradU[0][0] = stretch*cos(angle) - 1.0;
        gradU[0][thicknessIdx] = sin(angle);
        gradU[thicknessIdx][0] = -stretch*sin(angle);
        gradU[thicknessIdx][thicknessIdx] = cos(angle) - 1.0;
        return gradU;
    }

private:
    Scalar thicknessCoordinate_(const GlobalPosition& pos) const
    { return pos[thicknessIdx] - this->gridDiscretization().bBoxMin()[thicknessIdx] - 0.5*thickness_; }

    Tensor neoHookeanStress_(const Tensor& F) const
    {
        using std::log;
        auto FinvT = F;
        FinvT.invert();
        FinvT = getTransposed(FinvT);
        return mu_*F + (lambda_*log(F.determinant()) - mu_)*FinvT;
    }

    static constexpr Scalar eps_ = 1e-7;
    Scalar mu_, lambda_, thickness_;
    Scalar curvature_ = 0.0;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
