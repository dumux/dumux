// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Two-layer hollow cylinder (plane strain, 2D) or sphere (3D) with isotropic growth in each layer
 *
 * The layers \f$ R_i \leq R < R_c \f$ and \f$ R_c < R \leq R_o \f$ grow isotropically with
 * the stretches \f$ \vartheta_\mathrm{inner} \f$ and \f$ \vartheta_\mathrm{outer} \f$ and are
 * compressible neo-Hookean with their own moduli. Both surfaces are traction free. A growth
 * mismatch between the layers leaves residual stress, so the solution depends on the pull-back
 * \f$ J_\mathrm{g} \mathbf{P}_\mathrm{e} \mathbf{F}_\mathrm{g}^{-T} \f$ of the stress. A quarter
 * (2D) or an eighth (3D) of the body is computed with symmetry conditions and compared with the
 * radially symmetric solution of Dumux::LayeredShellRadialReference.
 */
#ifndef DUMUX_TEST_HYPERELASTIC_GROWTH_LAYERED_SHELL_PROBLEM_HH
#define DUMUX_TEST_HYPERELASTIC_GROWTH_LAYERED_SHELL_PROBLEM_HH

#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <string>
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

#include "radialreference.hh"

namespace Dumux {

template<class TypeTag>
class LayeredShellProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using ConstraintInfo = Dumux::DirichletConstraintInfo<numEq>;
    using GridIndex = typename IndexTraits<GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, PrimaryVariables, GridIndex>;

    static constexpr int dim = GridView::dimension;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;
    using Reference = LayeredShellRadialReference<Scalar>;
    using Layer = typename Reference::Layer;

public:
    LayeredShellProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , radii_{getParam<Scalar>("Problem.InnerRadius"),
             getParam<Scalar>("Problem.InterfaceRadius"),
             getParam<Scalar>("Problem.OuterRadius")}
    {
        for (int i = 0; i < 2; ++i)
        {
            const std::string group = i == 0 ? "Problem.Inner." : "Problem.Outer.";
            layers_[i] = Layer{getParam<Scalar>(group + "GrowthStretch"),
                               getParam<Scalar>(group + "ShearModulus"),
                               getParam<Scalar>(group + "LameLambda")};
        }
        reference_ = std::make_unique<Reference>(dim - 1, radii_, layers_);

        // symmetry planes x_i = 0
        std::map<GridIndex, ConstraintInfo> symmetry;
        auto elemDisc = localView(*gridGeometry);
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            elemDisc.bind(element);
            for (const auto& boundaryFace : boundaryFaces(elemDisc))
                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                {
                    const auto globalPos = ipData(elemDisc, localDof).global();
                    for (int i = 0; i < dim; ++i)
                        if (std::abs(globalPos[i]) < eps_)
                            symmetry[localDof.dofIndex()].set(i);
                }
        }
        for (const auto& [dofIdx, info] : symmetry)
            constraints_.push_back(DirichletConstraintData{info, PrimaryVariables(0.0), dofIdx});
    }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    //! Sets the growth stretches to the fraction s of their final values
    void setGrowthFraction(const Scalar s)
    { growthFraction_ = s; }

    template<class IpData>
    Tensor firstPiolaKirchhoffStressTensor(const Tensor& F, const ElementDiscretization& elemDisc, const IpData&) const
    {
        const auto& layer = layer_(elemDisc);
        const Scalar stretch = 1.0 + growthFraction_*(layer.growthStretch - 1.0);
        return Growth::firstPiolaKirchhoffStressTensor(F, Growth::isotropicTensor<dim>(stretch),
                                                       [&](const Tensor& Fe) { return neoHookeanStress_(Fe, layer); });
    }

    Dune::FieldVector<Scalar, dim> exactDisplacement(const GlobalPosition& pos) const
    {
        const Scalar R = pos.two_norm();
        auto u = pos;
        u *= reference_->radius(R)/R - 1.0;
        return u;
    }

    Tensor exactDisplacementGradient(const GlobalPosition& pos) const
    {
        const Scalar R = pos.two_norm();
        const Scalar radialStretch = reference_->radialStretch(R);
        const Scalar hoopStretch = reference_->radius(R)/R;
        Tensor gradU(0.0);
        for (int i = 0; i < dim; ++i)
        {
            gradU[i][i] = hoopStretch - 1.0;
            for (int j = 0; j < dim; ++j)
                gradU[i][j] += (radialStretch - hoopStretch)*pos[i]*pos[j]/(R*R);
        }
        return gradU;
    }

private:
    const Layer& layer_(const ElementDiscretization& elemDisc) const
    { return elemDisc.element().geometry().center().two_norm() < radii_[1] ? layers_[0] : layers_[1]; }

    Tensor neoHookeanStress_(const Tensor& F, const Layer& layer) const
    {
        using std::log;
        auto FinvT = F;
        FinvT.invert();
        FinvT = getTransposed(FinvT);
        return layer.shearModulus*F + (layer.lameLambda*log(F.determinant()) - layer.shearModulus)*FinvT;
    }

    static constexpr Scalar eps_ = 1e-8;
    std::array<Scalar, 3> radii_;
    std::array<Layer, 2> layers_;
    std::unique_ptr<Reference> reference_;
    Scalar growthFraction_ = 1.0;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
