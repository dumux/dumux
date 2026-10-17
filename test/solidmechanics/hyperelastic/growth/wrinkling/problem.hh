// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief A film growing on a substrate in plane strain
 *
 * The substrate \f$ 0 < y < H \f$ is clamped at \f$ y = 0 \f$, the film \f$ H < y < H + h \f$
 * grows along \f$ x \f$ with \f$ \mathbf{F}_\mathrm{g} = \mathbf{I} + (\vartheta - 1)\mathbf{e}_x \otimes \mathbf{e}_x \f$
 * and has a traction-free top. The sides \f$ x = 0 \f$ and \f$ x = L \f$ slide
 * (\f$ u_x = 0 \f$). Both layers are compressible neo-Hookean with their own moduli. The flat state,
 * substrate unstrained and film compressed along \f$ x \f$ and stretched along \f$ y \f$, is a
 * solution for every \f$ \vartheta \f$ and becomes unstable at a critical growth stretch, where the
 * film wrinkles.
 */
#ifndef DUMUX_TEST_HYPERELASTIC_GROWTH_WRINKLING_PROBLEM_HH
#define DUMUX_TEST_HYPERELASTIC_GROWTH_WRINKLING_PROBLEM_HH

#include <cmath>
#include <map>
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
class FilmSubstrateProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using ConstraintInfo = Dumux::DirichletConstraintInfo<numEq>;
    using GridIndex = typename IndexTraits<GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, PrimaryVariables, GridIndex>;

    static constexpr int dim = GridView::dimension;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;
    using Vector = Dune::FieldVector<Scalar, dim>;

public:
    FilmSubstrateProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , substrateDepth_(getParam<Scalar>("Problem.SubstrateDepth"))
    , film_{getParam<Scalar>("Problem.FilmShearModulus"), getParam<Scalar>("Problem.FilmLameLambda")}
    , substrate_{getParam<Scalar>("Problem.SubstrateShearModulus"), getParam<Scalar>("Problem.SubstrateLameLambda")}
    {
        const auto& lowerLeft = gridGeometry->bBoxMin();
        const auto& upperRight = gridGeometry->bBoxMax();
        std::map<GridIndex, ConstraintInfo> constrained;
        auto elemDisc = localView(*gridGeometry);
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            elemDisc.bind(element);
            for (const auto& boundaryFace : boundaryFaces(elemDisc))
                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                {
                    const auto pos = ipData(elemDisc, localDof).global();
                    auto& info = constrained[localDof.dofIndex()];
                    if (pos[1] < lowerLeft[1] + eps_)
                        info.setAll();
                    if (pos[0] < lowerLeft[0] + eps_ || pos[0] > upperRight[0] - eps_)
                        info.set(0);
                }
        }
        for (const auto& [dofIdx, info] : constrained)
            if (info.isConstraintEquation(0) || info.isConstraintEquation(1))
                constraints_.push_back(DirichletConstraintData{info, PrimaryVariables(0.0), dofIdx});
    }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    void setGrowthStretch(const Scalar stretch)
    { growthStretch_ = stretch; }

    template<class IpData>
    Tensor firstPiolaKirchhoffStressTensor(const Tensor& F, const ElementDiscretization& elemDisc, const IpData&) const
    {
        const bool isFilm = elemDisc.element().geometry().center()[1] > substrateDepth_;
        const auto& layer = isFilm ? film_ : substrate_;
        Vector axis(0.0);
        axis[0] = 1.0;
        const auto Fg = Growth::transverselyIsotropicTensor(axis, isFilm ? growthStretch_ : Scalar(1.0), Scalar(1.0));
        return Growth::firstPiolaKirchhoffStressTensor(F, Fg, [&](const Tensor& Fe)
        {
            using std::log;
            auto FeInvT = Fe;
            FeInvT.invert();
            FeInvT = getTransposed(FeInvT);
            return layer.shearModulus*Fe + (layer.lameLambda*log(Fe.determinant()) - layer.shearModulus)*FeInvT;
        });
    }

private:
    struct Layer { Scalar shearModulus, lameLambda; };

    static constexpr Scalar eps_ = 1e-8;
    Scalar substrateDepth_;
    Layer film_, substrate_;
    Scalar growthStretch_ = 1.0;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
