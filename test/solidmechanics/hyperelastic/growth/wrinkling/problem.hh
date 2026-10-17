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
 *
 * Beyond the onset, the wrinkled branch is traced by prescribing the vertical displacement of the top
 * corner at \f$ x = 0 \f$ in addition (`enablePointConstraint`) and searching the growth stretch at
 * which its reaction vanishes.
 */
#ifndef DUMUX_TEST_HYPERELASTIC_GROWTH_WRINKLING_PROBLEM_HH
#define DUMUX_TEST_HYPERELASTIC_GROWTH_WRINKLING_PROBLEM_HH

#include <algorithm>
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

        // the top corner at x = 0
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            elemDisc.bind(element);
            for (const auto& localDof : localDofs(elemDisc))
            {
                const auto pos = ipData(elemDisc, localDof).global();
                if (pos[0] < lowerLeft[0] + eps_ && pos[1] > upperRight[1] - eps_)
                    pointDofIndex_ = localDof.dofIndex();
            }
        }
        pointConstraintInfo_ = constrained[pointDofIndex_];
        pointConstraintInfo_.set(1);
        constraintsWithPoint_ = constraints_;
        std::erase_if(constraintsWithPoint_, [&](const auto& c) { return c.dofIndex() == pointDofIndex_; });
        constraintsWithPoint_.push_back(DirichletConstraintData{pointConstraintInfo_, PrimaryVariables(0.0), pointDofIndex_});
    }

    const std::vector<DirichletConstraintData>& constraints() const
    { return pointConstraintEnabled_ ? constraintsWithPoint_ : constraints_; }

    //! Prescribes the vertical displacement of the top corner at x = 0 in addition to the boundary conditions
    void enablePointConstraint(bool enable)
    { pointConstraintEnabled_ = enable; }

    void setPointDisplacement(Scalar displacement)
    {
        PrimaryVariables values(0.0);
        values[1] = displacement;
        constraintsWithPoint_.back() = DirichletConstraintData{pointConstraintInfo_, values, pointDofIndex_};
    }

    GridIndex pointDofIndex() const
    { return pointDofIndex_; }

    //! The stretch of the film thickness in the flat state, where the film is free of stress along y
    Scalar filmThicknessStretch(Scalar growthStretch) const
    {
        const auto stressYY = [&](Scalar stretch)
        {
            Tensor F(0.0);
            F[0][0] = 1.0;
            F[1][1] = stretch;
            return layerStress_(F, film_, growthStretch)[1][1];
        };

        Scalar stretch = 1.0;
        for (int iter = 0; iter < 50; ++iter)
        {
            constexpr Scalar h = 1e-7;
            const Scalar f = stressYY(stretch);
            stretch -= f/((stressYY(stretch + h) - stressYY(stretch - h))/(2.0*h));
            if (std::abs(f) < 1e-14)
                break;
        }
        return stretch;
    }

    Scalar filmThickness() const
    { return this->gridDiscretization().bBoxMax()[1] - substrateDepth_; }

    void setGrowthStretch(const Scalar stretch)
    { growthStretch_ = stretch; }

    template<class IpData>
    Tensor firstPiolaKirchhoffStressTensor(const Tensor& F, const ElementDiscretization& elemDisc, const IpData&) const
    {
        const bool isFilm = elemDisc.element().geometry().center()[1] > substrateDepth_;
        return isFilm ? layerStress_(F, film_, growthStretch_) : layerStress_(F, substrate_, 1.0);
    }

private:
    struct Layer { Scalar shearModulus, lameLambda; };

    static Tensor layerStress_(const Tensor& F, const Layer& layer, Scalar growthStretch)
    {
        Vector axis(0.0);
        axis[0] = 1.0;
        const auto Fg = Growth::transverselyIsotropicTensor(axis, growthStretch, Scalar(1.0));
        return Growth::firstPiolaKirchhoffStressTensor(F, Fg, [&](const Tensor& Fe)
        {
            using std::log;
            auto FeInvT = Fe;
            FeInvT.invert();
            FeInvT = getTransposed(FeInvT);
            return layer.shearModulus*Fe + (layer.lameLambda*log(Fe.determinant()) - layer.shearModulus)*FeInvT;
        });
    }

    static constexpr Scalar eps_ = 1e-8;
    Scalar substrateDepth_;
    Layer film_, substrate_;
    Scalar growthStretch_ = 1.0;
    std::vector<DirichletConstraintData> constraints_, constraintsWithPoint_;
    ConstraintInfo pointConstraintInfo_;
    GridIndex pointDofIndex_ = 0;
    bool pointConstraintEnabled_ = false;
};

} // end namespace Dumux

#endif
