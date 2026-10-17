// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Growth bending strip in the mixed displacement-pressure formulation
 *
 * The compressible neo-Hookean law of Dumux::BendingStripProblem is split as
 * \f$ \mathbf{P}_\mathrm{e} = \mu(\mathbf{F}_\mathrm{e} - \mathbf{F}_\mathrm{e}^{-T})
 *     + p_\mathrm{e} J_\mathrm{e} \mathbf{F}_\mathrm{e}^{-T} \f$ with
 * \f$ p_\mathrm{e} = \lambda \ln J_\mathrm{e} / J_\mathrm{e} \f$. With growth, the first part is
 * pulled back with Dumux::Growth::firstPiolaKirchhoffStressTensor and the second becomes
 * \f$ p J \mathbf{F}^{-T} \f$ with \f$ p = \lambda \ln J_\mathrm{e} / J_\mathrm{e} \f$,
 * \f$ J_\mathrm{e} = J / J_\mathrm{g} \f$. The exact solution is stress free, so
 * \f$ J_\mathrm{e} = 1 \f$ and \f$ p = 0 \f$, while \f$ J = J_\mathrm{g} = 1 + \kappa T \f$
 * varies within the elements.
 */
#ifndef DUMUX_TEST_HYPERELASTIC_GROWTH_BENDING_STRIP_PROBLEM_MIXED_HH
#define DUMUX_TEST_HYPERELASTIC_GROWTH_BENDING_STRIP_PROBLEM_MIXED_HH

#include <cmath>
#include <memory>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/solidmechanics/hyperelastic/growth.hh>

namespace Dumux {

//! Geometry and growth of the bending strip, shared by the momentum and the pressure problem
template<class Scalar, int dim>
class BendingStripGrowth
{
    static constexpr int thicknessIdx = dim - 1;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;
    using Vector = Dune::FieldVector<Scalar, dim>;

public:
    BendingStripGrowth(const Vector& lowerLeft, const Vector& upperRight)
    : lowerLeft_(lowerLeft)
    , thickness_(upperRight[thicknessIdx] - lowerLeft[thicknessIdx])
    {}

    void setCurvature(const Scalar curvature)
    { curvature_ = curvature; }

    Tensor growthTensor(const Vector& pos) const
    {
        Vector axis(0.0);
        axis[0] = 1.0;
        return Growth::transverselyIsotropicTensor(axis, 1.0 + curvature_*thicknessCoordinate_(pos), Scalar(1.0));
    }

    Vector exactDisplacement(const Vector& pos) const
    {
        using std::sin; using std::cos;
        Vector u(0.0);
        if (curvature_ == 0.0)
            return u;

        const auto radius = 1.0/curvature_ + thicknessCoordinate_(pos);
        const auto angle = curvature_*pos[0];
        const auto midSurface = lowerLeft_[thicknessIdx] + 0.5*thickness_;
        u[0] = radius*sin(angle) - pos[0];
        u[thicknessIdx] = midSurface + radius*cos(angle) - 1.0/curvature_ - pos[thicknessIdx];
        return u;
    }

    Tensor exactDisplacementGradient(const Vector& pos) const
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
    Scalar thicknessCoordinate_(const Vector& pos) const
    { return pos[thicknessIdx] - lowerLeft_[thicknessIdx] - 0.5*thickness_; }

    Vector lowerLeft_;
    Scalar thickness_;
    Scalar curvature_ = 0.0;
};

template<class TypeTag>
class BendingStripMixedMomentumProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using DirichletConstraintData = Dumux::DirichletConstraintData<Dumux::DirichletConstraintInfo<numEq>, PrimaryVariables,
                                                                   typename IndexTraits<GridView>::GridIndex>;

    static constexpr int dim = GridView::dimension;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    BendingStripMixedMomentumProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                     std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , growth_(gridGeometry->bBoxMin(), gridGeometry->bBoxMax())
    , mu_(getParam<Scalar>("Problem.ShearModulus"))
    {
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

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    void setCurvature(const Scalar curvature)
    { growth_.setCurvature(curvature); }

    template<class GlobalPosition>
    Tensor isochoricFirstPiolaKirchhoffStressTensor(const Tensor& F, const ElementDiscretization&,
                                                    const GlobalPosition& ipGlobal) const
    {
        return Growth::firstPiolaKirchhoffStressTensor(F, growth_.growthTensor(ipGlobal), [&](const Tensor& Fe)
        {
            auto FeInv = Fe;
            FeInv.invert();
            return mu_*(Fe - getTransposed(FeInv));
        });
    }

    const BendingStripGrowth<Scalar, dim>& growth() const
    { return growth_; }

private:
    static constexpr Scalar eps_ = 1e-7;
    std::shared_ptr<CouplingManager> couplingManager_;
    BendingStripGrowth<Scalar, dim> growth_;
    Scalar mu_;
    std::vector<DirichletConstraintData> constraints_;
};

template<class TypeTag>
class BendingStripMixedPressureProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<PrimaryVariables::size()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    static constexpr int dim = GridView::dimension;

public:
    BendingStripMixedPressureProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                     std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , growth_(gridGeometry->bBoxMin(), gridGeometry->bBoxMax())
    , lambda_(getParam<Scalar>("Problem.LameLambda"))
    {}

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    void setCurvature(const Scalar curvature)
    { growth_.setCurvature(curvature); }

    //! The equilibrium pressure \f$ \lambda \ln J_\mathrm{e} / J_\mathrm{e} \f$, \f$ J_\mathrm{e} = J / J_\mathrm{g} \f$
    template<class GlobalPosition_>
    Scalar volumetricPressure(Scalar J, const ElementDiscretization&, const GlobalPosition_& ipGlobal) const
    {
        using std::log;
        const Scalar Je = J/growth_.growthTensor(ipGlobal).determinant();
        return lambda_*log(Je)/Je;
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    {
        BoundaryTypes values;
        values.setAllFluxBoundary();
        return values;
    }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    BendingStripGrowth<Scalar, dim> growth_;
    Scalar lambda_;
};

} // end namespace Dumux

#endif
