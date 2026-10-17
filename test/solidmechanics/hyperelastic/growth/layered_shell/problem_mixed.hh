// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Growing layered shell in the mixed displacement-pressure formulation
 *
 * The compressible neo-Hookean law of Dumux::LayeredShellProblem is split as
 * \f$ \mathbf{P}_\mathrm{e} = \mu(\mathbf{F}_\mathrm{e} - \mathbf{F}_\mathrm{e}^{-T})
 *     + p_\mathrm{e} J_\mathrm{e} \mathbf{F}_\mathrm{e}^{-T} \f$ with
 * \f$ p_\mathrm{e} = \lambda \ln J_\mathrm{e} / J_\mathrm{e} \f$. With growth, the first part is
 * pulled back with Dumux::Growth::firstPiolaKirchhoffStressTensor and the second becomes
 * \f$ p J \mathbf{F}^{-T} \f$ with \f$ p = \lambda \ln J_\mathrm{e} / J_\mathrm{e} \f$,
 * \f$ J_\mathrm{e} = J / J_\mathrm{g} \f$, so both formulations solve the same problem.
 */
#ifndef DUMUX_TEST_HYPERELASTIC_GROWTH_LAYERED_SHELL_PROBLEM_MIXED_HH
#define DUMUX_TEST_HYPERELASTIC_GROWTH_LAYERED_SHELL_PROBLEM_MIXED_HH

#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <string>
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

#include "radialreference.hh"

namespace Dumux {

template<class Scalar>
struct LayeredShellParameters
{
    using Reference = LayeredShellRadialReference<Scalar>;
    using Layer = typename Reference::Layer;

    LayeredShellParameters()
    : radii{getParam<Scalar>("Problem.InnerRadius"),
            getParam<Scalar>("Problem.InterfaceRadius"),
            getParam<Scalar>("Problem.OuterRadius")}
    {
        for (int i = 0; i < 2; ++i)
        {
            const std::string group = i == 0 ? "Problem.Inner." : "Problem.Outer.";
            layers[i] = Layer{getParam<Scalar>(group + "GrowthStretch"),
                              getParam<Scalar>(group + "ShearModulus"),
                              getParam<Scalar>(group + "LameLambda")};
        }
    }

    template<class Element>
    const Layer& layer(const Element& element) const
    { return element.geometry().center().two_norm() < radii[1] ? layers[0] : layers[1]; }

    std::array<Scalar, 3> radii;
    std::array<Layer, 2> layers;
};

template<class TypeTag>
class LayeredShellMixedMomentumProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using ConstraintInfo = Dumux::DirichletConstraintInfo<numEq>;
    using GridIndex = typename IndexTraits<GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, PrimaryVariables, GridIndex>;

    static constexpr int dim = GridView::dimension;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    LayeredShellMixedMomentumProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                     std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , reference_(dim - 1, params_.radii, params_.layers)
    {
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

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    //! Sets the growth stretches to the fraction s of their final values
    void setGrowthFraction(const Scalar s)
    { growthFraction_ = s; }

    template<class GlobalPosition_>
    Tensor isochoricFirstPiolaKirchhoffStressTensor(const Tensor& F, const ElementDiscretization& elemDisc,
                                                    const GlobalPosition_&) const
    {
        const auto& layer = params_.layer(elemDisc.element());
        const Scalar stretch = 1.0 + growthFraction_*(layer.growthStretch - 1.0);
        return Growth::firstPiolaKirchhoffStressTensor(F, Growth::isotropicTensor<dim>(stretch),
                                                       [&](const Tensor& Fe)
        {
            auto FeInv = Fe;
            FeInv.invert();
            return layer.shearModulus*(Fe - getTransposed(FeInv));
        });
    }

    Dune::FieldVector<Scalar, dim> exactDisplacement(const GlobalPosition& pos) const
    {
        const Scalar R = pos.two_norm();
        auto u = pos;
        u *= reference_.radius(R)/R - 1.0;
        return u;
    }

    Tensor exactDisplacementGradient(const GlobalPosition& pos) const
    {
        const Scalar R = pos.two_norm();
        const Scalar radialStretch = reference_.radialStretch(R);
        const Scalar hoopStretch = reference_.radius(R)/R;
        Tensor gradU(0.0);
        for (int i = 0; i < dim; ++i)
        {
            gradU[i][i] = hoopStretch - 1.0;
            for (int j = 0; j < dim; ++j)
                gradU[i][j] += (radialStretch - hoopStretch)*pos[i]*pos[j]/(R*R);
        }
        return gradU;
    }

    //! The exact pressure \f$ \lambda \ln J_\mathrm{e} / J_\mathrm{e} \f$
    Scalar exactPressure(const GlobalPosition& pos, Scalar lameLambda, Scalar growthStretch) const
    {
        using std::log; using std::pow;
        const Scalar R = pos.two_norm();
        const Scalar J = reference_.radialStretch(R)*pow(reference_.radius(R)/R, dim - 1);
        const Scalar Je = J/pow(growthStretch, dim);
        return lameLambda*log(Je)/Je;
    }

    const LayeredShellParameters<Scalar>& parameters() const
    { return params_; }

private:
    static constexpr Scalar eps_ = 1e-8;
    std::shared_ptr<CouplingManager> couplingManager_;
    LayeredShellParameters<Scalar> params_;
    LayeredShellRadialReference<Scalar> reference_;
    Scalar growthFraction_ = 1.0;
    std::vector<DirichletConstraintData> constraints_;
};

template<class TypeTag>
class LayeredShellMixedPressureProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
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
    LayeredShellMixedPressureProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                     std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    {}

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    void setGrowthFraction(const Scalar s)
    { growthFraction_ = s; }

    //! The equilibrium pressure \f$ \lambda \ln J_\mathrm{e} / J_\mathrm{e} \f$, \f$ J_\mathrm{e} = J / J_\mathrm{g} \f$
    template<class GlobalPosition_>
    Scalar volumetricPressure(Scalar J, const ElementDiscretization& elemDisc, const GlobalPosition_&) const
    {
        using std::log; using std::pow;
        const auto& layer = params_.layer(elemDisc.element());
        const Scalar stretch = 1.0 + growthFraction_*(layer.growthStretch - 1.0);
        const Scalar Je = J/pow(stretch, dim);
        return layer.lameLambda*log(Je)/Je;
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    {
        BoundaryTypes values;
        values.setAllFluxBoundary();
        return values;
    }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    LayeredShellParameters<Scalar> params_;
    Scalar growthFraction_ = 1.0;
};

} // end namespace Dumux

#endif
