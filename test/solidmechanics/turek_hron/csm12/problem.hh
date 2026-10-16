// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_DYNAMIC_HYPERELASTICITY_TEST_PROBLEM_HH
#define DUMUX_DYNAMIC_HYPERELASTICITY_TEST_PROBLEM_HH

#include <vector>

#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>

#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

// This test case corresponds to the CSM benchmark problem
// of a hyperelastic solid material under gravity given in
// Turek, S., Hron, J. (2006). "Proposal for Numerical Benchmarking of Fluid-Structure Interaction
// between an Elastic Object and Laminar Incompressible Flow."
// https://doi.org/10.1007/3-540-34596-5_15
template<class TypeTag>
class DynamicHyperelasticityProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using DirichletConstraintData = Dumux::DirichletConstraintData<Dumux::DirichletConstraintInfo<numEq>, PrimaryVariables,
                                                                   typename IndexTraits<GridView>::GridIndex>;
    static constexpr int dim = GridView::dimension;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

public:
    DynamicHyperelasticityProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , name_(getParam<std::string>("Problem.Name"))
    , gravity_(getParam<Scalar>("Problem.Gravity"))
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
                    if (isClamped_(globalPos) && !isConstrained[localDof.dofIndex()])
                    {
                        isConstrained[localDof.dofIndex()] = true;
                        typename DirichletConstraintData::ConstraintInfo info;
                        info.setAll();
                        constraints_.push_back(DirichletConstraintData{std::move(info), dirichletAtPos(globalPos), localDof.dofIndex()});
                    }
                }
        }
        makeConstraintsConsistent(*gridGeometry, constraints_);
    }

    const std::string& name() const
    { return name_; }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        // gravity forcing
        return {0.0, -this->spatialParams().solidDensity()*gravity_};
    }

    // Saint-Venant Kirchhoff material
    Tensor firstPiolaKirchhoffStressTensor(const Tensor& F) const
    {
        // material parameters
        const auto mu = this->spatialParams().shearModulus();
        const auto lambda = this->spatialParams().firstLameParameter();

        // Lagrangian Green strain E = 1/2*(F^T F - I)
        auto E = multiplyMatrices(transpose(F), F);
        E *= 0.5;
        Scalar trace = 0.0;
        for (int i = 0; i < dim; ++i)
        {
            E[i][i] -= 0.5;
            trace += E[i][i];
        }

        // 2nd Piola Kirchhoff stress tensor S = λtr(E)I + 2µE
        auto& S = E;
        S *= 2*mu;
        for (int i = 0; i < dim; ++i)
            S[i][i] += lambda*trace;

        // 1st Piola Kirchhoff stress tensor P = FS
        return multiplyMatrices(F, S);
    }

    // compute displacement at benchmark evaluation point A
    template<class SolutionVector, class GridVariables>
    GlobalPosition evalControlPointDisplacement(const GridVariables&, const SolutionVector& p) const
    {
        if (!controlPointElement_)
        {
            const auto& tree = this->gridDiscretization().boundingBoxTree();
            auto entities = intersectingEntities(controlPointA_, tree);
            if (entities.empty()) DUNE_THROW(Dune::Exception, "Control point A not found!");
            controlPointElement_ = entities[0];
        }

        const auto& gg = this->gridDiscretization();
        const auto element = gg.element(*controlPointElement_);
        const auto elemSol = elementSolution(element, p, gg);

        return evalSolution(element, element.geometry(), gg, elemSol, controlPointA_);
    }


private:
    bool isClamped_(const GlobalPosition& globalPos) const
    { return std::hypot(globalPos[0]-0.2, globalPos[1]-0.2) < 0.05 + eps_; }

    static constexpr Scalar eps_ = 1e-7;
    std::vector<DirichletConstraintData> constraints_;
    std::string name_;
    Scalar gravity_;
    GlobalPosition controlPointA_{{0.6, 0.2}};
    mutable std::optional<std::size_t> controlPointElement_;
};

} // end namespace Dumux

#endif
