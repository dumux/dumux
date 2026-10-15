// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTICITY_TEST_PROBLEM_HH
#define DUMUX_HYPERELASTICITY_TEST_PROBLEM_HH

#include <vector>

#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/math.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

// This test case is adapted from
// https://fenicsproject.org/olddocs/dolfin/1.6.0/python/demo/documented/hyperelasticity/python/documentation.html
// using a slightly different material law and finite volumes
template<class TypeTag>
class HyperelasticityProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using DirichletConstraintData = Dumux::DirichletConstraintData<Dumux::DirichletConstraintInfo<numEq>, PrimaryVariables,
                                                                   typename IndexTraits<GridView>::GridIndex>;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int dim = GridView::dimension;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    HyperelasticityProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
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
                    if (isDirichletBoundary_(globalPos) && !isConstrained[localDof.dofIndex()])
                    {
                        isConstrained[localDof.dofIndex()] = true;
                        typename DirichletConstraintData::ConstraintInfo info;
                        info.setAll();
                        constraints_.push_back(DirichletConstraintData{std::move(info), dirichletAtPos(globalPos), localDof.dofIndex()});
                    }
                }
        }
    }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        if (globalPos[0] > this->gridDiscretization().bBoxMax()[0] - eps_)
            return {0.0, 0.0, 0.0};

        PrimaryVariables values(0.0);
        const auto y = globalPos[1];
        const auto z = globalPos[2];
        values[0] = 0.0;
        values[1] = 0.5*(0.5 + (y-0.5)*std::cos(M_PI/3.0) - (z-0.5)*std::sin(M_PI/3.0) - y);
        values[2] = 0.5*(0.5 + (y-0.5)*std::sin(M_PI/3.0) + (z-0.5)*std::cos(M_PI/3.0) - z);
        return values;
    }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    { return {0.0, -0.5, 0.0}; }

    // Some Neo-Hookean material
    // from Simo and Armero, 1992 (https://doi.org/10.1002/nme.1620330705) Eq. 77/78
    // ψ = K(0.5(J*J - 1) - ln J) + 0.5µ (J^(-2/3) tr(C) - 3)
    // P = 2F∂ψ/∂C = µ J^(-2/3) (F - 1/3*tr(C)*F^-T) + K*(J*J - 1)F^-T
    Tensor firstPiolaKirchhoffStressTensor(Tensor F) const
    {
        // invariants
        const auto J = F.determinant();
        const auto Jm23 = std::pow(J, -2.0/3.0);
        auto trC = 0.0;
        for (int i = 0; i < dimWorld; ++i)
            for (int j = 0; j < dimWorld; ++j)
                trC += F[i][j]*F[i][j];

        // inverse transpose of deformation gradient
        const auto invFT = [&](){
            auto invFT = F;
            if (J != 0.0)
                invFT.invert();
            return transpose(invFT);
        }();

        // material parameters
        const auto mu = this->spatialParams().shearModulus();
        const auto K = this->spatialParams().bulkModulus();

        // assemble 1st Piola Kirchhoff stress tensor
        Tensor P(0.0);
        P.axpy(K*(J*J - 1), invFT);
        auto& FTerm = F;
        FTerm.axpy(-1.0/3.0*trC, invFT);
        P.axpy(mu*Jm23, FTerm);
        return P;
    }

private:
    bool isDirichletBoundary_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < this->gridDiscretization().bBoxMin()[0] + eps_
            || globalPos[0] > this->gridDiscretization().bBoxMax()[0] - eps_;
    }

    static constexpr Scalar eps_ = 1e-7;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
