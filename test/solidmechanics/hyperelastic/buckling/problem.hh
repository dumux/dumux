// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTICITY_TEST_PROBLEM_HH
#define DUMUX_HYPERELASTICITY_TEST_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>

namespace Dumux {

template<class TypeTag>
class HyperelasticityProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int dim = GridView::dimension;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    HyperelasticityProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , displacement_(0.0)
    {}

    void setCurrentDisplacement(const Scalar displacement)
    { displacement_ = displacement; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_ || globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        return NumEqVector(0.0);
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        if (globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
            return {0.0, 0.0, 0.0};
        else
            return {displacement_, 0.0, 0.0};
    }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    { return {0.0, 0.0, 0.0}; }

    // Some Neo-Hookean material
    // from Simo and Armero, 1992 (https://doi.org/10.1002/nme.1620330705) Eq. 77/78
    // ψ = K(0.5(J*J - 1) - ln J) + 0.5µ (J^(-2/3) tr(C) - 3)
    // P = 2F∂ψ/∂C = µ J^(-2/3) (F - 1/3*tr(C)*F^-T) + K*(J*J - 1)F^-T
    Tensor firstPiolaKirchhoffStressTensor(Tensor F, const GlobalPosition& globalPos) const
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
        const auto mu = this->spatialParams().shearModulus(globalPos);
        const auto K = this->spatialParams().bulkModulus(globalPos);

        // assemble 1st Piola Kirchhoff stress tensor
        Tensor P(0.0);
        P.axpy(K*(J*J - 1), invFT);
        auto& FTerm = F;
        FTerm.axpy(-1.0/3.0*trC, invFT);
        P.axpy(mu*Jm23, FTerm);
        return P;
    }

private:
    static constexpr Scalar eps_ = 1e-7;
    Scalar displacement_;
};

} // end namespace Dumux

#endif
