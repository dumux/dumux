// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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

// This test case is adapted from
// https://fenicsproject.org/olddocs/dolfin/1.6.0/python/demo/documented/hyperelasticity/python/documentation.html
// using a slightly different material law and finite volumes
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
    {
        detJ_.resize(this->gridGeometry().gridView().size(0), 0.0);
        magnSigma_.resize(this->gridGeometry().gridView().size(0), 0.0);
        mu_.resize(this->gridGeometry().gridView().size(0), 0.0);
        muFactor_ = getParam<double>("Problem.MuFactor");
        std::cout << "mu_s: " << this->spatialParams().shearModulus() << std::endl;
        std::cout << "mu: " << this->spatialParams().shearModulus()*muFactor_ << std::endl;
        std::cout << "wave length estimate: " << 2*M_PI*thickness_*std::cbrt(muFactor_/3.0) << std::endl;
        std::cout << "num waves estimate: " << 1.0/(thickness_*std::cbrt(muFactor_/3.0)) << std::endl;
    }


    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        const auto r = std::hypot(globalPos[0], globalPos[1]);
        if (r < 0.5 + 1e-2)
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
        return {0.0, 0.0, 0.0};
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        return {0.0, 0.0, 0.0};
    }

    void setGrowth(Scalar g)
    {
        std::cout << "g: " << g << std::endl;
        growth_ = g;
    }

    template<class GridVariables, class SolutionVector>
    void updateOutput(const GridVariables& gridVars, const SolutionVector& sol)
    {
        const auto& gg = gridVars.gridGeometry();
        auto fvGeometry = localView(gg);
        auto elemVolVars = localView(gridVars.curGridVolVars());
        for (const auto& element : elements(gg.gridView()))
        {
            const auto geometry = element.geometry();
            const auto center = geometry.center();
            const auto r = std::hypot(center[0], center[1]);
            fvGeometry.bindElement(element);
            elemVolVars.bindElement(element, fvGeometry, sol);

            using FeLocalBasis = typename GridGeometry::FeCache::FiniteElementType::Traits::LocalBasisType;
            using ShapeJacobian = typename FeLocalBasis::Traits::JacobianType;

            const auto F = [&]()
            {
                const auto ipLocal = geometry.local(center);
                const auto jacInvT = geometry.jacobianInverseTransposed(ipLocal);
                const auto& localBasis = fvGeometry.feLocalBasis();
                std::vector<ShapeJacobian> shapeJacobian;
                localBasis.evaluateJacobian(ipLocal, shapeJacobian);

                std::vector<GlobalPosition> gradN(fvGeometry.numScv());
                for (const auto& scv: scvs(fvGeometry))
                    jacInvT.mv(shapeJacobian[scv.localDofIndex()][0], gradN[scv.indexInElement()]);

                Dune::FieldMatrix<Scalar, dimWorld, dimWorld> F(0.0);
                for (const auto& scv : scvs(fvGeometry))
                {
                    const auto& volVars = elemVolVars[scv];
                    for (int dir = 0; dir < dimWorld; ++dir)
                        F[dir].axpy(volVars.displacement(dir), gradN[scv.indexInElement()]);
                }

                for (int dir = 0; dir < dimWorld; ++dir)
                    F[dir][dir] += 1;

                // compute elastic part of F
                if (r < 1.0-thickness_)
                    F *= 1.0/growth_;

                return F;
            }();

            const auto eIdx = gg.elementMapper().index(element);
            detJ_[eIdx] = F.determinant();

            mu_[eIdx] = this->spatialParams().shearModulus();
            if (r > 1.0-thickness_)
                mu_[eIdx] *= muFactor_;

            auto sigma = firstPiolaKirchhoffStressTensor(element, fvGeometry, F);
            sigma *= 1.0/detJ_[eIdx];
            sigma.rightmultiply(transpose(F));
            magnSigma_[eIdx] = sigma.infinity_norm();
        }
    }

    const std::vector<double>& detJ() const
    { return detJ_; }

    const std::vector<double>& magnSigma() const
    { return magnSigma_; }

    const std::vector<double>& mu() const
    { return mu_; }

    // Some Neo-Hookean material
    // from Simo and Armero, 1992 (https://doi.org/10.1002/nme.1620330705) Eq. 77/78
    // ψ = K(0.5(J*J - 1) - ln J) + 0.5µ (J^(-2/3) tr(C) - 3)
    // P = 2F∂ψ/∂C = µ J^(-2/3) (F - 1/3*tr(C)*F^-T) + K*(J*J - 1)F^-T
    Tensor firstPiolaKirchhoffStressTensor(const Element& element, const FVElementGeometry& fvGeometry, Tensor F) const
    {
        const auto center = element.geometry().center();
        const auto r = std::hypot(center[0], center[1]);

        // compute elastic part of F
        if (r < 1.0-thickness_)
            F *= 1.0/growth_;

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
        auto mu = this->spatialParams().shearModulus();
        auto K = this->spatialParams().bulkModulus();
        if (r > 1.0-thickness_)
        {
            mu *= muFactor_;
            K *= muFactor_;
        }

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
    Scalar growth_ = 1.0;
    Scalar thickness_ = 0.05;
    Scalar muFactor_;
    std::vector<Scalar> detJ_, magnSigma_, mu_;
};

} // end namespace Dumux

#endif
