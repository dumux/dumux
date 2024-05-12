// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTICITY_TEST_PROBLEM_HH
#define DUMUX_HYPERELASTICITY_TEST_PROBLEM_HH

#include <adpp/backward.hpp>
#include <adpp/backward/tensor.hpp>

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
    {}

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_ || globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
            values.setAllDirichlet();
        return values;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        return {0.1, 0.0, 0.0};
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
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
    Tensor firstPiolaKirchhoffStressTensor(Tensor FIn) const
    {
        static_assert(dimWorld == 3, "1st Piola-Kirchhoff stress tensor only implemented in 3d");

        using namespace adpp;
        using namespace adpp::backward;

        static constexpr var K;
        static constexpr var mu;
        static constexpr tensor F{adpp::shape<3, 3>};

        static constexpr auto trC = F.dot(F);
        static constexpr auto J = F.det();
        static constexpr auto psi
            = K*(cval<0.5>*(J*J - cval<1.0>) - log(J))
            + cval<0.5>*mu*(pow(J, cval<-2.0/3.0>)*trC - cval<3.0>);

        // TODO: bindings from tensors...
        const auto P = derivatives_of(psi, wrt(F.vars()), at(
            F = {
                FIn[0][0], FIn[0][1], FIn[0][2],
                FIn[1][0], FIn[1][1], FIn[1][2],
                FIn[2][0], FIn[2][1], FIn[2][2]
            },
            mu = this->spatialParams().shearModulus(),
            K = this->spatialParams().bulkModulus())
        );

        // TODO: add export_to function? Or allow for substitution of underlying array type?
        Tensor result;
        result[0][0] = P[F[md_index<0, 0>]];
        result[0][1] = P[F[md_index<0, 1>]];
        result[0][2] = P[F[md_index<0, 2>]];

        result[1][0] = P[F[md_index<1, 0>]];
        result[1][1] = P[F[md_index<1, 1>]];
        result[1][2] = P[F[md_index<1, 2>]];

        result[2][0] = P[F[md_index<2, 0>]];
        result[2][1] = P[F[md_index<2, 1>]];
        result[2][2] = P[F[md_index<2, 2>]];
        return result;
    }

private:
    static constexpr Scalar eps_ = 1e-7;
};

} // end namespace Dumux

#endif
