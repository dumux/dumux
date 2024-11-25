// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HYPERELASTICITY_TEST_PROBLEM_HH
#define DUMUX_HYPERELASTICITY_TEST_PROBLEM_HH

#include <xpress/xp.hpp>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>


// add compatibility traits for Dune::FieldMatrix (TODO: where to put such thing!? in deps folder as well?)
#ifndef DOXYGEN
namespace xp {

template<typename T, int rows, int cols>
struct shape_of<Dune::FieldMatrix<T, rows, cols>> { using type = xp::md_shape<rows, cols>; };

}  // namespace xp
#endif  // DOXYGEN

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

        using real = xp::dtype::real;
        using shape = xp::md_shape<3, 3>;

        const xp::var<real> K;
        const xp::var<real> mu;
        const xp::tensor<shape, real> F;

        const auto trC = F*F;
        const auto J = det(F);
        const auto psi = K*(xp::val<0.5>*(J*J - xp::val<1.0>) - log(J))
                        + xp::val<0.5>*mu*(pow(J, xp::val<-2.0/3.0>)*trC - xp::val<3.0>);

        const auto J_value = value_of(J, at(F = FIn));
        derivative_of(psi, wrt(F), at(
            mu = this->spatialParams().shearModulus(),
            K = this->spatialParams().bulkModulus(),
            F = FIn,
            // precompute a few expressions to speed up computations
            J = J_value,
            trC = value_of(trC, at(F = FIn)),
            J*J = value_of(J*J, at(J = J_value))
        )).export_to(FIn);

        return FIn;
    }

private:
    static constexpr Scalar eps_ = 1e-7;
};

} // end namespace Dumux

#endif
