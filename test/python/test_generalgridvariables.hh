// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_PYTHON_GENERAL_GRIDVARIABLES_HH
#define DUMUX_TEST_PYTHON_GENERAL_GRIDVARIABLES_HH

#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

/*!
 * \brief A homogeneous deformation imposed on the whole boundary of a hyperelastic body.
 *
 * Without a body force, u = (F - I)x solves the problem for any deformation gradient F: the stress
 * is the same everywhere, so it has no divergence. A scheme with linear shape functions represents
 * u exactly and balances the constant stress on every control volume, so the discrete solution is
 * u itself (the patch test).
 */
template<class TypeTag>
class HomogeneousDeformationProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Geometry::GlobalCoordinate;
    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = PrimaryVariables::size();
    using DirichletConstraintData = Dumux::DirichletConstraintData<
        DirichletConstraintInfo<numEq>, PrimaryVariables, typename IndexTraits<GridView>::GridIndex>;

public:
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

    HomogeneousDeformationProblem(std::shared_ptr<const GridGeometry> gridGeometry, const Tensor& F)
    : ParentType(gridGeometry)
    , gridGeometry_(gridGeometry)
    , F_(F)
    {
        appendBoundaryDirichletConstraints(
            constraints_, *gridGeometry,
            [] (const auto&, const auto&, const auto&)
            {
                typename DirichletConstraintData::ConstraintInfo info;
                info.setAll();
                return info;
            },
            [&] (const auto& elemDisc, const auto& localDof)
            { return displacement(ipData(elemDisc, localDof).global()); }
        );
    }

    std::shared_ptr<const GridGeometry> gridGeometryPtr() const
    { return gridGeometry_; }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    PrimaryVariables displacement(const GlobalPosition& x) const
    {
        PrimaryVariables u(0.0);
        F_.mv(x, u);
        u -= x;
        return u;
    }

    //! Saint Venant-Kirchhoff, P = F S with S = lambda tr(E) I + 2 mu E and E = (F^T F - I)/2
    Tensor firstPiolaKirchhoffStressTensor(const Tensor& F) const
    {
        Tensor E(0.0);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                for (int k = 0; k < dim; ++k)
                    E[i][j] += 0.5*F[k][i]*F[k][j];
        for (int i = 0; i < dim; ++i)
            E[i][i] -= 0.5;

        Scalar trE = 0.0;
        for (int i = 0; i < dim; ++i)
            trE += E[i][i];

        Tensor S = E;
        S *= 2.0*this->spatialParams().shearModulus();
        for (int i = 0; i < dim; ++i)
            S[i][i] += this->spatialParams().firstLameParameter()*trE;

        return F*S;
    }

private:
    std::shared_ptr<const GridGeometry> gridGeometry_;
    Tensor F_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

namespace Dumux::Python {

template<class Problem, class... options>
void registerHomogeneousDeformationProblem(pybind11::handle scope,
                                           pybind11::class_<Problem, options...> cls)
{
    using GridGeometry = typename Problem::GridGeometry;
    using Tensor = typename Problem::Tensor;
    cls.def(pybind11::init([](std::shared_ptr<const GridGeometry> gridGeometry,
                              const std::vector<std::vector<double>>& F){
        Tensor tensor;
        if (F.size() != Tensor::rows)
            DUNE_THROW(Dune::InvalidStateException, "F needs one row per dimension");
        for (std::size_t i = 0; i < Tensor::rows; ++i)
        {
            if (F[i].size() != Tensor::cols)
                DUNE_THROW(Dune::InvalidStateException, "F needs one column per dimension");
            for (std::size_t j = 0; j < Tensor::cols; ++j)
                tensor[i][j] = F[i][j];
        }
        return std::make_shared<Problem>(gridGeometry, tensor);
    }), pybind11::keep_alive<1, 2>());
    cls.def("gridGeometry", [](Problem& self){
        return std::const_pointer_cast<GridGeometry>(self.gridGeometryPtr());
    });
}

} // end namespace Dumux::Python

#endif
