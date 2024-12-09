// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_MULTIDOMAIN_FLUIDSTRUCTURE_STRUCTURE_MESHMOTION_PROBLEM_HH
#define DUMUX_TEST_MULTIDOMAIN_FLUIDSTRUCTURE_STRUCTURE_MESHMOTION_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>

namespace Dumux {

// This test case corresponds to the CSM benchmark problem
// of a hyperelastic solid material under gravity given in
// Turek, S., Hron, J. (2006). "Proposal for Numerical Benchmarking of Fluid-Structure Interaction
// between an Elastic Object and Laminar Incompressible Flow."
// https://doi.org/10.1007/3-540-34596-5_15
template<class TypeTag>
class StructureProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    static constexpr int dim = GridView::dimension;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
public:
    StructureProblem(
        std::shared_ptr<const GridGeometry> gridGeometry,
        std::shared_ptr<const CouplingManager> couplingManager
    )
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , gravity_(getParam<Scalar>("Problem.Gravity"))
    {}

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        const auto r = std::hypot(globalPos[0]-0.2, globalPos[1]-0.2);
        if (r < 0.05 + eps_)
        {
            values.setDirichlet(0);
            values.setDirichlet(1);
        }
        else
            values.setAllNeumann();
        return values;
    }

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

    // the following methods are needed to solve structural dynamics
    // with the Newmark-beta time integration scheme

    // we use the Newmark scheme for time integration
    void setNewmarkScheme(std::shared_ptr<const Experimental::NewmarkBeta<Scalar, SolutionVector>> newmark)
    { newmark_ = std::move(newmark); }

    // the effective density of the solid material
    Scalar solidDensity(const Element&, const SubControlVolume&) const
    { return this->spatialParams().solidDensity(); }

    // the Newmark scheme is used for time integration and this
    // computes the acceleration at the current time step for us
    auto acceleration(const Element& element,
                      const SubControlVolume& scv,
                      const Scalar dt,
                      const PrimaryVariables& d) const
    {
        const auto dofIndex = scv.dofIndex();
        return newmark_->acceleration(dofIndex, dt, d);
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

private:
    static constexpr Scalar eps_ = 1e-7;
    std::shared_ptr<const Experimental::NewmarkBeta<Scalar, SolutionVector>> newmark_;
    std::shared_ptr<const CouplingManager> couplingManager_;
    Scalar gravity_;
};

} // end namespace Dumux

#endif
