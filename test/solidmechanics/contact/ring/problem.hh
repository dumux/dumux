// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief A thin elastic ring (plane strain) collapsing under external pressure
 *
 * The ring of midline radius \f$ R \f$ and wall thickness \f$ t \f$ is loaded on its outer surface
 * by the follower pressure \f$ p (1 - \epsilon \cos 2\varphi) \f$, whose small imperfection
 * \f$ \epsilon \f$ selects the elliptic mode with the minor axis along \f$ y \f$. The nodes on the
 * coordinate axes are held on them, which removes the rigid motions without restricting the doubly
 * symmetric deformation. The material is neo-Hookean,
 * \f$ \mathbf{P} = \mu(\mathbf{F} - \mathbf{F}^{-T}) + \lambda \ln J\, \mathbf{F}^{-T} \f$, with the
 * Lame parameters of the Young's modulus and Poisson ratio.
 */
#ifndef DUMUX_TEST_CONTACT_RING_PROBLEM_HH
#define DUMUX_TEST_CONTACT_RING_PROBLEM_HH

#include <cmath>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/boundarytypes_.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/math.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

template<class TypeTag>
class RingCollapseProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<PrimaryVariables::size()>;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using GridIndex = typename IndexTraits<GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<Dumux::DirichletConstraintInfo<numEq>, PrimaryVariables, GridIndex>;
    static constexpr int dim = GridView::dimension;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

public:
    RingCollapseProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , radius_(getParam<Scalar>("Problem.Radius"))
    , thickness_(getParam<Scalar>("Problem.Thickness"))
    , youngsModulus_(getParam<Scalar>("Problem.YoungsModulus"))
    , poissonRatio_(getParam<Scalar>("Problem.PoissonRatio"))
    , imperfection_(getParam<Scalar>("Problem.Imperfection"))
    {
        const Scalar eps = 1e-8*radius_;
        std::vector<bool> visited(gridGeometry->numDofs(), false);
        auto elemDisc = localView(*gridGeometry);
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            elemDisc.bind(element);
            for (const auto& localDof : localDofs(elemDisc))
            {
                const GridIndex dofIdx = localDof.dofIndex();
                if (visited[dofIdx])
                    continue;
                visited[dofIdx] = true;
                const auto pos = ipData(elemDisc, localDof).global();
                const bool onYAxis = std::abs(pos[0]) < eps;
                const bool onXAxis = std::abs(pos[1]) < eps;
                if (!onXAxis && !onYAxis)
                    continue;
                typename DirichletConstraintData::ConstraintInfo info;
                if (onYAxis)
                    info.set(0);
                if (onXAxis)
                    info.set(1);
                constraints_.push_back(DirichletConstraintData{std::move(info), PrimaryVariables(0.0), dofIdx});
            }
        }
    }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    void setPressure(Scalar pressure)
    { pressure_ = pressure; }

    //! The bending stiffness per unit length of the wall in plane strain
    Scalar bendingStiffness() const
    { return youngsModulus_*thickness_*thickness_*thickness_/(12.0*(1.0 - poissonRatio_*poissonRatio_)); }

    Scalar radius() const
    { return radius_; }

    Scalar thickness() const
    { return thickness_; }

    Tensor firstPiolaKirchhoffStressTensor(const Tensor& F) const
    {
        const Scalar mu = youngsModulus_/(2.0*(1.0 + poissonRatio_));
        const Scalar lambda = youngsModulus_*poissonRatio_/((1.0 + poissonRatio_)*(1.0 - 2.0*poissonRatio_));
        auto FinvT = F;
        FinvT.invert();
        FinvT = getTransposed(FinvT);
        auto P = F;
        P -= FinvT;
        P *= mu;
        P.axpy(lambda*std::log(F.determinant()), FinvT);
        return P;
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    {
        BoundaryTypes values;
        values.setAllFluxBoundary();
        return values;
    }

    //! The pressure on the deformed outer surface, \f$ p\, \mathrm{cof}(\mathbf{F}) \mathbf{N} \f$ per reference area
    template<class ElementVariables, class FaceIpData>
    NumEqVector boundaryFlux(const ElementDiscretization& elemDisc,
                             const ElementVariables& elemVars,
                             const FaceIpData& ipData) const
    {
        const auto& pos = ipData.global();
        if (pos.two_norm() < radius_)
            return NumEqVector(0.0);

        const auto& ipCache = cache(elemVars, ipData);
        Tensor F(0.0);
        for (const auto& localDof : localDofs(elemDisc))
            for (int i = 0; i < dim; ++i)
                F[i].axpy(elemVars[localDof].displacement(i), ipCache.gradN(localDof.index()));
        for (int i = 0; i < dim; ++i)
            F[i][i] += 1.0;

        const Tensor cofactor{{F[1][1], -F[1][0]}, {-F[0][1], F[0][0]}};
        const Scalar cos2phi = (pos[0]*pos[0] - pos[1]*pos[1])/pos.two_norm2();
        NumEqVector flux(0.0);
        cofactor.umv(ipData.unitOuterNormal(), flux);
        flux *= pressure_*(1.0 - imperfection_*cos2phi);
        return flux;
    }

private:
    Scalar radius_, thickness_, youngsModulus_, poissonRatio_, imperfection_;
    Scalar pressure_ = 0.0;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
