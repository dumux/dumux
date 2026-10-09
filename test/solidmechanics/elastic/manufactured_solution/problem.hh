// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsTests
 * \brief Definition of a test problem for the linear elastic model with a manufactured solution.
 *
 * On the unit square, \f$ \mathbf{u} = ((x - x^2)\sin(2\pi y), \sin(2\pi x)\sin(2\pi y)) \f$ solves
 * \f$ -\nabla\cdot\boldsymbol{\sigma} = \mathbf{f} \f$. The displacement vanishes on the boundary; on the
 * left edge the traction in \f$ x \f$ direction and on the upper edge the traction in \f$ y \f$ direction
 * of the exact solution are prescribed instead.
 */
#ifndef DUMUX_ELASTICPROBLEM_HH
#define DUMUX_ELASTICPROBLEM_HH

#include <cmath>
#include <vector>

#include <dune/common/fmatrix.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/math.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

/*!
 * \ingroup GeomechanicsTests
 * \brief Problem definition for the deformation of an elastic body.
 */
template<class TypeTag>
class ElasticProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<PrimaryVariables::size()>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using ConstraintInfo = Dumux::DirichletConstraintInfo<numEq>;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, PrimaryVariables,
                                                                   typename IndexTraits<GridView>::GridIndex>;

    static constexpr Scalar pi = M_PI;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GradU = Dune::FieldMatrix<Scalar, dim, dimWorld>;

public:
    ElasticProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    { appendDirichletConstraints_(); }

    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

    //! The traction is prescribed in x direction on the left and in y direction on the upper edge
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] < eps_ && globalPos[1] > eps_ && globalPos[1] < 1.0 - eps_)
            values.setFluxBoundary(Indices::momentum(/*x-dir*/0));
        if (globalPos[1] > 1.0 - eps_ && globalPos[0] > eps_ && globalPos[0] < 1.0 - eps_)
            values.setFluxBoundary(Indices::momentum(/*y-dir*/1));
        return values;
    }

    //! The traction of the exact solution, \f$ \mathbf{g} = -\boldsymbol{\sigma}\mathbf{n} \f$
    template<class ElementVariables, class FaceIpData>
    NumEqVector boundaryFlux(const ElementDiscretization&, const ElementVariables&, const FaceIpData& ipData) const
    {
        const auto& globalPos = ipData.global();
        GradU gradU = exactGradient(globalPos);
        GradU epsilon;
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dimWorld; ++j)
                epsilon[i][j] = 0.5*(gradU[i][j] + gradU[j][i]);

        const auto& lameParams = this->spatialParams().lameParamsAtPos(globalPos);
        GradU sigma(0.0);
        const auto traceEpsilon = trace(epsilon);
        for (int i = 0; i < dim; ++i)
        {
            sigma[i][i] = lameParams.lambda()*traceEpsilon;
            for (int j = 0; j < dimWorld; ++j)
                sigma[i][j] += 2.0*lameParams.mu()*epsilon[i][j];
        }

        NumEqVector values;
        sigma.mv(ipData.unitOuterNormal(), values);
        values *= -1.0;
        return values;
    }

    //! The body force \f$ \mathbf{f} = -\nabla\cdot\boldsymbol{\sigma} \f$ of the exact solution
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector source = divSigma_(globalPos);
        source *= -1.0;
        return source;
    }

    /*!
     * \brief Evaluates the exact displacement to this problem at a given position.
     */
    PrimaryVariables exactSolution(const GlobalPosition& globalPos) const
    {
        using std::sin;

        const auto x = globalPos[0];
        const auto y = globalPos[1];

        PrimaryVariables exact(0.0);
        exact[Indices::momentum(/*x-dir*/0)] = (x-x*x)*sin(2*pi*y);
        exact[Indices::momentum(/*y-dir*/1)] = sin(2*pi*x)*sin(2*pi*y);
        return exact;
    }

    /*!
     * \brief Evaluates the exact displacement gradient to this problem at a given position.
     */
    GradU exactGradient(const GlobalPosition& globalPos) const
    {
        using std::sin;
        using std::cos;

        const auto x = globalPos[0];
        const auto y = globalPos[1];

        static constexpr int xIdx = Indices::momentum(/*x-dir*/0);
        static constexpr int yIdx = Indices::momentum(/*y-dir*/1);

        GradU exactGrad(0.0);
        exactGrad[xIdx][xIdx] = (1-2*x)*sin(2*pi*y);
        exactGrad[xIdx][yIdx] = (x - x*x)*2*pi*cos(2*pi*y);
        exactGrad[yIdx][xIdx] = 2*pi*cos(2*pi*x)*sin(2*pi*y);
        exactGrad[yIdx][yIdx] = 2*pi*sin(2*pi*x)*cos(2*pi*y);
        return exactGrad;
    }

private:
    //! The displacement components without prescribed traction are constrained to the exact solution
    void appendDirichletConstraints_()
    {
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);

            for (const auto& boundaryFace : boundaryFaces(elemDisc))
            {
                const auto bcTypes = this->boundaryTypes(elemDisc, boundaryFace);
                if (bcTypes.hasOnlyFluxBoundary())
                    continue;

                ConstraintInfo info;
                for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
                    if (!bcTypes.isFluxBoundary(eqIdx))
                        info.set(eqIdx);

                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                {
                    const auto& globalPos = ipData(elemDisc, localDof).global();
                    constraints_.push_back(
                        DirichletConstraintData{info, exactSolution(globalPos), localDof.dofIndex()}
                    );
                }
            }
        }
    }

    //! The divergence of the stress of the exact solution
    PrimaryVariables divSigma_(const GlobalPosition& globalPos) const
    {
        using std::sin;
        using std::cos;

        const auto x = globalPos[0];
        const auto y = globalPos[1];

        // the lame parameters (we know they only depend on the position here)
        const auto& lameParams = this->spatialParams().lameParamsAtPos(globalPos);
        const auto lambda = lameParams.lambda();
        const auto mu = lameParams.mu();

        // precalculated products
        const Scalar pi_2 = 2.0*pi;
        const Scalar pi_2_square = pi_2*pi_2;
        const Scalar cos_2pix = cos(pi_2*x);
        const Scalar sin_2pix = sin(pi_2*x);
        const Scalar cos_2piy = cos(pi_2*y);
        const Scalar sin_2piy = sin(pi_2*y);

        const Scalar dE11_dx = -2.0*sin_2piy;
        const Scalar dE22_dx = pi_2_square*cos_2pix*cos_2piy;
        const Scalar dE11_dy = pi_2*(1.0-2.0*x)*cos_2piy;
        const Scalar dE22_dy = -1.0*pi_2_square*sin_2pix*sin_2piy;
        const Scalar dE12_dy = 0.5*pi_2_square*(cos_2pix*cos_2piy - (x-x*x)*sin_2piy);
        const Scalar dE21_dx = 0.5*((1.0-2*x)*pi_2*cos_2piy - pi_2_square*sin_2pix*sin_2piy);

        // compute exact divergence of sigma
        PrimaryVariables divSigma(0.0);
        divSigma[Indices::momentum(/*x-dir*/0)] = lambda*(dE11_dx + dE22_dx) + 2*mu*(dE11_dx + dE12_dy);
        divSigma[Indices::momentum(/*y-dir*/1)] = lambda*(dE11_dy + dE22_dy) + 2*mu*(dE21_dx + dE22_dy);
        return divSigma;
    }

    static constexpr Scalar eps_ = 3e-6;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
