// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Approximation of the pressure Schur complement of the transient Stokes system
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_CHANNEL_TRANSIENT_PRESSURE_OPERATOR_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_CHANNEL_TRANSIENT_PRESSURE_OPERATOR_HH

#include <memory>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include <dumux/assembly/jacobianpattern.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief Approximates the pressure Schur complement of the transient Stokes system on the cells of the
 *        mass balance, for use with StokesSolver::setPressureMatrix and StokesSolver::setPressureDiagonal
 *
 * For a small time step the velocity block is dominated by its storage term. Eliminating the velocity with
 * that term alone gives the Poisson operator \f$ \Delta t \, \nabla\cdot(\rho^{-1} \nabla) \f$, discretized
 * with two-point fluxes. Boundaries contribute nothing, and a relative shift of the diagonal removes the
 * singularity of the constant pressure. The viscous part of the velocity block adds a mass matrix scaled by
 * the inverse viscosity; the inverse of that diagonal is the second part of the approximation
 * (Cahouet & Chabard, Int. J. Numer. Meth. Fluids 8 (1988) 869-895).
 */
template<class MassGridGeometry>
class TransientPressureOperator
{
    using Scalar = double;
    using Extrusion = Extrusion_t<MassGridGeometry>;

public:
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, 1, 1>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<Scalar, 1>>;

    TransientPressureOperator(std::shared_ptr<const MassGridGeometry> gridGeometry,
                              const Scalar density, const Scalar viscosity,
                              const Scalar regularization = 1e-8)
    : gridGeometry_(std::move(gridGeometry))
    , density_(density)
    , regularization_(regularization)
    , matrix_(std::make_shared<Matrix>())
    , viscousDiagonal_(std::make_shared<Vector>(gridGeometry_->numDofs()))
    {
        matrix_->setBuildMode(Matrix::random);
        getJacobianPattern<true>(*gridGeometry_).exportIdx(*matrix_);

        auto fvGeometry = localView(*gridGeometry_);
        for (const auto& element : elements(gridGeometry_->gridView()))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
                (*viscousDiagonal_)[scv.dofIndex()] = viscosity/Extrusion::volume(fvGeometry, scv);
        }
    }

    //! Assemble the operator for the given time step size
    void update(const Scalar timeStepSize)
    {
        auto& matrix = *matrix_;
        matrix = 0.0;

        auto fvGeometry = localView(*gridGeometry_);
        for (const auto& element : elements(gridGeometry_->gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                    continue;

                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                const auto distance = (outsideScv.center() - insideScv.center()).two_norm();
                const auto coefficient = timeStepSize*Extrusion::area(fvGeometry, scvf)/(distance*density_);

                matrix[insideScv.dofIndex()][insideScv.dofIndex()] += coefficient;
                matrix[insideScv.dofIndex()][outsideScv.dofIndex()] -= coefficient;
            }
        }

        for (std::size_t i = 0; i < matrix.N(); ++i)
            matrix[i][i] *= 1.0 + regularization_;
    }

    std::shared_ptr<const Matrix> matrix() const
    { return matrix_; }

    std::shared_ptr<const Vector> viscousDiagonal() const
    { return viscousDiagonal_; }

private:
    std::shared_ptr<const MassGridGeometry> gridGeometry_;
    Scalar density_;
    Scalar regularization_;
    std::shared_ptr<Matrix> matrix_;
    std::shared_ptr<Vector> viscousDiagonal_;
};

} // end namespace Dumux

#endif
