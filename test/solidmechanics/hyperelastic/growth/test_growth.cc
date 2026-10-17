// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test of the growth tensors and of the first Piola-Kirchhoff stress with multiplicative growth
 *
 * The stress is compared with central differences of the strain energy
 * \f$ W(\mathbf{F}) = J_\mathrm{g}\, \psi(\mathbf{F}\mathbf{F}_\mathrm{g}^{-1}) \f$,
 * and a deformation \f$ \mathbf{F} = \mathbf{R}\mathbf{F}_\mathrm{g} \f$ with a rotation
 * \f$ \mathbf{R} \f$ has to be stress free.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <random>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/solidmechanics/hyperelastic/growth.hh>

namespace Dumux {

template<int dim>
using Tensor = Dune::FieldMatrix<double, dim, dim>;

template<int dim>
using Vector = Dune::FieldVector<double, dim>;

constexpr double mu = 1.3;
constexpr double lambda = 4.1;

// compressible neo-Hookean strain energy and its first Piola-Kirchhoff stress
template<int dim>
double strainEnergy(const Tensor<dim>& F)
{
    using std::log;
    const double lnJ = log(F.determinant());
    return 0.5*mu*(F.frobenius_norm2() - dim) - mu*lnJ + 0.5*lambda*lnJ*lnJ;
}

template<int dim>
Tensor<dim> stress(const Tensor<dim>& F)
{
    using std::log;
    auto FinvT = F;
    FinvT.invert();
    FinvT = getTransposed(FinvT);
    return mu*F + (lambda*log(F.determinant()) - mu)*FinvT;
}

template<int dim>
Tensor<dim> growthStressByFiniteDifferences(const Tensor<dim>& F, const Tensor<dim>& Fg)
{
    auto FgInv = Fg;
    FgInv.invert();
    const auto energy = [&](const Tensor<dim>& G) { return Fg.determinant()*strainEnergy<dim>(G*FgInv); };

    constexpr double h = 1e-6;
    Tensor<dim> P(0.0);
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
        {
            auto Fp = F, Fm = F;
            Fp[i][j] += h;
            Fm[i][j] -= h;
            P[i][j] = (energy(Fp) - energy(Fm))/(2.0*h);
        }
    return P;
}

template<int dim>
Tensor<dim> randomRotation(std::mt19937& gen)
{
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    Tensor<dim> Q;
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            Q[i][j] = dist(gen);

    // Gram-Schmidt on the columns
    for (int k = 0; k < dim; ++k)
    {
        for (int l = 0; l < k; ++l)
        {
            double proj = 0.0;
            for (int i = 0; i < dim; ++i)
                proj += Q[i][k]*Q[i][l];
            for (int i = 0; i < dim; ++i)
                Q[i][k] -= proj*Q[i][l];
        }
        double norm = 0.0;
        for (int i = 0; i < dim; ++i)
            norm += Q[i][k]*Q[i][k];
        for (int i = 0; i < dim; ++i)
            Q[i][k] /= std::sqrt(norm);
    }

    if (Q.determinant() < 0.0)
        for (int i = 0; i < dim; ++i)
            Q[i][0] *= -1.0;
    return Q;
}

template<int dim>
Tensor<dim> randomDeformationGradient(std::mt19937& gen)
{
    std::uniform_real_distribution<double> dist(-0.2, 0.2);
    Tensor<dim> F;
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            F[i][j] = (i == j ? 1.0 : 0.0) + dist(gen);
    return F;
}

template<int dim>
void checkEigenpair(const Tensor<dim>& Fg, const Vector<dim>& v, double stretch, const std::string& name)
{
    Vector<dim> Fgv;
    Fg.mv(v, Fgv);
    Fgv.axpy(-stretch, v);
    if (Fgv.two_norm() > 1e-14)
        DUNE_THROW(Dune::Exception, name << ": growth tensor does not stretch " << v << " by " << stretch);
}

template<int dim>
void checkSymmetric(const Tensor<dim>& Fg, const std::string& name)
{
    if ((Fg - getTransposed(Fg)).infinity_norm() > 1e-15)
        DUNE_THROW(Dune::Exception, name << ": growth tensor is not symmetric");
}

template<int dim>
void checkStress(const Tensor<dim>& Fg, std::mt19937& gen, const std::string& name)
{
    const auto F = randomDeformationGradient<dim>(gen);
    const auto P = Growth::firstPiolaKirchhoffStressTensor(F, Fg, stress<dim>);
    const auto PRef = growthStressByFiniteDifferences<dim>(F, Fg);
    const auto error = (P - PRef).infinity_norm()/PRef.infinity_norm();
    std::cout << dim << "D " << name << ": relative deviation from the energy derivative " << error << std::endl;
    if (error > 1e-7)
        DUNE_THROW(Dune::Exception, name << ": stress deviates from the derivative of the strain energy by " << error);

    const auto R = randomRotation<dim>(gen);
    const auto PRotated = Growth::firstPiolaKirchhoffStressTensor(R*Fg, Fg, stress<dim>);
    if (PRotated.infinity_norm() > 1e-13*mu)
        DUNE_THROW(Dune::Exception, name << ": rotated growth deformation is not stress free, |P| = " << PRotated.infinity_norm());
}

template<int dim>
void runTests()
{
    std::mt19937 gen(42);

    const auto isotropic = Growth::isotropicTensor<dim>(1.3);
    checkSymmetric(isotropic, "isotropic");
    for (int i = 0; i < dim; ++i)
    {
        Vector<dim> e(0.0);
        e[i] = 1.0;
        checkEigenpair(isotropic, e, 1.3, "isotropic");
    }
    checkStress(isotropic, gen, "isotropic");

    const auto Q = randomRotation<dim>(gen);
    Vector<dim> n;
    for (int i = 0; i < dim; ++i)
        n[i] = Q[i][0];
    const auto transverse = Growth::transverselyIsotropicTensor(n, 1.2, 0.8);
    checkSymmetric(transverse, "transversely isotropic");
    checkEigenpair(transverse, n, 1.2, "transversely isotropic");
    for (int k = 1; k < dim; ++k)
    {
        Vector<dim> t;
        for (int i = 0; i < dim; ++i)
            t[i] = Q[i][k];
        checkEigenpair(transverse, t, 0.8, "transversely isotropic");
    }
    checkStress(transverse, gen, "transversely isotropic");

    Vector<dim> stretches;
    for (int k = 0; k < dim; ++k)
        stretches[k] = 0.9 + 0.2*k;
    const auto orthotropic = Growth::orthotropicTensor(Q, stretches);
    checkSymmetric(orthotropic, "orthotropic");
    for (int k = 0; k < dim; ++k)
    {
        Vector<dim> a;
        for (int i = 0; i < dim; ++i)
            a[i] = Q[i][k];
        checkEigenpair(orthotropic, a, stretches[k], "orthotropic");
    }
    checkStress(orthotropic, gen, "orthotropic");

    // growth tensors need not be symmetric, e.g. a growth by shear
    auto shear = Growth::isotropicTensor<dim>(1.0);
    shear[0][1] = 0.4;
    checkStress(shear, gen, "shear");
}

} // end namespace Dumux

int main()
{
    Dumux::runTests<2>();
    Dumux::runTests<3>();
    return 0;
}
