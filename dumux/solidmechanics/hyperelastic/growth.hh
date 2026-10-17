// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Hyperelastic
 * \brief Multiplicative decomposition of the deformation gradient into growth and elastic parts
 *
 * Growth, swelling, shrinkage and thermal expansion change the stress-free shape of the material
 * locally. Following @cite Rodriguez1994, the deformation gradient is decomposed as
 * \f[ \mathbf{F} = \mathbf{F}_\mathrm{e} \mathbf{F}_\mathrm{g}, \f]
 * where the growth tensor \f$ \mathbf{F}_\mathrm{g} \f$ maps the reference configuration to a
 * stress-free grown configuration, which need not be compatible, and the elastic part
 * \f$ \mathbf{F}_\mathrm{e} = \mathbf{F}\mathbf{F}_\mathrm{g}^{-1} \f$ restores compatibility.
 * With the strain energy density \f$ \psi(\mathbf{F}_\mathrm{e}) \f$ per unit volume of the grown
 * configuration, the strain energy per unit reference volume is
 * \f$ W(\mathbf{F}) = J_\mathrm{g}\, \psi(\mathbf{F}\mathbf{F}_\mathrm{g}^{-1}) \f$ with
 * \f$ J_\mathrm{g} = \det \mathbf{F}_\mathrm{g} \f$, and the first Piola-Kirchhoff stress is
 * \f[
 *     \mathbf{P} = \frac{\partial W}{\partial \mathbf{F}}
 *                = J_\mathrm{g}\, \mathbf{P}_\mathrm{e}(\mathbf{F}_\mathrm{e})\, \mathbf{F}_\mathrm{g}^{-T},
 *     \quad \mathbf{P}_\mathrm{e} = \frac{\partial \psi}{\partial \mathbf{F}_\mathrm{e}}.
 * \f]
 * A problem of the hyperelastic models implements its stress hook with
 * Dumux::Growth::firstPiolaKirchhoffStressTensor, passing the growth tensor at the integration
 * point and the elastic constitutive law.
 */
#ifndef DUMUX_SOLIDMECHANICS_HYPERELASTIC_GROWTH_HH
#define DUMUX_SOLIDMECHANICS_HYPERELASTIC_GROWTH_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>

namespace Dumux::Growth {

/*!
 * \ingroup Hyperelastic
 * \brief The first Piola-Kirchhoff stress \f$ J_\mathrm{g} \mathbf{P}_\mathrm{e}(\mathbf{F}\mathbf{F}_\mathrm{g}^{-1}) \mathbf{F}_\mathrm{g}^{-T} \f$
 *        of a material with growth tensor \f$ \mathbf{F}_\mathrm{g} \f$
 * \param F the deformation gradient
 * \param Fg the growth tensor
 * \param elasticStress the elastic constitutive law, a callable returning the first Piola-Kirchhoff
 *        stress \f$ \mathbf{P}_\mathrm{e} \f$ for the elastic part \f$ \mathbf{F}_\mathrm{e} \f$
 */
template<class Scalar, int dim, class ElasticStress>
Dune::FieldMatrix<Scalar, dim, dim> firstPiolaKirchhoffStressTensor(const Dune::FieldMatrix<Scalar, dim, dim>& F,
                                                                    const Dune::FieldMatrix<Scalar, dim, dim>& Fg,
                                                                    const ElasticStress& elasticStress)
{
    auto FgInv = Fg;
    FgInv.invert();

    Dune::FieldMatrix<Scalar, dim, dim> P = elasticStress(F*FgInv);
    P.rightmultiply(getTransposed(FgInv));
    P *= Fg.determinant();
    return P;
}

/*!
 * \ingroup Hyperelastic
 * \brief Isotropic growth \f$ \mathbf{F}_\mathrm{g} = \vartheta \mathbf{I} \f$
 * \param stretch the growth stretch \f$ \vartheta \f$ in every direction
 */
template<int dim, class Scalar>
Dune::FieldMatrix<Scalar, dim, dim> isotropicTensor(const Scalar stretch)
{
    Dune::FieldMatrix<Scalar, dim, dim> Fg(0.0);
    for (int i = 0; i < dim; ++i)
        Fg[i][i] = stretch;
    return Fg;
}

/*!
 * \ingroup Hyperelastic
 * \brief Transversely isotropic growth
 *        \f$ \mathbf{F}_\mathrm{g} = \vartheta_n \mathbf{n}\otimes\mathbf{n}
 *        + \vartheta_t (\mathbf{I} - \mathbf{n}\otimes\mathbf{n}) \f$
 *
 * Fibre growth along \f$ \mathbf{n} \f$ is \f$ \vartheta_t = 1 \f$, area growth of surfaces
 * normal to \f$ \mathbf{n} \f$ (e.g. tangential growth of a layer) is \f$ \vartheta_n = 1 \f$.
 * \param normal the unit vector \f$ \mathbf{n} \f$
 * \param normalStretch the growth stretch \f$ \vartheta_n \f$ along \f$ \mathbf{n} \f$
 * \param tangentialStretch the growth stretch \f$ \vartheta_t \f$ in all directions normal to \f$ \mathbf{n} \f$
 */
template<class Scalar, int dim>
Dune::FieldMatrix<Scalar, dim, dim> transverselyIsotropicTensor(const Dune::FieldVector<Scalar, dim>& normal,
                                                                const Scalar normalStretch,
                                                                const Scalar tangentialStretch)
{
    auto Fg = isotropicTensor<dim>(tangentialStretch);
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            Fg[i][j] += (normalStretch - tangentialStretch)*normal[i]*normal[j];
    return Fg;
}

/*!
 * \ingroup Hyperelastic
 * \brief Orthotropic growth
 *        \f$ \mathbf{F}_\mathrm{g} = \sum_i \vartheta_i \mathbf{a}_i\otimes\mathbf{a}_i \f$
 * \param frame the matrix whose columns are the orthonormal material directions \f$ \mathbf{a}_i \f$
 * \param stretches the growth stretches \f$ \vartheta_i \f$ along the material directions
 */
template<class Scalar, int dim>
Dune::FieldMatrix<Scalar, dim, dim> orthotropicTensor(const Dune::FieldMatrix<Scalar, dim, dim>& frame,
                                                      const Dune::FieldVector<Scalar, dim>& stretches)
{
    Dune::FieldMatrix<Scalar, dim, dim> Fg(0.0);
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            for (int k = 0; k < dim; ++k)
                Fg[i][j] += stretches[k]*frame[i][k]*frame[j][k];
    return Fg;
}

} // end namespace Dumux::Growth

#endif
