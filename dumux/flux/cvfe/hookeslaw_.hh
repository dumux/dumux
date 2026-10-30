// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEFlux
 * \brief Hooke's law for control-volume finite element schemes, evaluated at an interpolation point
 */
#ifndef DUMUX_FLUX_CVFE_HOOKES_LAW__HH
#define DUMUX_FLUX_CVFE_HOOKES_LAW__HH

#include <dune/common/fmatrix.hh>

#include <dumux/common/concepts/ipdata_.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup CVFEFlux
 * \brief Hooke's law for control-volume finite element schemes, evaluated at an interpolation point
 *
 * The stress of a linear elastic isotropic solid is
 * \f$ \boldsymbol{\sigma} = \lambda \operatorname{tr}(\boldsymbol{\varepsilon}) \mathbf{I} + 2\mu \boldsymbol{\varepsilon} \f$
 * with the strain \f$ \boldsymbol{\varepsilon} = \frac{1}{2}(\nabla \mathbf{u} + \nabla \mathbf{u}^T) \f$ of the
 * interpolated displacement. The spatial parameters provide the Lamé parameters,
 * `lameParams(elemDisc, ipData)`.
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \tparam GridDiscretization the grid discretization
 */
template<class Scalar, class GridDiscretization>
class CVFEHookesLawAtIp
{
    using GridView = typename GridDiscretization::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementDiscretization = typename GridDiscretization::LocalView;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static_assert(dim == dimWorld, "Hooke's law is not implemented for network or surface grids");

public:
    using StressTensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;

    /*!
     * \brief Returns the stress tensor at an interpolation point
     * \note Looks up the cache of the interpolation point, see the overload below.
     */
    template<class Problem, class ElementVariables, Concept::IpData IpData>
    static StressTensor stressTensor(const Problem& problem,
                                     const Element& element,
                                     const ElementDiscretization& elemDisc,
                                     const ElementVariables& elemVars,
                                     const IpData& ipData)
    { return stressTensor(problem, element, elemDisc, elemVars, ipData, cache(elemVars, ipData)); }

    //! Returns the stress tensor at an interpolation point
    template<class Problem, class ElementVariables, Concept::IpData IpData, class IpCache>
    static StressTensor stressTensor(const Problem& problem,
                                     const Element& element,
                                     const ElementDiscretization& elemDisc,
                                     const ElementVariables& elemVars,
                                     const IpData& ipData,
                                     const IpCache& ipCache)
    {
        const auto gradU = displacementGradient(elemDisc, elemVars, ipCache);
        const auto& lame = problem.spatialParams().lameParams(elemDisc, ipData);

        Scalar trace = 0.0;
        for (int i = 0; i < dim; ++i)
            trace += gradU[i][i];

        StressTensor sigma(0.0);
        for (int i = 0; i < dim; ++i)
        {
            for (int j = 0; j < dimWorld; ++j)
                sigma[i][j] = lame.mu()*(gradU[i][j] + gradU[j][i]);
            sigma[i][i] += lame.lambda()*trace;
        }

        return sigma;
    }

    //! Returns the gradient of the interpolated displacement at an interpolation point, one row per direction
    template<class ElementVariables, class IpCache>
    static StressTensor displacementGradient(const ElementDiscretization& elemDisc,
                                             const ElementVariables& elemVars,
                                             const IpCache& ipCache)
    {
        StressTensor gradU(0.0);
        for (const auto& localDof : localDofs(elemDisc))
            for (int dir = 0; dir < dim; ++dir)
                gradU[dir].axpy(elemVars[localDof].displacement(dir), ipCache.gradN(localDof.index()));
        return gradU;
    }
};

} // end namespace Dumux::Experimental

#endif
