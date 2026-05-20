// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEFlux
 * \brief Hooke's law for CVFE schemes. Integrates the stress tensor over the
 *        sub-control volume face using the scheme's quadrature rule, making
 *        it generic for all CVFE methods (Box, PQ1Bubble, PQ2, ...).
 */
#ifndef DUMUX_FLUX_CVFE_HOOKES_LAW_HH
#define DUMUX_FLUX_CVFE_HOOKES_LAW_HH

#include <dune/common/fmatrix.hh>

#include <dumux/common/math.hh>
#include <dumux/flux/hookeslaw.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/cvfe/localdof.hh>

namespace Dumux {

/*!
 * \ingroup CVFEFlux
 * \brief Hooke's law for all CVFE schemes.
 *
 * The stress integral over an scvf is computed by summing over all quadrature
 * points provided by the scheme's quadrature rule. For single-QP schemes (Box,
 * PQ1Bubble) this reduces to sigma(midpoint) * n * area, identical to the
 * former Box-only implementation. For multi-QP schemes (PQ2) it integrates
 * the linearly-varying stress exactly.
 */
template<class ScalarType, class GridGeometry>
class HookesLaw<ScalarType, GridGeometry, typename GridGeometry::DiscretizationMethod>
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static_assert(dim == dimWorld, "HookesLaw not implemented for network/surface grids");

public:
    using Scalar = ScalarType;
    using StressTensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;
    using ForceVector = typename StressTensor::row_type;

    using DiscretizationMethod = typename GridGeometry::DiscretizationMethod;
    static constexpr DiscretizationMethod discMethod{};

    /*!
     * \brief Returns the integrated force (in Newton) on a sub-control volume face.
     *        Sums contributions over all quadrature points of the face.
     */
    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static ForceVector force(const Problem& problem,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf,
                             const ElementFluxVarsCache& elemFluxVarCache)
    {
        ForceVector scvfForce(0.0);
        for (const auto& qpData : CVFE::quadratureRule(fvGeometry, scvf))
        {
            const auto sigma = stressTensor(problem, element, fvGeometry, elemVolVars,
                                            elemFluxVarCache[qpData.ipData()]);
            ForceVector contribution(0.0);
            sigma.mv(scvf.unitOuterNormal(), contribution);
            contribution *= qpData.weight();
            scvfForce += contribution;
        }
        return scvfForce;
    }

    //! assembles the stress tensor at a given integration point
    template<class Problem, class ElementVolumeVariables, class FluxVarsCache>
    static StressTensor stressTensor(const Problem& problem,
                                     const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables& elemVolVars,
                                     const FluxVarsCache& fluxVarCache)
    {
        const auto& lameParams = problem.spatialParams().lameParams(element, fvGeometry, elemVolVars, fluxVarCache);

        StressTensor gradU(0.0);
        for (int dir = 0; dir < dim; ++dir)
            for (const auto& localDof : localDofs(fvGeometry))
                gradU[dir].axpy(elemVolVars[localDof.index()].displacement(dir), fluxVarCache.gradN(localDof.index()));

        StressTensor epsilon;
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dimWorld; ++j)
                epsilon[i][j] = 0.5*(gradU[i][j] + gradU[j][i]);

        StressTensor sigma(0.0);
        const auto traceEpsilon = trace(epsilon);
        for (int i = 0; i < dim; ++i)
        {
            sigma[i][i] = lameParams.lambda()*traceEpsilon;
            for (int j = 0; j < dimWorld; ++j)
                sigma[i][j] += 2.0*lameParams.mu()*epsilon[i][j];
        }

        return sigma;
    }
};

} // end namespace Dumux

#endif
