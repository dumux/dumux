// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxFlux
 * \brief Specialization of the effective stress law for the box scheme. This computes the stress
 *        tensor and surface forces resulting from mechanical deformation and the pore pressure.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_EFFECTIVE_STRESS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_EFFECTIVE_STRESS_LAW_HH

#include <dumux/flux/effectivestresslaw.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

namespace Dumux {

/*!
 * \ingroup BoxFlux
 * \brief Effective stress law for box scheme
 * \tparam StressType type used for the computation of
 *         purely mechanical stresses (i.e. material law)
 * \tparam GridGeometry the finite volume grid geometry
 */
template<class StressType, class GridGeometry>
class EffectiveStressLaw<StressType, GridGeometry, typename GridGeometry::DiscretizationMethod>
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static_assert(dim == dimWorld, "EffectiveStressLaw not implemented for network/surface grids");
    static_assert(StressType::discMethod == DiscretizationMethods::box, "The provided stress type must be specialized for the box scheme");

public:
    //! export the type used for scalar values
    using Scalar = typename StressType::Scalar;
    //! export the type used for the stress tensor
    using StressTensor = typename StressType::StressTensor;
    //! export the type used for force vectors
    using ForceVector = typename StressType::ForceVector;
    //! state the discretization method this implementation belongs to

    using DiscretizationMethod = DiscretizationMethods::Box;
    // state the discretization method this implementation belongs to
    static constexpr DiscretizationMethod discMethod{};

    /*!
     * \brief Computes the force (in Newton) acting on a sub-control volume face.
     */
    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static ForceVector force(const Problem& problem,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf,
                             const ElementFluxVarsCache& elemFluxVarCache)
    {
        const auto sigma = stressTensor(problem, element, fvGeometry, elemVolVars, elemFluxVarCache[scvf]);

        ForceVector scvfForce(0.0);
        sigma.mv(scvf.unitOuterNormal(), scvfForce);
        scvfForce *= Extrusion::area(fvGeometry, scvf);

        return scvfForce;
    }

    //! assembles the (total) stress tensor of the porous medium at a given integration point
    template<class Problem, class ElementVolumeVariables, class FluxVarsCache>
    static StressTensor stressTensor(const Problem& problem,
                                     const Element& element,
                                     const FVElementGeometry& fvGeometry,
                                     const ElementVolumeVariables& elemVolVars,
                                     const FluxVarsCache& fluxVarsCache)
    {
        // compute the purely mechanical stress
        auto sigma = StressType::stressTensor(problem, element, fvGeometry, elemVolVars, fluxVarsCache);

        // obtain biot coefficient and effective pore pressure
        const auto biotCoeff = problem.spatialParams().biotCoefficient(element, fvGeometry, elemVolVars, fluxVarsCache);
        const auto effPress = problem.spatialParams().effectivePorePressure(element, fvGeometry, elemVolVars, fluxVarsCache);

        // subtract pore pressure from the diagonal entries
        const auto bcp = biotCoeff*effPress;
        for (int i = 0; i < dim; ++i)
            sigma[i][i] -= bcp;

        return sigma;
    }

    //! assembles the (effective) stress tensor of the solid skeleton at a given integration point
    template<class Problem, class ElementVolumeVariables, class FluxVarsCache>
    static StressTensor effectiveStressTensor(const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const ElementVolumeVariables& elemVolVars,
                                              const FluxVarsCache& fluxVarsCache)
    { return StressType::stressTensor(problem, element, fvGeometry, elemVolVars, fluxVarsCache); }
};

} // end namespace Dumux

#endif
