// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxFlux
 * \brief This file contains the data which is required to calculate
 *        energy fluxes due to molecular diffusion with Fourier's law.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_FOURIERS_LAW_HH

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/facetensoraverage.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class DiscretizationMethod>
class FouriersLawImplementation;

/*!
 * \ingroup BoxFlux
 * \brief Specialization of Fourier's Law for the box method.
 */
template <class TypeTag>
class FouriersLawImplementation<TypeTag, DiscretizationMethods::Box>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    /*!
     * \brief Returns the heat flux within the porous medium
     *        (in J/s) across the given sub-control volume face.
     * \note This law assumes thermal equilibrium between the fluid
     *       and solid phases, and uses an effective thermal conductivity
     *       for the overall aggregate.
     */
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        // get inside and outside diffusion tensors and calculate the harmonic mean
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];

        // effective diffusion tensors
        auto insideLambda = insideVolVars.effectiveThermalConductivity();
        auto outsideLambda = outsideVolVars.effectiveThermalConductivity();

        // scale by extrusion factor
        insideLambda *= insideVolVars.extrusionFactor();
        outsideLambda *= outsideVolVars.extrusionFactor();

        // the resulting averaged diffusion tensor
        const auto lambda = faceTensorAverage(insideLambda, outsideLambda, scvf.unitOuterNormal());

        // evaluate gradTemp at integration point
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // compute the temperature gradient with the shape functions
        Dune::FieldVector<Scalar, GridView::dimensionworld> gradTemp(0.0);
        for (auto&& scv : scvs(fvGeometry))
            gradTemp.axpy(elemVolVars[scv].temperature(), fluxVarsCache.gradN(scv.indexInElement()));

        // compute the heat conduction flux
        return -1.0*vtmv(scvf.unitOuterNormal(), lambda, gradTemp)*Extrusion::area(fvGeometry, scvf);
    }
};

} // end namespace Dumux

#endif
