// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup StaggeredFlux
 * \brief Specialization of Fourier's Law for the staggered free flow method.
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FOURIERS_LAW_HH

#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/fluxvariablescaching.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, class DiscretizationMethod>
class FouriersLawImplementation;

/*!
 * \ingroup StaggeredFlux
 * \brief Specialization of Fourier's Law for the staggered free flow method.
 */
template <class TypeTag>
class FouriersLawImplementation<TypeTag, DiscretizationMethods::Staggered>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

public:
    using DiscretizationMethod = DiscretizationMethods::Staggered;
    // state the discretization method this implementation belongs to
    static constexpr DiscretizationMethod discMethod{};

    //! state the type for the corresponding cache
    //! We don't cache anything for this law
    using Cache = FluxVariablesCaching::EmptyDiffusionCache;

    /*!
     * \brief Returns the heat flux within the porous medium
     *        (in J/s) across the given sub-control volume face.
     * \note This law assumes thermal equilibrium between the fluid
     *       and solid phases, and uses an effective thermal conductivity
     *       for the overall aggregate.
     */
    template<class Problem>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace &scvf)
    {
        Scalar flux(0.0);

        // conductive energy flux is zero for outflow boundary conditions
        if (scvf.boundary() && problem.boundaryTypes(element, scvf).isOutflow(Indices::energyEqIdx))
            return flux;

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        const Scalar insideTemperature = insideVolVars.temperature();
        const Scalar outsideTemperature = outsideVolVars.temperature();

        const Scalar insideLambda = insideVolVars.effectiveThermalConductivity() * insideVolVars.extrusionFactor();
        const Scalar insideDistance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();

        if (scvf.boundary())
        {
            flux = insideLambda * (insideTemperature - outsideTemperature) / insideDistance;
        }
        else
        {
            const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
            const Scalar outsideLambda = outsideVolVars.effectiveThermalConductivity() * outsideVolVars.extrusionFactor();
            const Scalar outsideDistance = (outsideScv.dofPosition() - scvf.ipGlobal()).two_norm();
            const Scalar avgLambda = harmonicMean(insideLambda, outsideLambda, insideDistance, outsideDistance);

            flux = avgLambda * (insideTemperature - outsideTemperature) / (insideDistance + outsideDistance);
        }

        flux *= Extrusion::area(fvGeometry, scvf);
        return flux;
    }
};

} // end namespace Dumux

#endif
