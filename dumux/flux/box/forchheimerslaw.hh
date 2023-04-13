// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxFlux
 * \brief Forchheimers's law for the box method
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FORCHHEIMERS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_FORCHHEIMERS_LAW_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/typetraits/typetraits.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/cvfe/darcyslaw.hh>
#include <dumux/flux/facetensoraverage.hh>

namespace Dumux {

// forward declarations
template<class TypeTag, class ForchheimerVelocity, class DiscretizationMethod>
class ForchheimersLawImplementation;

/*!
 * \ingroup BoxFlux
 * \brief Forchheimer's law for box scheme
 * \note Forchheimer's law is specialized for network and surface grids (i.e. if grid dim < dimWorld)
 * \tparam Scalar the scalar type for scalar physical quantities
 * \tparam GridGeometry the grid geometry
 * \tparam ForchheimerVelocity class for the calculation of the Forchheimer velocity
 * \tparam isNetwork whether we are computing on a network grid embedded in a higher world dimension
 */
template<class Scalar, class GridGeometry, class ForchheimerVelocity>
class BoxForchheimersLaw;

/*!
 * \ingroup BoxFlux
 * \brief Forchheimer's law for box scheme
 */
template <class TypeTag, class ForchheimerVelocity>
class ForchheimersLawImplementation<TypeTag, ForchheimerVelocity, DiscretizationMethods::Box>
: public BoxForchheimersLaw<GetPropType<TypeTag, Properties::Scalar>,
                            GetPropType<TypeTag, Properties::GridGeometry>,
                            ForchheimerVelocity>
{};

/*!
 * \ingroup BoxFlux
 * \brief Specialization of the BoxForchheimersLaw
 */
template<class ScalarType, class GridGeometry, class ForchheimerVelocity>
class BoxForchheimersLaw
{
    using ThisType = BoxForchheimersLaw<ScalarType, GridGeometry, ForchheimerVelocity>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using DimWorldVector = typename ForchheimerVelocity::DimWorldVector;
    using DimWorldMatrix = typename ForchheimerVelocity::DimWorldMatrix;

    using DarcysLaw = CVFEDarcysLaw<ScalarType, GridGeometry>;

public:
    //! state the scalar type of the law
    using Scalar = ScalarType;

    using DiscretizationMethod = DiscretizationMethods::Box;
    //! state the discretization method this implementation belongs to
    static constexpr DiscretizationMethod discMethod{};

    /*! \brief Compute the advective flux of a phase across
    *          the given sub-control volume face using the Forchheimer equation.
    *
    *          The flux is given in N*m, and can be converted
    *          into a volume flux (m^3/s) or mass flux (kg/s) by applying an upwind scheme
    *          for the mobility or the product of density and mobility, respectively.
    */
    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];

        auto insideK = insideVolVars.permeability();
        auto outsideK = outsideVolVars.permeability();

        // scale with correct extrusion factor
        insideK *= insideVolVars.extrusionFactor();
        outsideK *= outsideVolVars.extrusionFactor();

        const auto K = faceTensorAverage(insideK, outsideK, scvf.unitOuterNormal());
        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");

        const auto& shapeValues = fluxVarCache.shapeValues();

        // evaluate gradP - rho*g at integration point
        DimWorldVector gradP(0.0);
        Scalar rho(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];

            if (enableGravity)
                rho += volVars.density(phaseIdx)*shapeValues[scv.indexInElement()][0];

            // the global shape function gradient
            gradP.axpy(volVars.pressure(phaseIdx), fluxVarCache.gradN(scv.indexInElement()));
        }

        if (enableGravity)
            gradP.axpy(-rho, problem.spatialParams().gravity(scvf.center()));

        DimWorldVector darcyVelocity = mv(K, gradP);
        darcyVelocity *= -1;

        auto upwindTerm = [phaseIdx](const auto& volVars){ return volVars.mobility(phaseIdx); };
        darcyVelocity *= ForchheimerVelocity::UpwindScheme::multiplier(elemVolVars, scvf, upwindTerm, darcyVelocity*scvf.unitOuterNormal(), phaseIdx);

        const auto velocity = ForchheimerVelocity::velocity(fvGeometry,
                                                            elemVolVars,
                                                            scvf,
                                                            phaseIdx,
                                                            calculateHarmonicMeanSqrtPermeability(problem, elemVolVars, scvf),
                                                            darcyVelocity);

        Scalar flux = velocity * scvf.unitOuterNormal();
        flux *= Extrusion::area(fvGeometry, scvf);
        return flux;
    }

    //! The flux variables cache has to be bound to an element prior to flux calculations
    //! During the binding, the transmissibility will be computed and stored using the method below.
    template<class Problem, class ElementVolumeVariables>
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const SubControlVolumeFace& scvf)
    {
        return DarcysLaw::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf);
    }

    /*! \brief Returns the harmonic mean of \f$\sqrt{K_0}\f$ and \f$\sqrt{K_1}\f$.
     *
     * This is a specialization for scalar-valued permeabilities which returns a tensor with identical diagonal entries.
     */
    template<class Problem, class ElementVolumeVariables>
    static DimWorldMatrix calculateHarmonicMeanSqrtPermeability(const Problem& problem,
                                                                const ElementVolumeVariables& elemVolVars,
                                                                const SubControlVolumeFace& scvf)
    {
        return ForchheimerVelocity::calculateHarmonicMeanSqrtPermeability(problem, elemVolVars, scvf);
    }
};

} // end namespace Dumux

#endif
