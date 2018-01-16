// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup CCTpfaDiscretization
 * \brief Darcy's law for cell-centered finite volume schemes with two-point flux approximation
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_DARCYS_LAW_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

namespace Dumux
{
// forward declarations
template<class TypeTag, DiscretizationMethods discMethod>
class DarcysLawImplementation;

template<class TypeTag, bool isNetwork>
class CCTpfaDarcysLaw;

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Darcy's law for cell-centered finite volume schemes with two-point flux approximation
 * \note Darcy's law is speialized for network and surface grids (i.e. if grid dim < dimWorld)
 */
template <class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethods::CCTpfa>
: public CCTpfaDarcysLaw<TypeTag, (GET_PROP_TYPE(TypeTag, Grid)::dimension < GET_PROP_TYPE(TypeTag, Grid)::dimensionworld) >
{};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Class that fills the cache corresponding to tpfa Darcy's Law
 */
template<class TypeTag>
class TpfaDarcysLawCacheFiller
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

public:
    //! Function to fill a TpfaDarcysLawCache of a given scvf
    //! This interface has to be met by any advection-related cache filler class
    template<class FluxVariablesCacheFiller>
    static void fill(FluxVariablesCache& scvfFluxVarsCache,
                     const Problem& problem,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolumeFace& scvf,
                     const FluxVariablesCacheFiller& fluxVarsCacheFiller)
    {
        scvfFluxVarsCache.updateAdvection(problem, element, fvGeometry, elemVolVars, scvf);
    }
};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief The cache corresponding to tpfa Darcy's Law
 */
template<class TypeTag>
class TpfaDarcysLawCache
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

public:
    using Filler = TpfaDarcysLawCacheFiller<TypeTag>;

    void updateAdvection(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace &scvf)
    {
        tij_ = AdvectionType::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf);
    }

    const Scalar& advectionTij() const
    { return tij_; }

private:
    Scalar tij_;
};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Specialization of the CCTpfaDarcysLaw grids where dim=dimWorld
 */
template<class TypeTag>
class CCTpfaDarcysLaw<TypeTag, /*isNetwork*/ false>
{
    using Implementation = DarcysLawImplementation<TypeTag, DiscretizationMethods::CCTpfa>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

  public:
    //! state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCTpfa;

    //! state the type for the corresponding cache
    using Cache = TpfaDarcysLawCache<TypeTag>;

    //! Compute the advective flux
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        static const bool enableGravity = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Problem.EnableGravity");

        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        if (enableGravity)
        {
            // do averaging for the density over all neighboring elements
            const auto rho = scvf.boundary() ? outsideVolVars.density(phaseIdx)
                                             : (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;

            // Obtain inside and outside pressures
            const auto pInside = insideVolVars.pressure(phaseIdx);
            const auto pOutside = outsideVolVars.pressure(phaseIdx);

            const auto& tij = fluxVarsCache.advectionTij();
            const auto& g = problem.gravityAtPos(scvf.ipGlobal());

            //! compute alpha := n^T*K*g
            const auto alpha_inside = vtmv(scvf.unitOuterNormal(), insideVolVars.permeability(), g)*insideVolVars.extrusionFactor();

            Scalar flux = tij*(pInside - pOutside) + rho*scvf.area()*alpha_inside;

            //! On interior faces we have to add K-weighted gravitational contributions
            if (!scvf.boundary())
            {
                const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                const auto outsideK = outsideVolVars.permeability();
                const auto outsideTi = computeTpfaTransmissibility(scvf, outsideScv, outsideK, outsideVolVars.extrusionFactor());
                const auto alpha_outside = vtmv(scvf.unitOuterNormal(), outsideK, g)*outsideVolVars.extrusionFactor();

                flux += rho*tij/outsideTi*(alpha_inside - alpha_outside);
            }

            return flux;
        }
        else
        {
            // Obtain inside and outside pressures
            const auto pInside = insideVolVars.pressure(phaseIdx);
            const auto pOutside = outsideVolVars.pressure(phaseIdx);

            // return flux
            return fluxVarsCache.advectionTij()*(pInside - pOutside);
        }
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibility will be computed and stored using the method below.
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const SubControlVolumeFace& scvf)
    {
        Scalar tij;

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        // check if we evaluate the permeability in the volume (for discontinuous fields, default)
        // or at the scvf center for analytical permeability fields (e.g. convergence studies)
        auto getPermeability = [&problem](const VolumeVariables& volVars,
                                          const GlobalPosition& scvfIpGlobal) -> typename SpatialParams::PermeabilityType
                               {
                                    if (GET_PROP_VALUE(TypeTag, EvaluatePermeabilityAtScvfIP))
                                    {
                                        // TODO: Make PermeabilityType a property!!
                                        // We do an implicit cast to permeability type in case EvaluatePermeabilityAtScvfIP is not true so that it compiles
                                        // If it is true make sure that no cast is necessary and the correct return type is specified in the spatial params
                                        static_assert(!GET_PROP_VALUE(TypeTag, EvaluatePermeabilityAtScvfIP)
                                                || std::is_same<std::decay_t<typename SpatialParams::PermeabilityType>, std::decay_t<decltype(problem.spatialParams().permeabilityAtPos(scvfIpGlobal))>>::value,
                                                "permeabilityAtPos doesn't return PermeabilityType stated in the spatial params!");
                                        return problem.spatialParams().permeabilityAtPos(scvfIpGlobal);
                                    }
                                    else
                                        return volVars.permeability();
                               };

        const Scalar ti = computeTpfaTransmissibility(scvf, insideScv, getPermeability(insideVolVars, scvf.ipGlobal()),
                                                      insideVolVars.extrusionFactor());

        // on the boundary (dirichlet) we only need ti
        if (scvf.boundary())
            tij = scvf.area()*ti;

        // otherwise we compute a tpfa harmonic mean
        else
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            // as we assemble fluxes from the neighbor to our element the outside index
            // refers to the scv of our element, so we use the scv method
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const Scalar tj = -1.0*computeTpfaTransmissibility(scvf, outsideScv, getPermeability(outsideVolVars, scvf.ipGlobal()),
                                                               outsideVolVars.extrusionFactor());

            // harmonic mean (check for division by zero!)
            // TODO: This could lead to problems!? Is there a better way to do this?
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = scvf.area()*(ti * tj)/(ti + tj);
        }

        return tij;
    }
};

/*!
 * \ingroup CCTpfaDiscretization
 * \brief Specialization of the CCTpfaDarcysLaw grids where dim < dimWorld (network/surface grids)
 */
template<class TypeTag>
class CCTpfaDarcysLaw<TypeTag, /*isNetwork*/ true>
{
    using Implementation = DarcysLawImplementation<TypeTag, DiscretizationMethods::CCTpfa>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    //! state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCTpfa;

    //! state the type for the corresponding cache
    using Cache = TpfaDarcysLawCache<TypeTag>;

    //! Compute the advective flux
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        static const bool gravity = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Problem.EnableGravity");

        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        if (gravity)
        {
            // do averaging for the density over all neighboring elements
            const auto rho = [&]()
            {
                // boundaries
                if (scvf.boundary())
                    return outsideVolVars.density(phaseIdx);

                // inner faces with two neighboring elements
                else if (scvf.numOutsideScvs() == 1)
                    return (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;

                // inner faces in networks (general case)
                else
                {
                    Scalar rho(insideVolVars.density(phaseIdx));
                    for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
                    {
                        const auto outsideScvIdx = scvf.outsideScvIdx(i);
                        const auto& outsideVolVars = elemVolVars[outsideScvIdx];
                        rho += outsideVolVars.density(phaseIdx);
                    }
                    return rho/(scvf.numOutsideScvs()+1);
                }
            }();

            const auto& tij = fluxVarsCache.advectionTij();
            const auto& g = problem.gravityAtPos(scvf.ipGlobal());

            // Obtain inside and outside pressures
            const auto pInside = insideVolVars.pressure(phaseIdx);
            const auto pOutside = [&]()
            {
                // Dirichlet boundaries and inner faces with one neighbor
                if (scvf.numOutsideScvs() == 1)
                    return outsideVolVars.pressure(phaseIdx);

                // inner faces in networks (general case)
                else
                {
                    Scalar sumTi(tij);
                    Scalar sumPTi(tij*pInside);

                    // add inside gravitational contribution
                    sumPTi += rho*scvf.area()
                              *insideVolVars.extrusionFactor()
                              *vtmv(scvf.unitOuterNormal(), insideVolVars.permeability(), g);

                    for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
                    {
                        const auto outsideScvIdx = scvf.outsideScvIdx(i);
                        const auto& flippedScvf = fvGeometry.flipScvf(scvf.index(), i);
                        const auto& outsideVolVars = elemVolVars[outsideScvIdx];
                        const auto& outsideFluxVarsCache = elemFluxVarsCache[flippedScvf];
                        sumTi += outsideFluxVarsCache.advectionTij();
                        sumPTi += outsideFluxVarsCache.advectionTij()*outsideVolVars.pressure(phaseIdx);

                        // add outside gravitational contribution
                        sumPTi += rho*scvf.area()
                                  *outsideVolVars.extrusionFactor()
                                  *vtmv(flippedScvf.unitOuterNormal(), outsideVolVars.permeability(), g);
                    }
                    return sumPTi/sumTi;
                }
            }();

            //! precompute alpha := n^T*K*g
            const auto alpha_inside = vtmv(scvf.unitOuterNormal(), insideVolVars.permeability(), g)*insideVolVars.extrusionFactor();

            Scalar flux = tij*(pInside - pOutside) + scvf.area()*rho*alpha_inside;

            //! On interior faces with one neighbor we have to add K-weighted gravitational contributions
            if (!scvf.boundary() && scvf.numOutsideScvs() == 1)
            {
                const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                const auto& outsideScvf = fvGeometry.flipScvf(scvf.index());
                const auto outsideK = outsideVolVars.permeability();
                const auto outsideTi = computeTpfaTransmissibility(outsideScvf, outsideScv, outsideK, outsideVolVars.extrusionFactor());
                const auto alpha_outside = vtmv(outsideScvf.unitOuterNormal(), outsideK, g)*outsideVolVars.extrusionFactor();

                flux -= rho*tij/outsideTi*(alpha_inside + alpha_outside);
            }

            return flux;
        }
        else
        {
            // Obtain inside and outside pressures
            const auto pInside = insideVolVars.pressure(phaseIdx);
            const auto pOutside = [&]()
            {
                // Dirichlet boundaries and inner faces with two neighboring elements
                if (scvf.numOutsideScvs() <= 1)
                    return outsideVolVars.pressure(phaseIdx);

                // inner faces in networks (general case)
                else
                {
                    const auto& insideFluxVarsCache = elemFluxVarsCache[scvf];
                    Scalar sumTi(insideFluxVarsCache.advectionTij());
                    Scalar sumPTi(insideFluxVarsCache.advectionTij()*pInside);

                    for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
                    {
                        const auto outsideScvIdx = scvf.outsideScvIdx(i);
                        const auto& flippedScvf = fvGeometry.flipScvf(scvf.index(), i);
                        const auto& outsideVolVars = elemVolVars[outsideScvIdx];
                        const auto& outsideFluxVarsCache = elemFluxVarsCache[flippedScvf];
                        sumTi += outsideFluxVarsCache.advectionTij();
                        sumPTi += outsideFluxVarsCache.advectionTij()*outsideVolVars.pressure(phaseIdx);
                    }
                    return sumPTi/sumTi;
                }
            }();

            // return flux
            return fluxVarsCache.advectionTij()*(pInside - pOutside);
        }
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibility will be computed and stored using the method below.
    static Scalar calculateTransmissibility(const Problem& problem,
                                            const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const SubControlVolumeFace& scvf)
    {
        Scalar tij;

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        // check if we evaluate the permeability in the volume (for discontinuous fields, default)
        // or at the scvf center for analytical permeability fields (e.g. convergence studies)
        auto getPermeability = [&problem](const VolumeVariables& volVars,
                                          const GlobalPosition& scvfIpGlobal) -> typename SpatialParams::PermeabilityType
                               {
                                    if (GET_PROP_VALUE(TypeTag, EvaluatePermeabilityAtScvfIP))
                                    {
                                        // TODO: Make PermeabilityType a property!!
                                        // We do an implicit cast to permeability type in case EvaluatePermeabilityAtScvfIP is not true so that it compiles
                                        // If it is true make sure that no cast is necessary and the correct return type is specified in the spatial params
                                        static_assert(!GET_PROP_VALUE(TypeTag, EvaluatePermeabilityAtScvfIP)
                                                || std::is_same<std::decay_t<typename SpatialParams::PermeabilityType>, std::decay_t<decltype(problem.spatialParams().permeabilityAtPos(scvfIpGlobal))>>::value,
                                                "permeabilityAtPos doesn't return PermeabilityType stated in the spatial params!");
                                        return problem.spatialParams().permeabilityAtPos(scvfIpGlobal);
                                    }
                                    else
                                        return volVars.permeability();
                               };

        const Scalar ti = computeTpfaTransmissibility(scvf, insideScv, getPermeability(insideVolVars, scvf.ipGlobal()),
                                                      insideVolVars.extrusionFactor());

        // for the boundary (dirichlet) or at branching points we only need ti
        if (scvf.boundary() || scvf.numOutsideScvs() > 1)
            tij = scvf.area()*ti;

        // otherwise we compute a tpfa harmonic mean
        else
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            // as we assemble fluxes from the neighbor to our element the outside index
            // refers to the scv of our element, so we use the scv method
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const Scalar tj = computeTpfaTransmissibility(fvGeometry.flipScvf(scvf.index()), outsideScv, getPermeability(outsideVolVars, scvf.ipGlobal()),
                                                          outsideVolVars.extrusionFactor());

            // harmonic mean (check for division by zero!)
            // TODO: This could lead to problems!? Is there a better way to do this?
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = scvf.area()*(ti * tj)/(ti + tj);
        }

        return tij;
    }
};

} // end namespace Dumux

#endif
