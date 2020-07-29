// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup CCTpfaFlux
 * \brief Darcy's law for cell-centered finite volume schemes with two-point flux approximation
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_DARCYS_LAW_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

namespace Dumux {

// forward declarations
template<class TypeTag, DiscretizationMethod discMethod>
class DarcysLawImplementation;

/*!
 * \ingroup CCTpfaFlux
 * \brief Darcy's law for cell-centered finite volume schemes with two-point flux approximation
 * \note Darcy's law is specialized for network and surface grids (i.e. if grid dim < dimWorld)
 * \tparam Scalar the scalar type for scalar physical quantities
 * \tparam GridGeometry the grid geometry
 * \tparam isNetwork whether we are computing on a network grid embedded in a higher world dimension
 */
template<class Scalar, class GridGeometry, bool isNetwork>
class CCTpfaDarcysLaw;

/*!
 * \ingroup CCTpfaFlux
 * \brief Darcy's law for cell-centered finite volume schemes with two-point flux approximation
 * \note Darcy's law is specialized for network and surface grids (i.e. if grid dim < dimWorld)
 */
template <class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethod::cctpfa>
: public CCTpfaDarcysLaw<GetPropType<TypeTag, Properties::Scalar>,
                         GetPropType<TypeTag, Properties::GridGeometry>,
                         (GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension < GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimensionworld)>
{};

/*!
 * \ingroup CCTpfaFlux
 * \brief Class that fills the cache corresponding to tpfa Darcy's Law
 */
template<class GridGeometry>
class TpfaDarcysLawCacheFiller
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

public:
    //! Function to fill a TpfaDarcysLawCache of a given scvf
    //! This interface has to be met by any advection-related cache filler class
    //! TODO: Probably get cache type out of the filler
    template<class FluxVariablesCache, class Problem, class ElementVolumeVariables, class FluxVariablesCacheFiller>
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
 * \ingroup CCTpfaFlux
 * \brief The cache corresponding to tpfa Darcy's Law
 */
template<class AdvectionType, class GridGeometry>
class TpfaDarcysLawCache
{
    using Scalar = typename AdvectionType::Scalar;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

public:
    using Filler = TpfaDarcysLawCacheFiller<GridGeometry>;

    template<class Problem, class ElementVolumeVariables>
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
 * \ingroup CCTpfaFlux
 * \brief Specialization of the CCTpfaDarcysLaw grids where dim=dimWorld
 */
template<class ScalarType, class GridGeometry>
class CCTpfaDarcysLaw<ScalarType, GridGeometry, /*isNetwork*/ false>
{
    using ThisType = CCTpfaDarcysLaw<ScalarType, GridGeometry, /*isNetwork*/ false>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

  public:
    //! state the scalar type of the law
    using Scalar = ScalarType;

    //! state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;

    //! state the type for the corresponding cache
    using Cache = TpfaDarcysLawCache<ThisType, GridGeometry>;

    /*!
     * \brief Returns the advective flux of a fluid phase
     *        across the given sub-control volume face.
     * \note This assembles the term
     *       \f$-|\sigma| \mathbf{n}^T \mathbf{K} \left( \nabla p - \rho \mathbf{g} \right)\f$,
     *       where \f$|\sigma|\f$ is the area of the face and \f$\mathbf{n}\f$ is the outer
     *       normal vector. Thus, the flux is given in N*m, and can be converted
     *       into a volume flux (m^3/s) or mass flux (kg/s) by applying an upwind scheme
     *       for the mobility or the product of density and mobility, respectively.
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
        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");

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
            const auto& g = problem.spatialParams().gravity(scvf.ipGlobal());

            //! compute alpha := n^T*K*g
            const auto alpha_inside = vtmv(scvf.unitOuterNormal(), insideVolVars.permeability(), g)*insideVolVars.extrusionFactor();

            Scalar flux = tij*(pInside - pOutside) + rho*Extrusion::area(scvf)*alpha_inside;

            //! On interior faces we have to add K-weighted gravitational contributions
            if (!scvf.boundary())
            {
                const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                const auto outsideK = outsideVolVars.permeability();
                const auto outsideTi = fvGeometry.gridGeometry().isPeriodic()
                    ? computeTpfaTransmissibility(fvGeometry.flipScvf(scvf.index()), outsideScv, outsideK, outsideVolVars.extrusionFactor())
                    : -1.0*computeTpfaTransmissibility(scvf, outsideScv, outsideK, outsideVolVars.extrusionFactor());
                const auto alpha_outside = vtmv(scvf.unitOuterNormal(), outsideK, g)*outsideVolVars.extrusionFactor();

                flux -= rho*tij/outsideTi*(alpha_inside - alpha_outside);
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
    template<class Problem, class ElementVolumeVariables>
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

        const Scalar ti = computeTpfaTransmissibility(scvf, insideScv,
                                                      getPermeability_(problem, insideVolVars, scvf.ipGlobal()),
                                                      insideVolVars.extrusionFactor());

        // on the boundary (dirichlet) we only need ti
        if (scvf.boundary())
            tij = Extrusion::area(scvf)*ti;

        // otherwise we compute a tpfa harmonic mean
        else
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            // as we assemble fluxes from the neighbor to our element the outside index
            // refers to the scv of our element, so we use the scv method
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const Scalar tj = fvGeometry.gridGeometry().isPeriodic()
                ? computeTpfaTransmissibility(fvGeometry.flipScvf(scvf.index()), outsideScv, getPermeability_(problem, outsideVolVars, scvf.ipGlobal()), outsideVolVars.extrusionFactor())
                : -1.0*computeTpfaTransmissibility(scvf, outsideScv, getPermeability_(problem, outsideVolVars, scvf.ipGlobal()), outsideVolVars.extrusionFactor());

            // harmonic mean (check for division by zero!)
            // TODO: This could lead to problems!? Is there a better way to do this?
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = Extrusion::area(scvf)*(ti * tj)/(ti + tj);
        }

        return tij;
    }

private:
    template<class Problem, class VolumeVariables,
             std::enable_if_t<!Problem::SpatialParams::evaluatePermeabilityAtScvfIP(), int> = 0>
    static decltype(auto) getPermeability_(const Problem& problem,
                                           const VolumeVariables& volVars,
                                           const GlobalPosition& scvfIpGlobal)
    { return volVars.permeability(); }

    template<class Problem, class VolumeVariables,
             std::enable_if_t<Problem::SpatialParams::evaluatePermeabilityAtScvfIP(), int> = 0>
    static decltype(auto) getPermeability_(const Problem& problem,
                                           const VolumeVariables& volVars,
                                           const GlobalPosition& scvfIpGlobal)
    { return problem.spatialParams().permeabilityAtPos(scvfIpGlobal); }
};

/*!
 * \ingroup CCTpfaFlux
 * \brief Specialization of the CCTpfaDarcysLaw grids where dim < dimWorld (network/surface grids)
 */
template<class ScalarType, class GridGeometry>
class CCTpfaDarcysLaw<ScalarType, GridGeometry, /*isNetwork*/ true>
{
    using ThisType = CCTpfaDarcysLaw<ScalarType, GridGeometry, /*isNetwork*/ true>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! state the scalar type of the law
    using Scalar = ScalarType;

    //! state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;

    //! state the type for the corresponding cache
    using Cache = TpfaDarcysLawCache<ThisType, GridGeometry>;

    /*!
     * \brief Returns the advective flux of a fluid phase
     *        across the given sub-control volume face.
     * \note This assembles the term
     *       \f$-|\sigma| \mathbf{n}^T \mathbf{K} \left( \nabla p - \rho \mathbf{g} \right)\f$,
     *       where \f$|\sigma|\f$ is the area of the face and \f$\mathbf{n}\f$ is the outer
     *       normal vector. Thus, the flux is given in N*m, and can be converted
     *       into a volume flux (m^3/s) or mass flux (kg/s) by applying an upwind scheme
     *       for the mobility or the product of density and mobility, respectively.
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
        static const bool gravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");

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
            const auto& g = problem.spatialParams().gravity(scvf.ipGlobal());

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
                    sumPTi += rho*Extrusion::area(scvf)
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
                        sumPTi += rho*Extrusion::area(scvf)
                                  *outsideVolVars.extrusionFactor()
                                  *vtmv(flippedScvf.unitOuterNormal(), outsideVolVars.permeability(), g);
                    }
                    return sumPTi/sumTi;
                }
            }();

            //! precompute alpha := n^T*K*g
            const auto alpha_inside = vtmv(scvf.unitOuterNormal(), insideVolVars.permeability(), g)*insideVolVars.extrusionFactor();

            Scalar flux = tij*(pInside - pOutside) + Extrusion::area(scvf)*rho*alpha_inside;

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
    template<class Problem, class ElementVolumeVariables>
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

        const Scalar ti = computeTpfaTransmissibility(scvf, insideScv,
                                                      getPermeability_(problem, insideVolVars, scvf.ipGlobal()),
                                                      insideVolVars.extrusionFactor());

        // for the boundary (dirichlet) or at branching points we only need ti
        if (scvf.boundary() || scvf.numOutsideScvs() > 1)
            tij = Extrusion::area(scvf)*ti;

        // otherwise we compute a tpfa harmonic mean
        else
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            // as we assemble fluxes from the neighbor to our element the outside index
            // refers to the scv of our element, so we use the scv method
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const Scalar tj = computeTpfaTransmissibility(fvGeometry.flipScvf(scvf.index()), outsideScv,
                                                          getPermeability_(problem, outsideVolVars, scvf.ipGlobal()),
                                                          outsideVolVars.extrusionFactor());

            // harmonic mean (check for division by zero!)
            // TODO: This could lead to problems!? Is there a better way to do this?
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = Extrusion::area(scvf)*(ti * tj)/(ti + tj);
        }

        return tij;
    }

private:
    template<class Problem, class VolumeVariables,
             std::enable_if_t<!Problem::SpatialParams::evaluatePermeabilityAtScvfIP(), int> = 0>
    static decltype(auto) getPermeability_(const Problem& problem,
                                           const VolumeVariables& volVars,
                                           const GlobalPosition& scvfIpGlobal)
    { return volVars.permeability(); }

    template<class Problem, class VolumeVariables,
             std::enable_if_t<Problem::SpatialParams::evaluatePermeabilityAtScvfIP(), int> = 0>
    static decltype(auto) getPermeability_(const Problem& problem,
                                           const VolumeVariables& volVars,
                                           const GlobalPosition& scvfIpGlobal)
    { return problem.spatialParams().permeabilityAtPos(scvfIpGlobal); }
};

} // end namespace Dumux

#endif
