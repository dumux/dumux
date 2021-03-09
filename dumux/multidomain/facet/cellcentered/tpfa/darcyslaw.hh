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
 * \ingroup FacetCoupling
 * \copydoc Dumux::CCTpfaFacetCouplingDarcysLaw
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FACET_COUPLING_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FACET_COUPLING_DARCYS_LAW_HH

#include <array>
#include <cmath>

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/flux/cctpfa/darcyslaw.hh>

namespace Dumux {

//! Forward declaration of the implementation
template<class ScalarType, class GridGeometry, bool isNetwork>
class CCTpfaFacetCouplingDarcysLawImpl;

/*!
 * \ingroup FacetCoupling
 * \brief The cache corresponding to tpfa Darcy's Law with facet coupling
 * \note We distinguish between network and non-network grids here. Specializations
 *       for the two cases can be found below.
 */
template<class AdvectionType, class GridGeometry, bool isNetwork>
class CCTpfaFacetCouplingDarcysLawCache;

/*!
 * \ingroup FacetCoupling
 * \brief Darcy's law for cell-centered finite volume schemes with two-point flux approximation
 *        in the context of coupled models where the coupling occurs across the facets of the bulk
 *        domain elements with a lower-dimensional domain living on these facets.
 */
template<class ScalarType, class GridGeometry>
using CCTpfaFacetCouplingDarcysLaw =
      CCTpfaFacetCouplingDarcysLawImpl< ScalarType, GridGeometry, ( int(GridGeometry::GridView::dimension) <
                                                                      int(GridGeometry::GridView::dimensionworld) ) >;

/*!
 * \ingroup FacetCoupling
 * \brief Specialization of FacetCouplingTpfaDarcysLawCache for non-network grids.
 */
template<class AdvectionType, class GridGeometry>
class CCTpfaFacetCouplingDarcysLawCache<AdvectionType, GridGeometry, /*isNetwork*/false>
{
    using Scalar = typename AdvectionType::Scalar;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

public:
    //! export the corresponding filler class
    using Filler = TpfaDarcysLawCacheFiller<GridGeometry>;

    //! we store the transmissibilities associated with the interior
    //! cell, outside cell, and the fracture facet in an array. Access
    //! to this array should be done using the following indices:
    static constexpr int insideTijIdx = 0;
    static constexpr int outsideTijIdx = 1;
    static constexpr int facetTijIdx = 2;

    //! Export transmissibility storage type
    using AdvectionTransmissibilityContainer = std::array<Scalar, 3>;

    //! Export the type used for the gravity coefficients
    using GravityCoefficients = std::array<Scalar, 2>;

    //! update subject to a given problem
    template< class Problem, class ElementVolumeVariables >
    void updateAdvection(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace &scvf)
    {
        tij_ = AdvectionType::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf, g_);
    }

    //! We use the same name as in the TpfaDarcysLawCache so
    //! that this cache and the law implementation for non-coupled
    //! models can be reused here on facets that do not lie on an
    //! interior boundary, i.e. do not coincide with a facet element
    Scalar advectionTij() const
    { return tij_[insideTijIdx]; }

    //! returns the transmissibility associated with the inside cell
    Scalar advectionTijInside() const
    { return tij_[insideTijIdx]; }

    //! returns the transmissibility associated with the outside cell
    Scalar advectionTijOutside() const
    { return tij_[outsideTijIdx]; }

    //! returns the transmissibility associated with the outside cell
    Scalar advectionTijFacet() const
    { return tij_[facetTijIdx]; }

    //! return the coefficients for the computation of gravity at the scvf
    const GravityCoefficients& gravityCoefficients() const
    { return g_; }

private:
    std::array<Scalar, 3> tij_;
    GravityCoefficients g_;
};

/*!
 * \ingroup FacetCoupling
 * \brief Specialization of CCTpfaFacetCouplingDarcysLawImpl for dim=dimWorld
 */
template<class ScalarType, class GridGeometry>
class CCTpfaFacetCouplingDarcysLawImpl<ScalarType, GridGeometry, /*isNetwork*/false>
{
    using ThisType = CCTpfaFacetCouplingDarcysLawImpl<ScalarType, GridGeometry, false>;
    using TpfaDarcysLaw = CCTpfaDarcysLaw<ScalarType, GridGeometry, false>;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

  public:
    //! state the scalar type of the law
    using Scalar = ScalarType;
    //! export the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;
    //! export the type for the corresponding cache
    using Cache = CCTpfaFacetCouplingDarcysLawCache<ThisType, GridGeometry, /*isNetwork*/false>;
    //! export the type used to store transmissibilities
    using TijContainer = typename Cache::AdvectionTransmissibilityContainer;


    //! Compute the advective flux
    template< class Problem, class ElementVolumeVariables, class ElementFluxVarsCache >
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        if (!problem.couplingManager().isOnInteriorBoundary(element, scvf))
            return TpfaDarcysLaw::flux(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarsCache);

        // Obtain inside and fracture pressures
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);
        const auto pInside = insideVolVars.pressure(phaseIdx);
        const auto pFacet = facetVolVars.pressure(phaseIdx);

        // compute and return flux
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        Scalar flux = fluxVarsCache.advectionTijInside()*pInside + fluxVarsCache.advectionTijFacet()*pFacet;

        // maybe add gravitational acceleration
        static const Scalar gravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
        if (gravity)
        {
            // compute alpha := n^T*K*g and add to flux (use arithmetic mean for density)
            const auto& g = problem.spatialParams().gravity(scvf.ipGlobal());
            const auto rho = 0.5*(insideVolVars.density(phaseIdx) + facetVolVars.density(phaseIdx));
            const auto rhoTimesArea = rho*Extrusion::area(scvf);
            const auto alpha_inside = rhoTimesArea*insideVolVars.extrusionFactor()
                                      *vtmv(scvf.unitOuterNormal(), insideVolVars.permeability(), g);

            flux += alpha_inside;
            if (!scvf.boundary())
            {
                const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

                // add further gravitational contributions
                if ( problem.interiorBoundaryTypes(element, scvf).hasOnlyNeumann() )
                {
                    static const Scalar xi = getParamFromGroup<Scalar>(problem.paramGroup(), "FacetCoupling.Xi", 1.0);
                    const auto alpha_facet = rhoTimesArea*insideVolVars.extrusionFactor()
                                             *vtmv(scvf.unitOuterNormal(), facetVolVars.permeability(), g);
                    const auto alpha_outside = rhoTimesArea*outsideVolVars.extrusionFactor()
                                               *vtmv(scvf.unitOuterNormal(), outsideVolVars.permeability(), g);

                    flux -= fluxVarsCache.gravityCoefficients()[0]*(xi*alpha_inside - alpha_facet + (1.0 - xi)*alpha_outside);
                    flux += fluxVarsCache.gravityCoefficients()[1]*(xi*alpha_outside - alpha_facet + (1.0 - xi)*alpha_inside);
                }

                // add outside contribution
                flux += fluxVarsCache.advectionTijOutside()*outsideVolVars.pressure(phaseIdx);
            }

            return flux;
        }
        else
            return scvf.boundary() ? flux
                                   : flux + fluxVarsCache.advectionTijOutside()*elemVolVars[scvf.outsideScvIdx()].pressure(phaseIdx);
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibility will be computed and stored using the method below.
    template< class Problem, class ElementVolumeVariables >
    static TijContainer calculateTransmissibility(const Problem& problem,
                                                  const Element& element,
                                                  const FVElementGeometry& fvGeometry,
                                                  const ElementVolumeVariables& elemVolVars,
                                                  const SubControlVolumeFace& scvf)
    {
        typename Cache::GravityCoefficients g;
        return calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf, g);
    }

    // This overload additionally receives a container in which the coefficients required
    // for the computation of the gravitational acceleration ar the scvf are stored
    template< class Problem, class ElementVolumeVariables >
    static TijContainer calculateTransmissibility(const Problem& problem,
                                                  const Element& element,
                                                  const FVElementGeometry& fvGeometry,
                                                  const ElementVolumeVariables& elemVolVars,
                                                  const SubControlVolumeFace& scvf,
                                                  typename Cache::GravityCoefficients& g)
    {
        TijContainer tij;
        if (!problem.couplingManager().isCoupled(element, scvf))
        {
            //! use the standard darcy's law and only compute one transmissibility
            tij[Cache::insideTijIdx] = TpfaDarcysLaw::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf);
            return tij;
        }

        //! xi factor for the coupling conditions
        static const Scalar xi = getParamFromGroup<Scalar>(problem.paramGroup(), "FacetCoupling.Xi", 1.0);

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const auto wIn = Extrusion::area(scvf)
                         *computeTpfaTransmissibility(scvf, insideScv,
                                                      insideVolVars.permeability(),
                                                      insideVolVars.extrusionFactor());

        // proceed depending on the interior BC types used
        const auto iBcTypes = problem.interiorBoundaryTypes(element, scvf);

        // neumann-coupling
        if (iBcTypes.hasOnlyNeumann())
        {
            const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);
            const auto wFacet = 2.0*Extrusion::area(scvf)*insideVolVars.extrusionFactor()
                                   /facetVolVars.extrusionFactor()
                                   *vtmv(scvf.unitOuterNormal(), facetVolVars.permeability(), scvf.unitOuterNormal());

            // The fluxes across this face and the outside face can be expressed in matrix form:
            // \f$\mathbf{C} \bar{\mathbf{u}} + \mathbf{D} \mathbf{u} + \mathbf{E} \mathbf{u}_\gamma\f$,
            // where \f$\gamma$\f denotes the domain living on the facets and \f$\bar{\mathbf{u}}$\f are
            // intermediate face unknowns in the matrix domain. Equivalently, flux continuity reads:
            // \f$\mathbf{A} \bar{\mathbf{u}} = \mathbf{B} \mathbf{u} + \mathbf{M} \mathbf{u}_\gamma\f$.
            // Combining the two, we can eliminate the intermediate unknowns and compute the transmissibilities
            // that allow the description of the fluxes as functions of the cell and Dirichlet pressures only.
            if (!scvf.boundary())
            {
                const auto outsideScvIdx = scvf.outsideScvIdx();
                const auto& outsideVolVars = elemVolVars[outsideScvIdx];
                const auto wOut = -1.0*Extrusion::area(scvf)
                                  *computeTpfaTransmissibility(scvf, fvGeometry.scv(outsideScvIdx),
                                                               outsideVolVars.permeability(),
                                                               outsideVolVars.extrusionFactor());

                if ( !Dune::FloatCmp::eq(xi, 1.0, 1e-6) )
                {
                    // The gravity coefficients are the first row of the inverse of the A matrix in the local eq system
                    // multiplied with wIn. Note that we never compute the inverse but use an optimized implementation below.
                    // The A matrix has the following coefficients:
                    // A = | xi*wIn + wFacet, (xi - 1.0)*wOut  |  -> AInv = 1/detA | xi*wOut + wFacet, -(xi - 1.0)*wOut |
                    //     | wIn*(xi - 1.0) , xi*wOut + wFacet |                   | -wIn*(xi - 1.0) , xi*wIn + wFacet  |
                    const Scalar xiMinusOne = (xi - 1.0);
                    const Scalar a01 = xiMinusOne*wOut;
                    const Scalar a11 = xi*wOut + wFacet;
                    const Scalar detA = (xi*wIn + wFacet)*a11 - xiMinusOne*wIn*a01;
                    g[0] = wIn*a11/detA; g[1] = -wIn*a01/detA;

                    // optimized implementation: factorization obtained using sympy
                    const Scalar factor = wIn * wFacet / ( wIn * wOut * ( 2.0 * xi - 1.0 ) + wFacet * ( xi * ( wIn + wOut ) + wFacet ) );
                    tij[Cache::insideTijIdx]  = factor * ( wOut * xi + wFacet );
                    tij[Cache::outsideTijIdx] = factor * ( wOut * ( 1.0 - xi ) );
                    tij[Cache::facetTijIdx]   = factor * ( - wOut - wFacet );
                }
                else
                {
                    g[0] = wIn/(wIn+wFacet); g[1] = 0.0;
                    tij[Cache::insideTijIdx] = wFacet*g[0];
                    tij[Cache::facetTijIdx] = -tij[Cache::insideTijIdx];
                    tij[Cache::outsideTijIdx] = 0.0;
                }
            }
            else
            {
                // TODO: check for division by zero??
                tij[Cache::insideTijIdx] = wFacet*wIn/(wIn+wFacet);
                tij[Cache::facetTijIdx] = -tij[Cache::insideTijIdx];
                tij[Cache::outsideTijIdx] = 0.0;
            }
        }
        else if (iBcTypes.hasOnlyDirichlet())
        {
            tij[Cache::insideTijIdx] = wIn;
            tij[Cache::outsideTijIdx] = 0.0;
            tij[Cache::facetTijIdx] = -wIn;
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Interior boundary types other than pure Dirichlet or Neumann");

        return tij;
    }
};

/*!
 * \ingroup FacetCoupling
 * \brief Specialization of FacetCouplingTpfaDarcysLawCache for network grids
 */
template<class AdvectionType, class GridGeometry>
class CCTpfaFacetCouplingDarcysLawCache<AdvectionType, GridGeometry, /*isNetwork*/true>
{
    using Scalar = typename AdvectionType::Scalar;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

public:
    //! export the corresponding filler class
    using Filler = TpfaDarcysLawCacheFiller<GridGeometry>;

    //! we store the transmissibilities associated with the interior
    //! cell and the fracture facet in an array. Access to this array
    //! should be done using the following indices:
    static constexpr int insideTijIdx = 0;
    static constexpr int facetTijIdx = 1;

    //! Export transmissibility storage type
    using AdvectionTransmissibilityContainer = std::array<Scalar, 2>;

    //! Export the type used for the gravity coefficients
    using GravityCoefficients = std::array<Scalar, 1>;

    //! update subject to a given problem
    template< class Problem, class ElementVolumeVariables >
    void updateAdvection(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace &scvf)
    {
        tij_ = AdvectionType::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf, g_);
    }

    //! We use the same name as in the TpfaDarcysLawCache so
    //! that this cache and the law implementation for non-coupled
    //! models can be reused here on facets that do not lie on an
    //! interior boundary, i.e. do not coincide with a facet element
    Scalar advectionTij() const
    { return tij_[insideTijIdx]; }

    //! returns the transmissibility associated with the inside cell
    Scalar advectionTijInside() const
    { return tij_[insideTijIdx]; }

    //! returns the transmissibility associated with the outside cell
    Scalar advectionTijFacet() const
    { return tij_[facetTijIdx]; }

    //! return the coefficients for the computation of gravity at the scvf
    const GravityCoefficients& gravityCoefficients() const
    { return g_; }

private:
    AdvectionTransmissibilityContainer tij_;
    GravityCoefficients g_;
};

/*!
 * \ingroup FacetCoupling
 * \brief Specialization of CCTpfaFacetCouplingDarcysLawImpl for dim<dimWorld
 */
template<class ScalarType, class GridGeometry>
class CCTpfaFacetCouplingDarcysLawImpl<ScalarType, GridGeometry, /*isNetwork*/true>
{
    using ThisType = CCTpfaFacetCouplingDarcysLawImpl<ScalarType, GridGeometry, true>;
    using TpfaDarcysLaw = CCTpfaDarcysLaw<ScalarType, GridGeometry, true>;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

  public:
    //! state the scalar type of the law
    using Scalar = ScalarType;
    //! state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;
    //! state the type for the corresponding cache
    using Cache = CCTpfaFacetCouplingDarcysLawCache<ThisType, GridGeometry, /*isNetwork*/true>;
    //! export the type used to store transmissibilities
    using TijContainer = typename Cache::AdvectionTransmissibilityContainer;

    //! Compute the advective flux
    template< class Problem, class ElementVolumeVariables, class ElementFluxVarsCache >
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        if (!problem.couplingManager().isOnInteriorBoundary(element, scvf))
            return TpfaDarcysLaw::flux(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarsCache);

        // On surface grids only xi = 1.0 can be used, as the coupling condition
        // for xi != 1.0 does not generalize for surface grids where there can be
        // seveal neighbor meeting at a branching point.
        static const Scalar xi = getParamFromGroup<Scalar>(problem.paramGroup(), "FacetCoupling.Xi", 1.0);
        if (Dune::FloatCmp::ne(xi, 1.0, 1e-6))
            DUNE_THROW(Dune::InvalidStateException, "Xi != 1.0 cannot be used on surface grids");

        // Obtain inside and fracture pressures
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);
        const auto pInside = insideVolVars.pressure(phaseIdx);
        const auto pFacet = facetVolVars.pressure(phaseIdx);

        // compute and return flux
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        Scalar flux = fluxVarsCache.advectionTijInside()*pInside + fluxVarsCache.advectionTijFacet()*pFacet;

        static const Scalar gravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
        if (gravity)
        {
            // compute alpha := n^T*K*g and add to flux (use arithmetic mean for density)
            const auto& g = problem.spatialParams().gravity(scvf.ipGlobal());
            const auto rho = 0.5*(insideVolVars.density(phaseIdx) + facetVolVars.density(phaseIdx));
            const auto rhoTimesArea = rho*Extrusion::area(scvf);
            const auto alpha_inside = rhoTimesArea*insideVolVars.extrusionFactor()
                                      *vtmv(scvf.unitOuterNormal(), insideVolVars.permeability(), g);

            flux += alpha_inside;

            // maybe add further gravitational contributions
            if ( !scvf.boundary() && problem.interiorBoundaryTypes(element, scvf).hasOnlyNeumann() )
            {
                const auto alpha_facet = rhoTimesArea*insideVolVars.extrusionFactor()
                                         *vtmv(scvf.unitOuterNormal(), facetVolVars.permeability(), g);

                flux -= fluxVarsCache.gravityCoefficients()[0]*(alpha_inside - alpha_facet);
            }
        }

        return flux;
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibility will be computed and stored using the method below.
    template< class Problem, class ElementVolumeVariables >
    static TijContainer calculateTransmissibility(const Problem& problem,
                                                  const Element& element,
                                                  const FVElementGeometry& fvGeometry,
                                                  const ElementVolumeVariables& elemVolVars,
                                                  const SubControlVolumeFace& scvf)
    {
        typename Cache::GravityCoefficients g;
        return calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf, g);
    }

    // This overload additionally receives a container in which the coefficients required
    // for the computation of the gravitational acceleration ar the scvf are stored
    template< class Problem, class ElementVolumeVariables >
    static TijContainer calculateTransmissibility(const Problem& problem,
                                                  const Element& element,
                                                  const FVElementGeometry& fvGeometry,
                                                  const ElementVolumeVariables& elemVolVars,
                                                  const SubControlVolumeFace& scvf,
                                                  typename Cache::GravityCoefficients& g)
    {
        TijContainer tij;
        if (!problem.couplingManager().isCoupled(element, scvf))
        {
            //! use the standard darcy's law and only compute one transmissibility
            tij[Cache::insideTijIdx] = TpfaDarcysLaw::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf);
            return tij;
        }

        //! xi factor for the coupling conditions
        static const Scalar xi = getParamFromGroup<Scalar>(problem.paramGroup(), "FacetCoupling.Xi", 1.0);

        // On surface grids only xi = 1.0 can be used, as the coupling condition
        // for xi != 1.0 does not generalize for surface grids where the normal
        // vectors of the inside/outside elements have different orientations.
        if (Dune::FloatCmp::ne(xi, 1.0, 1e-6))
            DUNE_THROW(Dune::InvalidStateException, "Xi != 1.0 cannot be used on surface grids");

        const auto area = Extrusion::area(scvf);
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const auto wIn = area*computeTpfaTransmissibility(scvf, insideScv, insideVolVars.permeability(), insideVolVars.extrusionFactor());

        // proceed depending on the interior BC types used
        const auto iBcTypes = problem.interiorBoundaryTypes(element, scvf);

        // neumann-coupling
        if (iBcTypes.hasOnlyNeumann())
        {
            const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);

            // Here we use the square root of the facet extrusion factor
            // as an approximate average distance from scvf ip to facet center
            using std::sqrt;
            const auto wFacet = 2.0*area*insideVolVars.extrusionFactor()
                                        /sqrt(facetVolVars.extrusionFactor())
                                        *vtmv(scvf.unitOuterNormal(), facetVolVars.permeability(), scvf.unitOuterNormal());

            // TODO: check for division by zero??
            g[0] = wIn/(wIn+wFacet);
            tij[Cache::insideTijIdx] = wFacet*g[0];
            tij[Cache::facetTijIdx] = -tij[Cache::insideTijIdx];
        }
        else if (iBcTypes.hasOnlyDirichlet())
        {
            tij[Cache::insideTijIdx] = wIn;
            tij[Cache::facetTijIdx] = -wIn;
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Interior boundary types other than pure Dirichlet or Neumann");

        return tij;
    }
};

} // end namespace Dumux

#endif
