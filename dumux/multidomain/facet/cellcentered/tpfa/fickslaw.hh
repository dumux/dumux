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
 * \copydoc Dumux::CCTpfaFacetCouplingFicksLaw
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FACET_COUPLING_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FACET_COUPLING_FICKS_LAW_HH

#include <array>
#include <cmath>

#include <dune/common/float_cmp.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/flux/cctpfa/fickslaw.hh>

namespace Dumux {

//! Forward declaration of the implementation
template<class TypeTag,
         ReferenceSystemFormulation referenceSystem,
         bool isNetwork>
class CCTpfaFacetCouplingFicksLawImpl;

/*!
 * \ingroup FacetCoupling
 * \brief Fick's law for cell-centered finite volume schemes with two-point flux approximation
 *        in the context of coupled models where the coupling occurs across the facets of the bulk
 *        domain elements with a lower-dimensional domain defined on these facets.
 */
template<class TypeTag, ReferenceSystemFormulation referenceSystem =  ReferenceSystemFormulation::massAveraged>
using CCTpfaFacetCouplingFicksLaw =
      CCTpfaFacetCouplingFicksLawImpl< TypeTag, referenceSystem,
                                       ( int(GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimension)
                                          < int(GetPropType<TypeTag, Properties::GridGeometry>::GridView::dimensionworld) ) >;

/*!
 * \ingroup FacetCoupling
 * \brief Specialization of CCTpfaFacetCouplingFicksLawImpl for dim=dimWorld
 */
template<class TypeTag, ReferenceSystemFormulation referenceSystem>
class CCTpfaFacetCouplingFicksLawImpl<TypeTag, referenceSystem, /*isNetwork*/false>
: public FicksLawImplementation<TypeTag, DiscretizationMethod::cctpfa, referenceSystem>
{
    using Implementation = CCTpfaFacetCouplingFicksLawImpl<TypeTag, referenceSystem, false>;
    using ParentType = FicksLawImplementation<TypeTag, DiscretizationMethod::cctpfa, referenceSystem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using BalanceEqOpts = GetPropType<TypeTag, Properties::BalanceEqOpts>;

    static const int numPhases = ModelTraits::numFluidPhases();
    static const int numComponents = ModelTraits::numFluidComponents();

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;

    /*!
     * \brief The cache class used with this specialization of Fick's law.
     */
    class FacetCouplingFicksLawCache
    {
    public:
        //! export the corresponding filler class
        using Filler = typename ParentType::Cache::Filler;

        //! we store the transmissibilities associated with the interior
        //! cell, outside cell, and the fracture facet in an array. Access
        //! to this array should be done using the following indices:
        static constexpr int insideTijIdx = 0;
        static constexpr int outsideTijIdx = 1;
        static constexpr int facetTijIdx = 2;

        //! Export transmissibility storage type
        using DiffusionTransmissibilityContainer = std::array<Scalar, 3>;

        //! update subject to a given problem
        template< class Problem, class ElementVolumeVariables >
        void updateDiffusion(const Problem& problem,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace &scvf,
                             unsigned int phaseIdx,
                             unsigned int compIdx)
        {
            tij_[phaseIdx][compIdx] = Implementation::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, compIdx);
        }

        //! We use the same name as in the TpfaFicksLawCache so
        //! that this cache and the law implementation for non-coupled
        //! models can be reused here on facets that do not lie on an
        //! interior boundary, i.e. do not coincide with a facet element
        Scalar diffusionTij(unsigned int phaseIdx, unsigned int compIdx) const
        { return tij_[phaseIdx][compIdx][insideTijIdx]; }

        //! returns the transmissibility associated with the inside cell
        Scalar diffusionTijInside(unsigned int phaseIdx, unsigned int compIdx) const
        { return tij_[phaseIdx][compIdx][insideTijIdx]; }

        //! returns the transmissibility associated with the outside cell
        Scalar diffusionTijOutside(unsigned int phaseIdx, unsigned int compIdx) const
        { return tij_[phaseIdx][compIdx][outsideTijIdx]; }

        //! returns the transmissibility associated with the outside cell
        Scalar diffusionTijFacet(unsigned int phaseIdx, unsigned int compIdx) const
        { return tij_[phaseIdx][compIdx][facetTijIdx]; }

    private:
        std::array< std::array<DiffusionTransmissibilityContainer, numComponents>, numPhases> tij_;
    };

public:
    //! export the type for the corresponding cache
    using Cache = FacetCouplingFicksLawCache;

    //! state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;

    //! Return the reference system
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    //! Compute the diffusive fluxes
    template< class Problem, class ElementVolumeVariables, class ElementFluxVarsCache >
    static ComponentFluxVector flux(const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolumeFace& scvf,
                                    int phaseIdx,
                                    const ElementFluxVarsCache& elemFluxVarsCache)
    {
        if (!problem.couplingManager().isOnInteriorBoundary(element, scvf))
            return ParentType::flux(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarsCache);

        ComponentFluxVector componentFlux(0.0);
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            if(compIdx == FluidSystem::getMainComponent(phaseIdx))
                continue;

            // get inside/outside volume variables
            const auto& fluxVarsCache = elemFluxVarsCache[scvf];
            const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
            const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

            // the inside and outside mass/mole fractions fractions
            const Scalar xInside = massOrMoleFraction(insideVolVars, referenceSystem, phaseIdx, compIdx);
            const Scalar xFacet = massOrMoleFraction(facetVolVars, referenceSystem, phaseIdx, compIdx);
            const Scalar xOutside = massOrMoleFraction(outsideVolVars, referenceSystem, phaseIdx, compIdx);

            const Scalar rhoInside = massOrMolarDensity(insideVolVars, referenceSystem, phaseIdx);
            const Scalar rhoFacet = massOrMolarDensity(facetVolVars, referenceSystem, phaseIdx);
            const Scalar rho = 0.5*(rhoInside + rhoFacet);

            componentFlux[compIdx] = fluxVarsCache.diffusionTijInside(phaseIdx, compIdx)*xInside
                                     + fluxVarsCache.diffusionTijFacet(phaseIdx, compIdx)*xFacet;

            if (!scvf.boundary())
                componentFlux[compIdx] += fluxVarsCache.diffusionTijOutside(phaseIdx, compIdx)*xOutside;
            componentFlux[compIdx] *= rho;

            if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx) && !FluidSystem::isTracerFluidSystem())
                componentFlux[FluidSystem::getMainComponent(phaseIdx)] -= componentFlux[compIdx];
        }

        return componentFlux;
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibility will be computed and stored using the method below.
    template< class Problem, class ElementVolumeVariables >
    static typename Cache::DiffusionTransmissibilityContainer
    calculateTransmissibility(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolumeFace& scvf,
                              unsigned int phaseIdx, unsigned int compIdx)
    {
        typename Cache::DiffusionTransmissibilityContainer tij;
        if (!problem.couplingManager().isCoupled(element, scvf))
        {
            //! use the standard Fick's law and only compute one transmissibility
            tij[Cache::insideTijIdx] = ParentType::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, compIdx);
            return tij;
        }

        //! xi factor for the coupling conditions
        static const Scalar xi = getParamFromGroup<Scalar>(problem.paramGroup(), "FacetCoupling.Xi", 1.0);

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const auto wIn = Extrusion::area(scvf)
                         *computeTpfaTransmissibility(scvf, insideScv,
                                                      insideVolVars.effectiveDiffusionCoefficient(phaseIdx, phaseIdx, compIdx),
                                                      insideVolVars.extrusionFactor());

        // proceed depending on the interior BC types used
        const auto iBcTypes = problem.interiorBoundaryTypes(element, scvf);

        // neumann-coupling
        if (iBcTypes.hasOnlyNeumann())
        {
            const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);
            const auto wFacet = 2.0*Extrusion::area(scvf)*insideVolVars.extrusionFactor()
                                   /facetVolVars.extrusionFactor()
                                   *vtmv(scvf.unitOuterNormal(),
                                         facetVolVars.effectiveDiffusionCoefficient(phaseIdx, phaseIdx, compIdx),
                                         scvf.unitOuterNormal());

            // The fluxes across this face and the outside face can be expressed in matrix form:
            // \f$\mathbf{C} \bar{\mathbf{u}} + \mathbf{D} \mathbf{u} + \mathbf{E} \mathbf{u}_\gamma\f$,
            // where \f$\gamma$\f denotes the domain living on the facets and \f$\bar{\mathbf{u}}$\f are
            // intermediate face unknowns in the matrix domain. Equivalently, flux continuity reads:
            // \f$\mathbf{A} \bar{\mathbf{u}} = \mathbf{B} \mathbf{u} + \mathbf{M} \mathbf{u}_\gamma\f$.
            // Combining the two, we can eliminate the intermediate unknowns and compute the transmissibilities
            // that allow the description of the fluxes as functions of the cell and Dirichlet mass/mole fractions only.
            if (!scvf.boundary())
            {
                const auto outsideScvIdx = scvf.outsideScvIdx();
                const auto& outsideVolVars = elemVolVars[outsideScvIdx];
                const auto wOut = -1.0*Extrusion::area(scvf)
                                  *computeTpfaTransmissibility(scvf, fvGeometry.scv(outsideScvIdx),
                                                               outsideVolVars.effectiveDiffusionCoefficient(phaseIdx, phaseIdx, compIdx),
                                                               outsideVolVars.extrusionFactor());

                if ( !Dune::FloatCmp::eq(xi, 1.0, 1e-6) )
                {
                    // optimized implementation: factorization obtained using sympy
                    // see CCTpfaFacetCouplingDarcysLaw for more details
                    const Scalar factor = wIn * wFacet / ( wIn * wOut * ( 2.0 * xi - 1.0 ) + wFacet * ( xi * ( wIn + wOut ) + wFacet ) );
                    tij[Cache::insideTijIdx]  = factor * ( wOut * xi + wFacet );
                    tij[Cache::outsideTijIdx] = factor * ( wOut * ( 1.0 - xi ) );
                    tij[Cache::facetTijIdx]   = factor * ( - wOut - wFacet );
                }
                else
                {
                    tij[Cache::insideTijIdx] = wFacet*wIn/(wIn+wFacet);
                    tij[Cache::facetTijIdx] = -tij[Cache::insideTijIdx];
                    tij[Cache::outsideTijIdx] = 0.0;
                }
            }
            else
            {
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
 * \brief Specialization of CCTpfaFacetCouplingFicksLawImpl for dim<dimWorld
 */
template<class TypeTag, ReferenceSystemFormulation referenceSystem>
class CCTpfaFacetCouplingFicksLawImpl<TypeTag, referenceSystem, /*isNetwork*/true>
: public FicksLawImplementation<TypeTag, DiscretizationMethod::cctpfa, referenceSystem>
{
    using Implementation = CCTpfaFacetCouplingFicksLawImpl<TypeTag, referenceSystem, true>;
    using ParentType = FicksLawImplementation<TypeTag, DiscretizationMethod::cctpfa, referenceSystem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using BalanceEqOpts = GetPropType<TypeTag, Properties::BalanceEqOpts>;

    static const int numPhases = ModelTraits::numFluidPhases();
    static const int numComponents = ModelTraits::numFluidComponents();

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;

    /*!
     * \brief The cache class used with this specialization of Fick's law.
     */
    class FacetCouplingFicksLawCache
    {
    public:
        //! export the corresponding filler class
        using Filler = typename ParentType::Cache::Filler;

        //! we store the transmissibilities associated with the interior
        //! cell and the fracture facet in an array. Access to this array
        //! should be done using the following indices:
        static constexpr int insideTijIdx = 0;
        static constexpr int facetTijIdx = 1;

        //! Export transmissibility storage type
        using DiffusionTransmissibilityContainer = std::array<Scalar, 2>;

        //! update subject to a given problem
        template< class Problem, class ElementVolumeVariables >
        void updateDiffusion(const Problem& problem,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace &scvf,
                             unsigned int phaseIdx,
                             unsigned int compIdx)
        {
            tij_[phaseIdx][compIdx] = Implementation::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, compIdx);
        }

        //! We use the same name as in the TpfaFicksLawCache so
        //! that this cache and the law implementation for non-coupled
        //! models can be reused here on facets that do not lie on an
        //! interior boundary, i.e. do not coincide with a facet element
        Scalar diffusionTij(unsigned int phaseIdx, unsigned int compIdx) const
        { return tij_[phaseIdx][compIdx][insideTijIdx]; }

        //! returns the transmissibility associated with the inside cell
        Scalar diffusionTijInside(unsigned int phaseIdx, unsigned int compIdx) const
        { return tij_[phaseIdx][compIdx][insideTijIdx]; }

        //! returns the transmissibility associated with the outside cell
        Scalar diffusionTijFacet(unsigned int phaseIdx, unsigned int compIdx) const
        { return tij_[phaseIdx][compIdx][facetTijIdx]; }

    private:
        std::array< std::array<DiffusionTransmissibilityContainer, numComponents>, numPhases> tij_;
    };

public:
    //! export the type for the corresponding cache
    using Cache = FacetCouplingFicksLawCache;

    //! state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;

    //! Return the reference system
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    //! Compute the diffusive fluxes
    template< class Problem, class ElementVolumeVariables, class ElementFluxVarsCache >
    static ComponentFluxVector flux(const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolumeFace& scvf,
                                    int phaseIdx,
                                    const ElementFluxVarsCache& elemFluxVarsCache)
    {
        if (!problem.couplingManager().isOnInteriorBoundary(element, scvf))
            return ParentType::flux(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarsCache);

        ComponentFluxVector componentFlux(0.0);
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            if(compIdx == FluidSystem::getMainComponent(phaseIdx))
                continue;

            // get inside/outside volume variables
            const auto& fluxVarsCache = elemFluxVarsCache[scvf];
            const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
            const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);

            // the inside and outside mass/mole fractions fractions
            const Scalar xInside = massOrMoleFraction(insideVolVars, referenceSystem, phaseIdx, compIdx);
            const Scalar xFacet = massOrMoleFraction(facetVolVars, referenceSystem, phaseIdx, compIdx);

            const Scalar rhoInside = massOrMolarDensity(insideVolVars, referenceSystem, phaseIdx);
            const Scalar rhoFacet = massOrMolarDensity(facetVolVars, referenceSystem, phaseIdx);
            const Scalar rho = 0.5*(rhoInside + rhoFacet);

            componentFlux[compIdx] = rho*(fluxVarsCache.diffusionTijInside(phaseIdx, compIdx)*xInside
                                          + fluxVarsCache.diffusionTijFacet(phaseIdx, compIdx)*xFacet);

            if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx) && !FluidSystem::isTracerFluidSystem())
                componentFlux[FluidSystem::getMainComponent(phaseIdx)] -= componentFlux[compIdx];
        }

        return componentFlux;
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibility will be computed and stored using the method below.
    template< class Problem, class ElementVolumeVariables >
    static typename Cache::DiffusionTransmissibilityContainer
    calculateTransmissibility(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolumeFace& scvf,
                              unsigned int phaseIdx, unsigned int compIdx)
    {
        typename Cache::DiffusionTransmissibilityContainer tij;
        if (!problem.couplingManager().isCoupled(element, scvf))
        {
            //! use the standard Fick's law and only compute one transmissibility
            tij[Cache::insideTijIdx] = ParentType::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, compIdx);
            return tij;
        }

        //! xi factor for the coupling conditions
        static const Scalar xi = getParamFromGroup<Scalar>(problem.paramGroup(), "FacetCoupling.Xi", 1.0);

        // On surface grids only xi = 1.0 can be used, as the coupling condition
        // for xi != 1.0 does not generalize for surface grids where the normal
        // vectors of the inside/outside elements have different orientations.
        if (Dune::FloatCmp::ne(xi, 1.0, 1e-6))
            DUNE_THROW(Dune::InvalidStateException, "Xi != 1.0 cannot be used on surface grids");

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const auto wIn = Extrusion::area(scvf)
                         *computeTpfaTransmissibility(scvf, insideScv,
                                                      insideVolVars.effectiveDiffusionCoefficient(phaseIdx, phaseIdx, compIdx),
                                                      insideVolVars.extrusionFactor());

        // proceed depending on the interior BC types used
        const auto iBcTypes = problem.interiorBoundaryTypes(element, scvf);

        // neumann-coupling
        if (iBcTypes.hasOnlyNeumann())
        {
            // Here we use the square root of the facet extrusion factor
            // as an approximate average distance from scvf ip to facet center
            using std::sqrt;
            const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);
            const auto wFacet = 2.0*Extrusion::area(scvf)*insideVolVars.extrusionFactor()
                                   /sqrt(facetVolVars.extrusionFactor())
                                   *vtmv(scvf.unitOuterNormal(),
                                         facetVolVars.effectiveDiffusionCoefficient(phaseIdx, phaseIdx, compIdx),
                                         scvf.unitOuterNormal());

            tij[Cache::insideTijIdx] = wFacet*wIn/(wIn+wFacet);
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
