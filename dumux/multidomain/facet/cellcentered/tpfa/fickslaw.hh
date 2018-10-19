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
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \copydoc Dumux::CCTpfaFacetCouplingFicksLaw
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FACET_COUPLING_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FACET_COUPLING_FICKS_LAW_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/discretization/cellcentered/tpfa/fickslaw.hh>

namespace Dumux {

//! Forward declaration of the implementation
template<class TypeTag, bool isNetwork>
class CCTpfaFacetCouplingFicksLawImpl;

/*!
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief The cache corresponding to tpfa Fick's Law with facet coupling
 * \note We distinguish between network and non-network grids here. Specializations
 *       for the two cases can be found below.
 */
template<class TypeTag, bool isNetwork>
class CCTpfaFacetCouplingFicksLawCache;

/*!
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief Fick's law for cell-centered finite volume schemes with two-point flux approximation
 *        in the context of coupled models where the coupling occurs across the facets of the bulk
 *        domain elements with a lower-dimensional domain living on these facets.
 *
 * \tparam TypeTag the problem type tag
 */
template<class TypeTag>
using CCTpfaFacetCouplingFicksLaw =
      CCTpfaFacetCouplingFicksLawImpl< TypeTag, ( int(GET_PROP_TYPE(TypeTag, GridView)::dimension) <
                                                  int(GET_PROP_TYPE(TypeTag, GridView)::dimensionworld) ) >;

/*!
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief Specialization of the FacetCouplingTpfaFicksLawCache for non-network grids.
 */
template<class TypeTag>
class CCTpfaFacetCouplingFicksLawCache<TypeTag, /*isNetwork*/false>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;

    static constexpr int numPhases = GET_PROP_TYPE(TypeTag, ModelTraits)::numPhases();
    static constexpr int numComponents = GET_PROP_TYPE(TypeTag, ModelTraits)::numComponents();

    // the standard tpfa fick's law implementation
    using DiffusionType = typename GET_PROP_TYPE(TypeTag, MolecularDiffusionType);
    using TpfaFicksLaw = FicksLawImplementation<TypeTag, DiscretizationMethod::cctpfa>;

public:
    //! export the corresponding filler class (use standard tpfa one)
    using Filler = typename TpfaFicksLaw::Cache::Filler;

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
                         unsigned int phaseIdx, unsigned int compIdx)
    {
        tij_[phaseIdx][compIdx] = DiffusionType::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, compIdx);
    }

    //! We use the same name as in the TpfaFicksLawCache so
    //! that this cache and the law implementation for non-coupled
    //! models can be reused here on facets that do not lie on an
    //! interior boundary, i.e. do not coincide with a fracture
    Scalar diffusionTij(unsigned int phaseIdx, unsigned int compIdx) const
    { return tij_[phaseIdx][compIdx][insideTijIdx]; }

    //! returns the transmissibility associated with the inside cell
    Scalar diffusionTijInside(unsigned int phaseIdx, unsigned int compIdx) const
    { return tij_[phaseIdx][compIdx][insideTijIdx]; }

    //! returns the transmissibility associated with the outside cell
    Scalar diffusionTijOutside(unsigned int phaseIdx, unsigned int compIdx) const
    {return tij_[phaseIdx][compIdx][outsideTijIdx]; }

    //! returns the transmissibility associated with the outside cell
    Scalar diffusionTijFacet(unsigned int phaseIdx, unsigned int compIdx) const
    {return tij_[phaseIdx][compIdx][facetTijIdx]; }

private:
    std::array< std::array< std::array<Scalar, 3>, numComponents>, numPhases > tij_;
};

/*!
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief Specialization of the CCTpfaFicksLaw grids where dim=dimWorld
 */
template<class TypeTag>
class CCTpfaFacetCouplingFicksLawImpl<TypeTag, /*isNetwork*/false>
{
    using ThisType = CCTpfaFacetCouplingFicksLawImpl<TypeTag, false>;
    using TpfaFicksLaw = FicksLawImplementation<TypeTag, DiscretizationMethod::cctpfa>;

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;

    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using IndexType = typename GridView::IndexSet::IndexType;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BalanceEqOpts = typename GET_PROP_TYPE(TypeTag, BalanceEqOpts);

    using ScalarType = typename GET_PROP_TYPE(TypeTag, Scalar);
    static constexpr int numPhases = GET_PROP_TYPE(TypeTag, ModelTraits)::numPhases();
    static constexpr int numComponents = GET_PROP_TYPE(TypeTag, ModelTraits)::numComponents();
    using ComponentFluxVector = Dune::FieldVector<ScalarType, numComponents>;

  public:
    //! state the scalar type of the law
    using Scalar = ScalarType;
    //! export the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;
    //! export the type for the corresponding cache
    using Cache = CCTpfaFacetCouplingFicksLawCache<TypeTag, /*isNetwork*/false>;
    //! export the type used to store transmissibilities
    using TijContainer = typename Cache::DiffusionTransmissibilityContainer;


    //! Compute the diffusive flux
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
            return TpfaFicksLaw::flux(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarsCache);

        ComponentFluxVector componentFlux(0.0);
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            if(compIdx == FluidSystem::getMainComponent(phaseIdx))
                continue;

            // Obtain inside and fracture pressures
            const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
            const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);
            const auto xInside = insideVolVars.moleFraction(phaseIdx, compIdx);
            const auto xFacet = facetVolVars.moleFraction(phaseIdx, compIdx);

            // compute and return flux
            const auto& fluxVarsCache = elemFluxVarsCache[scvf];
            Scalar flux = fluxVarsCache.diffusionTijInside(phaseIdx, compIdx)*xInside
                          + fluxVarsCache.diffusionTijFacet(phaseIdx, compIdx)*xFacet;

            if (!scvf.boundary())
                flux += fluxVarsCache.diffusionTijOutside(phaseIdx, compIdx)
                        *elemVolVars[scvf.outsideScvIdx()].moleFraction(phaseIdx, compIdx);

            // for the density, use arithmetic average
            const auto rho = 0.5*(insideVolVars.molarDensity(phaseIdx) + facetVolVars.molarDensity(phaseIdx));

            componentFlux[compIdx] = rho*flux;
            if (BalanceEqOpts::mainComponentIsBalanced(phaseIdx) && !FluidSystem::isTracerFluidSystem())
                componentFlux[FluidSystem::getMainComponent(phaseIdx)] -= componentFlux[compIdx];
        }

        return componentFlux;
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibility will be computed and stored using the method below.
    template< class Problem, class ElementVolumeVariables >
    static TijContainer calculateTransmissibility(const Problem& problem,
                                                  const Element& element,
                                                  const FVElementGeometry& fvGeometry,
                                                  const ElementVolumeVariables& elemVolVars,
                                                  const SubControlVolumeFace& scvf,
                                                  unsigned int phaseIdx, unsigned int compIdx)
    {
        TijContainer tij;
        if (!problem.couplingManager().isCoupled(element, scvf))
        {
            //! use the standard darcy's law and only compute one transmissibility
            tij[Cache::insideTijIdx] = TpfaFicksLaw::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, compIdx);
            return tij;
        }

        //! xi factor for the coupling conditions
        static const Scalar xi = getParamFromGroup<Scalar>(problem.paramGroup(), "FacetCoupling.Xi", 1.0);
        static const Scalar oneMinusXi = 1.0 - xi;

        using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const auto insideD = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(),
                                                                insideVolVars.saturation(phaseIdx),
                                                                insideVolVars.diffusionCoefficient(phaseIdx, compIdx));
        const auto wIn = scvf.area()*computeTpfaTransmissibility(scvf, insideScv, insideD, insideVolVars.extrusionFactor());

        // proceed depending on the interior BC types used
        const auto iBcTypes = problem.interiorBoundaryTypes(element, scvf);

        // neumann-coupling
        if (iBcTypes.hasOnlyNeumann())
        {
            const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);
            // TODO How can we tell this class which effective diffusivity law is used
            //      in the facet domain? Here, we simply assume it uses the same law
            //      as this domain. But, we cannot make it an additional template parameter
            //      because that leads to a compiler error for models that do not specify
            //      an effective diffusivity law, e.g. models that do not consider diffusion.
            const auto facetD = EffDiffModel::effectiveDiffusivity(facetVolVars.porosity(),
                                                                   facetVolVars.saturation(phaseIdx),
                                                                   facetVolVars.diffusionCoefficient(phaseIdx, compIdx));
            const auto wFacet = 2.0*scvf.area()*insideVolVars.extrusionFactor()
                                   /facetVolVars.extrusionFactor()
                                   *vtmv(scvf.unitOuterNormal(), facetD, scvf.unitOuterNormal());

            // The fluxes across this face and the outside face can be expressed in matrix form:
            // \f$\mathbf{C} \bar{\mathbf{u}} + \mathbf{D} \mathbf{u} + \mathbf{E} \mathbf{u}_\gamma\f$,
            // where \f$\gamma$\f denotes the domain living on the facets and \f$\bar{\mathbf{u}}$\f are
            // intermediate face unknowns in the matrix domain. Equivalently, flux continuity reads:
            // \f$\mathbf{A} \bar{\mathbf{u}} = \mathbf{B} \mathbf{u} + \mathbf{M} \mathbf{u}_\gamma\f$.
            // Combining the two, we can eliminate the intermediate unknowns and compute the transmissibilities.
            if (!scvf.boundary() && xi != 1.0)
            {
                const auto outsideScvIdx = scvf.outsideScvIdx();
                const auto& outsideVolVars = elemVolVars[outsideScvIdx];
                const auto outsideD = EffDiffModel::effectiveDiffusivity(outsideVolVars.porosity(),
                                                                         outsideVolVars.saturation(phaseIdx),
                                                                         outsideVolVars.diffusionCoefficient(phaseIdx, compIdx));
                const auto wOut = -1.0*scvf.area()*computeTpfaTransmissibility(scvf,
                                                                               fvGeometry.scv(outsideScvIdx),
                                                                               outsideD,
                                                                               outsideVolVars.extrusionFactor());
                const Scalar xiWIn = xi*wIn;
                const Scalar xiWOut = xi*wOut;
                const Scalar oneMinusXiWIn = oneMinusXi*wIn;
                const Scalar oneMinusXiWOut = oneMinusXi*wOut;

                // assemble matrices
                Dune::FieldMatrix<Scalar, 2, 2> A, B;
                A[0][0] = xiWIn+wFacet;
                A[0][1] = oneMinusXiWOut;
                A[1][0] = oneMinusXiWIn;
                A[1][1] = xiWOut-wFacet;

                B[0][0] = xiWIn;
                B[0][1] = oneMinusXiWOut;
                B[1][0] = oneMinusXiWIn;
                B[1][1] = xiWOut;

                // tij = C(A^-1)B
                const Scalar detA = A[0][0]*A[1][1] - A[1][0]*A[0][1];
                tij[Cache::insideTijIdx] = xiWIn - xiWIn*(A[1][1]*B[0][0] - A[0][1]*B[1][0])/detA;
                tij[Cache::outsideTijIdx] = -xiWIn*(A[1][1]*B[0][1] - A[0][1]*B[1][1])/detA;
                tij[Cache::facetTijIdx] = -xiWIn*wFacet*(A[1][1] + A[0][1])/detA;
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

} // end namespace Dumux

#endif
