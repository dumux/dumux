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
 * \copydoc Dumux::CCTpfaFacetCouplingFouriersLaw
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FACET_COUPLING_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FACET_COUPLING_FOURIERS_LAW_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/discretization/cellcentered/tpfa/fourierslaw.hh>

namespace Dumux {

//! Forward declaration of the implementation
template<class TypeTag, bool isNetwork>
class CCTpfaFacetCouplingFouriersLawImpl;

/*!
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief The cache corresponding to tpfa Fourier's Law with facet coupling
 * \note We distinguish between network and non-network grids here. Specializations
 *       for the two cases can be found below.
 */
template<class TypeTag, bool isNetwork>
class CCTpfaFacetCouplingFouriersLawCache;

/*!
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief Fourier's law for cell-centered finite volume schemes with two-point flux approximation
 *        in the context of coupled models where the coupling occurs across the facets of the bulk
 *        domain elements with a lower-dimensional domain living on these facets.
 *
 * \tparam TypeTag the problem type tag
 */
template<class TypeTag>
using CCTpfaFacetCouplingFouriersLaw =
      CCTpfaFacetCouplingFouriersLawImpl< TypeTag, ( int(GET_PROP_TYPE(TypeTag, GridView)::dimension) <
                                                     int(GET_PROP_TYPE(TypeTag, GridView)::dimensionworld) ) >;

/*!
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief Specialization of the FacetCouplingTpfaFouriersLawCache for non-network grids.
 */
template<class TypeTag>
class CCTpfaFacetCouplingFouriersLawCache<TypeTag, /*isNetwork*/false>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;

    // the standard tpfa fourier's law implementation
    using HeatConductionType = typename GET_PROP_TYPE(TypeTag, HeatConductionType);
    using TpfaFouriersLaw = FouriersLawImplementation<TypeTag, DiscretizationMethod::cctpfa>;

public:
    //! export the corresponding filler class (use standard tpfa one)
    using Filler = typename TpfaFouriersLaw::Cache::Filler;

    //! we store the transmissibilities associated with the interior
    //! cell, outside cell, and the fracture facet in an array. Access
    //! to this array should be done using the following indices:
    static constexpr int insideTijIdx = 0;
    static constexpr int outsideTijIdx = 1;
    static constexpr int facetTijIdx = 2;

    //! Export transmissibility storage type
    using HeatConductionTransmissibilityContainer = std::array<Scalar, 3>;

    //! update subject to a given problem
    template< class Problem, class ElementVolumeVariables >
    void updateHeatConduction(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolumeFace &scvf)
    {
        tij_ = HeatConductionType::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf);
    }

    //! We use the same name as in the TpfaFicksLawCache so
    //! that this cache and the law implementation for non-coupled
    //! models can be reused here on facets that do not lie on an
    //! interior boundary, i.e. do not coincide with a fracture
    Scalar heatConductionTij() const
    { return tij_[insideTijIdx]; }

    //! returns the transmissibility associated with the inside cell
    Scalar heatConductionTijInside() const
    { return tij_[insideTijIdx]; }

    //! returns the transmissibility associated with the outside cell
    Scalar heatConductionTijOutside() const
    {return tij_[outsideTijIdx]; }

    //! returns the transmissibility associated with the outside cell
    Scalar heatConductionTijFacet() const
    {return tij_[facetTijIdx]; }

private:
    HeatConductionTransmissibilityContainer tij_;
};

/*!
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief Specialization of the CCTpfaFouriersLaw grids where dim=dimWorld
 */
template<class TypeTag>
class CCTpfaFacetCouplingFouriersLawImpl<TypeTag, /*isNetwork*/false>
{
    using ThisType = CCTpfaFacetCouplingFouriersLawImpl<TypeTag, false>;
    using TpfaFouriersLaw = FouriersLawImplementation<TypeTag, DiscretizationMethod::cctpfa>;

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;

    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using IndexType = typename GridView::IndexSet::IndexType;

  public:
    //! state the scalar type of the law
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    //! export the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;
    //! export the type for the corresponding cache
    using Cache = CCTpfaFacetCouplingFouriersLawCache<TypeTag, /*isNetwork*/false>;
    //! export the type used to store transmissibilities
    using TijContainer = typename Cache::HeatConductionTransmissibilityContainer;


    //! Compute the diffusive flux
    template< class Problem, class ElementVolumeVariables, class ElementFluxVarsCache >
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        if (!problem.couplingManager().isOnInteriorBoundary(element, scvf))
            return TpfaFouriersLaw::flux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        // Obtain inside and fracture pressures
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);
        const auto tInside = insideVolVars.temperature();
        const auto tFacet = facetVolVars.temperature();

        // compute and return flux
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        Scalar flux = fluxVarsCache.heatConductionTijInside()*tInside
                      + fluxVarsCache.heatConductionTijFacet()*tFacet;

        if (!scvf.boundary())
            flux += fluxVarsCache.heatConductionTijOutside()*elemVolVars[scvf.outsideScvIdx()].temperature();

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
        TijContainer tij;
        if (!problem.couplingManager().isCoupled(element, scvf))
        {
            //! use the standard darcy's law and only compute one transmissibility
            tij[Cache::insideTijIdx] = TpfaFouriersLaw::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf);
            return tij;
        }

        //! xi factor for the coupling conditions
        static const Scalar xi = getParamFromGroup<Scalar>(problem.paramGroup(), "FacetCoupling.Xi", 1.0);
        static const Scalar oneMinusXi = 1.0 - xi;

        using EffThermCondModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const auto insideLambda = EffThermCondModel::effectiveThermalConductivity(insideVolVars,
                                                                                  problem.spatialParams(),
                                                                                  element,
                                                                                  fvGeometry,
                                                                                  insideScv);
        const auto wIn = scvf.area()*computeTpfaTransmissibility(scvf, insideScv, insideLambda, insideVolVars.extrusionFactor());

        // proceed depending on the interior BC types used
        const auto iBcTypes = problem.interiorBoundaryTypes(element, scvf);

        // neumann-coupling
        if (iBcTypes.hasOnlyNeumann())
        {
            const auto& facetData = problem.couplingManager().getLowDimCouplingData(element, scvf);

            // TODO How can we tell this class which effective thermal conductivity law is used
            //      in the facet domain? Here, we simply assume it uses the same law
            //      as this domain. But, we cannot make it an additional template parameter
            //      because that leads to a compiler error for models that do not specify
            //      an effective diffusivity law, e.g. models that do not consider diffusion.
            const auto facetLambda = EffThermCondModel::effectiveThermalConductivity(facetData.volVars(),
                                                                                     facetData.problem().spatialParams(),
                                                                                     facetData.element(),
                                                                                     facetData.fvGeometry(),
                                                                                     facetData.scv());

            const auto wFacet = 2.0*scvf.area()*insideVolVars.extrusionFactor()
                                   /facetData.volVars().extrusionFactor()
                                   *vtmv(scvf.unitOuterNormal(), facetLambda, scvf.unitOuterNormal());

            // The fluxes across this face and the outside face can be expressed in matrix form:
            // \f$\mathbf{C} \bar{\mathbf{u}} + \mathbf{D} \mathbf{u} + \mathbf{E} \mathbf{u}_\gamma\f$,
            // where \f$\gamma$\f denotes the domain living on the facets and \f$\bar{\mathbf{u}}$\f are
            // intermediate face unknowns in the matrix domain. Equivalently, flux continuity reads:
            // \f$\mathbf{A} \bar{\mathbf{u}} = \mathbf{B} \mathbf{u} + \mathbf{M} \mathbf{u}_\gamma\f$.
            // Combining the two, we can eliminate the intermediate unknowns and compute the transmissibilities.
            if (!scvf.boundary() && xi != 1.0)
            {
                const auto outsideScvIdx = scvf.outsideScvIdx();
                const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
                const auto& outsideVolVars = elemVolVars[outsideScvIdx];
                const auto outsideElement = problem.fvGridGeometry().element(outsideScvIdx);
                const auto outsideLambda = EffThermCondModel::effectiveThermalConductivity(outsideVolVars,
                                                                                           problem.spatialParams(),
                                                                                           outsideElement,
                                                                                           fvGeometry,
                                                                                           outsideScv);
                const auto wOut = -1.0*scvf.area()*computeTpfaTransmissibility(scvf,
                                                                               outsideScv,
                                                                               outsideLambda,
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
