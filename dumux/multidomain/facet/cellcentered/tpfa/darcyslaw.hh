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
 * \copydoc Dumux::CCTpfaFacetCouplingDarcysLaw
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FACET_COUPLING_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FACET_COUPLING_DARCYS_LAW_HH

#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/discretization/cellcentered/tpfa/darcyslaw.hh>

namespace Dumux {

//! Forward declaration of the implementation
template<class TypeTag, bool isNetwork>
class CCTpfaFacetCouplingDarcysLawImpl;

/*!
 * \ingroup MixedDimensionFacet
 * \brief The cache corresponding to tpfa Darcy's Law with facet coupling
 * \note We distinguish between network and non-network grids here. Specializations
 *       for the two cases can be found below.
 */
template<class TypeTag, bool isNetwork>
class CCTpfaFacetCouplingDarcysLawCache;

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Darcy's law for cell-centered finite volume schemes with two-point flux approximation
 *        in the context of coupled models where the coupling occurs across the facets of the bulk
 *        domain elements with a lower-dimensional domain living on these facets.
 */
template<class TypeTag>
using CCTpfaFacetCouplingDarcysLaw =
      CCTpfaFacetCouplingDarcysLawImpl< TypeTag, ( int(GET_PROP_TYPE(TypeTag, GridView)::dimension) <
                                                   int(GET_PROP_TYPE(TypeTag, GridView)::dimensionworld) ) >;

/*!
 * \ingroup MixedDimensionFacet
 * \brief Specialization of the FacetCouplingTpfaDarcysLawCache for non-network grids.
 */
template<class TypeTag>
class CCTpfaFacetCouplingDarcysLawCache<TypeTag, /*isNetwork*/false>
{
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using Scalar = typename GridVariables::Scalar;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;

public:
    //! export the corresponding filler class
    using Filler = TpfaDarcysLawCacheFiller<FVGridGeometry>;

    //! we store the transmissibilities associated with the interior
    //! cell, outside cell, and the fracture facet in an array. Access
    //! to this array should be done using the following indices:
    static constexpr int insideTijIdx = 0;
    static constexpr int outsideTijIdx = 1;
    static constexpr int facetTijIdx = 2;

    //! Export transmissibility storage type
    using AdvectionTransmissibilityContainer = std::array<Scalar, 3>;

    //! update subject to a given problem
    template< class Problem >
    void updateAdvection(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace &scvf)
    {
        tij_ = AdvectionType::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf);
    }

    //! We use the same name as in the TpfaDarcysLawCache so
    //! that this cache and the law implementation for non-coupled
    //! models can be reused here on facets that do not lie on an
    //! interior boundary, i.e. do not coincide with a fracture
    Scalar advectionTij() const { return tij_[insideTijIdx]; }

    //! returns the transmissibility associated with the inside cell
    Scalar advectionTijInside() const { return tij_[insideTijIdx]; }

    //! returns the transmissibility associated with the outside cell
    Scalar advectionTijOutside() const {return tij_[outsideTijIdx]; }

    //! returns the transmissibility associated with the outside cell
    Scalar advectionTijFacet() const {return tij_[facetTijIdx]; }

private:
    std::array<Scalar, 3> tij_;
};

/*!
 * \ingroup MixedDimensionFacet
 * \brief Specialization of the CCTpfaDarcysLaw grids where dim=dimWorld
 */
template<class TypeTag>
class CCTpfaFacetCouplingDarcysLawImpl<TypeTag, /*isNetwork*/false>
{
    using TpfaDarcysLaw = DarcysLawImplementation<TypeTag, DiscretizationMethod::cctpfa>;

    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using Scalar = typename GridVariables::Scalar;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVarsCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using FluxVariablesCache = typename ElementFluxVarsCache::FluxVariablesCache;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;

    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using IndexType = typename GridView::IndexSet::IndexType;

  public:
    //! export the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;
    //! export the type for the corresponding cache
    using Cache = CCTpfaFacetCouplingDarcysLawCache<TypeTag, /*isNetwork*/false>;
    //! export the type used to store transmissibilities
    using TijContainer = typename Cache::AdvectionTransmissibilityContainer;

    //! Compute the advective flux
    template< class Problem >
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        static const Scalar gravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
        if (gravity)
            DUNE_THROW(Dune::NotImplemented, "gravity for darcys law with facet coupling");

        if (!problem.couplingManager().isCoupled(element, scvf))
            return TpfaDarcysLaw::flux(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarsCache);

        // Obtain inside and fracture pressures
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto pInside = insideVolVars.pressure(phaseIdx);
        const auto pFacet = problem.couplingManager().getLowDimVolVars(element, scvf).pressure(phaseIdx);

        // return flux
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        return scvf.boundary() ? fluxVarsCache.advectionTijInside()*pInside + fluxVarsCache.advectionTijFacet()*pFacet
                               : fluxVarsCache.advectionTijInside()*pInside
                                 + fluxVarsCache.advectionTijOutside()*elemVolVars[scvf.outsideScvIdx()].pressure(phaseIdx)
                                 + fluxVarsCache.advectionTijFacet()*pFacet;
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibility will be computed and stored using the method below.
    template< class Problem >
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
            tij[FluxVariablesCache::insideTijIdx] = TpfaDarcysLaw::calculateTransmissibility(problem,
                                                                                             element,
                                                                                             fvGeometry,
                                                                                             elemVolVars,
                                                                                             scvf);
            return tij;
        }

        //! xi factor for the coupling conditions
        static const Scalar xi = getParamFromGroup<Scalar>(problem.paramGroup(), "FacetCoupling.Xi", 1.0);
        static const Scalar oneMinusXi = 1.0 - xi;

        const auto area = scvf.area();
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
            const auto wFacet = 2.0*area*insideVolVars.extrusionFactor()
                                        /facetVolVars.extrusionFactor()
                                        *vtmv(scvf.unitOuterNormal(), facetVolVars.permeability(), scvf.unitOuterNormal());

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
                const auto wOut = -1.0*area*computeTpfaTransmissibility(scvf,
                                                                        fvGeometry.scv(outsideScvIdx),
                                                                        outsideVolVars.permeability(),
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
                tij[FluxVariablesCache::insideTijIdx] = xiWIn - xiWIn*(A[1][1]*B[0][0] - A[0][1]*B[1][0])/detA;
                tij[FluxVariablesCache::outsideTijIdx] = -xiWIn*(A[1][1]*B[0][1] - A[0][1]*B[1][1])/detA;
                tij[FluxVariablesCache::facetTijIdx] = -xiWIn*wFacet*(A[1][1] + A[0][1])/detA;
            }
            else
            {
                // TODO: check for division by zero??
                tij[FluxVariablesCache::insideTijIdx] = wFacet*wIn/(wIn+wFacet);
                tij[FluxVariablesCache::facetTijIdx] = -tij[FluxVariablesCache::insideTijIdx];
                tij[FluxVariablesCache::outsideTijIdx] = 0.0;
            }
        }
        else if (iBcTypes.hasOnlyDirichlet())
        {
            tij[FluxVariablesCache::insideTijIdx] = wIn;
            tij[FluxVariablesCache::outsideTijIdx] = 0.0;
            tij[FluxVariablesCache::facetTijIdx] = -wIn;
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Interior boundary types other than pure Dirichlet or Neumann");

        return tij;
    }
};

/*!
 * \ingroup MixedDimensionFacet
 * \brief Specialization of the FacetCouplingTpfaDarcysLawCache for network grids
 */
template<class TypeTag>
class CCTpfaFacetCouplingDarcysLawCache<TypeTag, /*isNetwork*/true>
{
    using AdvectionType = typename GET_PROP_TYPE(TypeTag, AdvectionType);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using Scalar = typename GridVariables::Scalar;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;

public:
    //! export the corresponding filler class
    using Filler = TpfaDarcysLawCacheFiller<FVGridGeometry>;

    //! we store the transmissibilities associated with the interior
    //! cell, outside cell, and the fracture facet in an array. Access
    //! to this array should be done using the following indices:
    static constexpr int insideTijIdx = 0;
    static constexpr int facetTijIdx = 1;
    static constexpr int firstOutsideTijIdx = 2;

    //! Export transmissibility storage type
    using AdvectionTransmissibilityContainer = std::vector<Scalar>;

    //! update subject to a given problem
    template< class Problem >
    void updateAdvection(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace &scvf)
    {
        tij_ = AdvectionType::calculateTransmissibility(problem, element, fvGeometry, elemVolVars, scvf);
    }

    //! We use the same name as in the TpfaDarcysLawCache so
    //! that this cache and the law implementation for non-coupled
    //! models can be reused here on facets that do not lie on an
    //! interior boundary, i.e. do not coincide with a fracture
    Scalar advectionTij() const { return tij_[insideTijIdx]; }

    //! returns the transmissibility associated with the inside cell
    Scalar advectionTijInside() const { return tij_[insideTijIdx]; }

    //! returns the transmissibility associated with the outside cell
    Scalar advectionTijOutside(unsigned int idxInOutside) const {return tij_[firstOutsideTijIdx+idxInOutside]; }

    //! returns the transmissibility associated with the outside cell
    Scalar advectionTijFacet() const {return tij_[facetTijIdx]; }

private:
    std::vector<Scalar> tij_;
};

/*!
 * \ingroup MixedDimensionFacet
 * \brief Specialization of the CCTpfaDarcysLaw grids where dim<dimWorld
 */
template<class TypeTag>
class CCTpfaFacetCouplingDarcysLawImpl<TypeTag, /*isNetwork*/true>
{
    using TpfaDarcysLaw = DarcysLawImplementation<TypeTag, DiscretizationMethod::cctpfa>;

    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using Scalar = typename GridVariables::Scalar;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVarsCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using FluxVariablesCache = typename ElementFluxVarsCache::FluxVariablesCache;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;

    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using IndexType = typename GridView::IndexSet::IndexType;

  public:
    //! state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;
    //! state the type for the corresponding cache
    using Cache = CCTpfaFacetCouplingDarcysLawCache<TypeTag, /*isNetwork*/true>;
    //! export the type used to store transmissibilities
    using TijContainer = typename Cache::AdvectionTransmissibilityContainer;

    //! Compute the advective flux
    template< class Problem >
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        static const Scalar gravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
        if (gravity)
            DUNE_THROW(Dune::NotImplemented, "gravity for darcys law with facet coupling");

        if (!problem.couplingManager().isCoupled(element, scvf))
            return TpfaDarcysLaw::flux(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, elemFluxVarsCache);

        // Obtain inside and fracture pressures
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto pInside = insideVolVars.pressure(phaseIdx);
        const auto pFacet = problem.couplingManager().getLowDimVolVars(element, scvf).pressure(phaseIdx);

        // return flux
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        if (scvf.boundary())
            return fluxVarsCache.advectionTijInside()*pInside + fluxVarsCache.advectionTijFacet()*pFacet;
        else
        {
            Scalar flux = fluxVarsCache.advectionTijInside()*pInside + fluxVarsCache.advectionTijFacet()*pFacet;
            for (unsigned int idxInOutside = 0; idxInOutside < scvf.numOutsideScvs(); ++idxInOutside)
                flux += fluxVarsCache.advectionTijOutside(idxInOutside)
                        *elemVolVars[scvf.outsideScvIdx(idxInOutside)].pressure(phaseIdx);
            return flux;
        }
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibility will be computed and stored using the method below.
    template< class Problem >
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
            tij.resize(1);
            tij[FluxVariablesCache::insideTijIdx] = TpfaDarcysLaw::calculateTransmissibility(problem,
                                                                                             element,
                                                                                             fvGeometry,
                                                                                             elemVolVars,
                                                                                             scvf);
            return tij;
        }

        //! xi factor for the coupling conditions
        static const Scalar xi = getParamFromGroup<Scalar>(problem.paramGroup(), "FacetCoupling.Xi", 1.0);
        static const Scalar oneMinusXi = 1.0 - xi;

        const auto area = scvf.area();
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        const auto wIn = area*computeTpfaTransmissibility(scvf, insideScv, insideVolVars.permeability(), insideVolVars.extrusionFactor());

        // resize transmissibility container
        const auto numOutsideScvs = scvf.numOutsideScvs();
        tij.resize(2+numOutsideScvs);
        std::fill(tij.begin(), tij.end(), 0.0);

        // proceed depending on the interior BC types used
        const auto iBcTypes = problem.interiorBoundaryTypes(element, scvf);

        // neumann-coupling
        if (iBcTypes.hasOnlyNeumann())
        {
            const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);
            using std::sqrt;
            // Here we use the square root of the facet extrusion factor
            // as an approximate average distance from scvf ip to facet center
            const auto wFacet = 2.0*area*insideVolVars.extrusionFactor()
                                        /sqrt(facetVolVars.extrusionFactor())
                                        *vtmv(scvf.unitOuterNormal(), facetVolVars.permeability(), scvf.unitOuterNormal());

            // The fluxes across this face and the outside face can be expressed in matrix form:
            // \f$\mathbf{C} \bar{\mathbf{u}} + \mathbf{D} \mathbf{u} + \mathbf{E} \mathbf{u}_\gamma\f$,
            // where \f$\gamma$\f denotes the domain living on the facets and \f$\bar{\mathbf{u}}$\f are
            // intermediate face unknowns in the matrix domain. Equivalently, flux continuity reads:
            // \f$\mathbf{A} \bar{\mathbf{u}} = \mathbf{B} \mathbf{u} + \mathbf{M} \mathbf{u}_\gamma\f$.
            // Combining the two, we can eliminate the intermediate unknowns and compute the transmissibilities.
            if (!scvf.boundary() && xi != 1.0)
            {
                // assemble matrices
                const Scalar xiWIn = xi*wIn;
                const Scalar oneMinusXiWIn = oneMinusXi*wIn;
                const auto numDofs = numOutsideScvs+1;

                Dune::DynamicMatrix<Scalar> A(numDofs, numDofs, 0.0);
                Dune::DynamicMatrix<Scalar> B(numDofs, numDofs, 0.0);
                Dune::DynamicVector<Scalar> M(numDofs, 0.0);

                A[0][0] = xiWIn+wFacet;
                B[0][0] = xiWIn;
                M[0] = wFacet;

                for (unsigned int i = 0; i < numOutsideScvs; ++i)
                {
                    const auto outsideScvIdx = scvf.outsideScvIdx(i);
                    const auto& outsideVolVars = elemVolVars[outsideScvIdx];
                    const auto& flipScvf = fvGeometry.flipScvf(scvf.index(), i);
                    const auto wOut = area*computeTpfaTransmissibility(flipScvf,
                                                                       fvGeometry.scv(outsideScvIdx),
                                                                       outsideVolVars.permeability(),
                                                                       outsideVolVars.extrusionFactor());
                    const auto wFacetOut = 2.0*area*insideVolVars.extrusionFactor()
                                                   /sqrt(facetVolVars.extrusionFactor())
                                                   *vtmv(flipScvf.unitOuterNormal(),
                                                         facetVolVars.permeability(),
                                                         flipScvf.unitOuterNormal());
                    // assemble local system matrices
                    const auto xiWOut = xi*wOut;
                    const auto oneMinusXiWOut = oneMinusXi*wOut;
                    const auto curDofIdx = i+1;

                    M[curDofIdx] = wFacetOut;
                    A[curDofIdx][0] += oneMinusXiWIn;
                    B[curDofIdx][0] += oneMinusXiWIn;

                    for (unsigned int otherDofIdx = 0; otherDofIdx < numDofs; ++otherDofIdx)
                    {
                        if (otherDofIdx == curDofIdx)
                        {
                            A[curDofIdx][curDofIdx] += xiWOut+wFacetOut;
                            B[curDofIdx][curDofIdx] += xiWOut;
                        }
                        else
                        {
                            A[otherDofIdx][curDofIdx] += oneMinusXiWOut;
                            B[otherDofIdx][curDofIdx] += oneMinusXiWOut;
                        }
                    }
                }

                A.invert();

                // compute transmissibilities
                for (unsigned int i = 0; i < numDofs; ++i)
                {
                    tij[FluxVariablesCache::insideTijIdx] -= A[0][i]*B[i][0];
                    tij[FluxVariablesCache::facetTijIdx] -= A[0][i]*M[i];
                    for (unsigned int idxInOutside = 0; idxInOutside < numOutsideScvs; ++idxInOutside)
                        tij[FluxVariablesCache::firstOutsideTijIdx+idxInOutside] -= A[0][i]*B[i][idxInOutside+1];
                }
                std::for_each(tij.begin(), tij.end(), [xiWIn] (auto& t) { t *= xiWIn; });
                tij[FluxVariablesCache::insideTijIdx] += xiWIn;
            }
            else
            {
                // TODO: check for division by zero??
                tij[FluxVariablesCache::insideTijIdx] = wFacet*wIn/(wIn+wFacet);
                tij[FluxVariablesCache::facetTijIdx] = -tij[FluxVariablesCache::insideTijIdx];
            }
        }
        else if (iBcTypes.hasOnlyDirichlet())
        {
            tij[FluxVariablesCache::insideTijIdx] = wIn;
            tij[FluxVariablesCache::facetTijIdx] = -wIn;
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Interior boundary types other than pure Dirichlet or Neumann");

        return tij;
    }
};

} // end namespace Dumux

#endif
