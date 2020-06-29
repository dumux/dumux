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
 * \ingroup DarcyDarcyCoupling
 * \brief Coupling manager for equal-dimension boundary coupling
 */

#ifndef DUMUX_MULTIDOMAIN_DARCYDARCY_BOUNDARY_COUPLINGMANAGER_HH
#define DUMUX_MULTIDOMAIN_DARCYDARCY_BOUNDARY_COUPLINGMANAGER_HH

#include <iostream>
#include <vector>

#include <dune/common/indices.hh>
#include <dune/common/exceptions.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/boundary/darcydarcy/couplingmapper.hh>

namespace Dumux {

/*!
 * \ingroup DarcyDarcyCoupling
 * \brief Coupling manager for equal-dimension boundary coupling of darcy models
 */
template<class MDTraits>
class DarcyDarcyBoundaryCouplingManager
: public CouplingManager<MDTraits>
{
    using ParentType = CouplingManager<MDTraits>;

    using Scalar = typename MDTraits::Scalar;
    using SolutionVector = typename MDTraits::SolutionVector;

    template<std::size_t i> using SubDomainTypeTag = typename MDTraits::template SubDomain<i>::TypeTag;
    template<std::size_t i> using Problem = GetPropType<SubDomainTypeTag<i>, Properties::Problem>;
    template<std::size_t i> using PrimaryVariables = GetPropType<SubDomainTypeTag<i>, Properties::PrimaryVariables>;
    template<std::size_t i> using NumEqVector = GetPropType<SubDomainTypeTag<i>, Properties::NumEqVector>;
    template<std::size_t i> using ElementVolumeVariables = typename GetPropType<SubDomainTypeTag<i>, Properties::GridVolumeVariables>::LocalView;
    template<std::size_t i> using VolumeVariables = typename GetPropType<SubDomainTypeTag<i>, Properties::GridVolumeVariables>::VolumeVariables;
    template<std::size_t i> using GridGeometry = typename MDTraits::template SubDomain<i>::GridGeometry;
    template<std::size_t i> using FVElementGeometry = typename GridGeometry<i>::LocalView;
    template<std::size_t i> using SubControlVolumeFace = typename GridGeometry<i>::SubControlVolumeFace;
    template<std::size_t i> using SubControlVolume = typename GridGeometry<i>::SubControlVolume;
    template<std::size_t i> using GridView = typename GridGeometry<i>::GridView;
    template<std::size_t i> using Element = typename GridView<i>::template Codim<0>::Entity;

    template<std::size_t i>
    static constexpr auto domainIdx()
    { return typename MDTraits::template SubDomain<i>::Index{}; }

    template<std::size_t i>
    static constexpr bool isCCTpfa()
    { return GridGeometry<i>::discMethod == DiscretizationMethod::cctpfa; }

    using CouplingStencil = std::vector<std::size_t>;
public:
    //! export traits
    using MultiDomainTraits = MDTraits;
    //! export stencil types
    using CouplingStencils = std::unordered_map<std::size_t, CouplingStencil>;

    /*!
     * \brief Methods to be accessed by main
     */
    // \{

    void init(std::shared_ptr<Problem<domainIdx<0>()>> problem0,
              std::shared_ptr<Problem<domainIdx<1>()>> problem1,
              const SolutionVector& curSol)
    {
        static_assert(isCCTpfa<0>() && isCCTpfa<1>(), "Only cctpfa implemented!");

        this->setSubProblems(std::make_tuple(problem0, problem1));
        this->updateSolution(curSol);
        couplingMapper_.update(*this);
    }

    // \}

    /*!
     * \brief Methods to be accessed by the assembly
     */
    // \{

    /*!
     * \brief returns an iteratable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the element residual of the given element of domain i
     *
     * \param domainI the domain index of domain i
     * \param element the coupled element of domain í
     * \param domainJ the domain index of domain j
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil& couplingStencil(Dune::index_constant<i> domainI,
                                           const Element<i>& element,
                                           Dune::index_constant<j> domainJ) const
    {
        static_assert(i != j, "A domain cannot be coupled to itself!");
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(element);
        return couplingMapper_.couplingStencil(domainI, eIdx, domainJ);
    }

    // \}

    /*!
     * \brief If the boundary entity is on a coupling boundary
     * \param domainI the domain index for which to compute the flux
     * \param scvf the sub control volume face
     */
    template<std::size_t i>
    bool isCoupled(Dune::index_constant<i> domainI,
                   const SubControlVolumeFace<i>& scvf) const
    {
        return couplingMapper_.isCoupled(domainI, scvf);
    }

    /*!
     * \brief Evaluate the boundary conditions for a coupled Neumann boundary segment.
     *
     * This is the method for the case where the Neumann condition is
     * potentially solution dependent
     *
     * \param domainI the domain index for which to compute the flux
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scvf The sub control volume face
     * \param phaseIdx the phase for which to compute the flux
     *
     * Negative values mean influx.
     * E.g. for the mass balance that would the mass flux in \f$ [ kg / (m^2 \cdot s)] \f$.
     */
    template<std::size_t i>
    Scalar advectiveFluxCoupling(Dune::index_constant<i> domainI,
                                 const Element<i>& element,
                                 const FVElementGeometry<i>& fvGeometry,
                                 const ElementVolumeVariables<i>& elemVolVars,
                                 const SubControlVolumeFace<i>& scvf,
                                 int phaseIdx = 0) const
    {
        Scalar flux = 0.0;
        static const bool enableGravity = getParamFromGroup<bool>(this->problem(domainI).paramGroup(), "Problem.EnableGravity");
        constexpr auto otherDomainIdx = domainIdx<1-i>();

        const auto& outsideElement = this->problem(otherDomainIdx).gridGeometry().element(couplingMapper_.outsideElementIndex(domainI, scvf));
        auto fvGeometryOutside = localView(this->problem(otherDomainIdx).gridGeometry());
        fvGeometryOutside.bindElement(outsideElement);

        const auto& flipScvf = fvGeometryOutside.scvf(couplingMapper_.flipScvfIndex(domainI, scvf));
        const auto& outsideScv = fvGeometryOutside.scv(flipScvf.insideScvIdx());
        const auto outsideVolVars = volVars(otherDomainIdx, outsideElement, outsideScv);

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        using Extrusion = typename GridGeometry<i>::Extrusion;
        const auto ti = computeTpfaTransmissibility(scvf, insideScv, insideVolVars.permeability(), insideVolVars.extrusionFactor());
        const auto tj = computeTpfaTransmissibility(flipScvf, outsideScv, outsideVolVars.permeability(), outsideVolVars.extrusionFactor());
        Scalar tij = 0.0;
        if (ti*tj > 0.0)
            tij = Extrusion::area(scvf)*(ti * tj)/(ti + tj);

        if (enableGravity)
        {
            // do averaging for the density over all neighboring elements
            const auto rho = (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;
            const auto& g = this->problem(domainI).spatialParams().gravity(scvf.ipGlobal());

            //! compute alpha := n^T*K*g
            const auto alpha_inside = vtmv(scvf.unitOuterNormal(), insideVolVars.permeability(), g)*insideVolVars.extrusionFactor();
            const auto alpha_outside = vtmv(scvf.unitOuterNormal(), outsideVolVars.permeability(), g)*outsideVolVars.extrusionFactor();

            // Obtain inside and outside pressures
            const auto pInside = insideVolVars.pressure(0);
            const auto pOutside = outsideVolVars.pressure(0);

            flux = tij*(pInside - pOutside) + rho*Extrusion::area(scvf)*alpha_inside - rho*tij/tj*(alpha_inside - alpha_outside);

        }
        else
        {
            // Obtain inside and outside pressures
            const auto pInside = insideVolVars.pressure(phaseIdx);
            const auto pOutside = outsideVolVars.pressure(phaseIdx);

            // return flux
            flux = tij*(pInside - pOutside);
        }

        // upwind scheme
        static const Scalar upwindWeight = getParam<Scalar>("Flux.UpwindWeight");
        auto upwindTerm = [phaseIdx](const auto& volVars){ return volVars.density(phaseIdx)*volVars.mobility(phaseIdx); };
        if (std::signbit(flux)) // if sign of flux is negative
            flux *= (upwindWeight*upwindTerm(outsideVolVars)
                     + (1.0 - upwindWeight)*upwindTerm(insideVolVars));
        else
            flux *= (upwindWeight*upwindTerm(insideVolVars)
                     + (1.0 - upwindWeight)*upwindTerm(outsideVolVars));

        // make this a area-specific flux
        flux /= Extrusion::area(scvf)*insideVolVars.extrusionFactor();
        return flux;
    }

    //! Return the volume variables of domain i for a given element and scv
    template<std::size_t i>
    VolumeVariables<i> volVars(Dune::index_constant<i> domainI,
                               const Element<i>& element,
                               const SubControlVolume<i>& scv) const
    {
        VolumeVariables<i> volVars;
        const auto elemSol = elementSolution(element, this->curSol()[domainI], this->problem(domainI).gridGeometry());
        volVars.update(elemSol, this->problem(domainI), element, scv);
        return volVars;
    }

private:
    DarcyDarcyBoundaryCouplingMapper<MDTraits> couplingMapper_;
};

} // end namespace Dumux

#endif
