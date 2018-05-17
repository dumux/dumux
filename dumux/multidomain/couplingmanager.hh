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
 * \brief The interface of the coupling manager for multi domain problems
 */

#ifndef DUMUX_MULTIDOMAIN_COUPLING_MANAGER_HH
#define DUMUX_MULTIDOMAIN_COUPLING_MANAGER_HH

#include <dune/common/exceptions.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup MultiDomain
 * \brief The interface of the coupling manager for multi domain problems
 */
template<class Traits>
class CouplingManager
{
public:
    //! default type used for coupling element stencils
    template<std::size_t domainI, std::size_t domainJ>
    using CouplingStencilType = std::vector< std::size_t >;

    /*!
     * \brief The coupling element stencil, i.e. which elements of domain J
     *        are coupled to the given element of domain I
     */
    template<class Element, std::size_t i, std::size_t j>
    const CouplingStencilType<i, j>& couplingStencil(const Element& element,
                                                     Dune::index_constant<i> domainI,
                                                     Dune::index_constant<j> domainJ) const
    { DUNE_THROW(Dune::NotImplemented, "Coupling manager does not implement couplingElementStencil() function"); }

    //! Return an empty list of additional dof dependencies per default
    template<std::size_t id>
    std::vector<std::vector<std::size_t>> additionalDofDependencies(Dune::index_constant<id>)
    { return std::vector<std::vector<std::size_t>>(); }

    //! Return an empty list of additional dof dependencies per default
    template<std::size_t id, class Element>
    std::vector<std::size_t> getAdditionalDofDependencies(Dune::index_constant<id>, const Element& element) const
    { return std::vector<std::size_t>(); }

    //! Return an empty list of additional dof dependencies per default
    template<std::size_t id, class Element>
    std::vector<std::size_t> getAdditionalDofDependenciesInverse(Dune::index_constant<id>, const Element& element) const
    { return std::vector<std::size_t>(); }

    /*!
     * \brief Prepares all data on the other domains necessary for the assembly of an element of domain i
     */
    template<class Element, std::size_t i, class Assembler>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element& element, const Assembler& assembler)
    { DUNE_THROW(Dune::NotImplemented, "Coupling manager does not implement bindCouplingContext() function"); }

    /*!
     * \brief Update the context for a derivative i->j
     */
    template<std::size_t i, class IndexTypeJ, class PrimaryVariablesJ, class Assembler>
    void updateCouplingContext(Dune::index_constant<i> domainI, Dune::index_constant<i> domainJ,
                               IndexTypeJ globalJ, const PrimaryVariablesJ& priVarsJ, unsigned int pvIdxJ,
                               const Assembler& assembler)
    { DUNE_THROW(Dune::NotImplemented, "Coupling manager does not implement updateCouplingContext() function"); }

    /*!
     * \brief Update the local views of domain i after the coupling context has been updated.
     *        This is only necessary if some of these containers depend on quantities of the other domain
     *        stored in the coupling context here in the coupling manager.
     */
    template<std::size_t i, class Element, class FVElementGeometry,
             class VolumeVariablesContainer, class FluxVariablesCacheContainer>
    void updateSelf(Dune::index_constant<i> domainI,
                    const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const VolumeVariablesContainer& elemVolVars,
                    const FluxVariablesCacheContainer& elemFluxVarsCache)
    {  DUNE_THROW(Dune::NotImplemented, "Coupling manager does not implement updateSelf() function"); }
};

} //end namespace Dumux

#endif
