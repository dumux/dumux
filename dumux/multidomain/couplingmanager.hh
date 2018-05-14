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

    //! default type for encapsulating both the element-local and grid index of a dof
    template<std::size_t domainI, std::size_t domainJ>
    struct DofData
    {
        std::size_t index;
        unsigned int localIndex;
    };

    //! default type for storing data of all dofs within a coupled element
    template<std::size_t domainI, std::size_t domainJ>
    using CoupledElementDofData = std::vector< DofData<domainI, domainJ> >;

    /*!
     * \brief The coupling element stencil, i.e. which elements of domain J
     *        are coupled to the given element of domain I
     */
    template<class Element, std::size_t i, std::size_t j>
    const CouplingStencilType<i, j>& couplingElementStencil(const Element& element,
                                                            Dune::index_constant<i> domainI,
                                                            Dune::index_constant<j> domainJ) const
    { DUNE_THROW(Dune::NotImplemented, "Coupling manager does not implement couplingElementStencil() function"); }

    /*!
     * \brief returns data on all dofs inside an element of domain j
     *        that is coupled to an element of domain i with the given index
     */
    template<class ElementI, std::size_t i, class IndexTypeJ, std::size_t j>
    const CoupledElementDofData<i, j>& coupledElementDofData(Dune::index_constant<i> domainI,
                                                             const ElementI& elementI,
                                                             Dune::index_constant<j> domainJ,
                                                             IndexTypeJ globalJ) const
    { DUNE_THROW(Dune::NotImplemented, "Coupling manager does not implement coupledElementDofData() function"); }

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
    template<std::size_t i, class Element, class ElementSolution, class Assembler>
    void updateCouplingContext(Dune::index_constant<i> domainI, Dune::index_constant<i> domainJ,
                               const Element& element, const ElementSolution& elemSol, std::size_t pvIdx,
                               const Assembler& assembler)
    { DUNE_THROW(Dune::NotImplemented, "Coupling manager does not implement updateCouplingContext() function"); }

    /*!
     * \brief Update the local views of domain i after the coupling context has been updated.
     *        This is only necessary if some of these containers depend on quantities of the other domain
     *        stored in the coupling context here in the coupling manager.
     */
    template<class Element, class FVElementGeometry,
             class ElementVolumeVariables, class ElementFluxVariablesCache, class GridVariables>
    void updateSelf(const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const ElementVolumeVariables& elemVolVars,
                    const ElementFluxVariablesCache& elemFluxVarsCache,
                    const GridVariables& gridVariables)
    {}
};

} //end namespace Dumux

#endif
