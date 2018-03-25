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

namespace Dumux {

/*!
 * \file
 * \ingroup MultiDomain
 * \brief The interface of the coupling manager for multi domain problems
 */
template<class Traits, class Imp = void>
class CouplingManager
{

public:
    using StencilType = std::vector<std::size_t>;
    using Implementation = Imp;

    /*!
     * \brief The coupling stencil, i.e. which dofs of domain i are in the stencil of the element
     */
    template<class Element, std::size_t i, std::size_t j>
    const StencilType& couplingStencil(const Element& element,
                                       Dune::index_constant<i> domainI,
                                       Dune::index_constant<j> domainJ) const
    {
        return asImp_().couplingStencil(element, domainI, domainJ);
    }

    // Return an empty list of additional dof dependencies per default
    template<std::size_t id>
    std::vector<std::vector<std::size_t>> additionalDofDependicies(Dune::index_constant<id>) const
    {
        return std::vector<std::vector<std::size_t>>();
    }

    /*!
     * \brief The coupling stencil, i.e. which dofs of domain i are in the stencil of the element
     */
    template<class Element, std::size_t i, class Assembler>
    void bindCouplingContext(Dune::index_constant<i> domainI, const Element& element, const Assembler& assembler)
    {
        return asImp_().bindCouplingContext(element, domainI, assembler);
    }

    /*!
     * \brief Update the context for a derivative i->j
     */
    template<std::size_t i, class Element, class PrimaryVariables, class Assembler>
    void updateCouplingContext(Dune::index_constant<i> domainI, Dune::index_constant<i> domainJ,
                               const Element& element, const PrimaryVariables& priVars, const Assembler& assembler)
    {
        return asImp_().updateCouplingContext(domainI, domainJ, element, priVars, assembler);
    }

    /*!
     * \brief The coupling stencil, i.e. which dofs of domain i are in the stencil of the element
     */
    template<class Element, class FVElementGeometry,
             class ElementVolumeVariables, class ElementFluxVariablesCache, class GridVariables>
    void updateSelf(const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const ElementVolumeVariables& elemVolVars,
                    const ElementFluxVariablesCache& elemFluxVarsCache,
                    const GridVariables& gridVariables)
    {
        // return asImp_().updateSelf(element, priVars);
    }

private:
    Implementation &asImp_()
    { return *static_cast<Imp*>(this); }

    const Implementation &asImp_() const
    { return *static_cast<const Imp*>(this); }
};

} //end namespace Dumux

#endif
