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
 * \ingroup MultiDomain
 * \copydoc Dumux::IdentityIndexMap
 */
#ifndef DUMUX_MULTIDOMAIN_IDENTITY_INDEX_MAP_HH
#define DUMUX_MULTIDOMAIN_IDENTITY_INDEX_MAP_HH

namespace Dumux {

/*!
 * \file
 * \ingroup MultiDomain
 * \brief Identity index map between two sub-domains.
 * \note In some MultiDomain applications, the coupling occurs between different
 *       processes, but on the same grid. In this context, index maps can be used
 *       to map the element index of a sub-domain to the index of another the
 *       overlapping element in another sub-domain. This mapping is trivial, if
 *       the same grid manager and the same ordering is chosen. However, one could
 *       want to use a different grid manager and/or different ordering in one of
 *       the sub-domains, in which case the element indices have to be mapped.
 *       This implementation here is the trivial identity index map in case the
 *       orderings are identical.
 */
template< class SubDomainId1, class IndexType1,
          class SubDomainId2, class IndexType2 >
struct IdentityIndexMap
{
    /*!
     * \brief Maps the index of an element of sub-domain 1 to the
     *        corresponding index of the overlapping element in sub-domain 2.
     */
    static constexpr IndexType1 map(SubDomainId2 id, IndexType2 idx)
    { return idx; }

    /*!
    * \brief Maps the index of an element of sub-domain 2 to the
    *        corresponding index of the overlapping element in sub-domain 1.
     */
    static constexpr IndexType2 map(SubDomainId1 id, IndexType1 idx)
    { return idx; }
};

} //end namespace Dumux

#endif
