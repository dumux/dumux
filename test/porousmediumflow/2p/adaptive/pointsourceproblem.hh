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
 * \ingroup TwoPTests
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_LENSPROBLEM_POINTSOURCE_HH
#define DUMUX_LENSPROBLEM_POINTSOURCE_HH

#include "problem.hh"

namespace Dumux {

/*!
 * \ingroup TwoPTests
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */
template <class TypeTag >
class PointSourceTestProblem : public TwoPTestProblemAdaptive<TypeTag>
{
    using ParentType = TwoPTestProblemAdaptive<TypeTag>;
    using PointSource =  GetPropType<TypeTag, Properties::PointSource>;

public:
    //! Use parent's constructor
    using ParentType::ParentType;

    /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of PointSources that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in untis
     * \f$ [ \textnormal{unit of conserved quantity} / s ] \f$.
     * Positive values mean that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    {
        // inject 2 kg/s of nonwetting phase at position (1, 1);
        pointSources.push_back(PointSource({0.502, 3.02}, {0, 0.1}));
    }
};

} // end namespace Dumux

#endif
