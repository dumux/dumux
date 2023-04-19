// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
