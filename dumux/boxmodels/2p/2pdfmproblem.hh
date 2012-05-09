// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *                 2012 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Base class for all problems which use the two-phase box model
 */
#ifndef DUMUX_2PDFM_PROBLEM_HH
#define DUMUX_2PDFM_PROBLEM_HH

#include "2pdfmproperties.hh"

#include <dumux/common/deprecated.hh>
#include <dumux/boxmodels/common/porousmediaboxproblem.hh>

namespace Dumux
{
/*!
 * \ingroup BoxBaseProblems
 * \ingroup TwoPBoxModel
 * \brief Base class for all problems which use the two-phase box model
 */
template<class TypeTag>
class TwoPDFMProblem : public PorousMediaBoxProblem<TypeTag>
{
    typedef PorousMediaBoxProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     * \param verbose Turn verbosity on or off
     */
    DUMUX_DEPRECATED_MSG("use PorousMediaBoxProblem instead")
    TwoPDFMProblem(TimeManager &timeManager,
                const GridView &gridView,
                bool verbose = true)
        : ParentType(timeManager, gridView)
    {}

    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     * \param spatialParameters The spatial parameters object
     * \param verbose Turn verbosity on or off
     */
    TwoPDFMProblem(TimeManager &timeManager,
                const GridView &gridView,
                SpatialParameters &spatialParameters,
                bool verbose = true)
        : ParentType(timeManager, gridView)
    {
        this->newSpatialParams_ = false;
	delete this->spatialParameters_;
	this->spatialParameters_ = &spatialParameters;
    }
};

}

#endif
