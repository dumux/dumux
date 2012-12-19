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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \brief Base class for all problems which use the two-phase DFM box model
 */
#ifndef DUMUX_BOXMODELS_2PDFM_PROBLEM_HH
#define DUMUX_BOXMODELS_2PDFM_PROBLEM_HH

#include <dumux/implicit/common/implicitporousmediaproblem.hh>
#include "2pdfmproperties.hh"

namespace Dumux
{
/*!
 * \ingroup BoxBaseProblems
 * \ingroup TwoPBoxModel
 * \brief Base class for all problems which use the two-phase box model
 */
template<class TypeTag>
class TwoPDFMProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParams) SpatialParams;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     * \param verbose Turn verbosity on or off
     */
    DUNE_DEPRECATED_MSG("use ImplicitPorousMediaProblem instead")
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
     * \param spatialParams The spatial parameters object
     * \param verbose Turn verbosity on or off
     */
    TwoPDFMProblem(TimeManager &timeManager,
                const GridView &gridView,
                SpatialParams &spatialParams,
                bool verbose = true)
        : ParentType(timeManager, gridView)
    {
        this->newSpatialParams_ = false;
    delete this->spatialParams_;
    this->spatialParams_ = &spatialParams;
    }
};
} // end namespace

#endif // DUMUX_BOXMODELS_2PDFM_PROBLEM_HH
