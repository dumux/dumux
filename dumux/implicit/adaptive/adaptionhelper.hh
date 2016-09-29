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
#ifndef DUMUX_IMPLICIT_ADAPTIONHELPER_HH
#define DUMUX_IMPLICIT_ADAPTIONHELPER_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/propertysystem.hh>

/**
 * @file
 * @brief  Base class holding the variables for implicit models.
 */

namespace Dumux
{
namespace Properties
{
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(ImplicitIsBox);
NEW_PROP_TAG(Problem);
}

/*!
 * \brief Base class holding the variables for implicit models.
 */
template<class TypeTag>
class ImplicitAdaptionHelper
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::Grid Grid;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Traits::template Codim<dofCodim>::Entity DofEntity;

public:
    //! Constructs an adaptive helper object
    /**
     * In addition to providing a storage object for cell-centered Methods, this class provides
     * mapping functionality to adapt the grid.
     *
     *  @param gridView a DUNE gridview object corresponding to diffusion and transport equation
     */
    ImplicitAdaptionHelper(Problem& problem)
    {}

    /*!
     * Store primary variables
     *
     * To reconstruct the solution in father elements, problem properties might
     * need to be accessed.
     * From upper level on downwards, the old solution is stored into an container
     * object, before the grid is adapted. Father elements hold averaged information
     * from the son cells for the case of the sons being coarsened.
     *
     * @param problem The current problem
     */
    void storePrimVars(Problem& problem)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a storePrimVars method.");
    }


    /*!
     * Reconstruct missing primary variables (where elements are created/deleted)
     *
     * To reconstruct the solution in father elements, problem properties might
     * need to be accessed.
     * Starting from the lowest level, the old solution is mapped on the new grid:
     * Where coarsened, new cells get information from old father element.
     * Where refined, a new solution is reconstructed from the old father cell,
     * and then a new son is created. That is then stored into the general data
     * structure (CellData).
     *
     * @param problem The current problem
     */
    void reconstructPrimVars(Problem& problem)
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide "
                   "a reconstructPrimVars method.");
    }

protected:

    int dofIndex(const Problem& problem, const DofEntity& entity) const
    {
                return problem.model().dofMapper().index(entity);
    }

    int dofIndex(const Problem& problem, const Element& element, const int scvIdx) const
    {
                return problem.model().dofMapper().subIndex(element, scvIdx, dofCodim);
    }

    int elementIndex(const Problem& problem, const Element& element) const
    {
                return problem.elementMapper().index(element);
    }

};
}
#endif
