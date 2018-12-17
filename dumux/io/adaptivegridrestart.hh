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
 * \ingroup InputOutput
 * \brief Provides a restart functionality for adaptive grids
 */
#ifndef DUMUX_ADAPTIVEGRIDRESTART_HH
#define DUMUX_ADAPTIVEGRIDRESTART_HH


#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dune/grid/common/backuprestore.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Indices denoting the different grid types.
 */

//! \cond \private
template<class Grid>
struct GridRestartCheck
{
    static const bool allowRestart = false;
};

// the specializations for grid managers that support restart
#if HAVE_DUNE_ALUGRID
template<int dim, int dimworld, Dune::ALUGridElementType elType, Dune::ALUGridRefinementType refinementType>
struct GridRestartCheck<Dune::ALUGrid<dim, dimworld, elType, refinementType> >
{
    static const bool allowRestart = true;
};
#endif
//! \endcond


/*!
 * \brief Default class for restart functionality for non-adaptive grids
 */
template <class Grid, bool allowGridRestart = GridRestartCheck<Grid>::allowRestart >
class AdaptiveGridRestart
{
public:
    /*!
     * \brief Write the grid to a file.
     */
    template<class Problem>
    static void serializeGrid(Problem& problem)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Adaptive restart functionality currently only works for dune-ALUGrid.");
    }

    /*!
     * \brief Restart the grid from the file.
     */
    template<class Problem>
    static void restartGrid(Problem& problem)
    {
        DUNE_THROW(Dune::NotImplemented,
                   "Adaptive restart functionality currently only works for dune-ALUGrid.");
    }
};

/*!
 * \brief Provides a restart functionality for adaptive grids
 */
template <class Grid>
class AdaptiveGridRestart<Grid, true>
{
public:
    /*!
     * \brief Write the grid to a file.
     */
    template<class Problem>
    static void serializeGrid(Problem& problem)
    {
        std::string gridName = restartGridFileName_(problem);
#if HAVE_DUNE_ALUGRID
        Dune::BackupRestoreFacility<Grid>::backup(problem.grid(), gridName);
#else
        problem.grid().template writeGrid(gridName, problem.timeManager().time()
                                                    + problem.timeManager().timeStepSize());
#endif
    }

    /*!
     * \brief Restart the grid from the file.
     */
    template<class Problem>
    static void restartGrid(Problem& problem)
    {}

private:
    //! \brief Return the restart file name.
    template<class Problem>
    static const std::string restartGridFileName_(Problem& problem)
    {
        int rank = problem.gridView().comm().rank();
        std::ostringstream oss;
        try {
            std::string name = getParam<std::string>("Problem.Name");
            oss << name;
        }
        catch (ParameterException &e)
        {
            std::cerr << e.what() << std::endl;
            std::cerr << "Taking name from problem.name(): " << problem.name() << std::endl;
            std::cerr << "Be sure to provide a parameter Problem.Name if you want to restart." << std::endl;
            oss << problem.name();
        }
        oss << "_time=" << problem.timeManager().time() + problem.timeManager().timeStepSize()
            << "_rank=" << rank << ".grs";
        return oss.str();
    }
};


} // namespace Dumux
#endif
