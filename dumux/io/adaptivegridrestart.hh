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
 * \brief Provides a restart functionality for adaptive grids
 */
#ifndef DUMUX_ADAPTIVEGRIDRESTART_HH
#define DUMUX_ADAPTIVEGRIDRESTART_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>
#if HAVE_ALUGRID
#include <dune/grid/alugrid/2d/alugrid.hh>
#include <dune/grid/alugrid/3d/alugrid.hh>
#endif
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#if HAVE_ALBERTA
#include <dune/grid/albertagrid/agrid.hh>
#endif
#include <dune/grid/utility/grapedataioformattypes.hh>

namespace Dumux
{
/*!
 * \brief Indices denoting the different grid types.
 */

//! \cond \private
template<class Grid, int dim>
struct GridRestartCheck
{
    static const bool allowRestart = false;
};

template<int dim>
struct GridRestartCheck<Dune::YaspGrid<dim>, dim>
{
    static const bool allowRestart = false;
};

template<int dim>
struct GridRestartCheck<Dune::SGrid<dim, dim>, dim>
{
    static const bool allowRestart = false;
};

#if HAVE_ALUGRID
template<int dim>
struct GridRestartCheck<Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming>, dim>
{
    static const bool allowRestart = true;
};

template<int dim>
struct GridRestartCheck<Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>, dim>
{
    static const bool allowRestart = true;
};
#endif

#if HAVE_UG
template<int dim>
struct GridRestartCheck<Dune::UGGrid<dim>, dim>
{
    static const bool allowRestart = false;
};
#endif
#if HAVE_ALBERTA
template<int dim>
struct GridRestartCheck<Dune::AlbertaGrid< dim, dim>, dim>
{
    static const bool allowRestart = true;
};
#endif
//! \endcond


/*!
 * \brief Default class for restart functionality for non-adaptive grids
 */
template <class Grid, int dim, bool allowGridRestart = GridRestartCheck<Grid, dim>::allowRestart >
class AdaptiveGridRestart
{
public:
    /*!
     * \brief Write the grid to a file.
     */
    template<class Problem>
    static void serializeGrid(Problem& problem)
    {
    };
    /*!
     * \brief Restart the grid from the file.
     */
    template<class Problem>
    static void restartGrid(Problem& problem)
    {
    };
};

/*!
 * \brief Provides a restart functionality for adaptive grids
 */
template <class Grid, int dim>
class AdaptiveGridRestart<Grid, dim, true>
{
public:
    /*!
     * \brief Write the grid to a file.
     */
    template<class Problem>
    static void serializeGrid(Problem& problem)
    {
        std::string gridName = restartGridFileName_(problem);
        double time = problem.timeManager().time();
        problem.grid().template writeGrid<Dune::xdr> (gridName, time);
    };
    /*!
     * \brief Restart the grid from the file.
     */

    template<class Problem>
    static void restartGrid(Problem& problem)
    {
        std::string gridName = restartGridFileName_(problem);
        double time = problem.timeManager().time();
        problem.grid().template readGrid<Dune::xdr> (gridName, time);
    };

private:
    //! \brief Return the restart file name.
    template<class Problem>
    static const std::string restartGridFileName_(Problem& problem)
    {
        int rank = problem.gridView().comm().rank();
        std::ostringstream oss;
        oss << problem.name()<<"_time="<<problem.timeManager().time()<<"_rank="<<rank<<".grs";
        return oss.str();
    }
};


} // namespace Dumux
#endif
