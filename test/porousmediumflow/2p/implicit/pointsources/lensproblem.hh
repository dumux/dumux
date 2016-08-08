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
 *
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 */

#ifndef DUMUX_LENSPROBLEM_POINTSOURCE_HH
#define DUMUX_LENSPROBLEM_POINTSOURCE_HH

#include "../lensproblem.hh"

#if HAVE_UG
namespace Dumux
{
template <class TypeTag>
class LensProblemPointSource;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
NEW_TYPE_TAG(LensProblemPointSource, INHERITS_FROM(LensCCAdaptiveProblem));
NEW_PROP_TAG(BaseProblem);
SET_TYPE_PROP(LensProblemPointSource, BaseProblem, ImplicitPorousMediaProblem<TypeTag>);
SET_BOOL_PROP(LensProblemPointSource, AdaptiveGrid, true);
SET_TYPE_PROP(LensProblemPointSource, AdaptionIndicator, TwoPImplicitGridAdaptIndicator<TypeTag>);
SET_TYPE_PROP(LensProblemPointSource, AdaptionInitializationIndicator, ImplicitGridAdaptInitializationIndicator<TypeTag>);
SET_TYPE_PROP(LensProblemPointSource, Grid, Dune::UGGrid<2>);
SET_TYPE_PROP(LensProblemPointSource, Problem, LensProblemPointSource<TypeTag>);
}

/*!
 * \ingroup TwoPModel
 * \ingroup ImplicitTestProblems
 * \brief Soil contamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeablility which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability. Note that
 * this problem is discretized using only two dimensions, so from the
 * point of view of the two-phase model, the depth of the domain
 * implicitly is 1 m everywhere.
 *
 * On the top and the bottom of the domain neumann boundary conditions
 * are used, while dirichlet conditions apply on the left and right
 * boundaries.
 *
 * DNAPL is injected at the top boundary from 3m to 4m at a rate of
 * 0.04 kg/(s m^2), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * The dirichlet boundaries on the left boundary is the hydrostatic
 * pressure scaled by a factor of 1.125, while on the right side it is
 * just the hydrostatic pressure. The DNAPL saturation on both sides
 * is zero.
 *
 * This problem uses the \ref TwoPModel.
 *
 * This problem should typically be simulated until \f$t_{\text{end}}
 * \approx 20\,000\;s\f$ is reached. A good choice for the initial time step
 * size is \f$t_{\text{inital}} = 250\;s\f$.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box2p -parameterFile test_box2p.input</tt> or
 * <tt>./test_cc2p -parameterFile test_cc2p.input</tt>
 */
template <class TypeTag >
class LensProblemPointSource : public LensProblem<TypeTag>
{
    typedef LensProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PointSource) PointSource;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    LensProblemPointSource(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {}

    /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the \a values method of the point source
     * has to return the absolute mass rate in untis
     * \f$ [ \textnormal{unit of conserved quantity} / s ] \f$.
     * Positive values mean that mass is created, negative ones mean that it vanishes.
     */
    void addPointSources(std::vector<PointSource>& pointSources) const
    {
        // inject 2 kg/s of non-wetting phase at position (1, 1);
        pointSources.push_back(PointSource({0.502, 3.02}, {0, 0.1}));
    }
};

} //end namespace

#endif // HAVE_UG

#endif
