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
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 */
#ifndef DUMUX_DISCRETIZATION_VOLUME_VARIABLES_HH
#define DUMUX_DISCRETIZATION_VOLUME_VARIABLES_HH

#include <dumux/implicit/properties.hh>

#include <dumux/common/valgrind.hh>

namespace Dumux
{

/*!
 * \ingroup ImplicitVolumeVariables
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 */
template <class TypeTag>
class ImplicitVolumeVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

public:

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param priVars A vector containing the primary variables for the control volume
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param fvGeometry The finite volume geometry for the element
     * \param scvIdx Local index of the sub control volume which is inside the element
     * \param isOldSol Specifies whether this is the previous solution or the current one
     *
     * \todo Eliminate the 'isOldSol' parameter. This implies that the
     *       'pseudo-primary variables' must be somehow be stored
     *       inside the PrimaryVariables. (e.g. we need to know the
     *       phase state in the 2p2c model)
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const SubControlVolume &scv)
    {
        Valgrind::CheckDefined(priVars);
        priVars_ = priVars;
        extrusionFactor_ = problem.boxExtrusionFactor(element, scv);
    }

    /*!
     * \brief Return the vector of primary variables
     */
    const PrimaryVariables &priVars() const
    { return priVars_; }

    /*!
     * \brief Return a component of primary variable vector
     *
     * \param pvIdx The index of the primary variable of interest
     */
    Scalar priVar(const int pvIdx) const
    {
        return priVars_[pvIdx];
    }

    /*!
     * \brief Return how much the sub-control volume is extruded.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    /*!
     * \brief If running in valgrind this makes sure that all
     *        quantities in the volume variables are defined.
     */
    void checkDefined() const
    {
#if !defined NDEBUG && HAVE_VALGRIND
        Valgrind::CheckDefined(priVars_);
#endif
    }

protected:
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    PrimaryVariables priVars_;
    Scalar extrusionFactor_;
};

} // end namespace

#endif
