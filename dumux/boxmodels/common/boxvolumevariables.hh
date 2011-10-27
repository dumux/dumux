/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
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
 *
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 */
#ifndef DUMUX_BOX_VOLUME_VARIABLES_HH
#define DUMUX_BOX_VOLUME_VARIABLES_HH

#include "boxproperties.hh"

#include <dumux/common/valgrind.hh>

namespace Dumux
{

/*!
 * \ingroup BoxModel
 * \ingroup BoxVolumeVariables
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 */
template <class TypeTag>
class BoxVolumeVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VolumeVariables)) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;

public:
    // default constructor
    BoxVolumeVariables()
    { evalPoint_ = 0; };

    // copy constructor
    BoxVolumeVariables(const BoxVolumeVariables &v)
    {
        evalPoint_ = 0;
        primaryVars_ = v.primaryVars_;
        extrusionFactor_ = v.extrusionFactor_;
    };

    /*!
     * \brief Assignment operator
     */
    BoxVolumeVariables &operator=(const BoxVolumeVariables &v)
    {
        evalPoint_ = 0;
        primaryVars_ = v.primaryVars_;
        extrusionFactor_ = v.extrusionFactor_;

        return *this;
    };

    /*!
     * \brief Sets the evaluation point used by the local jacobian.
     *
     * The evaluation point is only used by semi-smooth models.
     */
    void setEvalPoint(const Implementation *ep)
    {
        evalPoint_ = ep;
        Valgrind::CheckDefined(evalPoint_);
    }

    /*!
     * \brief Returns the evaluation point used by the local jacobian.
     *
     * The evaluation point is only used by semi-smooth models.
     */
    const Implementation &evalPoint() const
    { return (evalPoint_ == 0)?asImp_():*evalPoint_; }

    /*!
     * \brief Set the volume variables which should be used as initial
     *        conditions for complex calculations.
     */
    void setHint(const Implementation *hint)
    {};

    /*!
     * \brief Update all quantities for a given control volume.
     *
     * \param priVars The primary variables for the control volume
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param elemGeom The finite volume geometry for the element
     * \param scvIdx Local index of the sub control volume which is inside the element
     * \param isOldSol Specifies whether this is the previous solution or the current onw
     *
     * \todo Eliminate the 'isOldSol' parameter. This implies that the
     *       'pseudo-primary variables' must be somehow be stored
     *       inside the PrimaryVariables. (e.g. we need to know the
     *       phase state in the 2p2c model)
     */
    void update(const PrimaryVariables &priVars,
                const Problem &problem,
                const Element &element,
                const FVElementGeometry &elemGeom,
                int scvIdx,
                bool isOldSol)
    {
        Valgrind::CheckDefined(priVars);
        primaryVars_ = priVars;
        extrusionFactor_ = problem.boxExtrusionFactor(element, elemGeom, scvIdx);
    }

    /*!
     * \brief Return the vector of primary variables
     */
    const PrimaryVariables &primaryVars() const
    { return primaryVars_; }

    /*!
     * \brief Return a component of primary variable vector
     *
     * \param pvIdx The index of the primary variable of interest
     */
    Scalar primaryVar(int pvIdx) const
    {
        return primaryVars_[pvIdx];
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
        Valgrind::CheckDefined(primaryVars_);
        Valgrind::CheckDefined(evalPoint_);
        if (evalPoint_ && evalPoint_ != this)
            evalPoint_->checkDefined();
#endif
    };

protected:
    const Implementation &asImp_() const
    { return *static_cast<const Implementation*>(this); }
    Implementation &asImp_()
    { return *static_cast<Implementation*>(this); }

    // the evaluation point of the local jacobian
    const Implementation *evalPoint_;

    PrimaryVariables primaryVars_;
    Scalar extrusionFactor_;
};

} // end namepace

#endif
