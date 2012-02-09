// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 *
 * \brief Base class for all problems which use the single-phase,
 *        two-component box model.
 */
#ifndef DUMUX_1P2C_PROBLEM_HH
#define DUMUX_1P2C_PROBLEM_HH

#include <dumux/boxmodels/common/boxproblem.hh>
#include "1p2cproperties.hh"

namespace Dumux
{
/*!
 * \ingroup OnePTwoCBoxModel
 * \ingroup BoxBaseProblems
 * \brief Base class for all problems which use the single-phase, two-component box model.
 *
 */
template<class TypeTag>
class OnePTwoCBoxProblem : public BoxProblem<TypeTag>
{
    typedef BoxProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dim> Vector;

public:
    /*!
     * \brief The constructor.
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    OnePTwoCBoxProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView),
          gravity_(0),
          spatialParams_(gridView)
    {
        if (GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dim-1]  = -9.81;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Return the temperature \f$\mathrm{[K]}\f$ within a control volume.
     *
     * This is the discretization specific interface for the box
     * method. By default it just calls temperature(pos).
     *
     * \param element The DUNE Codim<0> entity which intersects with
     *                the finite volume
     * \param fvGeom The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume inside the element
     */
    Scalar boxTemperature(const Element &element,
                          const FVElementGeometry fvGeom,
                          int scvIdx) const
    { return asImp_().temperatureAtPos(fvGeom.subContVol[scvIdx].global); }

    /*!
     * \brief Return the temperature \f$\mathrm{[K]}\f$ at a given global position.
     *
     * This is not specific to the discretization. By default it just
     * calls temperature().
     *
     * \param pos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos(const GlobalPosition &pos) const
    { return asImp_().temperature(); }

    /*!
     * \brief Return the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    { DUNE_THROW(Dune::NotImplemented, "temperature() method not implemented by the actual problem"); };

    /*!
     * \brief Return the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *.
     * This is the box discretization specific interface. By default
     * it just calls gravityAtPos().
     *
     * \param element The DUNE Codim<0> entity which intersects with
     *                the finite volume
     * \param fvGeom The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume inside the element
     */
    const Vector &boxGravity(const Element &element,
                                     const FVElementGeometry &fvGeom,
                                     int scvIdx) const
    { return asImp_().gravityAtPos(fvGeom.subContVol[scvIdx].global); }

    /*!
     * \brief Return the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This is the discretization independent interface. By default it
     * just calls gravity().
     *
     * \param pos Coordinate vector of the global position
     */
    const Vector &gravityAtPos(const GlobalPosition &pos) const
    { return asImp_().gravity(); }

    /*!
     * \brief Return the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This method is used for problems where the gravitational
     * acceleration does not depend on the spatial position. The
     * default behaviour is that if the <tt>EnableGravity</tt>
     * property is true, \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$ holds,
     * else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$.
     */
    const Vector &gravity() const
    { return gravity_; }

    /*!
     * \brief Return the spatial parameters object.
     */
    SpatialParameters &spatialParameters()
    { return spatialParams_; }

    /*!
     * \copydoc spatialParameters()
     */
    const SpatialParameters &spatialParameters() const
    { return spatialParams_; }

    // \}

private:
    //! Return the implementation of the problem (i.e. static polymorphism).
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }
    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    Vector gravity_;

    // spatial parameters
    SpatialParameters spatialParams_;
};

}

#endif
