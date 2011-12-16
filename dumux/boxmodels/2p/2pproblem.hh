// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \brief Base class for all problems which use the two-phase box model
 */
#ifndef DUMUX_2P_PROBLEM_HH
#define DUMUX_2P_PROBLEM_HH

#include "2pproperties.hh"

#include <dumux/boxmodels/common/boxproblem.hh>

namespace Dumux
{
/*!
 * \ingroup BoxBaseProblems
 * \ingroup TwoPBoxModel
 * \brief Base class for all problems which use the two-phase box model
 */
template<class TypeTag>
class TwoPProblem : public BoxProblem<TypeTag>
{
    typedef BoxProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    // material properties
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> Vector;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     * \param verbose Turn verbosity on or off
     */
    TwoPProblem(TimeManager &timeManager,
                const GridView &gridView,
                bool verbose = true)
        : ParentType(timeManager, gridView),
          gravity_(0)
    {
        newSpatialParams_ = true;
        spatialParameters_ = new SpatialParameters(gridView);

        if (GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dim-1]  = -9.81;
    }

    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     * \param spatialParameters The spatial parameters object
     * \param verbose Turn verbosity on or off
     */
    TwoPProblem(TimeManager &timeManager,
                const GridView &gridView,
                SpatialParameters &spatialParameters,
                bool verbose = true)
        : ParentType(timeManager, gridView), spatialParameters_(&spatialParameters),
        gravity_(0)
    {
        newSpatialParams_ = false;

        if (GET_PARAM(TypeTag, bool, EnableGravity))
            gravity_[dim-1]  = -9.81;
    }

    ~TwoPProblem()
    {
        if (newSpatialParams_)
            delete spatialParameters_;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ within a control volume.
     *
     * This is the discretization specific interface for the box
     * method. By default it just calls temperature(pos).
     *
     * \param element The DUNE Codim<0> enitiy which intersects with
     *                the finite volume.
     * \param fvGeom The finite volume geometry of the element.
     * \param scvIdx The local index of the sub control volume inside the element
     */
    Scalar boxTemperature(const Element &element,
                          const FVElementGeometry fvGeom,
                          int scvIdx) const
    { return asImp_().temperatureAtPos(fvGeom.subContVol[scvIdx].global); }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at a given global position.
     *
     * This is not specific to the discretization. By default it just
     * calls temperature().
     *
     * \param pos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos(const GlobalPosition &pos) const
    { return asImp_().temperature(); }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    { DUNE_THROW(Dune::NotImplemented, "temperature() method not implemented by the actual problem"); };

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This is the box discretization specific interface. By default
     * it just calls gravityAtPos().
     */
    const Vector &boxGravity(const Element &element,
                                     const FVElementGeometry &fvGeom,
                                     int scvIdx) const
    { return asImp_().gravityAtPos(fvGeom.subContVol[scvIdx].global); }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * This is discretization independent interface. By default it
     * just calls gravity().
     */
    const Vector &gravityAtPos(const GlobalPosition &pos) const
    { return asImp_().gravity(); }

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
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
     * \brief Returns the spatial parameters object.
     */
    SpatialParameters &spatialParameters()
    { return *spatialParameters_; }

    /*!
     * \brief Returns the spatial parameters object.
     */
    const SpatialParameters &spatialParameters() const
    { return *spatialParameters_; }

    // \}

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }
    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    Vector gravity_;

    // fluids and material properties
    SpatialParameters*  spatialParameters_;
    bool newSpatialParams_;
};

}

#endif
