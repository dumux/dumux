// $Id$
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
 * \brief Base class for all problems which use the box scheme
 */
#ifndef DUMUX_RICHARDS_PROBLEM_HH
#define DUMUX_RICHARDS_PROBLEM_HH

#include <dumux/boxmodels/common/boxproblem.hh>

namespace Dumux
{
/*!
 * \ingroup RichardsProblems
 * \brief Base class for all problems which use the two-phase box model
 *
 * For a description of the Richards model, see Dumux::RichardsModel
 */
template<class TypeTag>
class RichardsBoxProblem : public BoxProblem<TypeTag>
{
    typedef BoxProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Implementation;
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

    typedef typename GridView::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor.
     *
     * The overloaded class must allocate all data structures
     * required, but _must not_ do any calls to the model, the
     * jacobian assembler, etc inside the constructor.
     *
     * If the problem requires information from these, the
     * BoxProblem::init() method be overloaded.
     *
     * \param timeManager The TimeManager which keeps track of time
     * \param gridView The GridView used by the problem.
     */
    RichardsBoxProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView),
          gravity_(0), spatialParams_(gridView)
    {
        gravity_ = 0;
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity)))
            gravity_[dim-1]  = -9.81;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ within a control volume.
     *
     * This method MUST be overwritten by the actual problem.
     *
     * \param element The DUNE Codim<0> enitiy which intersects with
     *                the finite volume.
     * \param fvGeom The finite volume geometry of the element.
     * \param scvIdx The local index of the sub control volume inside the element
     */
    Scalar temperature(const Element &element,
                       const FVElementGeometry fvGeom,
                       int scvIdx) const
    { DUNE_THROW(Dune::NotImplemented, "temperature() method not implemented by the actual problem"); };

    /*!
     * \brief Returns the reference pressure \f$\mathrm{[Pa]}\f$ of the non-wetting
     *        phase within a control volume.
     *
     * This method MUST be overwritten by the actual problem.
     *
          * \param element The DUNE Codim<0> enitiy which intersects with
     *                the finite volume.
     * \param fvGeom The finite volume geometry of the element.
     * \param scvIdx The local index of the sub control volume inside the element
     */
    Scalar referencePressure(const Element &element,
                             const FVElementGeometry fvGeom,
                             int scvIdx) const
    { DUNE_THROW(Dune::NotImplemented, "referencePressure() method not implemented by the actual problem"); };

    /*!
     * \brief Returns the acceleration due to gravity \f$\mathrm{[m/s^2]}\f$.
     *
     * If the <tt>EnableGravity</tt> property is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     */
    const GlobalPosition &gravity() const
    { return gravity_; }

    /*!
     * \brief Returns the spatial parameters object.
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
    // the gravity vector
    GlobalPosition gravity_;

    // material properties
    SpatialParameters spatialParams_;
};

}

#endif
