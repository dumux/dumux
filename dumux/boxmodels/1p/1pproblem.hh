// $Id: 1pboxproblem.hh 3759 2010-06-21 16:59:10Z bernd $
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Base class for all problems which use the single-phase box
 *        model
 */
#ifndef DUMUX_1P_PROBLEM_HH
#define DUMUX_1P_PROBLEM_HH

#include <dumux/boxmodels/common/boxproblem.hh>

namespace Dumux
{
/*!
 * \ingroup OnePProblems
 * \brief Base class for all problems which use the single-phase box model
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class OnePBoxProblem : public BoxProblem<TypeTag>
{
    typedef BoxProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::Grid Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    // material properties
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;

    enum {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    OnePBoxProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView),
          gravity_(0)
    {
        spatialParameters_ = new SpatialParameters(gridView);

        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity)))
            gravity_[dim-1]  = -9.81;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This method MUST be overwritten by the actual problem.
     */
    Scalar temperature() const
    { return asImp_()->temperature(); };

    /*!
     * \brief Returns the acceleration due to gravity.
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
    { return *spatialParameters_; }

    /*!
     * \copydoc spatialParameters()
     */
    const SpatialParameters &spatialParameters() const
    { return *spatialParameters_; }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation *asImp_()
    { return static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation *asImp_() const
    { return static_cast<const Implementation *>(this); }

    SpatialParameters* spatialParameters_;

    GlobalPosition gravity_;
};

}

#endif
