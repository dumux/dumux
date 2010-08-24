// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
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
 * \brief Base class for all problems which use the box scheme
 */
#ifndef DUMUX_TRANSPORTPROBLEM_2P_HH
#define DUMUX_TRANSPORTPROBLEM_2P_HH

#include <dumux/decoupled/common/onemodelproblem.hh>
#include <dumux/decoupled/2p/variableclass2p.hh>
#include <dumux/material/fluidsystems/2p_system.hh>
#include <dumux/decoupled/2p/2pproperties.hh>


namespace Dumux
{
/*!
 * \ingroup Decoupled
 * \brief  Base class for all 2-phase problems which use an impes algorithm
 *
 * \todo Please doc me more!
 */
template<class TypeTag, class Implementation>
class TransportProblem2P : public OneModelProblem<TypeTag, Implementation>
{
    typedef OneModelProblem<TypeTag, Implementation> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::Grid Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    // material properties
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::ScalarSolution Solution;


    enum {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld>      GlobalPosition;

public:
    TransportProblem2P(const GridView &gridView)
        : ParentType(gridView),
        gravity_(0),spatialParameters_(gridView)
    {
        gravity_ = 0;
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity)))
            gravity_[dim - 1] = - 9.81;
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
    { return this->asImp_()->temperature(); };

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
    { return spatialParameters_; }

    /*!
     * \copydoc spatialParameters()
     */
    const SpatialParameters &spatialParameters() const
    { return spatialParameters_; }

    void timeIntegration()
    {
        // allocate temporary vectors for the updates
        Solution k1 = asImp_().variables().saturation();

        dt_ = 1e100;
        Scalar t = timeManager().time();

        // obtain the first update and the time step size
        model().update(t, dt_, k1);

        //make sure t_old + dt is not larger than tend
        dt_ = std::min(dt_*cFLFactor_, timeManager().episodeMaxTimeStepSize());
        timeManager().setTimeStepSize(dt_);

        // explicit Euler: Sat <- Sat + dt*N(Sat)
        asImp_().variables().saturation() += (k1 *= dt_);
    }

    // \}

private:
    GlobalPosition gravity_;

    // fluids and material properties
    SpatialParameters spatialParameters_;

    static const Scalar cFLFactor_= GET_PROP_VALUE(TypeTag, PTAG(CFLFactor));
};

}

#endif
