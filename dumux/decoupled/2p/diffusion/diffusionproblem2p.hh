// $Id: diffusionproblem2p.hh 3784 2010-06-24 13:43:57Z bernd $
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
#ifndef DUMUX_DIFFUSIONPROBLEM_2P_HH
#define DUMUX_DIFFUSIONPROBLEM_2P_HH

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
class DiffusionProblem2P: public OneModelProblem<TypeTag, Implementation>
{
    typedef OneModelProblem<TypeTag, Implementation> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GridView::Grid Grid;typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    // material properties
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;

    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    DiffusionProblem2P(const GridView &gridView, bool verbose = true) :
        ParentType(gridView, verbose), gravity_(0)
    {
        spatialParameters_ = new SpatialParameters(gridView);
        newSpatialParams_ = true;
        gravity_ = 0;
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity)))
            gravity_[dim - 1] = -9.81;
    }

    DiffusionProblem2P(const GridView &gridView, SpatialParameters &spatialParameters, bool verbose = true) :
        ParentType(gridView, verbose), gravity_(0), spatialParameters_(&spatialParameters)
    {
        newSpatialParams_ = false;
        gravity_ = 0;
        if (GET_PROP_VALUE(TypeTag, PTAG(EnableGravity)))
            gravity_[dim - 1] = -9.81;
    }

    virtual ~DiffusionProblem2P()
    {
        if (newSpatialParams_)
        {
            delete spatialParameters_;
        }
    }

    /*!
     * \name Problem parameters
     */
    // \{

    void timeIntegration()
    {
        // set the initial condition of the model
        ParentType::init();

        //end simulation -> no time dependent problem!
        this->timeManager().setFinished();

        return;
    }

    void serialize()
    {
        return;
    }
    void deserialize(double t)
    {
        return;
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This method MUST be overwritten by the actual problem.
     */
    Scalar temperature() const
    {
        return asImp_()->temperature();
    }
    ;

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * If the <tt>EnableGravity</tt> property is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     */
    const GlobalPosition &gravity() const
    {
        return gravity_;
    }

    /*!
     * \brief Returns the spatial parameters object.
     */
    SpatialParameters &spatialParameters()
    {
        return *spatialParameters_;
    }

    /*!
     * \copydoc spatialParameters()
     */
    const SpatialParameters &spatialParameters() const
    {
        return *spatialParameters_;
    }

    // \}

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation *asImp_()
    {
        return static_cast<Implementation *> (this);
    }

    //! \copydoc asImp_()
    const Implementation *asImp_() const
    {
        return static_cast<const Implementation *> (this);
    }

    GlobalPosition gravity_;

    // fluids and material properties
    SpatialParameters* spatialParameters_;
    bool newSpatialParams_;
};

}

#endif
