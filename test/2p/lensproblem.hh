/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
#ifndef DUNE_LENSPROBLEM_HH
#define DUNE_LENSPROBLEM_HH

#define USE_UG 1

#ifdef USE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#else
#include <dune/grid/yaspgrid.hh>
#endif

#include <dumux/material/matrixproperties.hh>
#include <dumux/material/twophaserelations.hh>
#include <dumux/material/phaseproperties/phaseproperties2p.hh>


#include <dumux/auxiliary/timemanager.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/io/restart.hh>

#include <dumux/boxmodels/2p/2pboxmodel.hh>

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>

#include "lenssoil.hh"

namespace Dune
{

template <class TypeTag>
class LensProblem;

//////////
// Specify the properties for the lens problem
//////////
namespace Properties
{
NEW_TYPE_TAG(LensProblem, INHERITS_FROM(BoxTwoP));

// Set the grid type
SET_PROP(LensProblem, Grid)
{
#if USE_UG
    typedef Dune::UGGrid<2> type;
#else // USE_UG
    typedef Dune::YaspGrid<2> type;
#endif
};

// Set the problem property
SET_PROP(LensProblem, Problem)
{
    typedef Dune::LensProblem<TTAG(LensProblem)> type;
};
}

/*!
 * \ingroup TwoPBoxProblems
 * \brief Soil decontamination problem where DNAPL infiltrates a fully
 *        water saturated medium.
 *
 * The domain is sized 6m times 4m and features a rectangular lens
 * with low permeablility which spans from (1 m , 2 m) to (4 m, 3 m)
 * and is surrounded by a medium with higher permability.
 * 
 * On the top and the bottom of the domain neumann boundary conditions
 * are used, while dirichlet conditions apply on the left and right
 * boundaries.
 *
 * DNAPL is injected at the top boundary from 3m to 4m at a rate of
 * 0.04 kg/(s m), the remaining neumann boundaries are no-flow
 * boundaries.
 *
 * The dirichlet boundaries on the left boundary is the hydrostatic
 * pressure scaled by a factor of 1.125, while on the right side it is
 * just the hydrostatic pressure. The DNAPL saturation on both sides
 * is zero.
 *
 * This problem uses the \ref TwoPBoxModel.
 *
 * This problem should typically simulated until \f$t_{\text{end}} =
 * 50\,000\;s\f$ is reached. A good choice for the initial time step size
 * is \f$t_{\text{inital}} = 1\,000\;s\f$
 */
template <class TypeTag = TTAG(LensProblem) >
class LensProblem : public BasicDomain<typename GET_PROP_TYPE(TypeTag, PTAG(Grid)),
                                          typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) >
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar))     Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView))   GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model))      Model;
    typedef typename GridView::Grid                           Grid;

    typedef BasicDomain<Grid, Scalar>    ParentType;
    typedef LensProblem<TypeTag>      ThisType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;
    enum {
        numEq       = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        pressureIdx   = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        pWIdx = pressureIdx,
        sNIdx = saturationIdx,

        // Grid and world dimension
        dim         = GridView::dimension,
        dimWorld    = GridView::dimensionworld,
    };

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::PrimaryVarVector        PrimaryVarVector;
    typedef typename SolutionTypes::BoundaryTypeVector      BoundaryTypeVector;

    typedef typename GridView::template Codim<0>::Entity        Element;
    typedef typename GridView::template Codim<dim>::Entity      Vertex;
    typedef typename GridView::IntersectionIterator             IntersectionIterator;
  
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef Dune::FieldVector<Scalar, dim>       LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld>  GlobalPosition;

    enum Episode {}; // the type of an episode of the simulation
    typedef Dune::TimeManager<Episode>      TimeManager;
    typedef Dune::VtkMultiWriter<GridView>  VtkMultiWriter;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod))      NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController))  NewtonController;

    // material properties
    typedef Water                                  WettingPhase;
    typedef DNAPL                                  NonwettingPhase;
    typedef Dune::LensSoil<Grid, Scalar>           Soil;
    typedef Dune::TwoPhaseRelations<Grid, Scalar>  MaterialLaw;

public:
    LensProblem(Grid *grid,
                   const GlobalPosition &outerLowerLeft,
                   const GlobalPosition &outerUpperRight,
                   const GlobalPosition &innerLowerLeft,
                   const GlobalPosition &innerUpperRight,
                   Scalar dtInitial,
                   Scalar tEnd)
        : ParentType(grid),

          outerLowerLeft_(outerLowerLeft),
          outerUpperRight_(outerUpperRight),

          soil_(outerLowerLeft, outerUpperRight, innerLowerLeft, innerUpperRight),

          materialLaw_(soil_, wPhase_, nPhase_),
          timeManager_(tEnd,
                       this->grid().comm().rank() == 0),
          model_(*this),
          newtonMethod_(model_),
          resultWriter_("lens")
    {
        timeManager_.setStepSize(dtInitial);

        gravity_ = 0;
        gravity_[dim - 1] = -9.81;
        
        wasRestarted_ = false;
    }


    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \brief Start the simulation procedure. 
     *
     * This method is usually called by the main() function and simply
     * uses \ref Dune::TimeManager::runSimulation() to do the actual
     * work.
     */
    bool simulate()
    {
        timeManager_.runSimulation(*this);
        return true;
    };


    /*!
     * \brief Called by the \ref Dune::TimeManager in order to
     *        initialize the problem.
     */
    void init()
    {
        // set the initial condition of the model
        model_.initial();

        if (!wasRestarted_) {
            // write the inital solution to disk
            writeCurrentResult_();
        }
    }

    /*!
     * \brief Called by \ref Dune::TimeManager in order to do a time
     *        integration on the model.
     *
     * \note \a timeStepSize and \a nextStepSize are references and may
     *       be modified by the timeIntegration(). On exit of this
     *       function \a timeStepSize must contain the step size
     *       actually used by the time integration for the current
     *       steo, and \a nextStepSize must contain a suggestion for the 
     *       next time step size.
     */
    void timeIntegration(Scalar &stepSize, Scalar &nextStepSize)
    {
        model_.update(stepSize,
                      nextStepSize,
                      newtonMethod_,
                      newtonCtl_);
    }

    /*!
     * \brief Called by \ref Dune::TimeManager whenever a solution for a
     *        timestep has been computed.
     *
     * This is used to do some janitorial tasks like writing the
     * current solution to disk.
     */
    void timestepDone()
    {
        if (this->grid().comm().rank() == 0)
            std::cout << "Writing result file for current time step\n";

        // write the current result to disk
        writeCurrentResult_();

        // write restart file after every five steps
        static int dummy = 0;
        ++dummy;
        if (dummy % 5 == 0)
            serialize();
    };

    /*!
     * \brief Returns the current time step size [seconds].
     */
    Scalar timeStepSize() const
    { return timeManager_.stepSize(); }

    /*!
     * \brief Sets the current time step size [seconds].
     */
    void setTimeStepSize(Scalar dt)
    { return timeManager_.setStepSize(dt); }

    // \}

    /*!
     * \name Problem parameters
     */
    // \{

    /*! 
     * \brief Returns numerical model used for the problem.
     *
     * The lens problem uses \ref Dune::TwoPBoxModel .
     */
    Model &model()
    {
        return model_;
    }

    /*! 
     * \copydoc model()
     */
    const Model &model() const
    {
        return model_;
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    {
        return 283.15; // -> 10°C
    };

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * For this problem, this means \f$\boldsymbol{g} = ( 0,\ -9.81)^T \f$
     */
    const GlobalPosition &gravity () const
    {
        return gravity_;
    }

    /*! 
     * \brief Fluid properties of the wetting phase.
     *
     * For the lens problem, the wetting phase is \ref Dune::Water .
     */
    const WettingPhase &wettingPhase() const
    { return wPhase_; }

    /*! 
     * \brief Fluid properties of the non-wetting phase.
     *
     * For the lens problem, the non-wetting phase is \ref Dune::DNAPL .
     */
    const NonwettingPhase &nonwettingPhase() const
    { return nPhase_; }

    /*! 
     * \brief Returns the soil properties object.
     *
     * The lens problem uses \ref Dune::LensSoil .
     */
    Soil &soil()
    {  return soil_; }

    /*! 
     * \copydoc soil()
     */
    const Soil &soil() const
    {  return soil_; }

    /*! 
     * \brief Returns the material laws, i.e. capillary pressure -
     *        saturation and relative permeability-saturation
     *        relations.
     *
     * The lens problem uses the standard \ref Dune::TwoPhaseRelations
     * with Van-Genuchten capillary pressure.
     */
    MaterialLaw &materialLaw ()
    { return materialLaw_; }
    
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*! 
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    void boundaryTypes(BoundaryTypeVector         &values,
                       const Element              &element,
                       const FVElementGeometry    &fvElemGeom,
                       const IntersectionIterator &isIt,
                       int                         scvIdx,
                       int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        //const LocalPosition &localPos
        //    = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);

        values = BoundaryConditions::neumann;

        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
            values = BoundaryConditions::dirichlet;
    }

    /*! 
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVarVector           &values,
                   const Element              &element,
                   const FVElementGeometry    &fvElemGeom,
                   const IntersectionIterator &isIt,
                   int                         scvIdx,
                   int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        
        Scalar densityW = wettingPhase().density();
        
        if (onLeftBoundary_(globalPos))
        {
            Scalar height = outerUpperRight_[1] - outerLowerLeft_[1];
            Scalar depth = outerUpperRight_[1] - globalPos[1];
            Scalar alpha = (1 + 0.5/height);

            // hydrostatic pressure scaled by alpha
            values[pWIdx] = - alpha*densityW*gravity_[1]*depth;
            values[sNIdx] = 0.0;
        }
        else if (onRightBoundary_(globalPos))
        {
            Scalar depth = outerUpperRight_[1] - globalPos[1];

            // hydrostatic pressure
            values[pWIdx] = -densityW*gravity_[1]*depth;
            values[sNIdx] = 0.0;
        }
        else
            values = 0.0;
    }

    /*! 
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVarVector           &values,
                 const Element              &element,
                 const FVElementGeometry    &fvElemGeom,
                 const IntersectionIterator &isIt,
                 int                         scvIdx,
                 int                         boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        //const LocalPosition &localPos
        //    = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);

        values = 0.0;
        if (onInlet_(globalPos)) {
            values[sNIdx] = -0.04; // kg / (m * s)
        }
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*! 
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVarVector        &values,
                const Element           &element,
                const FVElementGeometry &,
                int subControlVolumeIdx) const
    {
        values = Scalar(0.0);
    }

    /*! 
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVarVector        &values,
                 const Element           &element,
                 const FVElementGeometry &fvElemGeom,
                 int                      scvIdx) const
    {
        // no DNAPL, some random pressure
        values[pWIdx] = 0.0;
        values[sNIdx] = 0.0;
    }
    // \}

    /*!
     * \name Restart mechanism
     */
    // \{

    /*!
     * \brief This method writes the complete state of the problem
     *        to the harddisk.
     *
     * The file will start with the prefix <tt>lens</lens>, contains
     * the current time of the simulation clock in it's name and has
     * the prefix <tt>.drs</tt>. (DuMuX Restart File.) See \ref
     * Dune::Restart for details.
     */
    void serialize()
    {
        typedef Dune::Restart<GridView> Restarter;

        Restarter res;
        res.serializeBegin(this->gridView(),
                           "lens",
                           timeManager_.time());

        timeManager_.serialize(res);
        resultWriter_.serialize(res);
        model_.serialize(res);

        res.serializeEnd();
    }

    /*!
     * \brief This method restores the complete state of the problem
     *        from disk.
     *
     * It is the inverse of the \ref serialize() method.
     */
    void deserialize(double t)
    {
        typedef Dune::Restart<GridView> Restarter;

        Restarter res;
        res.deserializeBegin(this->gridView(), "lens", t);

        timeManager_.deserialize(res);
        resultWriter_.deserialize(res);
        model_.deserialize(res);

        res.deserializeEnd();

        wasRestarted_ = true;
    };

    // \}

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < outerLowerLeft_[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > outerUpperRight_[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < outerLowerLeft_[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > outerUpperRight_[1] - eps_;
    }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar width = outerUpperRight_[0] - outerLowerLeft_[0];
        Scalar lambda = (outerUpperRight_[0] - globalPos[0])/width;
        return onUpperBoundary_(globalPos) && 0.5 < lambda  && lambda < 2.0/3.0;
    }
    
    // write the fields current solution into an VTK output file.
    void writeCurrentResult_()
    {
        resultWriter_.beginTimestep(timeManager_.time(),
                                    ParentType::grid().leafView());

        model_.addVtkFields(resultWriter_);

        resultWriter_.endTimestep();
    }


    static const Scalar eps_ = 3e-6;
    GlobalPosition  gravity_;

    GlobalPosition outerLowerLeft_;
    GlobalPosition outerUpperRight_;

    // fluids and material properties
    WettingPhase    wPhase_;
    NonwettingPhase nPhase_;
    Soil            soil_;
    MaterialLaw     materialLaw_;

    TimeManager     timeManager_;

    Model            model_;
    NewtonMethod     newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter  resultWriter_;

    bool wasRestarted_;
};
} //end namespace

#endif
