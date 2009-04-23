//$Id$
#ifndef DUNE_NEW_ELASTICMATRIXPROBLEM_HH
#define DUNE_NEW_ELASTICMATRIXPROBLEM_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include<iostream>
#include<iomanip>

#include<dune/grid/common/grid.hh>

#include<dumux/material/property_baseclasses.hh>

#include <dumux/material/matrixproperties.hh>
#include "elasticsoil.hh"

#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/auxiliary/timemanager.hh>

#include <dune/common/timer.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/uggrid.hh>
#include <dune/istl/io.hh>

#include<dumux/new_models/linearelasticity/elasticboxmodel.hh>
#include<dumux/new_models/linearelasticity/elasticnewtoncontroller.hh>

#include<dumux/nonlinear/new_newtonmethod.hh>

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>


/**
 * @file
 * @brief  Definition of a problem, where an elastic matrix is deformed
 * @author Melanie Darcis
 */

namespace Dune
{

    //! class that defines the parameters for the deformation of an elastic matrix
    /*! Problem definition for the deformation of an elastic matrix.
     * A solid displacement is posed on the left side of the box.
     * All remaining boundary faces are fixed.
     *
     *    Template parameters are:
     *
     *    - ScalarT  Floating point type used for scalars
     */
    template<class ScalarT>
    class ElasticProblem : public BasicDomain<Dune::UGGrid<3>,
                                                   ScalarT>
    {
        typedef Dune::UGGrid<3>                     Grid;
        typedef BasicDomain<Grid, ScalarT>          ParentType;
        typedef ElasticProblem<ScalarT>             ThisType;
        typedef ElasticBoxModel<ThisType>           Model;

        typedef Dune::ElasticSoil<Grid, ScalarT>     Soil;

    public:
        // the domain traits of the domain
        typedef typename ParentType::DomainTraits   DomainTraits;
        // the traits of the BOX scheme
        typedef typename Model::BoxTraits           BoxTraits;
        // the traits of the Pw-Sn model
        typedef typename Model::ElasticTraits       ElasticTraits;

    private:
        // some constants from the traits for convenience
        enum {
            numEq            = BoxTraits::numEq,
            uxIdx            = ElasticTraits::uxIdx,
            uyIdx            = ElasticTraits::uyIdx,
            uzIdx            = ElasticTraits::uzIdx,

            // Grid and world dimension
            dim      	= DomainTraits::dim,
            dimWorld 	= DomainTraits::dimWorld
        };

        // copy some types from the traits for convenience
        typedef typename DomainTraits::Scalar                     Scalar;
        typedef typename DomainTraits::Element                    Element;
        typedef typename DomainTraits::ElementIterator            ElementIterator;
        typedef typename DomainTraits::ReferenceElement           ReferenceElement;
        typedef typename DomainTraits::Vertex                     Vertex;
        typedef typename DomainTraits::VertexIterator             VertexIterator;
        typedef typename DomainTraits::IntersectionIterator       IntersectionIterator;
        typedef typename DomainTraits::LocalPosition              LocalPosition;
        typedef typename DomainTraits::GlobalPosition             GlobalPosition;

        typedef typename BoxTraits::FVElementGeometry             FVElementGeometry;
        typedef typename BoxTraits::SpatialFunction               SpatialFunction;
        typedef typename BoxTraits::SolutionVector                SolutionVector;
        typedef typename BoxTraits::BoundaryTypeVector            BoundaryTypeVector;

        typedef Dune::VtkMultiWriter<typename Grid::LeafGridView> VtkMultiWriter;

        enum Episode {}; // the type of an episode of the simulation
        typedef Dune::TimeManager<Episode>                  TimeManager;

        typedef typename Model::NewtonMethod                NewtonMethod;
        typedef ElasticNewtonController<NewtonMethod>      NewtonController;

    public:
        ElasticProblem(Grid *grid,
                            Scalar dtInitial,
                            Scalar tEnd)
            : ParentType(grid),
              timeManager_(tEnd, this->grid().comm().rank() == 0),
              model_(*this),
              newtonMethod_(model_),
              resultWriter_("newMatrix")
            {
                timeManager_.setStepSize(dtInitial);
        		eps_    = 1e-8;
                // specify the grid dimensions
             }

        ///////////////////////////////////
         // Strings pulled by the TimeManager during the course of the
         // simulation
         ///////////////////////////////////

         //! called by the time manager in order to create the initial
         //! solution
         void init()
         {
             // set the initial condition
             model_.initial();

                // write the inital solution to disk
                 writeCurrentResult_();
         }

        /*!
          * \brief Called by the TimeManager in order to get a time
          *        integration on the model.
          *
          * \note timeStepSize and nextStepSize are references and may
          *       be modified by the TimeIntegration. On exit of this
          *       function 'timeStepSize' must contain the step size
          *       actually used by the time integration for the current
          *       steo, and 'nextStepSize' must contain the suggested
          *       step size for the next time step.
          */
         void timeIntegration(Scalar &stepSize, Scalar &nextStepSize)
         {
             model_.update(stepSize,
                           nextStepSize,
                           newtonMethod_,
                           newtonCtl_);
         }

        //! called by the TimeManager whenever a solution for a
        //! timestep has been computed
        void timestepDone()
        {
            if (this->grid().comm().rank() == 0)
                std::cout << "Writing result file for current time step\n";

            // write the current result to disk
            writeCurrentResult_();
        };
        ///////////////////////////////////
        // End of simulation control stuff
        ///////////////////////////////////

        ///////////////////////////////////
        // Strings pulled by the TwoPTwoCBoxModel during the course of
        // the simulation (-> boundary conditions, initial conditions,
        // etc)
        ///////////////////////////////////
        //! Returns the current time step size in seconds
        Scalar timeStepSize() const
            { return timeManager_.stepSize(); }

        //! Set the time step size in seconds.
        void setTimeStepSize(Scalar dt)
            { return timeManager_.setStepSize(dt); }


        //! properties of the soil
        /*! properties of the soil
          \return    soil
        */
        const Soil &soil() const
            {  return soil_; }

        //! properties of the soil
        /*! properties of the soil
          \return    soil
        */
        Soil &soil()
            {  return soil_; }

        void boundaryTypes(BoundaryTypeVector         &values,
                           const Element              &element,
                           const FVElementGeometry    &fvElemGeom,
                           const IntersectionIterator &isIt,
                           int                         scvIdx,
                           int                         boundaryFaceIdx) const
            {
                const GlobalPosition &globalPos
                    = element.geometry().corner(scvIdx);
                const LocalPosition &localPos
                   = DomainTraits::referenceElement(element.geometry().type()).position(dim,scvIdx);

                if (globalPos[0] > 1.)
                    values = BoundaryConditions::neumann;
                else
                    values = BoundaryConditions::dirichlet;
            }

        /////////////////////////////
        // DIRICHLET boundaries
        /////////////////////////////
        void dirichlet(SolutionVector             &values,
                       const Element              &element,
                       const FVElementGeometry    &fvElemGeom,
                       const IntersectionIterator &isIt,
                       int                         scvIdx,
                       int                         boundaryFaceIdx) const
            {
                const GlobalPosition &globalPos
                    = element.geometry().corner(scvIdx);

                initial_(values, globalPos);

                values = 0.0;
             }

        /////////////////////////////
        // NEUMANN boundaries
        /////////////////////////////
        void neumann(SolutionVector             &values,
                     const Element              &element,
                     const FVElementGeometry    &fvElemGeom,
                     const IntersectionIterator &isIt,
                     int                         scvIdx,
                     int                         boundaryFaceIdx) const
            {
                const GlobalPosition &globalPos
                    = element.geometry().corner(scvIdx);
                values = 0.0;
                if(globalPos[0] > 10-eps_)
                {
                    values[0] = -1e-4;
                }
            }

        /////////////////////////////
        // sources and sinks
        /////////////////////////////
        void source(SolutionVector          &values,
                    const Element           &element,
                    const FVElementGeometry &fvElemGeom,
                    int                      scvIdx) const
        {
            values = Scalar(0.0);
        }

        //////////////////////////////

        /////////////////////////////
        // INITIAL values
        /////////////////////////////
        void initial(SolutionVector          &values,
                     const Element           &element,
                     const FVElementGeometry &fvElemGeom,
                     int                      scvIdx) const
            {
                const GlobalPosition &globalPos
                    = element.geometry().corner(scvIdx);
                initial_(values, globalPos);
            }

         bool simulate()
            {
                timeManager_.runSimulation(*this);
                return true;
            };


    private:
        // the internal method for the initial condition
        void initial_(SolutionVector         &values,
                      const GlobalPosition   &globalPos) const
        {
            values[uxIdx] = 0;
            values[uyIdx] = 0;
            values[uzIdx] = 0;
        }

        // write the fields current solution into an VTK output file.
        void writeCurrentResult_()
            {
                resultWriter_.beginTimestep(timeManager_.time(),
                                            ParentType::grid().leafView());

                model_.addVtkFields(resultWriter_);

                resultWriter_.endTimestep();
            }


        GlobalPosition outerLowerLeftBack_;
        GlobalPosition outerUpperRightFront_;
        GlobalPosition innerLowerLeftBack_;
        GlobalPosition innerUpperRightFront_;
        Scalar width_;
        Scalar height_;
        Scalar depth_;
        Scalar depthBOR_;
        Scalar eps_;
 //      Scalar densityW_, densityN_;
        GlobalPosition  gravity_;

        // fluids and material properties
        Soil            soil_;

        TimeManager     timeManager_;
        Scalar          initialTimeStepSize_;
        Scalar          endTime_;

        Model            model_;
        NewtonMethod     newtonMethod_;
        NewtonController newtonCtl_;

        VtkMultiWriter  resultWriter_;
    };
} //end namespace

#endif
