#ifndef DUNE_NEW_COUPLEPROBLEMPIPE_HH
#define DUNE_NEW_COUPLEPROBLEMPIPE_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include<iostream>
#include<iomanip>

#include<dune/grid/common/grid.hh>

#include<dumux/material/property_baseclasses.hh>
#include<dumux/material/relperm_pc_law.hh>

#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include <dumux/material/matrixproperties.hh>
#include <dumux/material/twophaserelations.hh>

#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/auxiliary/timemanager.hh>

#include <dune/common/timer.hh>

#include<dumux/new_models/1p/1pboxmodel.hh>


#include<dumux/nonlinear/new_newtonmethod.hh>
#include<dumux/nonlinear/new_newtoncontroller.hh>

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>

#include "pipecouplenewtoncontroller.hh"

// required for pipe
#include "pipe_headers.hh"

// timeloop options for pipe
void TimeloopOptsPipe( double& tstart, double& tend, double& max_dt, double& first_dt, double& CFL_factor, int& flag,
                       int& n_iter, double& max_def, int& modulo, int& stages )
{
    tstart      = 0.0;      // start time of simulation
    tend        = 0.1;    // end time of simulation
    max_dt      = 0.1;      // maximum time step
    first_dt    = 1;      // maximum first time step
    CFL_factor  = 1e0;     // safety-factor for Courant-Friedrich-Levy criterion
    flag        = 0;        // 0: no iteration
    // 1: iteration with n_iter steps
    // 2: iteration with max_def defect and timestep-refinement after n_iter iterations
    n_iter      = 50;
    max_def     = 1e-12;    // tolerance
    modulo      = 1;        // vtk - output is written for every modulo-th timestep

    stages      = 0;        // = 0: implicit Euler
    // > 0: number of stages for explicit RK method
}

/**
 * @file
 * @brief  Definition of a problem, where air is injected under a low permeable layer
 * @author Bernd Flemisch, Klaus Mosthaf
 */

namespace Dune
{

//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////--SOIL--//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

/** \todo Please doc me! */

template<class Grid, class ScalarT>
class OnePSoil: public Matrix2p<Grid,ScalarT>
{
public:
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef ScalarT Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {dim=Grid::dimension, dimWorld=Grid::dimensionworld};

    typedef Dune::FieldVector<CoordScalar,dim>      LocalPosition;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

    OnePSoil():Matrix2p<Grid,Scalar>()
    {
        Kout_ = 0.;
        for(int i = 0; i < dim; i++)
            Kout_[i][i] = 5e-10;
    }

    ~OnePSoil()
    {}

    const FieldMatrix<CoordScalar,dim,dim> &K (const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
    {
        return Kout_;
    }

    double porosity(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
    {
        return 0.4;
    }

    double Sr_w(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T) const
    {
        return 0.05;
    }

    double Sr_n(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T) const
    {
        return 0.0;
    }

    /* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
     * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
    double heatCap(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
    {
        return     790 /* spec. heat cap. of granite */
            * 2700 /* density of granite */
            * porosity(x, e, xi);
    }

    double heatCond(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double sat) const
    {
        static const double lWater = 0.6;
        static const double lGranite = 2.8;
        double poro = porosity(x, e, xi);
        double lsat = pow(lGranite, (1-poro)) * pow(lWater, poro);
        double ldry = pow(lGranite, (1-poro));
        return ldry + sqrt(sat) * (ldry - lsat);
    }

    std::vector<double> paramRelPerm(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T) const
    {
        // example for Brooks-Corey parameters
        std::vector<double> param(2);
        param[0] = 0; // pCMin
        param[1] = 1e5; // pCMax

        return param;
    }

    typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
    {
        return Matrix2p<Grid,Scalar>::linear;
    }

private:
    FieldMatrix<Scalar,dim,dim> Kout_;
};

//! class that defines the parameters of an air injection under a low permeable layer
/*! Problem definition of an air injection under a low permeable layer. Air enters the domain
 * at the right boundary and migrates upwards.
 * Problem was set up using the rect2d.dgf grid.
 *
 *    Template parameters are:
 *
 *    - ScalarT  Floating point type used for scalars
 */
template<class GridT, class ScalarT>
class NewCoupleProblemPipe : public BasicDomain<GridT,
                                                ScalarT>
{
    typedef GridT                                  Grid;
    typedef BasicDomain<Grid, ScalarT>             ParentType;
    typedef NewCoupleProblemPipe<GridT, ScalarT>   ThisType;
    typedef OnePBoxModel<ThisType>                 Model;

    typedef Dune::Water                            Fluid;
    typedef Dune::OnePSoil<Grid, ScalarT>          Soil;

public:
    // the domain traits of the domain
    typedef typename ParentType::DomainTraits   DomainTraits;
    // the traits of the BOX scheme
    typedef typename Model::BoxTraits           BoxTraits;
    // the traits of the 1 phase model
    typedef typename Model::OnePTraits          OnePTraits;

private:
    // some constants from the traits for convenience
    enum {
        numEq     = BoxTraits::numEq,
        pIdx      = OnePTraits::pIdx,

        // Grid and world dimension
        dim      = DomainTraits::dim,
        dimWorld = DomainTraits::dimWorld
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
    typedef Dune::TimeManager<Episode>                        TimeManager;

    typedef typename Model::NewtonMethod                      NewtonMethod;
    typedef Dune::PipeCoupleNewtonController<NewtonMethod>    NewtonController;

    //////////////////////////
    // required for pipeflow
    //////////////////////////

    typedef Dune::BlockVector< Dune::FieldVector<ScalarT,1> > PressureField;

    // set the initial condition
    typedef ICPressurePipe<GlobalPosition> InitialPressure;
    typedef ICVelocityPipe<GlobalPosition> InitialVelocity;

    // set the initial condition
    typedef sourceSinkPipe<GlobalPosition> SourceSinkTerm;

    // set friction coefficients
    typedef Lambda<Fluid> LambdaPipe;

    typedef LambdaLocal<GlobalPosition> LambdaPipeLocal;

    // make a mapper for codim DIM entities in the leaf grid
    typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid,PDimLayout> VertexMapper;



    //////////////////////////
    // required to extract the 1d pipe network
    //////////////////////////
    typedef VertexOnLine<Grid>         VertexType1;
    typedef VertexOutLine<Grid>        VertexType2;
    typedef std::vector<VertexType1>   VertexVectorOnLineType;
    typedef std::vector<VertexType2>   VertexVectorOutLineType;

    typedef std::map<unsigned, unsigned, GlobalNodeIdCompare> MapGlobalNodeIDtoPipeNodeOnOutlineIndexType;


    //////////////////////////
    // The pipe flow class
    //////////////////////////
    typedef ::PipeFlow<PressureBoundary,
                       VelocityBoundary,
                       InitialPressure,
                       InitialVelocity,
                       SourceSinkTerm,
                       PressureField,
                       LambdaPipe,
                       LambdaPipeLocal,
                       Grid,
                       VertexMapper,
                       MapGlobalNodeIDtoPipeNodeOnOutlineIndexType,
                       VertexVectorOnLineType,
                       VertexVectorOutLineType> PipeFlow;

public:
    NewCoupleProblemPipe(GridPtr<Grid> gridPtr,
                         Scalar dtInitial,
                         Scalar tEnd)
        : ParentType(&(*gridPtr)),
          gridPtr_(gridPtr),
          timeManager_(tEnd, this->grid().comm().rank() == 0),
          model_(*this),
          newtonMethod_(model_),
          newtonCtl_(1e-10, 20, 30, 100000),
          resultWriter_("new_couple_pipe"),
          vertexMapper_(*gridPtr_)
    {
        timeManager_.setStepSize(dtInitial);

        gravity_ = 0;
        //           gravity_[dim - 1] = -9.81;
        temperature_ = 283.15;

        // pipe parameters
        pipeRoughness_ = 0.3;
        pipeDiameter_ = 0.01; // in [m]

        // fluid exchange parameter
        alphaExchange_ = 5.0e-10 * 1.0e+2;

        pipeFlow_ = NULL;

        // isolate 1D grid network
        isolate1DNetwork_();
    };

    ~NewCoupleProblemPipe()
    {
        delete pipeFlow_;
    };

protected:
    void isolate1DNetwork_()
    {
        typedef std::map<unsigned, unsigned, GlobalNodeIdCompare> MapGlobalNodeIDtoPipeNodeOnOutlineIndexType;

        isolate(*gridPtr_,
                gridPtr_,
                vertexMapper_,
                mapGlobalNodeIDtoPipeNodeOnlineIndex_,
                mapGlobalNodeIDtoPipeNodeOutlineIndex_,
                vertexVectorOnLine_,
                vertexVectorOutLine_);
    }

public:
    ///////////////////////////////////
    // Strings pulled by the TimeManager during the course of the
    // simulation
    ///////////////////////////////////

    //! called by the time manager in order to create the initial
    //! solution
    void init()
    {
        // initialize the volume of the finite volume boxes
        boxVolumes_.resize(ParentType::numVertices());
        boxVolumes_ = Scalar(0.0);
        ElementIterator eIt = ParentType::elementBegin();
        const ElementIterator &eEndIt = ParentType::elementEnd();
        FVElementGeometry fvElemGeom;
        for (; eIt != eEndIt; ++eIt) {
            fvElemGeom.update(*eIt);
            for (int vertIdx = 0; vertIdx < fvElemGeom.numVertices; ++vertIdx) {
                boxVolumes_[ParentType::vertexIdx(*eIt, vertIdx)][0] += fvElemGeom.subContVol[vertIdx].volume;
            }
        }
        // set the initial condition
        model_.initial();

        // initialize the pipe data structures
        initPipe_();

        // write the inital solution to disk
        writeCurrentResult_();
    }

private:
    void initPipe_()
    {
        int numNodesPipeMesh = vertexVectorOnLine_.size();
        delete pipeFlow_;

        std::cout<<"number of pipe nodes: "<<numNodesPipeMesh<<std::endl;

        pipeFlow_ = new PipeFlow(numNodesPipeMesh,
                                 this->grid(),
                                 vertexMapper_,
                                 mapGlobalNodeIDtoPipeNodeOnlineIndex_,
                                 mapGlobalNodeIDtoPipeNodeOutlineIndex_,
                                 vertexVectorOnLine_,
                                 vertexVectorOutLine_,
                                 alphaExchange_,
                                 temperature_,
                                 fluid_.density(),
                                 fluid_.viscosity()/fluid_.density(),
                                 pipeRoughness_,
                                 pipeDiameter_,
                                 gravity_);
        pipeFlow_->pressure = Scalar(1.0e+5);
    }

public:
    /*!
     * \brief Update the pipe pressure using the current porous
     *        pressure field.
     */
    void updatePipeFlow()
    {
        pipeFlow_->pressurePorous = *model_.currentSolution();

        LocalPosition localPos(0);

        typedef typename MapGlobalNodeIDtoPipeNodeOnOutlineIndexType::iterator PipeMapIter;
        // calculate the mobility field
        VertexIterator vIt = ParentType::vertexBegin();
        const VertexIterator &vEndIt = ParentType::vertexEnd();
        for (; vIt != vEndIt; ++vIt) {
            int vIdx = ParentType::vertexIdx(*vIt);

            PipeMapIter pipeIt = mapGlobalNodeIDtoPipeNodeOnlineIndex_.find(vIdx);
            if (pipeIt == mapGlobalNodeIDtoPipeNodeOnlineIndex_.end())
                // vertex is not on a pipe
                continue;

            //            int pipeIdx = pipeIt->second; //TODO decrease mobility and porous media size to pipe size

            pipeFlow_->mobility[vIdx][0] = 1.0;
        };

        // recalculate the pressure inside the pipe network
        pipeFlow_->update();
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
        model_.update(stepSize, nextStepSize, newtonMethod_, newtonCtl_);
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

    /*! properties of the fluid
      \return    fluid
    */

    //! properties of the fluid
    const Fluid &fluid() const
    { return fluid_; }


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
        values = Dune::BoundaryConditions::neumann;
        switch (isIt->boundaryId()) {
            /*                case 1:
                              case 2:
                              case 3:
                              case 4:
                              values = Dune::BoundaryConditions::neumann;
                              break;
            */
        case 3:
        case 4:
            values = Dune::BoundaryConditions::dirichlet;
            break;
        }
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
        values[pIdx] = 2.0e+5;;

        switch (isIt->boundaryId()) {
        case 4:
            values[pIdx] = 2.0e+5;;
            break;
        }
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
        values = 0;
        switch (isIt->boundaryId()) {
        case 1:
        case 2:
        case 5:
        case 6:
            values[pIdx] = 0;
            break;
            /*                case 5:
                              values[pIdx] = -1.0;
                              break;
            */
        }
    }

    /////////////////////////////
    // sources and sinks
    /////////////////////////////
    void source(SolutionVector          &result,
                const Element           &element,
                const FVElementGeometry &fvElemGeom,
                int                      scvIdx) const
    {
        result = Scalar(0.0);

        int porousIdx = ParentType::vertexIdx(element, scvIdx);

        typedef typename MapGlobalNodeIDtoPipeNodeOnOutlineIndexType::const_iterator PipeMapIter;
        PipeMapIter pipeIt = mapGlobalNodeIDtoPipeNodeOnlineIndex_.find(porousIdx);
        if (pipeIt == mapGlobalNodeIDtoPipeNodeOnlineIndex_.end())
            // vertex is not on a pipe
            return;

        int pipeIdx = pipeIt->second;

        Scalar porousSol = (*model_.currentSolution())[porousIdx];
        Scalar pPipe = pipeFlow_->pressure[pipeIdx];

        double alphaEXCHANGE =
            fluid_.density() *
            pipeFlow_->mobility[porousIdx][0] *
            pipeFlow_->alphaExchange *
            (M_PI *
             pipeDiameter_ *
             vertexVectorOnLine_[pipeIdx].length(vertexVectorOnLine_) )
            / pipeDiameter_;

        // calculate total mass going into or out of the pipe
        result[pIdx] = alphaEXCHANGE * (pPipe - porousSol);

        // dvided with volume of the sub control volume because we
        // need a volume averaged quantity in the source term.
        result[pIdx] /= boxVolumes_[porousIdx][0];
    }

    //////////////////////////////

    /////////////////////////////
    // INITIAL values
    /////////////////////////////
    void initial(SolutionVector         &values,
                 const Element           &element,
                 const FVElementGeometry &fvElemGeom,
                 int                      scvIdx) const
    {
        values[pIdx] = 2.0e+5;
    }


    Scalar temperature() const
    {
        return temperature_; // 10°C
    };

    Scalar porosity(const Element &element, int localIdx) const
    {
        // TODO/HACK: porosity should be defined on the verts
        // as it is required on the verts!
        const LocalPosition &local =
            DomainTraits::referenceElement(element.type()).position(localIdx, dim);
        const GlobalPosition &globalPos = element.geometry().corner(localIdx);
        return soil().porosity(globalPos, *(ParentType::elementBegin()), local);
    };

    const GlobalPosition &gravity () const
    {
        return gravity_;
    }

    bool simulate()
    {
        timeManager_.runSimulation(*this);
        return true;
    };


private:
    // write the fields current solution into an VTK output file.
    void writeCurrentResult_()
    {
        resultWriter_.beginTimestep(timeManager_.time(),
                                    ParentType::grid().leafView());

        model_.addVtkFields(resultWriter_);
        pipeFlow_->addVtkFields(resultWriter_);

        resultWriter_.endTimestep();
    }


    GridPtr<Grid> gridPtr_;
    GlobalPosition  gravity_;

    // fluids and material properties
    Fluid           fluid_;
    Soil            soil_;
    Scalar          temperature_;

    TimeManager     timeManager_;

    Model            model_;
    NewtonMethod     newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter  resultWriter_;

    VertexMapper vertexMapper_;
    // required for pipe
    Scalar    pipeRoughness_;
    Scalar    pipeDiameter_;

    VertexVectorOnLineType vertexVectorOnLine_;
    VertexVectorOutLineType vertexVectorOutLine_;

    MapGlobalNodeIDtoPipeNodeOnOutlineIndexType mapGlobalNodeIDtoPipeNodeOnlineIndex_;
    MapGlobalNodeIDtoPipeNodeOnOutlineIndexType mapGlobalNodeIDtoPipeNodeOutlineIndex_;

    PipeFlow      *pipeFlow_;
    // total volumes of the FV boxes
    PressureField boxVolumes_;

    // the exchange parameter between pipe and porous flow
    Scalar alphaExchange_;
};
} //end namespace

#endif
