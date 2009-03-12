#ifndef DUNE_NEW_BENCHMARK3PROBLEM_HH
#define DUNE_NEW_BENCHMARK3PROBLEM_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include<iostream>
#include<iomanip>

#include<dune/grid/common/grid.hh>

#include<dumux/material/property_baseclasses.hh>
#include<dumux/material/relperm_pc_law.hh>

#include <dumux/material/phaseproperties/phaseproperties_brineco2.hh>
#include <dumux/material/matrixproperties.hh>
#include <dumux/material/twophaserelations.hh>
#include "benchmark3soil.hh"

#include <dumux/material/multicomponentrelations.hh>

#include <dumux/io/vtkmultiwriter.hh>
#include <dumux/auxiliary/timemanager.hh>
#include <dumux/io/restart.hh>

#include <dune/common/timer.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/istl/io.hh>

#include<dumux/new_models/2p2cni/2p2cniboxmodel.hh>
#include<dumux/new_models/2p2cni/2p2cninewtoncontroller.hh>

#include<dumux/timedisc/new_impliciteulerstep.hh>
#include<dumux/nonlinear/new_newtonmethod.hh>

#include <dumux/auxiliary/timemanager.hh>
#include <dumux/auxiliary/basicdomain.hh>

#define ISOTHERMAL 0

/**
 * @file
 * @brief  Definition of a problem, where co2 is injected into a deep saline aquifer (Johansen Formation, Norway)
 * @author Bernd Flemisch, Klaus Mosthaf, Melanie Darcis
 */

namespace Dune
{
//! class that defines the parameters of a co2 injection into a deep saline formation
/*! Problem definition of a co2 injection into a deep saline formation. co2 enters the domain
 * at the bottom center and migrates upwards.
 * Problem was set up using the benchmark3.dgf grid.
 *
 *    Template parameters are:
 *
 *    - ScalarT  Floating point type used for scalars
 */
template<class ScalarT>
class NewBenchmark3Problem : public BasicDomain<Dune::UGGrid<3>,
                                              ScalarT>
{
    typedef Dune::UGGrid<3>              Grid;
    typedef BasicDomain<Grid, ScalarT>     ParentType;
    typedef NewBenchmark3Problem<ScalarT>    ThisType;
    typedef CollectiveCommunication<ThisType> CollectiveCommunication;
#if !ISOTHERMAL
    typedef TwoPTwoCNIBoxModel<ThisType>   Model;
#else
    typedef TwoPTwoCBoxModel<ThisType>   Model;
#endif

    typedef Dune::GridPtr<Grid>                    GridPointer;

    typedef Dune::Liq_BrineCO2                     WettingPhase;
    typedef Dune::Gas_BrineCO2                     NonwettingPhase;

    typedef Dune::Benchmark3Soil<Grid, ScalarT>    Soil;
    typedef Dune::TwoPhaseRelations<Grid, ScalarT> MaterialLaw;
    typedef Dune::CBrineCO2                        Multicomp;

public:
    // the domain traits of the domain
    typedef typename ParentType::DomainTraits   DomainTraits;
    // the traits of the BOX scheme
    typedef typename Model::BoxTraits           BoxTraits;
    // the traits of the Pw-Sn model
#if !ISOTHERMAL
    typedef typename Model::TwoPTwoCNITraits    TwoPTwoCNITraits;
#else
    typedef typename Model::TwoPTwoCTraits      TwoPTwoCNITraits;
#endif

private:
    // some constants from the traits for convenience
    enum {
        numEq          = BoxTraits::numEq,
        pressureIdx    = TwoPTwoCNITraits::pressureIdx,
        switchIdx      = TwoPTwoCNITraits::switchIdx,
#if !ISOTHERMAL
        temperatureIdx = TwoPTwoCNITraits::temperatureIdx,
#endif

        // Phase State
        WPhaseOnly = TwoPTwoCNITraits::wPhaseOnly,
        NPhaseOnly = TwoPTwoCNITraits::nPhaseOnly,
        BothPhases = TwoPTwoCNITraits::bothPhases,

        // Grid and world dimension
        dim  = DomainTraits::dim,
        dimWorld = DomainTraits::dimWorld,

        pWsN        = 0,
        pNsW        = 1
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
    typedef TwoPTwoCNINewtonController<NewtonMethod>    NewtonController;

public:
    NewBenchmark3Problem(Grid *grid,
                       Scalar dtInitial,
                       Scalar tEnd)
        : ParentType(grid),
          soil_(this->grid(), true, "properties_johansen.dat"),
          materialLaw_(soil_, wPhase_, nPhase_),
          multicomp_(wPhase_, nPhase_),
          timeManager_(tEnd, this->grid().comm().rank() == 0),
          model_(*this),
          newtonMethod_(model_),
          resultWriter_("new_brineco2")
    {
        timeManager_.setStepSize(dtInitial);

        // specify the grid dimensions
        eps_    = 1e-8;

        depthBOR_ = 3238.20;

        gravity_[0] = 0;
        gravity_[1] = 0;
        gravity_[2] = -9.81;
        // choose primary variables
        formulation_ = pWsN;

        wasRestarted_ = false;
    }

    ///////////////////////////////////
    // Strings pulled by the TimeManager during the course of the
    // simulation
    ///////////////////////////////////

    //! called by the time manager in order to create the initial
    //! solution
    void init()
    {
        // set the episode length and initial time step size
        // timeManager_.setStepSize(initialTimeStepSize_);

        // read permeability and porosity from data file
        soil().readSoilProperties();
        soil().setSoilProperties();

        // set the initial condition
        model_.initial();

        // write the inital solution to disk
        if (!wasRestarted_)
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
        model_.update(stepSize, nextStepSize, newtonMethod_, newtonCtl_);
    }

    //! called by the TimeManager whenever a solution for a
    //! timestep has been computed
    void timestepDone()
    {
        if (this->grid().comm().rank() == 0)
            std::cout << "Writing result file for current time step\n";

        // write restart file after every five steps
         static int dummy = 0;
         ++dummy;
         if (dummy % 5 == 0)
             serialize();
        // write the current result to disk
        writeCurrentResult_();

        // update the domain with the current solution
        //                updateDomain_();
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


    //! properties of the wetting (liquid) phase
    /*! properties of the wetting (liquid) phase
      \return    wetting phase
    */
    const WettingPhase &wettingPhase() const
    { return wPhase_; }

    //! properties of the nonwetting (liquid) phase
    /*! properties of the nonwetting (liquid) phase
      \return    nonwetting phase
    */
    const NonwettingPhase &nonwettingPhase() const
    { return nPhase_; }


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

    //! object for multicomponent calculations
    /*! object for multicomponent calculations including mass fractions,
     * mole fractions and some basic laws
     \return    multicomponent object
    */
    MultiComp &multicomp ()
    //        const MultiComp &multicomp () const
    {
        return multicomp_;
    }

    //! object for definition of material law
    /*! object for definition of material law (e.g. Brooks-Corey, Van Genuchten, ...)
      \return    material law
    */
    MaterialLaw &materialLaw ()
    //        const MaterialLaw &materialLaw () const
    {
        return materialLaw_;
    }

    void boundaryTypes(BoundaryTypeVector         &values,
                       const Element              &element,
                       const FVElementGeometry    &fvElemGeom,
                       const IntersectionIterator &isIt,
                       int                         scvIdx,
                       int                         boundaryFaceIdx) const
    {
        // if(isIt->boundaryId() == 1) Boundary Ids only work for AluGrid !!!
        FieldVector<Scalar,dim-1> normalPos(0);
        const LocalPosition normal = isIt->unitOuterNormal(normalPos);
        Scalar normalZComp, absNormalZComp;
        normalZComp = normal[dim-1];
        absNormalZComp = std::fabs(normalZComp);

        if(absNormalZComp > 0.2)
           values = BoundaryConditions::neumann;
        else
           values = BoundaryConditions::dirichlet;


#if !ISOTHERMAL
        values[temperatureIdx] = BoundaryConditions::dirichlet;
#endif
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
            = fvElemGeom.subContVol[scvIdx].global; // WHY is dirichlet given for IP?
        //                const LocalPosition &localPos
        //                    = fvElemGeom.boundaryFace[boundaryFaceIdx].ipLocal;

        initial_(values, globalPos);
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
            = fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;
        //                const LocalPosition &localPos
        //                    = fvElemGeom.boundaryFace[boundaryFaceIdx].ipLocal;

        values = 0;



    }

    /////////////////////////////
    // sources and sinks
    /////////////////////////////
    void source(SolutionVector          &values,
                const Element           &element,
                const FVElementGeometry &fvElemGeom,
                int                      scvIdx) const
    {
         const GlobalPosition &globalPos = fvElemGeom.subContVol[scvIdx].global;
         Scalar m;
         values = Scalar(0.0);

         m = 2.573325096E-5;
         if (globalPos[0] > 5399. && globalPos[0] < 5508. && globalPos[1] > 3285. && globalPos[1] < 3393. && globalPos[2] < 286.)
             values[1] = m; //heat flux calculated with actual enthalpy values in
                            // 2p2cniboxmodel.hh
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
            = fvElemGeom.subContVol[scvIdx].global;
        /*                const LocalPosition &localPos
                          = fvElemGeom.subContVol[scvIdx].local;
        */

        initial_(values, globalPos);
    }

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(SolutionVector       &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar densityB = 1037.24; // Brine density
        values[pressureIdx] = 1.013e5 + (depthBOR_ - globalPos[2])*densityB*9.81;
        values[switchIdx] = 0.0;
#if !ISOTHERMAL
        values[temperatureIdx] = 283.15 + (depthBOR_ - globalPos[2])*0.03;
#endif
    }

public:
    Scalar porosity(const Element &element, int localIdx) const
    {
        // TODO/HACK: porosity should be defined on the verts
        // as it is required on the verts!
        const LocalPosition &local =
            DomainTraits::referenceElement(element.type()).position(localIdx, dim);
        const GlobalPosition &globalPos = element.geometry().corner(localIdx);
        return soil().porosity(globalPos, element, local, localIdx);
    };

    FieldMatrix<Scalar,dim,dim> &K(const Element &element, int localIdx)
    {
        const LocalPosition &local =
            DomainTraits::referenceElement(element.type()).position(localIdx, dim);
        const GlobalPosition &globalPos = element.geometry().corner(localIdx);
        return soil().K(globalPos, element, local, localIdx);
    }

    Scalar pC(Scalar satW, int globalIdx, const GlobalPosition &globalPos)
    {
        // TODO/HACK: porosity should be defined on the verts
        // as it is required on the verts!
        const LocalPosition &local =
            DomainTraits::referenceElement(ParentType::elementBegin()->type()).position(0, dim);
        return materialLaw().pC(satW, globalPos,*(ParentType::elementBegin()), local);
    };


    int initialPhaseState(const Vertex       &vert,
                          int              &globalIdx,
                          const GlobalPosition &globalPos) const
    {
        return WPhaseOnly;
    }


    const GlobalPosition &gravity () const
    {
        return gravity_;
    }

    double depthBOR () const
    {
        return depthBOR_;
    }

    bool simulate()
    {
        timeManager_.runSimulation(*this);
        return true;
    };

    int formulation () const
    {
        return formulation_;
    }

    void serialize()
    {
        typedef Dune::Restart<Grid> Restarter;

        Restarter res;
        res.serializeBegin(this->grid(),
                           "neweasy",
                           timeManager_.time());

        timeManager_.serialize(res);
        resultWriter_.serialize(res);
        model_.serialize(res);

        res.serializeEnd();
    }

    void deserialize(double t)
    {
        typedef Dune::Restart<Grid> Restarter;

        Restarter res;
        res.deserializeBegin(this->grid(), "neweasy", t);

        timeManager_.deserialize(res);
        resultWriter_.deserialize(res);
        model_.deserialize(res);

        res.deserializeEnd();

        wasRestarted_ = true;
    };

private:
    // write the fields current solution into an VTK output file.
    void writeCurrentResult_()
    {
        resultWriter_.beginTimestep(timeManager_.time(),
                                    ParentType::grid().leafView());

        model_.addVtkFields(resultWriter_);

        Dune::FieldVector<Scalar, 4> mass;
        model_.calculateMass(mass);

        if(collectiveCom_.rank() == 0){
        std::cout<< "Mass[kg] (nC, nC in nP, wC, wC in wP):"
        << mass[0] <<", "<<mass[1]<<", "<< mass[2]<<", "<< mass[3] <<". ";
        }

        resultWriter_.endTimestep();
    }


    Scalar width_;
    Scalar height_;
    Scalar depthBOR_;
    Scalar eps_;
    int formulation_;
    GlobalPosition  gravity_;

    // fluids and material properties
    WettingPhase    wPhase_;
    NonwettingPhase nPhase_;
    Soil            soil_;
    MaterialLaw     materialLaw_;
    Multicomp       multicomp_;

    TimeManager     timeManager_;

    Model            model_;
    NewtonMethod     newtonMethod_;
    NewtonController newtonCtl_;

    VtkMultiWriter  resultWriter_;
    CollectiveCommunication collectiveCom_;

    bool wasRestarted_;
};
} //end namespace

#endif
