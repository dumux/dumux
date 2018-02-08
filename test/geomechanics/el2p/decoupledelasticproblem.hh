// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/**
 * \file
 * \brief Definition of a problem, for the two-phase flow linear elasticity problem:
 * Problem definition for the deformation of an elastic solid.
 */
#ifndef DUMUX_EL2P_TESTPROBLEM_HH
#define DUMUX_EL2P_TESTPROBLEM_HH

#include <dune/pdelab/finiteelementmap/qkfem.hh>

#include <dumux/material/fluidsystems/brineco2.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/geomechanics/el2p/decoupled/model.hh>
#include <dumux/geomechanics/el2p/amgbackend.hh>

#include "co2tables.hh"
#include "decoupledelasticspatialparams.hh"

namespace Dumux
{
template<class TypeTag>
class El2P_TestProblem;


// initial conditions for momentum balance equation
template<class TypeTag, int dim>
class InitialDisplacement;

namespace Properties {
NEW_TYPE_TAG(El2P_TestProblem, INHERITS_FROM(BoxModel, BoxElasticTwoP, El2PSpatialParams));
NEW_PROP_TAG(InitialDisplacement); //!< The initial displacement function

// Set the grid type
SET_TYPE_PROP(El2P_TestProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);


SET_PROP(El2P_TestProblem, PressureFEM)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    typedef Dune::PDELab::QkLocalFiniteElementMap<GridView,Scalar,Scalar,1>  type;
};

SET_PROP(El2P_TestProblem, DisplacementFEM)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

public:
    typedef Dune::PDELab::QkLocalFiniteElementMap<GridView,Scalar,Scalar,1>  type;
};

// Set the problem property
SET_TYPE_PROP(El2P_TestProblem, Problem, El2P_TestProblem<TypeTag>);

// Set fluid configuration
SET_PROP(El2P_TestProblem, FluidSystem)
{
public:
    typedef BrineCO2FluidSystem<TypeTag> type;
};

// Set the CO2 table to be used; in this case not the the default table
SET_TYPE_PROP(El2P_TestProblem, CO2Table, CO2Tables);
// Set the salinity mass fraction of the brine in the reservoir
SET_SCALAR_PROP(El2P_TestProblem, ProblemSalinity, 0);

// Set the soil properties
SET_TYPE_PROP(El2P_TestProblem, SpatialParams, El2PSpatialParams<TypeTag>);

// Set the initial displacement function
SET_PROP(El2P_TestProblem, InitialDisplacement)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
public:
    typedef InitialDisplacement<TypeTag, dim> type;
};


SET_SCALAR_PROP(El2P_TestProblem, NewtonMaxRelativeShift, 1e-4);
SET_SCALAR_PROP(El2P_TestProblem, NewtonMaxSteps, 30);
// SET_SCALAR_PROP(El2P_TestProblem, NewtonUseLineSearch, true);

// use the algebraic multigrid
// SET_TYPE_PROP(El2P_TestProblem, LinearSolver, El2PAMGBackend<TypeTag>);
SET_TYPE_PROP(El2P_TestProblem, LinearSolver, SuperLUBackend<TypeTag> );

// central differences to calculate the jacobian by default
SET_INT_PROP(El2P_TestProblem, ImplicitNumericDifferenceMethod, 0);

// write the stress and displacement output according to rock mechanics
// sign convention (compressive stresses > 0)
SET_BOOL_PROP(El2P_TestProblem, VtkRockMechanicsSignConvention, true);
}

/*!
 * \ingroup ElTwoPBoxProblems
 *
 * \brief Problem definition for a two-phase flow process
 * in an elastic deformable matrix.
 *
 * This problem simulates an injection of CO2 into the center of a cube with 1000 m x 1000 m x 1000 m
 * dimension. The bottom boundary of this cube is in 2000 m depth. The initialization period is 1e6 s,
 * the real injection period is 1e6 s, the initial timestep is 10 s.
 * Apart from the pressure and the saturation distribution this problems solves for the changes in
 * solid displacement (ux, uy, uz [m]) due to injection. Based on the solid displacement vector
 * the injection-induced changes in the strain and stress tensors are evaluated.
 * Further the porosity and permeability are functions of the solid displacement.
 *
 * During an initialization period of length tInit [s] the pressure field is initialized.
 *
 * After the initialization the real simulation starts and the pressure field from the initialization
 * period is applied as initial condition and for the definition of the lateral Dirichlet
 * boundary conditions. The solid  displacement field is set to zero and the CO2 injection is started.
 */
template<class TypeTag = TTAG(El2P_TestProblem)>
class El2P_TestProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, LocalFEMSpace) LocalFEMSpace;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };
    enum {
        // indices of the primary variables
            uxIdx = Indices::uxIdx,
            uyIdx = Indices::uyIdx,
            uzIdx = Indices::uzIdx,
    };
    enum {
        // indices of the equations
            momentumXEqIdx = Indices::momentumXEqIdx,
            momentumYEqIdx = Indices::momentumYEqIdx,
            momentumZEqIdx = Indices::momentumZEqIdx,
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(VertexMapper)) VertexMapper;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     * \param tInitEnd End of initialization period
     */
    El2P_TestProblem(TimeManager &timeManager,
                    const GridView &gridView)
        : ParentType(timeManager, gridView),
        gridView_(gridView)
    {

         // resize the pressure field vector with the number of vertices
        pInit_.resize(gridView.size(dim));
        // fill the pressure field vector with zeros
        std::fill( pInit_.begin(), pInit_.end(), 0.0 );
        tInitEnd_ = GET_RUNTIME_PARAM(TypeTag, Scalar,TimeManager.TInitEnd);
        this->timeManager().startNextEpisode(tInitEnd_);

        // variable which determines if output should be written (initially set to false)
        output_ = false;
        // define if current run is initialization run
        // (initially set to true, will be set to false if initialization is over)
        initializationRun_ = true;
        coupled_ = false;

        depthBOR_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Injection.DepthBOR);

        brineDensity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Fluid.BrineDensity);

        this->spatialParams().setEpisode(this->timeManager().episodeIndex());
        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLengthInit);
    }

    void init()
    {
        if (this->timeManager().time() < 1e-8)
        {
            // set the initial approximated hydrostatic pressure distribution
            // based on an averaged brine density
            // or based on a pressure polynomial
            this->initializePressure();
            // output is written
//             this->setOutput(true);

            initializationRun_ = false;
        }

        ParentType::init();
    }

    // note: pInit is < 0 (just due to geomechanics sign convention applied here)
    // initialize the pressure field for initialization run
    // first an approximate hydrostatic pressure field is calculated based on an
    // averaged density. Then the model runs for the initialization period and
    // calculates the real hydrostatic pressure distribution based on the real
    // density distribution. The calculated pressure field is than applied for
    // initialization of the actual model run and for the pressure Dirichlet boundary values.

    void initializePressure()
    {
        for(const auto& vertex : vertices(gridView_))
        {
            int vIdxGlobal = this->vertexMapper().index(vertex);
            GlobalPosition globalPos = vertex.geometry().corner(0);

            // initial approximate pressure distribution at start of initialization run
            pInit_[vIdxGlobal] = -(1.0e5 + (depthBOR_ - globalPos[dimWorld-1]) * brineDensity_ * 9.81);
        }
    }

    // allows to change the coupled_ variable which defines if geomechanical feedback on flow is taken
    // into account
    void setCoupled(bool coupled)
    {
        coupled_ = coupled;
    }

    // returns the coupled_ variable which defines if geomechanical feedback on flow is taken
    // into account
    const bool coupled() const
    {
        return coupled_;
    }

    // allows to change the output_ variable which defines if output is written
    void setOutput(bool output)
    {
        output_ = output;
    }

//     // note: pInit is < 0 (just due to geomechanics sign convention applied here)
//     // function which is called after the initialization run in
//     // order to fill the pressure field vector pInit_ with the
//     // pressure result of the initialization
//     void setPressure()
//     {
//         initializationRun_ = false; // initialization run is now finished
//
//         this->setInitializationRun(initializationRun_);
//         std::cout<<"El2P_TestProblem: initialized pressure field copied to pInit_"<<std::endl;
//         for(const auto& vertex : vertices(gridView_))
//         {
//             int vIdxGlobal = this->vertexMapper().index(vertex);
//             pInit_[vIdxGlobal] = -this->model().curSol().base()[vIdxGlobal*2][0];
//         }
//     }

    // returns the initializationRun_ variable which defines if this is an initialization run
    bool initializationRun()
    {
        return initializationRun_;
    }

    // allows to set the initializationRun_ variable which defines if this is an initialization run
    void setInitializationRun(bool initializationRun)
    {
        initializationRun_ = initializationRun;
    }

    // function which returns an in-situ stress field that needs to be provided
    // for the principal stress calculation
    GlobalPosition initialStress(const GlobalPosition globalPos, const int dofIdxGlobal) const
    {
      GlobalPosition stress;
      Scalar porosity, rockDensity, gravity;
      gravity = -this->gravity()[dimWorld-1];
      porosity = this->spatialParams().porosity(globalPos);
      rockDensity = this->spatialParams().rockDensity(globalPos);

      // initial total stress field here assumed to be isotropic, lithostatic
      stress[0] = 0.6 * ( brineDensity_ * porosity * gravity * (depthBOR_ - globalPos[dimWorld-1])
                  + (1 - porosity) * rockDensity * gravity * (depthBOR_ - globalPos[dimWorld-1]) );
      if(dimWorld >=2)
      stress[1] = 1.0 * ( brineDensity_ * porosity * gravity * (depthBOR_ - globalPos[dimWorld-1])
                  + (1 - porosity) * rockDensity * gravity * (depthBOR_ - globalPos[dimWorld-1]) );
      if(dimWorld == 3)
      stress[2] = 1.0 * ( brineDensity_ * porosity * gravity * (depthBOR_ - globalPos[dimWorld-1])
                  + (1 - porosity) * rockDensity * gravity * (depthBOR_ - globalPos[dimWorld-1]) );

      return stress;
    }


    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    {
        return "el2pRutqvist";
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius at the ground surface
     * and a geothermal gradient of 0.03 K/m.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        Scalar T;
        T = 283.15 + (depthBOR_ - globalPos[dim-1]) * 0.025;

        return T;
    };


    // returns the bottom of reservoir value (depth in m)
    const Scalar depthBOR() const
    {
        return depthBOR_;
    }

    // note: pInit is < 0 (just due to geomechanics sign convention applied here)
    // function which returns the initialized pressure at an arbitrary location within the element
    // called from finite element method (el2plocaloperator.hh) and evaluated at Gauss points
    Scalar pInit(const GlobalPosition& globalPos, const GlobalPosition& localPos, const Element& element) const
    {
        Scalar pValue = 0.0;

        typename El2P_TestProblem<TypeTag>::LocalFEMSpace feMap(this->gridView());
        const typename LocalFEMSpace::Traits::FiniteElementType
        &localFiniteElement = feMap.find(element.geometry().type());
        typedef Dune::FieldVector<CoordScalar, 1> ShapeValue;
        std::vector<ShapeValue> shapeVal;
        localFiniteElement.localBasis().evaluateFunction(localPos, shapeVal);

        for (int i = 0; i < element.subEntities(dim); i++)
        {
            int vIdxGlobal = this->vertexMapper().subIndex(element, i, dim);
            pValue += pInit_[vIdxGlobal] * shapeVal[i];
        }

        return pValue;
    }

//     // note: pInit is < 0
//     // function which returns initial pressure distribution
//     std::vector<Scalar> pInit()
//     {
//         return pInit_;
//     }

    //function which returns brine density
    const Scalar brineDensity() const
        {
         return brineDensity_;
        }

    // returns true if the current solution should be written to
    // disk (i.e. as a VTK file)
    // during initialization no output is written
    // during actual simulation output is written initially and
    // at episode/simulation end
    bool shouldWriteOutput()
    {
        return output_;
    }

    // returns true if the current solution should be written to
    // disk (i.e. as a drs file)
    bool shouldWriteRestartFile() const
    {
        return output_;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The global position
     */
void boundaryTypesAtPos(BoundaryTypes &values, const GlobalPosition& globalPos) const
    {
        values.setAllNeumann();

        // The solid displacement on the left is fixed in x and y.
        if(globalPos[0] < eps_)
        {
            values.setDirichlet(Indices::u(0));
        }

        // The solid displacement at the bottom is fixed in y.
        // The pressure is set to the initial pressure.
        if(globalPos[1] < eps_)
        {
            values.setDirichlet(Indices::u(1));
//             values.setDirichlet(pressureIdx, contiWEqIdx);
//             values.setDirichlet(saturationIdx, contiNEqIdx);
        }

        // The pressure on top is set to the initial pressure.
        // The solid displacement is fixed in y
        if(globalPos[1] > this->bBoxMax()[1]-eps_)
        {
//             values.setDirichlet(pressureIdx, contiWEqIdx);
//             values.setDirichlet(saturationIdx, contiNEqIdx);
        }

//         // The pressure on the front is fixed in y
//         if(globalPos[1] < eps_)
//         {
//             values.setDirichlet(Indices::u(1));
//         }
//
//         // The pressure on the back is fixed in y
//         if(globalPos[1] > this->bBoxMax()[1]-eps_)
//         {
//             values.setDirichlet(Indices::u(1));
//         }
     }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex representing the "half volume on the boundary"
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        dirichletAtPos(values, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * This function is called directly from dumux/geomechanics/el2p/localoperator.hh
     * If it is renamed to dirichletAtPos it should be adjusted there as well.
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
//         if(initializationRun_ == true)
//         {
            values = 0.0;
//         }
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations [kg / (m^2 *s )]
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * This function is called directly from dumux/geomechanics/el2p/localoperator.hh
     * If it is renamed to neumannAtPos it should be adjusted there as well.
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
//         values[pressureIdx] = 0;
//         values[saturationIdx] = 0;
        values = 0.0;
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
     * For this method, the \a priVars parameter stores the rate momentum
     * is generated or annihilate per volume
     * unit. Positive values mean that momentum is created, negative ones
     * mean that it vanishes.
     */
    void sourceAtPos(PrimaryVariables &priVars,
                     const GlobalPosition &globalPos) const
    {
        priVars = Scalar(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param values The initial values for the primary variables
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    void preTimeStep()
    {
//         this->spatialParams().setEpisode(this->timeManager().episodeIndex());
    }

    void postTimeStep()
    {
//         this->spatialParams().setEpisode(this->timeManager().episodeIndex());
    }

    void advanceTimeLevel()
    {
//         // function in dumux/implicit/model:
//         // make the current solution the previous one.
//         uPrev_ = uCur_;
//         if (isBox)
//             prevHints_ = curHints_;
//
//         updatePrevHints();
        // is replaced by an empty function here to avoid uPrev_ during iterations
    }

        /*!
     * \brief Define length of next episode
     */
    void episodeEnd()
    {
//         Scalar oldTimeStep = this->timeManager().timeStepSize();
//         //calls suggestTimeStepSize function, which returns the new suggested TimeStepSize for no failure
//         //and the failureTimeStepSize for anyFailure = true
//         double newTimeStepSize = this->newtonController().suggestTimeStepSize(oldTimeStep);
//         std::cout << "elastic: newTimeStepSize is " << newTimeStepSize << "\n";
        if( (this->timeManager().time() + this->timeManager().timeStepSize() > tInitEnd_ + eps_)
            && ( this->timeManager().episodeIndex() < 2) )
        {
            episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLengthMainSimulation);
        }

        if( this->timeManager().episodeIndex() >= 2)
        {
            episodeLength_ = episodeLength_ * 2.0;
        }
        if( this->timeManager().episodeIndex() == 10)
            episodeLength_ = 389.0;
        if( this->timeManager().episodeIndex() == 12)
            episodeLength_ = 122.0;
        if( this->timeManager().episodeIndex() == 16)
            episodeLength_ = 92.0;
        if( this->timeManager().episodeIndex() == 21)
            episodeLength_ = 840;
        if( this->timeManager().episodeIndex() == 23)
            episodeLength_ = 1820;
        if( this->timeManager().episodeIndex() == 24)
            episodeLength_ = 100;
    }

    /*!
       * \brief Get the non-wetting saturation of an element in the transport problem.
    */
    Scalar getpw(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return pwVector_[eIdx][scvIdx];
    }

    /*!
     * \brief Get the effective pressure vector of the transport problem.
     */
    std::vector<std::vector<Scalar>> &setpw()
    { return pwVector_; }

    /*!
       * \brief Get the capillary pressure of an element in the transport problem.
    */

    Scalar getpc(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return pcVector_[eIdx][scvIdx];
    }

    /*!
     * \brief Get the effective pressure vector of the transport problem.
     */
    std::vector<std::vector<Scalar>> &setpc()
    { return pcVector_; }

    /*!
       * \brief Get the non-wetting saturation of an element in the transport problem.
    */
    Scalar getpn(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return pnVector_[eIdx][scvIdx];
    }

    /*!
     * \brief Get the effective pressure vector of the transport problem.
     */
    std::vector<std::vector<Scalar>> &setpn()
    { return pnVector_; }

    /*!
       * \brief Get the non-wetting saturation of an element in the transport problem.
    */
    Scalar getSw(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return SwVector_[eIdx][scvIdx];
    }

    /*!
     * \brief Get the effective pressure vector of the transport problem.
     */
    std::vector<std::vector<Scalar>> &setSw()
    { return SwVector_; }

        /*!
       * \brief Get the non-wetting saturation of an element in the transport problem.
    */
    Scalar getSn(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return SnVector_[eIdx][scvIdx];
    }

    /*!
     * \brief Get the non-wetting saturation vector of the transport problem.
     */
    std::vector<std::vector<Scalar>> &setSn()
    { return SnVector_; }
        /*!
       * \brief Get the non-wetting density of an element in the transport problem.
    */
    Scalar getRhon(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return rhonVector_[eIdx][scvIdx];
    }

    /*!
     * \brief Get the non-wetting density vector of the transport problem.
     */
    std::vector<std::vector<Scalar>> &setRhon()
    { return rhonVector_; }

   /*!
    * \brief Get the wetting density of an element in the transport problem.
    */
    Scalar getRhow(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return rhowVector_[eIdx][scvIdx];
    }

    /*!
     * \brief Get the wetting density vector of the transport problem.
     */
    std::vector<std::vector<Scalar>> &setRhow()
    { return rhowVector_; }

   /*!
    * \brief Get the effective Porosity of an element in the transport problem.
    */
    Scalar getEffPorosity(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        if(this->timeManager().time() + this->timeManager().timeStepSize() < tInitEnd_ + eps_)
        {
            return this->spatialParams().porosity(element, fvGeometry, scvIdx);
        }
        else
        {
            int eIdx = this->model().elementMapper().index(element);
            return effPorosityVector_[eIdx][scvIdx];
        }
    }

    /*!
     * \brief Set the effective permeability vector of the transport problem.
     */
    std::vector<std::vector<Scalar>> &setEffPorosity()
    { return effPorosityVector_; }

    /*!
       * \brief Get the effective Porosity of an element for the last timestep.
    */
    Scalar getEffPorosityOldTimestep(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return effPorosityVectorOldTimestep_[eIdx][scvIdx];
    }

    /*!
     * \brief Set the effective Porosity of an element for the last timestep.
     */
    std::vector<std::vector<Scalar>> &setEffPorosityOldTimestep()
    { return effPorosityVectorOldTimestep_; }

    Scalar getEffPermeability(const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        int eIdx = this->model().elementMapper().index(element);
        return effPermeabilityVector_[eIdx][scvIdx];
    }

    /*!
     * \brief Get the effective Permeability vector of the transport problem.
     */
    std::vector<std::vector<Scalar>> &setEffPermeability()
    { return effPermeabilityVector_; }

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &priVars,
                  const GlobalPosition &globalPos) const
    {
        priVars = 0.0; // initial condition for the solid displacement
    }

    static constexpr Scalar eps_ = 3e-6;
    Scalar depthBOR_, brineDensity_;
    Scalar episodeLength_;
    GridView gridView_;

    std::vector<Scalar> effPressureVector_;
    std::vector<Scalar> pInit_;

    std::vector<std::vector<Scalar>> pwVector_, pnVector_, pcVector_, SwVector_, SnVector_, rhonVector_, rhowVector_;
    std::vector<std::vector<Scalar>> effPorosityVector_, effPermeabilityVector_;
    std::vector<std::vector<Scalar>> effPorosityVectorOldIteration_;
    std::vector<std::vector<Scalar>> effPorosityVectorOldTimestep_;
    Scalar tInitEnd_;
public:
    bool initializationRun_, coupled_, output_;
};
} //end namespace

#endif
