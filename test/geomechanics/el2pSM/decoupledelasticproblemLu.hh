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

#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/geomechanics/el2p/decoupled/model.hh>
#include <dumux/geomechanics/el2p/amgbackend.hh>

#include "decoupled2pspatialparamsLu.hh"

namespace Dumux
{
template<class TypeTag>
class El2P_TestProblem;


// initial conditions for momentum balance equation
template<class TypeTag, int dim>
class InitialDisplacement;

namespace Properties {
NEW_TYPE_TAG(El2P_TestProblem, INHERITS_FROM(BoxModel, BoxElasticTwoP, TwoPSpatialParams));
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
    typedef H2OAirFluidSystem<TypeTag> type;
};

// Set the soil properties
SET_TYPE_PROP(El2P_TestProblem, SpatialParams, TwoPSpatialParams<TypeTag>);

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


SET_SCALAR_PROP(El2P_TestProblem, NewtonMaxRelativeShift, 1e-6);
SET_SCALAR_PROP(El2P_TestProblem, NewtonMaxSteps, 30);
// SET_SCALAR_PROP(El2P_TestProblem, NewtonUseLineSearch, true);

// use the algebraic multigrid
SET_TYPE_PROP(El2P_TestProblem, LinearSolver, El2PAMGBackend<TypeTag>);
// SET_TYPE_PROP(El2P_TestProblem, LinearSolver, SuperLUBackend<TypeTag> );

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
 * Apart from the pressure and the saturation distribution this problems solves for the changes in
 * solid displacement (ux, uy, uz [m]) due to injection. Based on the solid displacement vector
 * the injection-induced changes in the strain and stress tensors are evaluated.
 * Further the porosity and permeability are functions of the solid displacement.
 *
 * During an initialization period of length tInit [s] the pressure field is initialized.
 *
 * After the initialization the real simulation starts and the pressure field from the initialization
 * period is applied as initial condition and for the definition of the lateral Dirichlet
 * boundary conditions.
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

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;//khodam baraie pc
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;//khodam baraie pc

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

//         depthBOR_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Injection.DepthBOR);

//         waterDensity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Fluid.waterDensity);

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

    Scalar zMax(const GlobalPosition &globalPos) const
    {
      if (globalPos[0]<6.3+eps_)
          return 0;
      else if (globalPos[0]>6.3+eps_ && globalPos[0]<23.6+eps_)
          return (-10.-0.)/(23.6-6.3)*(globalPos[0])+ 3.64 +eps_;
      else if (globalPos[0]>23.6-eps_)
          return -10;
    }

        Scalar zMin(const GlobalPosition &globalPos) const
    {
          return -40.0;
    }

    Scalar depth(const GlobalPosition &globalPos) const
    {
        return  zMax(globalPos) - globalPos[1];//distance from the surface
    }

//     Scalar waterTable = -15;

    Scalar wTdepth(const GlobalPosition &globalPos) const
    {
       return (-15. - globalPos[1]);
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
//             pInit_[vIdxGlobal] = -(1.0e5 + (depthBOR_ - globalPos[dimWorld-1]) * waterDensity_ * 9.81);
            pInit_[vIdxGlobal] = -(1.0e5 + wTdepth(globalPos) * waterDensity_ * 9.81);
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
      porosity = this->spatialParams().porosityAtPos(globalPos);
      rockDensity = this->spatialParams().rockDensity(globalPos);

      //khodam az problem
       const auto& materialLawParams = this->spatialParams().materialLawParams(globalPos);
          const Scalar swr = materialLawParams.swr();
          const Scalar snr = materialLawParams.snr();


          const Scalar meterUeberGW = globalPos[1] +15;
          const Scalar pc = std::max(0.0, 9.81*1000.0*meterUeberGW);

          Scalar n_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.n);//khodam
          Scalar m_ = 1.0 - (1.0 / n_);//khodam
          Scalar alpha_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.alpha);//khodam

          //const Scalar sw = std::min(1.0-snr, std::max(swr, invertPcGW_(pc, materialLawParams)));// khodam: numerical solution

//           const Scalar sw = pow(pow(alpha_*pc, n_) + 1, -m_);//Khodam: numerical solution from vanGenuchten.hh
          const Scalar sw = std::min(1.0-snr, std::max(swr, pow(pow(alpha_*pc, n_) + 1, -m_)));//khodam: use the numerical solution for sw

      // initial total stress field here assumed to be isotropic, lithostatic

       stress[0] = 0.5 * sw * ( waterDensity_ * porosity * gravity * std::abs((zMax(globalPos) - globalPos[dimWorld-1])) + (1 - porosity) * rockDensity * gravity * std::abs((zMax(globalPos) - globalPos[dimWorld-1])) ); //khodam: (it is assumed that generated stress in X direction as a result of gravity, is 50% of the stress in Y direction

       if(dimWorld >=2)
       stress[1] =  sw * ( waterDensity_ * porosity * gravity * std::abs((zMax(globalPos) - globalPos[dimWorld-1])) + (1 - porosity) * rockDensity * gravity * std::abs((zMax(globalPos) - globalPos[dimWorld-1])) );

       if(dimWorld == 3)
       stress[2] = sw * ( waterDensity_ * porosity * gravity * std::abs((zMax(globalPos) - globalPos[dimWorld-1])) + (1 - porosity) * rockDensity * gravity * std::abs((zMax(globalPos) - globalPos[dimWorld-1])) );
      return stress;
    }

// //Numerical solution to calculate Pc
//         static Scalar invertPcGW_(const Scalar pcIn,
//                               const MaterialLawParams &pcParams)
//     {
//         Scalar lower(0.0);
//         Scalar upper(1.0);
//         const unsigned int maxIterations = 25;
//         const Scalar bisLimit = 1.0;
//
//         Scalar sw, pcGW;
//         for (unsigned int k = 1; k <= maxIterations; k++)
//         {
//             sw = 0.5*(upper + lower);
//             pcGW = MaterialLaw::pc(pcParams, sw);
//             const Scalar delta = std::abs(pcGW - pcIn);
//             if (delta < bisLimit)
//                 return sw;
//
//             if (k == maxIterations)
//                 return sw;
//
//             if (pcGW > pcIn)
//                 lower = sw;
//             else
//                 upper = sw;
//         }
//         return sw;
//     }


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
        //T = 298.15 /*+ (depthBOR_ - globalPos[dim-1]) * 0.025*/;
        T = 283.15;

        return T;
    };


    // returns the bottom of reservoir value (depth in m)
//     const Scalar depthBOR() const
//     {
//         return depthBOR_;
//     }

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
         return waterDensity_;
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

        // The solid displacement normal to the lateral boundaries is fixed.
        if(onInlet_(globalPos))
        {
            if (initializationRun_ == true)
            {
                values.setDirichlet(momentumXEqIdx);
                values.setDirichlet(momentumYEqIdx);
            }

          }

        if (leftBoundaryO_(globalPos)|| rightBoundaryO_(globalPos)){
            values.setDirichlet(momentumXEqIdx);//this part is aplied all the time

            if (initializationRun_ == true)
            {
                values.setDirichlet(momentumYEqIdx);
            }
          }

          if (leftBoundaryO_(globalPos)|| rightBoundaryO_(globalPos)){
            values.setDirichlet(momentumXEqIdx);//this part is aplied all the time
            //values.setDirichlet(momentumYEqIdx);
            if (initializationRun_ == true)
            {
                values.setDirichlet(momentumYEqIdx);
            }
          }

        if (leftBoundaryU_(globalPos)){
            values.setDirichlet(momentumXEqIdx);//this part is aplied all the time

            if (initializationRun_ == true)
            {
                values.setDirichlet(momentumYEqIdx);
            }
          }

        if ( rightBoundaryU_(globalPos)){
            values.setDirichlet(momentumXEqIdx);//this part is aplied all the time

            if (initializationRun_ == true)
            {
                values.setDirichlet(momentumYEqIdx);
            }
         }

        if (onLowerBoundary_(globalPos)){

             values.setDirichlet(momentumYEqIdx);

            if(initializationRun_ == true)
            {
                values.setDirichlet(momentumXEqIdx);
           }
        }
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
        this->spatialParams().setEpisode(this->timeManager().episodeIndex());
        this->timeManager().startNextEpisode(episodeLength_);
        std::cout << "el: episodeLength_ is " << episodeLength_ << "\n";
        this->timeManager().setTimeStepSize(episodeLength_);
        std::cout << "el: TimeStepSize_ " << this->timeManager().timeStepSize() << "\n";
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
        std::cout << "Episode control is excerted by the main problem coupling the subproblems" << std::endl;
    }

    void setEpisodeLength(Scalar episodeLength)
    {episodeLength_ = episodeLength; }

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
        Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();
        if((tInitEnd_ - time) < eps_)
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
    Scalar waterDensity_ = 1000;
    Scalar episodeLength_;
    GridView gridView_;

    std::vector<Scalar> effPressureVector_;
    std::vector<Scalar> pInit_;

    std::vector<std::vector<Scalar>> pwVector_, pnVector_, pcVector_, SwVector_, SnVector_, rhonVector_, rhowVector_;
    std::vector<std::vector<Scalar>> effPorosityVector_, effPermeabilityVector_;
    std::vector<std::vector<Scalar>> effPorosityVectorOldIteration_;
    std::vector<std::vector<Scalar>> effPorosityVectorOldTimestep_;
    Scalar tInitEnd_;
    bool onInlet_(const GlobalPosition &globalPos) const //is that small part of the upper boundary
    {
        return (globalPos[1]>zMax(globalPos)-eps_);
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return  (globalPos[1] < zMin(globalPos) + eps_);
    }

    bool leftBoundaryO_(const GlobalPosition &globalPos) const //is that small part of the upper boundary
    {
        return (globalPos[0] < -32. + eps_ && globalPos[1] > -15. - eps_);
    }

    bool leftBoundaryU_(const GlobalPosition &globalPos) const //is that small part of the upper boundary
    {
        return (globalPos[0] < -32. + eps_ && globalPos[1] < -15. + eps_);
    }

    bool rightBoundaryO_(const GlobalPosition &globalPos) const //is that small part of the upper boundary
    {
        return (globalPos[0] > 35. - eps_ && globalPos[1] > -15. - eps_);
    }

    bool rightBoundaryU_(const GlobalPosition &globalPos) const //is that small part of the upper boundary
    {
        return (globalPos[0] > 35. - eps_ && globalPos[1] < -15. + eps_);
    }
public:
    bool initializationRun_, coupled_, output_;
};
} //end namespace

#endif
