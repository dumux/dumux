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
// OLD #include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/components/simpleh2orelativepressure.hh>


// #include <dumux/material/fluidsystems/brineco2.hh>
// OLD #include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidsystems/h2oairrelativepressure.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/geomechanics/el2p/model.hh>
#include <dumux/geomechanics/el2p/amgbackend.hh>
#include <dumux/geomechanics/el2p/elementvolumevariables.hh>//Sh


#include "el2pspatialparams.hh"


#include <dune/common/version.hh>

namespace Dumux
{
template<class TypeTag>
class El2P_TestProblem;


// initial conditions for momentum balance equation
template<class TypeTag, int dim>
class InitialDisplacement;

// initial conditions for mass balance equations
template<class TypeTag>
class InitialPressSat;

namespace Properties {
NEW_TYPE_TAG(El2P_TestProblem, INHERITS_FROM(BoxModel, BoxElasticTwoP, El2PSpatialParams));
NEW_PROP_TAG(InitialDisplacement); //!< The initial displacement function
NEW_PROP_TAG(InitialPressSat); //!< The initial pressure and saturation function

// Set the grid type
// SET_TYPE_PROP(El2P_TestProblem, Grid, Dune::YaspGrid<3>);
SET_TYPE_PROP(El2P_TestProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);
// SET_TYPE_PROP(El2P_TestProblem, Grid, Dune :: UGGrid<3>);

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

//// Set fluid configuration
//SET_PROP(El2P_TestProblem, FluidSystem)
//{
////     typedef BrineCO2FluidSystem<TypeTag> type;
       //typedef Dumux::H2OAirRelPressFluidSystem<TypeTag> type;
//
//};

// Set the wetting phase
SET_PROP(El2P_TestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, SimpleH2ORelPress<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(El2P_TestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::GasPhase<Scalar, AirRelPress<Scalar> > type;
};


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

// Set the initial pressure and saturation function
SET_PROP(El2P_TestProblem, InitialPressSat)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
public:
    typedef InitialPressSat<TypeTag> type;
};

SET_SCALAR_PROP(El2P_TestProblem, NewtonMaxRelativeShift, 1e-1);//!!!Holger: this is added to adjust the first relative shift to make the initial step convergable, then we activete the Restart in input file and reduce this relative shift

// use the algebraic multigrid
// SET_TYPE_PROP(El2P_TestProblem, LinearSolver, El2PAMGBackend<TypeTag>);//didnt converge for fine mesh
SET_TYPE_PROP(El2P_TestProblem, LinearSolver, SuperLUBackend<TypeTag>);
// SET_TYPE_PROP(El2P_TestProblem, LinearSolver, UMFPackBackend<TypeTag> );//didnt converge for fine mesh

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
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;//khodam
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;//khodam
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;//khodam

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    enum {
        // indices of the primary variables
            pressureIdx = Indices::pwIdx,
            saturationIdx = Indices::snIdx,
            uxIdx = Indices::uxIdx,
            uyIdx = Indices::uyIdx,
            uzIdx = Indices::uzIdx

    };
    enum {
        // indices of the equations+
            contiWEqIdx = Indices::contiWEqIdx,
            contiNEqIdx = Indices::contiNEqIdx
    };
    enum {
        // indices of the phases
            wPhaseIdx = Indices::wPhaseIdx,
            nPhaseIdx = Indices::nPhaseIdx
    };


    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::BlockVector<GlobalPosition> InitialStressField;

    typedef typename GET_PROP_TYPE(TypeTag, LocalFEMSpace) LocalFEMSpace;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;

//     typedef typename GET_PROP_TYPE(TypeTag, CO2Table) CO2Table;
//     typedef Dumux::CO2<Scalar, CO2Table> CO2;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    El2P_TestProblem(TimeManager &timeManager,
                    const GridView &gridView)
        : ParentType(timeManager, gridView),
        gridView_(gridView)
    {
        GridCreator::grid().globalRefine(GET_RUNTIME_PARAM(TypeTag, Scalar,Grid.Refine));

        std::cout << "El2P_TestProblem: Initializing the fluid system for the el2p model\n";

        // initialize the tables of the fluid system
         //FluidSystem::init(/*Tmin=*/273,
                           ///*Tmax=*/300,
                           ///*nT=*/5,
                           ///*pmin=*/-1e5,
                           ///*pmax=*/1e6,
                           ///*np=*/10);

        // resize the pressure field vector with the number of vertices
        pInit_.resize(gridView.size(dim));

        // fill the pressure field vector with zeros
        std::fill( pInit_.begin(), pInit_.end(), 0.0 );

        // variable which determines if output should be written (initially set to false)
        output_ = false;
        // define if current run is initialization run
        // (initially set to true, will be set to false if initialization is over)
        initializationRun_ = true;
        // defines if feedback from geomechanics on flow is taken into account or not
        // (usually the coupling is switched off for the initialization run)
        coupled_ = false;
        // set initial episode length equal to length of initialization period
        Scalar tInitEnd = GET_RUNTIME_PARAM(TypeTag, Scalar,TimeManager.TInitEnd);
        this->timeManager().startNextEpisode(tInitEnd);
        // transfer the episode index to spatial parameters
        // (during intialization episode hydraulic different parameters might be applied)
        this->spatialParams().setEpisode(this->timeManager().episodeIndex());

//         depthBOR_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Injection.DepthBOR);
        episodeLength_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TimeManager.EpisodeLength);

        dt_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TimeManager.DtInitial);
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
            this->setOutput(true);
        }

        ParentType::init();
    }

         //instead of bBoxMax which was constant in the whole domain, we define the max elevation or the ground surface at each point
    Scalar zMax(const GlobalPosition &globalPos) const
    {
      if (globalPos[0]<6.3+eps_)
          return 15.;
      else if (globalPos[0]>6.3-eps_ && globalPos[0]<23.6+eps_)
          return (5.-15.)/(23.6-6.3)*(globalPos[0])+ (322./17.3);
      else if (globalPos[0]>23.6-eps_)
          return 5.;
    }

        Scalar zMin(const GlobalPosition &globalPos) const
    {
          return -25.0;
    }

    Scalar depth(const GlobalPosition &globalPos) const
    {
        return  zMax(globalPos) - globalPos[1];//distance from the surface
    }

//     Scalar waterTable = -15;

    Scalar wTdepth(const GlobalPosition &globalPos) const
    {
       return (0. - globalPos[1]);
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
            // pInit_[vIdxGlobal] = -(1.0e5 + wTdepth(globalPos) * waterDensity_ * 9.81);
            pInit_[vIdxGlobal] = -wTdepth(globalPos) * waterDensity_ * 9.81;
        }
    }

    // allows to change the coupled_ variable which defines if geomechanical feedback on flow is taken
    // into account
    void setCoupled(bool coupled)
    {
        coupled_ = coupled;
        std::cout << "coupled_ is now set to " << coupled << std::endl ;
    }

    // returns the coupled_ variable which defines if geomechanical feedback on flow is taken
    // into account
    bool coupled() const
    {
        return coupled_;
    }

    // allows to change the output_ variable which defines if output is written
    void setOutput(bool output)
    {
        output_ = output;
    }

    // note: pInit is < 0 (just due to geomechanics sign convention applied here)
    // function which is called after the initialization run in
    // order to fill the pressure field vector pInit_ with the
    // pressure result of the initialization
    void setPressure()
    {
        initializationRun_ = false; // initialization run is now finished

        this->setInitializationRun(initializationRun_);
        std::cout<<"El2P_TestProblem: initialized pressure field copied to pInit_"<<std::endl;
        for(const auto& vertex : vertices(gridView_))
        {
            int vIdxGlobal = this->vertexMapper().index(vertex);
            pInit_[vIdxGlobal] = -this->model().curSol().base()[vIdxGlobal*2][0];
        }
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
      porosity = this->spatialParams().porosity(globalPos);
      rockDensity = this->spatialParams().rockDensity(globalPos);

      //khodam az problem
       const auto& materialLawParams = this->spatialParams().materialLawParams(globalPos);
          const Scalar swr = materialLawParams.swr();
          const Scalar snr = materialLawParams.snr();


          const Scalar meterUeberGW = globalPos[1] +0;
          const Scalar pc = std::max(0.0, 9.81*1000.0*meterUeberGW);

          Scalar n_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.n);//khodam
          Scalar m_ = 1.0 - (1.0 / n_);//khodam
          Scalar alpha_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.alpha);//khodam

          const Scalar sw = std::min(1.0-snr, std::max(swr, pow(pow(alpha_*pc, n_) + 1, -m_)));//khodam: use the numerical solution for sw

      // initial total stress field here assumed to be isotropic, lithostatic
          //Asked Martin: Pref_ shouldnt be added because this is just stress and Pref is considered in Peff later on.


//         Scalar nu_ = GET_RUNTIME_PARAM(TypeTag, Scalar, FailureParameters.nu);
//         Scalar k0_ = nu_ / (1 - nu_);//khodam sigma h/sigma v = nu/(1-nu) and nu shouldnt be more than 0.5

//        stress[0] = k0_ * sw * ( waterDensity_ * porosity * gravity * std::abs((zMax(globalPos) - globalPos[dimWorld-1])) + (1 - porosity) * rockDensity * gravity * std::abs((zMax(globalPos) - globalPos[dimWorld-1])) ); //khodam: (it is assumed that generated stress in X direction as a result of gravity, is 50% of the stress in Y direction

       stress[0] = 0.5 *( sw * ( waterDensity_ * porosity * gravity * std::abs((zMax(globalPos) - globalPos[dimWorld-1]))) + (1 - porosity) * rockDensity * gravity * std::abs((zMax(globalPos) - globalPos[dimWorld-1])) ); //khodam: (it is assumed that generated stress in X direction as a result of gravity, is 50% of the stress in Y direction

       if(dimWorld >=2)
       stress[1] =  sw * ( waterDensity_ * porosity * gravity * std::abs((zMax(globalPos) - globalPos[dimWorld-1]))) + (1 - porosity) * rockDensity * gravity * std::abs((zMax(globalPos) - globalPos[dimWorld-1]));

       if(dimWorld == 3)
       stress[2] = sw * ( waterDensity_ * porosity * gravity * std::abs((zMax(globalPos) - globalPos[dimWorld-1]))) + (1 - porosity) * rockDensity * gravity * std::abs((zMax(globalPos) - globalPos[dimWorld-1]));


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
    std::string name() const
    {
        return "el2p.Comsol.test";
    }

    /*!
     * \brief Returns the temperature wcithin the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius at the ground surface
     * and a geothermal gradient of 0.03 K/m.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        Scalar T;
        T = 283.15; //+ (depthBOR_ - globalPos[dimWorld-1]) * 0.03;

        return T;
    };


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

    // note: pInit is < 0
    // function which returns initial pressure distribution
    std::vector<Scalar> pInit()
    {
        return pInit_;
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
      \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The center of the finite volume which ought to be set.
     *
     *               This function is called directly from dumux/geomechanics/el2p/localoperator.hh
     *               If it is renamed to boundaryTypesAtPos it should be adjusted there as well.
     */
    void boundaryTypesAtPos(BoundaryTypes &values, const GlobalPosition& globalPos) const
    {
        values.setAllNeumann();

        // The solid displacement normal to the lateral boundaries is fixed.
        if(onInlet_(globalPos))
        {
//             if (initializationRun_ == true)
//             {
// //                 values.setDirichlet(pressureIdx, contiWEqIdx);
// //                 values.setDirichlet(saturationIdx, contiNEqIdx);
//                 values.setDirichlet(uxIdx);
//                 values.setDirichlet(uyIdx);
//             }

          }


        if (leftBoundaryO_(globalPos)|| rightBoundaryO_(globalPos)){

            values.setDirichlet(uxIdx);//this part is aplied all the time

//             if (initializationRun_ == true)
//             {
//                 values.setDirichlet(uyIdx);
//             }
          }

        if (leftBoundaryU_(globalPos)){
            values.setDirichlet(uxIdx);//this part is aplied all the time
            //values.setDirichlet(momentumYEqIdx);
            if (initializationRun_ == true)
            {
                values.setDirichlet(pressureIdx, contiWEqIdx);
                values.setDirichlet(saturationIdx, contiNEqIdx);
//                 values.setDirichlet(uyIdx);
            }
          }


        if ( rightBoundaryU_(globalPos)){
            values.setDirichlet(pressureIdx, contiWEqIdx);
            values.setDirichlet(saturationIdx, contiNEqIdx);
            values.setDirichlet(uxIdx);//this part is aplied all the time
            if (initializationRun_ == true)
            {
//                 values.setDirichlet(uyIdx);
            }
         }


        if (onLowerBoundary_(globalPos)){ //fixed boundary

             values.setDirichlet(uxIdx);
             values.setDirichlet(uyIdx);

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
        //values[0] = -pInit_[this->vertexMapper().index(vertex)];
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
        values = 0.0;


//         if (onInlet_(globalPos)|| rightBoundaryO_(globalPos)){
//         if (onInlet_(globalPos)){
//
//           const auto& materialLawParams = this->spatialParams().materialLawParams(globalPos);
//           const Scalar swr = materialLawParams.swr();
//           const Scalar snr = materialLawParams.snr();
//
//
//           const Scalar meterUeberGW = globalPos[1] +0;
//           const Scalar pc = std::max(0.0, 9.81*1000.0*meterUeberGW);
//           const Scalar sw = std::min(1.0-snr, std::max(swr, invertPcGW_(pc, materialLawParams)));
//           values[pressureIdx] = (1.0e5 + 1000. * 9.81 * wTdepth(globalPos));
//           values[saturationIdx] = 1.-sw;
//         } else
//         {
          // OLD values[pressureIdx] = (1.0e5 + 1000. * 9.81 * wTdepth(globalPos));
          values[pressureIdx] = 1000. * 9.81 * wTdepth(globalPos);
          //std::cout<< "pressure=" <<values[pressureIdx] << std::endl;
          values[saturationIdx] = 0.0;
//         }
    }

//     static Scalar invertPcGW_(const Scalar pcIn,
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
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * \param values The neumann values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^2 \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * This function is called directly from dumux/geomechanics/el2p/localoperator.hh
     * If it is renamed to neumannAtPos it should be adjusted there as well.
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values[uxIdx] = 0.0;
        values[uyIdx] = 0.0;
    }

     void solDependentNeumann(PrimaryVariables &values,
                      const Element &element,
                      const FVElementGeometry &fvGeometry,
                      const Intersection &intersection,
                      const int scvIdx,
                      const int boundaryFaceIdx,
                      const ElementVolumeVariables &elemVolVars) const
    {

        values = 0.0;
        GlobalPosition globalPos = fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;


        Scalar avgRain_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.avgRain);//[m/s]

        Scalar rb_ = 2.0E-4; //[1/s] =K/B conductance
        Scalar ks_ = 1.39e-6;

        // OLD const Scalar pW = elemVolVars[scvIdx].pressure(pressureIdx);
        // OLD const Scalar pN = elemVolVars[scvIdx].pressure(nPhaseIdx);
        const Scalar pW = 1e5 + elemVolVars[scvIdx].pressure(pressureIdx);
        const Scalar pN = 1e5 + elemVolVars[scvIdx].pressure(nPhaseIdx);
        const Scalar satW = elemVolVars[scvIdx].saturation(wPhaseIdx);
        const Scalar satN = elemVolVars[scvIdx].saturation(nPhaseIdx);

        if (onInlet_(globalPos))
        {
            if(initializationRun_ == false)
            {
               if (satW > 1. - eps_){
                   //std::cout << "saturation above 1 \n" ;
                    values[contiWEqIdx] = rb_ * pW/(9.81); // [kg/(m2*s)]
//                     values[contiNEqIdx] = 0.0;

                    if (pN>1.e5) values[contiNEqIdx] = satN * (pN-1.e5) * rb_;
               }
                else {
                   //std::cout << "saturation below 1, infiltration \n" ;
                    values[contiWEqIdx] = -avgRain_*1000.;
                    if (pN>1.e5) { values[contiNEqIdx] = satN * (pN-1.e5) * rb_; }
                }
            }
            else
            {
                values = 0.0;
            }
        }

     }

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param values The source values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param scvIdx The local vertex index
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables &values,
            const Element &element,
            const FVElementGeometry &fvGeometry,
            int scvIdx) const
    {
        const GlobalPosition globalPos = element.geometry().corner(scvIdx);

        source(values, globalPos);
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param values The source values for the conservation equations in units of
     *                 \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = 0.0;

    }
    // \}

    /*!
     * \brief Transfer episode index to spatial parameters in order to apply different
     * hydraulic parameters during pressure initialization
     */
    void preTimeStep()
    {
        this->spatialParams().setEpisode(this->timeManager().episodeIndex());
    }

    /*!
     * \brief Write mass balance information for both fluid phases
     */
    void postTimeStep()
    {
        PrimaryVariables mass;
        this->model().globalStorage(mass);
        double time = this->timeManager().time()+this->timeManager().timeStepSize();

        if(time>1.0)
        this->newtonController().setMaxRelativeShift(1.e-5);//khodam from 1e-5

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "TIME, MASS NPhase (kg), MASS WPhase (kg): \n"
            <<"mass: "
            <<time<< " , "
            <<mass[1] << " , "
            <<mass[0]
            <<"\n"
            <<"***************************************" <<std::endl;
        }


    }

    /*!
     * \brief Define length of next episode
     */
    void episodeEnd()
    {
        this->timeManager().startNextEpisode(episodeLength_);
        // At the end of the initializationRun
        if (this->timeManager().time() == GET_RUNTIME_PARAM(TypeTag, Scalar,TimeManager.TInitEnd))
        {
            this->timeManager().setTimeStepSize(dt_);

            this->setCoupled(true);
            // pressure field resulting from the initialization period is applied for the initial
            // and the Dirichlet boundary conditions
            this->setPressure();
            // output is written
            this->setOutput(true);
        }
    }

private:
    static constexpr Scalar eps_ = 3e-6;
    Scalar depthBOR_;
    static constexpr Scalar waterDensity_ = 1000;
    Scalar episodeLength_;// = GET_RUNTIME_PARAM(TypeTag, Scalar, TimeManager.EpisodeLength);

    std::vector<Scalar> pInit_;
    GridView gridView_;
    Scalar dt_;

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        return (globalPos[1]>zMax(globalPos)-eps_);
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return  (globalPos[1] < zMin(globalPos) + eps_);
    }

    bool leftBoundaryO_(const GlobalPosition &globalPos) const
    {
        return (globalPos[0] < -32. + eps_ && globalPos[1] > 0. - eps_);
    }

    bool leftBoundaryU_(const GlobalPosition &globalPos) const
    {
        return (globalPos[0] < -32. + eps_ && globalPos[1] < 0. + eps_);
    }

    bool rightBoundaryO_(const GlobalPosition &globalPos) const
    {
        return (globalPos[0] > 35. - eps_ && globalPos[1] > 0. - eps_);
    }

    bool rightBoundaryU_(const GlobalPosition &globalPos) const
    {
        return (globalPos[0] > 35. - eps_ && globalPos[1] < 0. + eps_);
    }
public:
    bool initializationRun_, coupled_, output_;
    InitialStressField initialStressField_;
};


/*!
 * \ingroup ElTwoPBoxProblems
 *
 * \brief Initial conditions for momentum balance equation.
 *
 * Set initial conditions for solution of momentum balance equation
 * i.e. initialize solid displacement
 * This function is called from dumux/geomechanics/el2p/model.hh
 *  * This function is called from dumux/geomechanics/el2p/model.hh.
 *
 * The primary variables are initialized two times:
 * 1. before the initialization run.
 * 2. at the start of the actual simulation the solid displacement values which have
 *    changed during initialization of the pressure field are set to zero again.
 *
 */
template<class TypeTag, int dim>
class InitialDisplacement :
public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<typename GET_PROP_TYPE(TypeTag, GridView),typename GET_PROP_TYPE(TypeTag, Scalar),dim>,
    InitialDisplacement<TypeTag,dim> >
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::PDELab::AnalyticGridFunctionTraits<typename GET_PROP_TYPE(TypeTag, GridView),typename GET_PROP_TYPE(TypeTag, Scalar),dim> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, InitialDisplacement<TypeTag,dim> > BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;//khodam
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;//khodam

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    InitialDisplacement(const GridView & gridView) : BaseT(gridView) {}

    /*!
     * \brief Evaluate initial conditions for the momentum balance equation.
     *
     * \param position The position of the vertex
     * \param values The initial solid displacement vector at the vertex
     */
    inline void evaluateGlobal(const DomainType & position, RangeType & values) const
    {
        values = 0;
    }
};

/*!
 * \ingroup ElTwoPBoxProblems
 *
 * \brief Initial conditions for mass balance equations.
 *
 * Set initial conditions for solution of the mass balance equations
 * i.e. initialize wetting phase pressure and nonwetting phase saturation
 *
 * This function is called from dumux/geomechanics/el2p/model.hh.
 * The primary variables are initialized two times:
 * 1. before the initialization run.
 * 2. at the start of the actual simulation applying pressure field
 *    calculated during initialization
 *
 */
template<class TypeTag>
class InitialPressSat :
public Dune::PDELab::AnalyticGridFunctionBase<
    Dune::PDELab::AnalyticGridFunctionTraits<typename GET_PROP_TYPE(TypeTag, GridView),typename GET_PROP_TYPE(TypeTag, Scalar),2>,
    InitialPressSat<TypeTag> >
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::PDELab::AnalyticGridFunctionTraits<GridView,Scalar,2> Traits;
    typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, InitialPressSat<TypeTag>> BaseT;

    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;


    enum {
        // indices of the primary variables
        pressureIdx = Indices::pwIdx,
        saturationIdx = Indices::snIdx,

        dimWorld = GridView::dimensionworld
    };

    typedef typename Dune::MultipleCodimMultipleGeomTypeMapper<GridView,
                        Dune::MCMGVertexLayout> VertexMapper;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    InitialPressSat(const GridView & gridView)
    : BaseT(gridView)
    , gridView_(gridView)
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
    , vertexMapper_(gridView, Dune::mcmgVertexLayout())
#else
    , vertexMapper_(gridView)
#endif
    {
        // resize the pressure field vector with the number of vertices
        pInit_.resize(gridView.size(GridView::dimension));
        // fill the pressure field vector with zeros
        std::fill(pInit_.begin(), pInit_.end(), 0.0);
    }

    /*!
     * \brief Evaluate initial conditions for the mass balance equations.
     *
     * \param position The position of the vertex
     * \param values The initial pressure and saturation values at the vertex
     *
     * This function applies the pressure field pInit_ which is defined
     * in the problem.
     */

 //from: /temp/moradiC/dumux/dumux/material/fluidmatrixinteractions/2p    : vanGenuchten.hh

    inline void evaluateGlobal(const DomainType & position, RangeType & values) const
    {
            bool valueSet;
            valueSet = false;
            using std::pow;
            using std::min;
            using std::max;

            // loop over all vertices
            for (const auto& vertex : vertices(gridView_))
            {
                // get global index of current vertex
                int vIdxGlobal = vertexMapper_.index(vertex);
                Dune::FieldVector<double, dimWorld> globalPos =
                                (vertex).geometry().corner(0);

                // compare coordinates of current vertex with position coordinates
                if (globalPos[0] >= position[0] - eps_ && globalPos[0] <= position[0] + eps_
                                && globalPos[1] >= position[1] - eps_ && globalPos[1]
                                <= position[1] + eps_ && globalPos[dimWorld-1] >= position[dimWorld-1] - eps_
                                && globalPos[dimWorld-1] <= position[dimWorld-1] + eps_)
                {
                    // if coordinates are identical write the pressure value for this
                    // vertex (with index vIdxGlobal) into the values vector
                    values[pressureIdx] = pInit_[vIdxGlobal];
                    // the value of this vertex is set
                    valueSet = true;
                    //define saturation from Pc khodam

                    Scalar n_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.n);
                    Scalar m_ = 1.0 - (1.0 / n_);
                    Scalar alpha_ = GET_RUNTIME_PARAM(TypeTag, Scalar, TransportParameters.alpha);

                    const Scalar meterUeberGW_ = globalPos[1] + 0;
                    const Scalar pc = std::max(0.0, 9.81*1000.0*meterUeberGW_);//Pc>0 is unsaturated and Pc<=0 is saturated
//                     const Scalar sw = std::min(1.0-snr, std::max(swr, invertPcGW_(pc, materialLawParams)));
                    const Scalar sw = pow(pow(alpha_*pc, n_) + 1, -m_);

                    // initialize saturation values
//                     values[pressureIdx] = (1.0e5 + 1000. * 9.81 * (0. - globalPos[1]));
                    values[saturationIdx] = 1 - sw;

                }

            }

           // checkfin if the pressure value for this vertex has been initialized
            if (valueSet == false)
            {
                std::cout << " pressure value not initialized correctly "
                                << std::endl;
            }
        }


    /*!
     * \brief Fill the vector pInit_ for initialization
     *
     * \param pInit The pressure field vector defined in the problem class
     *
     * This function is called from dumux/geomechanics/el2p/model.hh.
     */
        void setPressure(std::vector<Scalar> pInit)
        {
            std::cout << "InitialPressSat: setPressure function called" << std::endl;
            for (const auto& vertex : vertices(gridView_))
            {
                int vIdxGlobal = vertexMapper_.index(vertex);
                pInit_[vIdxGlobal] = -pInit[vIdxGlobal];
            }
        }

private:
    static constexpr Scalar eps_ = 3e-6;
    Scalar depthBOR_;
    std::vector<Scalar> pInit_;
    GridView gridView_;
    VertexMapper vertexMapper_;
};

} //end namespace

#endif
                    //from: /temp/moradiC/dumux/dumux/material/fluidmatrixinteractions/2p    : vanGenuchten.hh
