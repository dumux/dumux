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
/*!
 * \file
 *
 * \brief Porous-medium subproblem with coupling at the top boundary.
 */
#ifndef DUMUX_2P2CNISUB_PROBLEM_HH
#define DUMUX_2P2CNISUB_PROBLEM_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dumux/implicit/common/implicitporousmediaproblem.hh>
#include <dumux/implicit/2p2cni/2p2cnimodel.hh>
#include <dumux/multidomain/couplinglocalresiduals/2p2cnicouplinglocalresidual.hh>
#include <dumux/multidomain/common/subdomainpropertydefaults.hh>
#include <dumux/multidomain/common/multidomainlocaloperator.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dune/pdelab/finiteelementmap/conformingconstraints.hh>

#include "2cnistokes2p2cniproblem.hh"

namespace Dumux
{
template <class TypeTag>
class TwoPTwoCNISubProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoPTwoCNISubProblem, 
	INHERITS_FROM(BoxTwoPTwoCNI, SubDomain, TwoCNIStokesTwoPTwoCNISpatialParams));

// Set the problem property
SET_TYPE_PROP(TwoPTwoCNISubProblem, Problem, TwoPTwoCNISubProblem<TTAG(TwoPTwoCNISubProblem)>);

//! Use the 2p2cni local jacobian operator for the 2p2cniCoupling model
SET_TYPE_PROP(TwoPTwoCNISubProblem, LocalResidual, TwoPTwoCNICouplingLocalResidual<TypeTag>);

// choose pn and Sw as primary variables
SET_INT_PROP(TwoPTwoCNISubProblem, Formulation, TwoPTwoCFormulation::pnsw);

// the gas component balance (air) is replaced by the total mass balance
SET_PROP(TwoPTwoCNISubProblem, ReplaceCompEqIdx)
{
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    static const int value = Indices::contiNEqIdx;
};

// Set the model parameter group for the TypeTag
//SET_PROP(TwoPTwoCNISubProblem, ModelParameterGroup)
//{ static const char *value() { return "PM"; }; };

// set the constraints for multidomain / pdelab
SET_TYPE_PROP(TwoPTwoCNISubProblem, Constraints,
        Dune::PDELab::NonoverlappingConformingDirichletConstraints);

// set the grid operator
SET_PROP(TwoPTwoCNISubProblem, GridOperator)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ConstraintsTrafo) ConstraintsTrafo;
    typedef typename GET_PROP_TYPE(TypeTag, GridFunctionSpace) GridFunctionSpace;
    typedef Dumux::PDELab::MultiDomainLocalOperator<TypeTag> LocalOperator;

    enum{numEq = GET_PROP_VALUE(TypeTag, NumEq)};

public:
    typedef Dune::PDELab::GridOperator<GridFunctionSpace,
            GridFunctionSpace, LocalOperator,
            Dune::PDELab::ISTLBCRSMatrixBackend<numEq, numEq>,
            Scalar, Scalar, Scalar,
            ConstraintsTrafo,
            ConstraintsTrafo,
            true> type;
};

// the fluidsystem is set in the coupled problem
SET_PROP(TwoPTwoCNISubProblem, FluidSystem)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag) CoupledTypeTag;
    typedef typename GET_PROP_TYPE(CoupledTypeTag, FluidSystem) FluidSystem;
public:
    typedef FluidSystem type;
};


//! Somerton is used as model to compute the effective thermal heat conductivity
SET_PROP(TwoPTwoCNISubProblem, ThermalConductivityModel)
{ private :
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
  public:
    typedef ThermalConductivitySomerton<Scalar> type;
};

// use formulation based on mass fractions
SET_BOOL_PROP(TwoPTwoCNISubProblem, UseMoles, false);

// enable/disable velocity output
SET_BOOL_PROP(TwoPTwoCNISubProblem, VtkAddVelocity, true);

// Enable gravity
SET_BOOL_PROP(TwoPTwoCNISubProblem, ProblemEnableGravity, true);
}


template <class TypeTag = TTAG(TwoPTwoCNISubProblem) >
class TwoPTwoCNISubProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GridView::Grid Grid;

    typedef TwoPTwoCNISubProblem<TypeTag> ThisType;
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    // the type tag of the coupled problem
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag) CoupledTypeTag;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { // the equation indices
        contiTotalMassIdx = Indices::contiNEqIdx,
        contiWEqIdx = Indices::contiWEqIdx,
        energyEqIdx = Indices::energyEqIdx
    };
    enum { // the indices of the primary variables
		pNIdx = Indices::pressureIdx,
        sWIdx = Indices::switchIdx,
        temperatureIdx = Indices::temperatureIdx
    };
    enum { // the indices for the phase presence
        wCompIdx = Indices::wCompIdx,
        nCompIdx = Indices::nCompIdx
    };
	enum { // the indices for the phase presence
        wPhaseOnly = Indices::wPhaseOnly,
        nPhaseOnly = Indices::nPhaseOnly,
        bothPhases = Indices::bothPhases
    };
    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };
    enum { // grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
//    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    // only required for radiation
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    //////////////

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    /*!
     * \brief docme
     *
     * \param timeManager The TimeManager which is used by the simulation
     * \param gridView The simulation's idea about physical space
     *
     */
public:
    TwoPTwoCNISubProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        try
        {
            Scalar noDarcyX = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, NoDarcyX);
            Scalar xMin = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, XMin);

            bboxMin_[0] = std::max(xMin,noDarcyX);
    		bboxMax_[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, XMax);

    		bboxMin_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, YMin);
    		bboxMax_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePos);

    		runUpDistanceX_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, RunUpDistanceX); // first part of the interface without coupling
            xMaterialInterface_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams, MaterialInterfaceX);
    		initializationTime_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, InitTime);

            refTemperature_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, RefTemperaturePM);
            refPressure_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, RefPressurePM);
            initialSw1_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, InitialSw1);
            initialSw2_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, PorousMedium, InitialSw2);

//            porosity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, SpatialParams.Medium, Porosity2);

            freqMassOutput_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqMassOutput);

            storageLastTimestep_ = Scalar(0);
            lastMassOutputTime_ = Scalar(0);

            outfile.open("evaporation.out");
            outfile << "time; massChange; evaporationRate; energyChange; massH2O" << std::endl;
        }
        catch (Dumux::ParameterException &e) {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }
        catch (...) {
            std::cerr << "Unknown exception thrown!\n";
            exit(1);
        }
    }

    ~TwoPTwoCNISubProblem()
    {
        outfile.close();
    }

    /*!
     * \name Problem parameters
     */
    // \{

    //! \copydoc Dumux::ImplicitProblem::name()
    const std::string &name() const
    { return GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Vtk, NamePM); }

    //! \copydoc Dumux::ImplicitProblem::init()
    void init()
    {
        // set the initial condition of the model
        ParentType::init();

        this->model().globalStorage(storageLastTimestep_);

        const int blModel = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, FreeFlow, UseBoundaryLayerModel);
        std::cout << "Using boundary layer model " << blModel << std::endl;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    //! \copydoc Dumux::ImplicitProblem::boundaryTypes()
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();
        Scalar time = this->timeManager().time();

        values.setAllNeumann();

        if (onLowerBoundary_(globalPos))
        {
//        	values.setAllDirichlet();
//        	values.setDirichlet(pressureIdx, contiTotalMassIdx);
        	values.setDirichlet(temperatureIdx, energyEqIdx);
        }

        if (onUpperBoundary_(globalPos))
        {
            if (time > initializationTime_)
            {
                if (globalPos[0] > runUpDistanceX_ - eps_)
                {
                    values.setAllCouplingInflow();
                }
                else
                    values.setAllNeumann();
            }
            else
            {
//                values.setAllDirichlet();
                values.setAllNeumann();
            }
        }
    }

    //! \copydoc Dumux::ImplicitProblem::dirichlet()
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        initial_(values, globalPos);
    }

    //! \copydoc Dumux::ImplicitProblem::neumann()
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &is,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        values = 0.;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{
    //! \copydoc Dumux::ImplicitProblem::solDependentSource()
    void solDependentSource(PrimaryVariables &values,
            		 const Element &element,
            		 const FVElementGeometry &fvGeometry,
            		 const int scvIdx,
            		 const ElementVolumeVariables &elemVolVars) const
    {
        values = 0;
//        const Scalar time = this->timeManager().time();
//        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);
//
//		const Scalar temperatureOut = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefTemperature);
//
////		const Scalar S = variationSun_(0);                       //Solar irradiance (Yamanaka 1998)
//		const Scalar sigma = 5.671e-8;                             //Boltzmann-Constant
//		const Scalar e_a = atmosphericEmissivity();                //atmosphere emissivity (Brutsaert 1982)
//		const Scalar e_s = surfaceEmissivity_(elemVolVars, scvIdx);
//		const Scalar a = albedo_(elemVolVars, scvIdx);
//		const Scalar T_a = 298.15;//variationT_a_(6);                       //atmosphere Temperature (Yamanaka 1998)
//		const Scalar T_s = elemVolVars[scvIdx].temperature();      //surface Temperature
//
//		const Scalar vertexPosition = 0.249248;//1.0; //GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InterfacePos);
//
//          //For d depth with b cells!!
////		Scalar radiation = S*(1.0-a)+
//		Scalar radiation = sigma*e_s*(e_a*std::pow(T_a,4.0)-(std::pow(T_s,4.0)));  //Equation from Novak 2010
//
//		if (time > initializationTime_
//				&& globalPos[1] > (vertexPosition - 0.00002)
//				&& globalPos[1] < (vertexPosition + 0.00002))
//		  values[energyEqIdx] += radiation/cellHeight_();
    }


    //! \copydoc Dumux::ImplicitProblem::initial()
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const int scvIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);

        values = 0.;

        initial_(values, globalPos);
    }

    /*!
     * \brief Return the initial phase state inside a control volume.
     *
     * \param vert The vertex
     * \param globalIdx The index of the global vertex
     * \param globalPos The global position
     */
    int initialPhasePresence(const Vertex &vert,
                             const int &globalIdx,
                             const GlobalPosition &globalPos) const
    {
        return bothPhases;
    }


//    /*!
//     * \brief docme
//     */
//    Scalar globalRadiation()
//    {
//        FVElementGeometry fvGeometry;
//        ElementVolumeVariables elemVolVars;
//        Scalar source = 0;
//        PrimaryVariables values(0);
//
//        ElementIterator elemIt = this->gridView().template begin<0>();
//        const ElementIterator elemEndIt = this->gridView().template end<0>();
//        for (elemIt = this->gridView().template begin<0>(); elemIt != elemEndIt; ++elemIt)
//        {
//            if(elemIt->partitionType() == Dune::InteriorEntity)
//            {
//                fvGeometry.update(this->gridView(), *elemIt);
//                elemVolVars.update(*this, *elemIt, fvGeometry, false);
//
//                for (int scvIdx = 0; scvIdx < elemIt->template count<dim>(); ++scvIdx)
//
//                {
//                    boxSDSource(values,
//                                 *elemIt,
//                                 fvGeometry,
//                                 scvIdx,
//                                 elemVolVars);
//                    Valgrind::CheckDefined(values);
//
//                    // thickness of the domain, for 2D usually 1m
//                    Scalar extrusionFactor =
//                        elemVolVars[scvIdx].extrusionFactor();
//                    values[energyEqIdx] *= fvGeometry.subContVol[scvIdx].volume * extrusionFactor;
//
//                    source += values[energyEqIdx];
//                }
//            }
//        }
//
//        if (this->gridView().comm().size() > 1)
//            source = this->gridView().comm().sum(source);
//
//        return source;
//    }

    //! \copydoc Dumux::ImplicitProblem::postTimeStep()
    void postTimeStep()
    {
        // Calculate masses
        PrimaryVariables storage;

        this->model().globalStorage(storage);
        const Scalar time = this->timeManager().time() +  this->timeManager().timeStepSize();

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0)
        {
            if (this->timeManager().timeStepIndex() % freqMassOutput_ == 0
                    || this->timeManager().episodeWillBeOver())
            {
                PrimaryVariables evaporationRate(0.);
                evaporationRate = storageLastTimestep_ - storage;

                assert(time - lastMassOutputTime_ != 0);
                evaporationRate /= (time - lastMassOutputTime_);
                // 2d: interface length has to be accounted for
                // in order to obtain kg/m²s
                evaporationRate /= (bboxMax_[0]-bboxMin_[0]);

                std::cout << "TotalMass: " << storage[contiTotalMassIdx]
                          << " MassH2O: " << storage[contiWEqIdx]
                          << " IntEnergy: " << storage[energyEqIdx]
                          << " EvaporationRate: " << evaporationRate[contiWEqIdx]
                          << " Time: " << time <<std::endl;
                if (this->timeManager().time() != 0.)
                    outfile << time <<
							"; " << evaporationRate[contiTotalMassIdx] <<
                            "; " << evaporationRate[contiWEqIdx] <<
                            "; " << evaporationRate[energyEqIdx] <<
                            "; " << storage[contiWEqIdx] <<
                            std::endl;

                storageLastTimestep_ = storage;
                lastMassOutputTime_ = time;
            }
        }
    }

    /*!
     * \brief Determine if we are on a corner of the grid
     *
     * \param globalPos The global position
     *
     */
    bool isCornerPoint(const GlobalPosition &globalPos)
    {
        return ((onLeftBoundary_(globalPos) && onLowerBoundary_(globalPos)) ||
            (onLeftBoundary_(globalPos) && onUpperBoundary_(globalPos)) ||
            (onRightBoundary_(globalPos) && onLowerBoundary_(globalPos)) ||
            (onRightBoundary_(globalPos) && onUpperBoundary_(globalPos)));
    }

    // required in case of mortar coupling
    // otherwise it should return false
    /*!
     * \brief docme
     *
     * \param global Pos The global position
     *
     */
    bool isInterfaceCornerPoint(const GlobalPosition &globalPos) const
    { return false; }

    // \}

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)

    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        //TODO: call density from fluidsystem
//        const MaterialLawParams &materialParams =
//            this->spatialParams().materialLawParams(globalPos);
//        Scalar pC = MaterialLaw::pC(materialParams, initialSw1_);

        values[pNIdx] = refPressure_
                + 1000.*this->gravity()[1]*(globalPos[1]-bboxMax_[1]);
//                + pC;

        if (globalPos[0] < xMaterialInterface_)
		     values[sWIdx] = initialSw1_;
	     else
	         values[sWIdx] = initialSw2_;
        values[temperatureIdx] = refTemperature_;
        //        if (onUpperBoundary_(globalPos) && onLeftBoundary_(globalPos))
//        {
//            values[sWIdx] = 0.99e-4;
//            values[temperatureIdx] = 283.55;
//        }

//    	values[pNIdx] = 1e5*(2.0 - globalPos[0]);
//    	if (onLeftBoundary_(globalPos))
//    		values[sWIdx] = 0.9e-4;
//    	else
//    		values[sWIdx] = 1e-4;
//        values[temperatureIdx] = 283.15;
    }

    /// RADIATION stuff
//    const Scalar albedo_(const ElementVolumeVariables &elemVolVars, const int scvIdx) const                                     //albedo
//    {
//  	  if (elemVolVars[scvIdx].saturation(wPhaseIdx)*porosity_ > 0.3)
//  		  return 0.075;
//  	  if (0.04 <= elemVolVars[scvIdx].saturation(wPhaseIdx)*porosity_ &&
//  			elemVolVars[scvIdx].saturation(wPhaseIdx)*porosity_ <= 0.3)
//  		  return 0.1846-0.3654*elemVolVars[scvIdx].saturation(wPhaseIdx)*porosity_;
//  	  if (elemVolVars[scvIdx].saturation(wPhaseIdx)*porosity_ > 0.4)
//  		  return 0.17;
//	  return 0.17;
//    }
//
//    const Scalar surfaceEmissivity_(const ElementVolumeVariables &elemVolVars, const int scvIdx) const                                    //surface emissivity
//    {
//  	  if (elemVolVars[scvIdx].saturation(wPhaseIdx)*porosity_ < 0.3)
//  			return 0.93+0.1333*elemVolVars[scvIdx].saturation(wPhaseIdx)*porosity_;
//  	  if (elemVolVars[scvIdx].saturation(wPhaseIdx)*porosity_ >= 0.3)
//  		    return 0.97;
//  	  return 0.;
//    }
//
//    //to simulate the Solar irradiance for a number of days (definable in the input file)
//    const Scalar variationSun_(const Scalar value) const
//    {
//	  const Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();
//	  const Scalar days = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, days);
//
//	  Scalar solar = 0;
//	  for (int n=0; n<days; ++n)
//	  {
//		  if (time > 21600 + n*86400 && time < 64800 + n*86400)
//			  solar = cos(2*M_PI*((time+43200)/(86400))) * value;
//	  }
//	  return solar;
//	}
//
//    //to simulate the atmospheric Temperature for a number of days (definable in the input file)
//    const Scalar variationT_a_(const Scalar value) const
//	{
//    	const Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();
//    	const Scalar T = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, FreeFlow, RefTemperature);
//    	{
//			  return T+cos(2*M_PI*((time+36000)/(86400))) * value;
//    	}
//    	std::cerr << "error with T_a at time " << time << std::endl;
//	}
//
//    //to calculate the atmospheric Emissivity depending on the Temperature near the Surface
//    const Scalar atmosphericEmissivity() const
//    {
//    	const Scalar T_a = variationT_a_(6);
//
//      	return 1.24*std::pow((((Dumux::H2O<Scalar>::vaporPressure(T_a)/100)/T_a)),(1.0/7.0));
//    }
//
//    const Scalar cellHeight_() const
//    {
//    	return 0.0004575;
//    }
    /////////

// TODO: what should be used as reference mass fraction??
//	void updateFluidStateForBC_(FluidState& fluidState) const
//	{
//		for (unsigned phaseIdx=0; phaseIdx<numPhases; ++phaseIdx)
//		{
//			fluidState.setTemperature(refTemperature_);
//			fluidState.setPressure(phaseIdx, refPressure_);
//
//			Scalar massFraction[numComponents];
//			massFraction[wCompIdx] = refMassfrac_;
//			massFraction[nCompIdx] = 1 - massFraction[wCompIdx];
//
//			// calculate average molar mass of the gas phase
//			Scalar M1 = FluidSystem::molarMass(wCompIdx);
//			Scalar M2 = FluidSystem::molarMass(nCompIdx);
//			Scalar X2 = massFraction[nCompIdx];
//			Scalar massToMoleDenominator = M2 + X2*(M1 - M2);
//
//			fluidState.setMoleFraction(phaseIdx, wCompIdx, massFraction[wCompIdx]*M2/massToMoleDenominator);
//			fluidState.setMoleFraction(phaseIdx, nCompIdx, massFraction[nCompIdx]*M1/massToMoleDenominator);
//		}
//	}


    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < bboxMin_[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > bboxMax_[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < bboxMin_[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > bboxMax_[1] - eps_; }

    bool onBoundary_(const GlobalPosition &globalPos) const
    {
    	return (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)
    			|| onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos));
    }

    static constexpr Scalar eps_ = 1e-8;
    GlobalPosition bboxMin_;
    GlobalPosition bboxMax_;
    Scalar xMaterialInterface_;

    int freqMassOutput_;

    PrimaryVariables storageLastTimestep_;
    Scalar lastMassOutputTime_;

    Scalar refTemperature_;
    Scalar refPressure_;
    Scalar initialSw1_;
    Scalar initialSw2_;

//    Scalar porosity_;
    Scalar runUpDistanceX_;
    Scalar initializationTime_;
    std::ofstream outfile;
};
} //end namespace

#endif
