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
 * \brief Definition of a problem, for the 1p2c problem:
 * Component transport of nitrogen dissolved in the water phase.
 */
#ifndef DUMUX_1P2C_CONVECTION_PROBLEM_HH
#define DUMUX_1P2C_CONVECTION_PROBLEM_HH

//#include <dune/grid/alugrid/2d/alugrid.hh>
#if HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#include <math.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/implicit/1p2c/1p2cmodel.hh>
#include <dumux/implicit/common/implicitporousmediaproblem.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidsystems/h2on2liquidphasefluidsystem.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include "1p2coutflowspatialparams.hh"

#define NONISOTHERMAL 1

namespace Dumux
{

template <class TypeTag>
class OnePTwoCConvectionProblem;

namespace Properties
{
#if NONISOTHERMAL
NEW_TYPE_TAG(OnePTwoCConvectionProblem, INHERITS_FROM(OnePTwoCNI));
NEW_TYPE_TAG(OnePTwoCConvectionBoxProblem, INHERITS_FROM(BoxModel, OnePTwoCConvectionProblem));
NEW_TYPE_TAG(OnePTwoCConvectionCCProblem, INHERITS_FROM(CCModel, OnePTwoCConvectionProblem));
#else
NEW_TYPE_TAG(OnePTwoCConvectionProblem, INHERITS_FROM(OnePTwoC));
NEW_TYPE_TAG(OnePTwoCConvectionBoxProblem, INHERITS_FROM(BoxModel, OnePTwoCConvectionProblem));
NEW_TYPE_TAG(OnePTwoCConvectionCCProblem, INHERITS_FROM(CCModel, OnePTwoCConvectionProblem));
#endif

// Set the grid type
SET_PROP(OnePTwoCConvectionProblem, Grid)
{
#if HAVE_UG
    typedef Dune::UGGrid<2> type;
#else
//    typedef Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming> type;
    typedef Dune::YaspGrid<2> type;
#endif
};

// Set the problem property
SET_PROP(OnePTwoCConvectionProblem, Problem)
{
    typedef Dumux::OnePTwoCConvectionProblem<TypeTag> type;
};

// Set fluid configuration
SET_PROP(OnePTwoCConvectionProblem, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::FluidSystems::H2ON2LiquidPhase<Scalar, true> type;
};

// Set the spatial parameters
SET_TYPE_PROP(OnePTwoCConvectionProblem,
              SpatialParams,
              Dumux::OnePTwoCOutflowSpatialParams<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(OnePTwoCConvectionProblem, UseMoles, true);

// Enable velocity output
SET_BOOL_PROP(OnePTwoCConvectionProblem, VtkAddVelocity, true);

// Disable gravity
SET_BOOL_PROP(OnePTwoCConvectionProblem, ProblemEnableGravity, false);
}


/*!
 * \ingroup OnePTwoCModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Definition of a problem, for the 1p2c problem:
 * Nitrogen is dissolved in the water phase and
 * is transported with the water flow from the left side to the right.
 *
 * The model domain is 1 m times 1 m with a discretization length of 0.05 m
 * and homogeneous soil properties (\f$ \mathrm{K=10e-10, \Phi=0.4, \tau=0.28}\f$).
 * Initially the domain is filled with pure water.
 *
 * At the left side, a Dirichlet condition defines a nitrogen mole fraction
 * of 0.3 mol/mol.
 * The water phase flows from the left side to the right due to the applied pressure
 * gradient of 1e5 Pa/m. The nitrogen is transported with the water flow
 * and leaves the domain at the right boundary
 * where an outflow boundary condition is applied.
 * 
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. The default setting for useMoles is true.
 *
 * This problem uses the \ref OnePTwoCModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box1p2c -parameterFile ./test_box1p2c.input</tt> or 
 * <tt>./test_cc1p2c -parameterFile ./test_cc1p2c.input</tt>
 */
template <class TypeTag>
class OnePTwoCConvectionProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel) ThermalConductivityModel;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef Dumux::H2O<Scalar> IapwsH2O;






    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // world dimension
        dimWorld = GridView::dimensionworld
    };

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dimWorld : 0 };

    enum {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,
#if NONISOTHERMAL
        temperatureIdx = Indices::temperatureIdx
#endif
    };
    enum {
        // index of the transport equation
        conti0EqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx,
#if NONISOTHERMAL
        energyEqIdx = Indices::energyEqIdx
#endif
    };


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    //! property that defines whether mole or mass fractions are used
        static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    OnePTwoCConvectionProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
        , eps_(1e-6)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, 
                                             std::string, 
                                             Problem, Name);
        outputInterval_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                int, Problem, OutputInterval);
        darcyVelocity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                Scalar, Problem, DarcyVelocity);

        temperatureHigh_ = 291.;
        temperatureLow_ = 290.;
        pressureHigh_ = 2e5;
        pressureLow_ = 1e5;




        //stateing in the console whether mole or mass fractions are used
        if(!useMoles)
        {
        	std::cout<<"problem uses mass-fractions"<<std::endl;
        }
        else
        {
        	std::cout<<"problem uses mole-fractions"<<std::endl;
        }
    }


    bool shouldWriteOutput() const
    {
        return
            this->timeManager().timeStepIndex() == 0 ||
            this->timeManager().timeStepIndex() % outputInterval_ == 0 ||
            this->timeManager().episodeWillBeOver() ||
            this->timeManager().willBeFinished();
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    void addOutputVtkFields()
        {
        //Here we calculate the analytical solution
            typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;
            unsigned numDofs = this->model().numDofs();

            //create required scalar fields
            ScalarField *temperatureExact = this->resultWriter().allocateManagedBuffer(numDofs);

            FVElementGeometry fvGeometry;
            VolumeVariables volVars;

            ElementIterator eIt = this->gridView().template begin<0>();
            ElementIterator eEndIt = this->gridView().template end<0>();
            fvGeometry.update(this->gridView(), *eIt);
            PrimaryVariables initialPriVars(0);
            GlobalPosition globalPos(0);
            initial_(initialPriVars, globalPos);

            //update the constant volume variables
            volVars.update(initialPriVars,
                           *this,
                           *eIt,
                           fvGeometry,
                           0,
                           false);

            Scalar porosity = this->spatialParams().porosity(*eIt, fvGeometry, 0);
            Scalar densityW = volVars.density();
            Scalar heatCapacityW = FluidSystem::heatCapacity(volVars.fluidState(), 0);
            Scalar effectiveHeatCapacityS = this->spatialParams().heatCapacity(*eIt, fvGeometry, 0);
            Scalar storageW =  densityW*heatCapacityW*porosity;
            Scalar storageTotal = storageW + effectiveHeatCapacityS;
            std::cout<<"storage: "<<storageTotal<<std::endl;
            Scalar time = std::max(this->timeManager().time() + this->timeManager().timeStepSize(), 1e-10);
            Scalar retardedFrontVelocity = darcyVelocity_*storageW/storageTotal/porosity;
            std::cout<<"retarded velocity: "<<retardedFrontVelocity<<std::endl;

            for (; eIt != eEndIt; ++eIt)
            {
                fvGeometry.update(this->gridView(), *eIt);
                for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                {
                    int globalIdx = this->model().dofMapper().map(*eIt, scvIdx, dofCodim);
                    if (isBox)
                        globalPos = eIt->geometry().corner(scvIdx);
                    else
                        globalPos = eIt->geometry().center();

                    (*temperatureExact)[globalIdx] = globalPos[0] < retardedFrontVelocity*time ? temperatureHigh_ : temperatureLow_;
                }
            }
            this->resultWriter().attachDofData(*temperatureExact, "temperatureExact", isBox);

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
        return name_.c_str();
    }

#if !NONISOTHERMAL
    /*!
     * \brief Returns the temperature within the domain [K].
     *
     * This problem assumes a temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; }; // in [K]

#endif

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
     * \param globalPos The position for which the bc type should be evaluated
     */
    void boundaryTypesAtPos(BoundaryTypes &values, 
                            const GlobalPosition &globalPos) const
    {
        if(globalPos[0] > this->bBoxMax()[0] - eps_)
        {
            values.setAllDirichlet();
        }
        else
        {
            values.setAllNeumann();
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position for which the bc type should be evaluated
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann
     *        boundary segment in dependency on the current solution.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local index of the sub-control volume
     * \param boundaryFaceIdx The index of the boundary face
     * \param elemVolVars All volume variables for the element
     *
     * This method is used for cases, when the Neumann condition depends on the
     * solution and requires some quantities that are specific to the fully-implicit method.
     * The \a values store the mass flux of each phase normal to the boundary.
     * Negative values indicate an inflow.
     */
     void solDependentNeumann(PrimaryVariables &values,
                      const Element &element,
                      const FVElementGeometry &fvGeometry,
                      const Intersection &intersection,
                      const int scvIdx,
                      const int boundaryFaceIdx,
                      const ElementVolumeVariables &elemVolVars) const
    {
        values = 0;
        GlobalPosition globalPos(0);
        if (isBox)
            globalPos = fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;
        else
            globalPos = fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;

        if(globalPos[0] < eps_)
        {
            values[pressureIdx] = -darcyVelocity_*elemVolVars[scvIdx].molarDensity()*elemVolVars[scvIdx].moleFraction(pressureIdx);
            values[massOrMoleFracIdx] = -darcyVelocity_*elemVolVars[scvIdx].molarDensity()*elemVolVars[scvIdx].moleFraction(massOrMoleFracIdx);
            values[temperatureIdx] = -darcyVelocity_*elemVolVars[scvIdx].density()
                                     *IapwsH2O::liquidEnthalpy(temperatureHigh_, elemVolVars[scvIdx].pressure());
//            std::cout<<"values[temperatureIdx]: "<<values[temperatureIdx]<<std::endl;
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
     * For this method, the \a priVars parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^3*s) or kg/(m^3*s))
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

    // \}

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &priVars,
                  const GlobalPosition &globalPos) const
    {
        priVars[pressureIdx] = pressureLow_; // initial condition for the pressure
        priVars[massOrMoleFracIdx] = 1e-10;  // initial condition for the N2 molefraction
#if NONISOTHERMAL
        priVars[temperatureIdx] = temperatureLow_;
#endif
    }

    Scalar temperatureHigh_;
    Scalar temperatureLow_;
    Scalar pressureHigh_;
    Scalar pressureLow_;
    Scalar darcyVelocity_;
    const Scalar eps_;
    std::string name_;
    int outputInterval_;
};

} //end namespace
#endif
