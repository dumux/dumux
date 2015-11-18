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
 * \brief Test for the OnePTwoCModel in combination with the NI model for a conduction problem:
 * The simulation domain is a tube where with an elevated temperature on the left hand side.
 *
 */
#ifndef DUMUX_1P2CNI_CONDUCTION_PROBLEM_HH
#define DUMUX_1P2CNI_CONDUCTION_PROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/implicit/1p2c/1p2cmodel.hh>
#include <dumux/implicit/common/implicitporousmediaproblem.hh>

#include <dumux/material/fluidsystems/h2on2liquidphasefluidsystem.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include "1p2cnispatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class OnePTwoCNIConductionProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePTwoCNIConductionProblem, INHERITS_FROM(OnePTwoCNI));
NEW_TYPE_TAG(OnePTwoCNIConductionBoxProblem, INHERITS_FROM(BoxModel, OnePTwoCNIConductionProblem));
NEW_TYPE_TAG(OnePTwoCNIConductionCCProblem, INHERITS_FROM(CCModel, OnePTwoCNIConductionProblem));



// Set the grid type
SET_TYPE_PROP(OnePTwoCNIConductionProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(OnePTwoCNIConductionProblem, Problem,
              Dumux::OnePTwoCNIConductionProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(OnePTwoCNIConductionProblem,
              FluidSystem,
              Dumux::FluidSystems::H2ON2LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar), true>);

// Set the spatial parameters
SET_TYPE_PROP(OnePTwoCNIConductionProblem,
              SpatialParams,
              Dumux::OnePTwoCNISpatialParams<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(OnePTwoCNIConductionProblem, UseMoles, true);

// Enable velocity output
SET_BOOL_PROP(OnePTwoCNIConductionProblem, VtkAddVelocity, true);

// Disable gravity
SET_BOOL_PROP(OnePTwoCNIConductionProblem, ProblemEnableGravity, false);
}


/*!
 * \ingroup OnePTwoCNIModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Test for the OnePTwoCModel in combination with the NI model for a conduction problem:
 * The simulation domain is a tube where with an elevated temperature on the left hand side.
 *
 * Initially the domain is fully saturated with water at a constant temperature.
 * On the left hand side there is a Dirichlet boundary condition with an increased temperature and on the right hand side
 * a Dirichlet boundary with constant pressure, saturation and temperature is applied.
 *
 * The results are compared to an analytical solution for a diffusion process:
  \f[
     T =T_{high} + (T_{init} - T_{high})erf \left(0.5\sqrt{\frac{x^2 S_{total}}{t \lambda_{eff}}}\right)
 \f]
 *
 * This problem uses the \ref OnePTwoCModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell: <br>
 * <tt>./test_box1p2cniconduction -ParameterFile ./test_box1p2cniconduction.input</tt> or <br>
 * <tt>./test_cc1p2cniconduction -ParameterFile ./test_cc1p2cniconduction.input</tt>
 */
template <class TypeTag>
class OnePTwoCNIConductionProblem : public ImplicitPorousMediaProblem<TypeTag>
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
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

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
        temperatureIdx = Indices::temperatureIdx
    };
    enum {
        // index of the transport equation
        conti0EqIdx = Indices::conti0EqIdx,
        transportEqIdx = Indices::transportEqIdx,
        energyEqIdx = Indices::energyEqIdx
    };


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    OnePTwoCNIConductionProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
        , eps_(1e-6)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);
        outputInterval_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                  int, Problem, OutputInterval);

        temperatureHigh_ = 300.;
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

            const auto firstElement = *this->gridView().template begin<0>();
            fvGeometry.update(this->gridView(), firstElement);
            PrimaryVariables initialPriVars(0);
            GlobalPosition globalPos(0);
            initial_(initialPriVars, globalPos);

            //update the constant volume variables
            volVars.update(initialPriVars,
                           *this,
                           firstElement,
                           fvGeometry,
                           0,
                           false);

            Scalar porosity = this->spatialParams().porosity(firstElement, fvGeometry, 0);
            Scalar densityW = volVars.density();
            Scalar heatCapacityW = FluidSystem::heatCapacity(volVars.fluidState(), 0);
            Scalar densityS = this->spatialParams().solidDensity(firstElement, fvGeometry, 0);
            Scalar heatCapacityS = this->spatialParams().solidHeatCapacity(firstElement, fvGeometry, 0);
            Scalar storage = densityW*heatCapacityW*porosity + densityS*heatCapacityS*(1 - porosity);
            Scalar effectiveThermalConductivity = ThermalConductivityModel::effectiveThermalConductivity(volVars, this->spatialParams(),
                                                                                                         firstElement, fvGeometry, 0);
            Scalar time = std::max(this->timeManager().time() + this->timeManager().timeStepSize(), 1e-10);


            for (const auto& element : Dune::elements(this->gridView()))
            {
                fvGeometry.update(this->gridView(), element);
                for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
                {
                    int dofIdxGlobal = this->model().dofMapper().subIndex(element, scvIdx, dofCodim);

                    if (isBox)
                        globalPos = element.geometry().corner(scvIdx);
                    else
                        globalPos = element.geometry().center();

                    (*temperatureExact)[dofIdxGlobal] = temperatureHigh_ + (initialPriVars[temperatureIdx] - temperatureHigh_)
                                                     *std::erf(0.5*std::sqrt(globalPos[0]*globalPos[0]*storage/time/effectiveThermalConductivity));
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
        if(globalPos[0] < eps_ || globalPos[0] > this->bBoxMax()[0] - eps_)
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

        // condition for the N2 molefraction at left boundary
        if (globalPos[0] < eps_)
        {
            values[temperatureIdx] = temperatureHigh_;
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     *
     * The units must be according to either using mole or mass fractions. (mole/(m^2*s) or kg/(m^2*s))
     */
    void neumann(PrimaryVariables &priVars,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &intersection,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        priVars = 0;
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
        priVars[pressureIdx] = 1e5; // initial condition for the pressure
        priVars[massOrMoleFracIdx] = 1e-5;  // initial condition for the N2 molefraction
        priVars[temperatureIdx] = 290.;
    }

    Scalar temperatureHigh_;
    const Scalar eps_;
    std::string name_;
    int outputInterval_;
};

} //end namespace
#endif // DUMUX_1P2CNI_CONDUCTION_PROBLEM_HH
