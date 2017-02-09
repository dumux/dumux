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
 *
 * \brief Test for the RichardsModel in combination with the NI model for a convection problem:
 * The simulation domain is a tube where water with an elevated temperature is injected
 * at a constant rate on the left hand side.
 */
#ifndef DUMUX_1PNI_CONVECTION_PROBLEM_HH
#define DUMUX_1PNI_CONVECTION_PROBLEM_HH

#include <math.h>

#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/porousmediumflow/richards/implicit/model.hh>
#include <dumux/porousmediumflow/richards/implicit/problem.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/material/fluidsystems/h2on2.hh>
#include "richardsnispatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class RichardsNIConvectionProblem;

namespace Properties
{
NEW_TYPE_TAG(RichardsNIConvectionProblem, INHERITS_FROM(RichardsNI, RichardsNISpatialParams));
NEW_TYPE_TAG(RichardsNIConvectionBoxProblem, INHERITS_FROM(BoxModel, RichardsNIConvectionProblem));
NEW_TYPE_TAG(RichardsNIConvectionCCProblem, INHERITS_FROM(CCModel, RichardsNIConvectionProblem));

// Set the grid type
SET_TYPE_PROP(RichardsNIConvectionProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(RichardsNIConvectionProblem, Problem,
              RichardsNIConvectionProblem<TypeTag>);

// Set the fluid system
SET_TYPE_PROP(RichardsNIConvectionProblem, FluidSystem, FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false>);

// Set the spatial parameters
SET_TYPE_PROP(RichardsNIConvectionProblem,
              SpatialParams,
              RichardsNISpatialParams<TypeTag>);
}


/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Test for the RichardsModel in combination with the NI model for a convection problem:
 * The simulation domain is a tube where water with an elevated temperature is injected
 * at a constant rate on the left hand side.
 *
 * Initially the domain is fully saturated with water at a constant temperature.
 * On the left hand side water is injected at a constant rate and on the right hand side
 * a Dirichlet boundary with constant pressure, saturation and temperature is applied.
 *
 * The results are compared to an analytical solution where a retarded front velocity is calculated as follows:
  \f[
     v_{Front}=\frac{q S_{water}}{\phi S_{total}}
 \f]
 *
 * The result of the analytical solution is written into the vtu files.
 * This problem uses the \ref RichardsModel and \ref NIModel model.
 *
 * To run the simulation execute the following line in shell: <br>
 * <tt>./test_boxrichardsniconvection -ParameterFile ./test_boxrichardsniconvection.input</tt> or <br>
 * <tt>./test_ccrichardsniconvection -ParameterFile ./test_ccrichardsniconvection.input</tt>
 */
template <class TypeTag>
class RichardsNIConvectionProblem : public RichardsProblem<TypeTag>
{
    typedef RichardsProblem<TypeTag> ParentType;
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
    typedef H2O<Scalar> IapwsH2O;

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
        pressureIdx = Indices::pwIdx,
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        temperatureIdx = Indices::temperatureIdx
    };
    enum {
        // index of the transport equation
        contiEqIdx = Indices::contiEqIdx,
        energyEqIdx = Indices::energyEqIdx
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    RichardsNIConvectionProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        //initialize fluid system
        FluidSystem::init();

        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        outputInterval_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Problem, OutputInterval);
        darcyVelocity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, DarcyVelocity);
        temperatureHigh_ = 291.;
        temperatureLow_ = 290.;
        pressureHigh_ = 2e5;
        pressureLow_ = 1e5;
    }


    bool shouldWriteOutput() const
    {
        return
            this->timeManager().timeStepIndex() == 0 ||
            this->timeManager().timeStepIndex() % outputInterval_ == 0 ||
            this->timeManager().episodeWillBeFinished() ||
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
        volVars.update(initialPriVars, *this, firstElement, fvGeometry, 0, false);

        Scalar porosity = this->spatialParams().porosity(firstElement, fvGeometry, 0);
        Scalar densityW = volVars.density(wPhaseIdx);
        Scalar heatCapacityW = FluidSystem::heatCapacity(volVars.fluidState(), 0);
        Scalar storageW =  densityW*heatCapacityW*porosity;
        Scalar densityS = this->spatialParams().solidDensity(firstElement, fvGeometry, 0);
        Scalar heatCapacityS = this->spatialParams().solidHeatCapacity(firstElement, fvGeometry, 0);
        Scalar storageTotal = storageW + densityS*heatCapacityS*(1 - porosity);
        std::cout<<"storage: "<<storageTotal<<std::endl;
        using std::max;
        Scalar time = max(this->timeManager().time() + this->timeManager().timeStepSize(), 1e-10);
        Scalar retardedFrontVelocity = darcyVelocity_*storageW/storageTotal/porosity;
        std::cout<<"retarded velocity: "<<retardedFrontVelocity<<std::endl;

        for (const auto& element : elements(this->gridView()))
        {
            fvGeometry.update(this->gridView(), element);
            for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
            {
                int globalIdx = this->model().dofMapper().subIndex(element, scvIdx, dofCodim);

                if (isBox)
                    globalPos = element.geometry().corner(scvIdx);
                else
                    globalPos = element.geometry().center();

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
            values[pressureIdx] = -darcyVelocity_*elemVolVars[scvIdx].density(wPhaseIdx);
            values[temperatureIdx] = -darcyVelocity_*elemVolVars[scvIdx].density(wPhaseIdx)
                                     *IapwsH2O::liquidEnthalpy(temperatureHigh_, elemVolVars[scvIdx].pressure(wPhaseIdx));
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
     */
    void sourceAtPos(PrimaryVariables &priVars,
                     const GlobalPosition &globalPos) const
    {
        priVars = Scalar(0.0);
    }

    /*!
     * \brief Returns the reference pressure [Pa] of the non-wetting
     *        fluid phase within a finite volume
     *
     * This problem assumes a constant reference pressure of 1 bar.
     *
     * \param element The DUNE Codim<0> entity which intersects with
     *                the finite volume in question
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The sub control volume index inside the finite
     *               volume geometry
     */
    Scalar referencePressure(const Element &element,
                             const FVElementGeometry &fvGeometry,
                             const int scvIdx) const
    { return 1e5; };

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
        priVars[temperatureIdx] = temperatureLow_;
    }

    Scalar temperatureHigh_;
    Scalar temperatureLow_;
    Scalar pressureHigh_;
    Scalar pressureLow_;
    Scalar darcyVelocity_;
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;
    int outputInterval_;
};

} //end namespace
#endif // DUMUX_RICHARDSNI_CONVECTION_PROBLEM_HH
