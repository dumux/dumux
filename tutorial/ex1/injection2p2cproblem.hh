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
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 */
#ifndef DUMUX_INJECTION_2P2C_PROBLEM_HH
#define DUMUX_INJECTION_2P2C_PROBLEM_HH

#include <dumux/porousmediumflow/2p2c/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/fluidsystems/h2on2.hh>

#include "injection2pspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class Injection2p2cProblem;

namespace Properties
{
NEW_TYPE_TAG(Injection2p2cProblem, INHERITS_FROM(TwoPTwoC, InjectionSpatialParams));
NEW_TYPE_TAG(Injection2p2cBoxProblem, INHERITS_FROM(BoxModel, Injection2p2cProblem));
NEW_TYPE_TAG(Injection2p2pcCCProblem, INHERITS_FROM(CCModel, Injection2p2cProblem));

// Set the grid type
SET_TYPE_PROP(Injection2p2cProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(Injection2p2cProblem, Problem, Injection2p2cProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(Injection2p2cProblem,
              FluidSystem,
              FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar), false /*useComplexRelations*/>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(Injection2p2cProblem, UseMoles, true);
}


/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief Problem where air is injected under a low permeable layer in a depth of 2700m.
 *
 * The domain is sized 60m times 40m and consists of two layers, a moderately
 * permeable one (\f$ K=10e-12\f$) and one with a lower permeablility (\f$ K=10e-13\f$)
 * in the rest of the domain.
 *
 * Nitrogen is injected into a water-filled aquifer through a well. First, we inject for one month.
 * Then, we continue simulating the development of the nitrogen plume for 10 years.This is realized with a solution-dependent Neumann boundary condition at the right boundary
 * (\f$ 7m<y<15m\f$). The aquifer is situated 2700m below sea level. The injected fluid phase migrates upwards due to buoyancy.
 * It accumulates and partially enters the lower permeable aquitard.
 *
 * The default setting for useMoles is true, i.e. each component is balaced in units of mole.
 * This problem uses the \ref TwoPTwoCModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./exercise1_2p2c -ParameterFile exercise1.input>
 */
template <class TypeTag>
class Injection2p2cProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // grid world dimension
    static constexpr auto dimWorld = GridView::dimensionworld;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    Injection2p2cProblem(TimeManager &timeManager,
                     const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
         // initialize the tables of the fluid system
        FluidSystem::init(/*tempMin=*/273.15,
                /*tempMax=*/423.15,
                /*numTemp=*/50,
                /*pMin=*/0.0,
                /*pMax=*/30e6,
                /*numP=*/300);

        // name of the problem and output file
        name_  = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
        // depth of the aquifer, units: m
        aquiferDepth_  = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, AquiferDepth);
        // the duration of the injection, units: second
        injectionDuration_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, InjectionDuration);

        // TODO: dumux-course-task
        // Get the specific inflow of 1e-4 kg/(s m^2) from the input file (totalAreaSpecificInflow_) here as it is done for the injectionDuration_.

    }

    /*!
     * \brief User defined output after the time integration
     *
     * Will be called diretly after the time integration.
     */
    void postTimeStep()
    {
        // Calculate storage terms
        PrimaryVariables storageW, storageN;
        this->model().globalPhaseStorage(storageW, Indices::wPhaseIdx);
        this->model().globalPhaseStorage(storageN, Indices::nPhaseIdx);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout<<"Storage: wetting=[" << storageW << "]"
                     << " nonwetting=[" << storageN << "]\n";
        }
    }


    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    std::string name() const
    { return name_+"-2p2c"; }

    /*!
     * \brief Returns the temperature \f$ K \f$
     */
    Scalar temperature() const
    {
        return 273.15 + 30; // [K]
    }

    /*!
     * \brief Returns the source term
     *
     * \param values Stores the source values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} / (m^\textrm{dim} \cdot s )] \f$
     * \param globalPos The global position
     */
    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
    {
        values = 0;
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment
     *
     * \param values Stores the value of the boundary type
     * \param globalPos The global position
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        if (globalPos[dimWorld-1] < eps_)
            values.setAllDirichlet();

        // and Neuman boundary conditions everywhere else
        // note that we don't differentiate between Neumann and Robin boundary types
        else
            values.setAllNeumann();
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet
     *        boundary segment
     *
     * \param values Stores the Dirichlet values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variable} ] \f$
     * \param globalPos The global position
     */
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
       initialAtPos(values, globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param globalPos The globalPosition of the boundary interface
     */
     void neumannAtPos(PrimaryVariables &values,
                       const GlobalPosition &globalPos) const
    {
         // initialize values to zero, i.e. no-flow Neumann boundary conditions
         values = 0;

        //if we are inside the injection zone set inflow Neumann boundary conditions
        if (this->timeManager().time() + this->timeManager().timeStepSize() < injectionDuration_
            && globalPos[1] < 15 + eps_ && globalPos[1] > 7 - eps_ && globalPos[0] > 0.9*this->bBoxMax()[0])
        {
            // TODO: dumux-course-task
            //instead of setting -1e-4 here directly use totalAreaSpecificInflow_ in the computation

            // set the Neumann values for the Nitrogen component balance
            // convert from units kg/(s*m^2) to mole/(s*m^2)
            values[Indices::contiNEqIdx] = -1e-4/FluidSystem::molarMass(FluidSystem::nCompIdx);
            values[Indices::contiWEqIdx] = 0.0;
        }
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    void initialAtPos(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        // get the water density at atmospheric conditions
        const Scalar densityW = FluidSystem::H2O::liquidDensity(temperature(), 1.0e5);

        // assume an intially hydrostatic liquid pressure profile
        // note: we subtract rho_w*g*h because g is defined negative
        const Scalar pw = 1.0e5 - densityW*this->gravity()[dimWorld-1]*(aquiferDepth_ - globalPos[dimWorld-1]);

        // initially we have some nitrogen dissolved
        // saturation mole fraction would be
        // moleFracLiquidN2 = (pw + pc + p_vap^sat)/henry;
        const Scalar moleFracLiquidN2 = pw*0.95/BinaryCoeff::H2O_N2::henry(temperature());

        // note that because we start with a single phase system the primary variables
        // are pl and x^w_N2. This will switch as soon after we start injecting to a two
        // phase system so the primary variables will be pl and Sn (non-wetting saturation).
        values[Indices::switchIdx] = moleFracLiquidN2;
        values[Indices::pressureIdx] = pw;
    }

    /*!
     * \brief Return the initial phase state inside a control volume.
     *
     * \param globalPos The global position
     */
    int initialPhasePresenceAtPos(const GlobalPosition &globalPos) const
    { return Indices::wPhaseOnly; }

    // \}

private:
    static constexpr Scalar eps_ = 1e-6;
    std::string name_; //! Problem name
    Scalar aquiferDepth_; //! Depth of the aquifer in m
    Scalar injectionDuration_; //! Duration of the injection in seconds
    //TODO: dumux-course-task
    //define the Scalar totalAreaSpecificInflow_ here

};
} //end namespace

#endif
