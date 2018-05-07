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
 * \brief Definition of a problem, where CO2 is injected in a reservoir.
 */
#ifndef DUMUX_HETEROGENEOUS_PROBLEM_HH
#define DUMUX_HETEROGENEOUS_PROBLEM_HH

#include <dumux/porousmediumflow/co2/implicit/model.hh>
#include <dumux/porousmediumflow/co2/implicit/volumevariables.hh>
#include <dumux/material/fluidsystems/brineco2.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/implicit/box/intersectiontovertexbc.hh>


#include "heterogeneousspatialparameters.hh"
#include "heterogeneousco2tables.hh"

namespace Dumux
{

template <class TypeTag>
class HeterogeneousProblem;

namespace Properties
{
NEW_TYPE_TAG(HeterogeneousProblem, INHERITS_FROM(TwoPTwoC, HeterogeneousSpatialParams));
NEW_TYPE_TAG(HeterogeneousBoxProblem, INHERITS_FROM(BoxModel, HeterogeneousProblem));
NEW_TYPE_TAG(HeterogeneousCCProblem, INHERITS_FROM(CCModel, HeterogeneousProblem));

// Set the grid type
#if HAVE_DUNE_ALUGRID
SET_TYPE_PROP(HeterogeneousProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);
#else
SET_TYPE_PROP(HeterogeneousProblem, Grid, Dune::YaspGrid<2>);
#endif

// Set the problem property
SET_TYPE_PROP(HeterogeneousProblem, Problem, HeterogeneousProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(HeterogeneousProblem, FluidSystem, BrineCO2FluidSystem<TypeTag>);

// Set the CO2 table to be used; in this case not the the default table
SET_TYPE_PROP(HeterogeneousProblem, CO2Table, HeterogeneousCO2Tables::CO2Tables);

// Set the salinity mass fraction of the brine in the reservoir
SET_SCALAR_PROP(HeterogeneousProblem, ProblemSalinity, 1e-1);

//! the CO2 Model and VolumeVariables properties
SET_TYPE_PROP(HeterogeneousProblem, Model, CO2Model<TypeTag>);
SET_TYPE_PROP(HeterogeneousProblem, VolumeVariables, CO2VolumeVariables<TypeTag>);

// Use Moles
SET_BOOL_PROP(HeterogeneousProblem, UseMoles, false);
}


/*!
 * \ingroup CO2Model
 * \ingroup ImplicitTestProblems
 * \brief Definition of a problem, where CO2 is injected in a reservoir.
 *
 * The domain is sized 200m times 100m and consists of four layers, a
 * permeable reservoir layer at the bottom, a barrier rock layer with reduced permeability, another reservoir layer
 * and at the top a barrier rock layer with a very low permeablility.
 *
 * CO2 is injected at the permeable bottom layer
 * from the left side. The domain is initially filled with brine.
 *
 * The grid is unstructered and permeability and porosity for the elements are read in from the grid file. The grid file
 * also contains so-called boundary ids which can be used assigned during the grid creation in order to differentiate
 * between different parts of the boundary.
 * These boundary ids can be imported into the problem where the boundary conditions can then be assigned accordingly.
 *
 * The model is able to use either mole or mass fractions. The property useMoles can be set to either true or false in the
 * problem file. Make sure that the according units are used in the problem setup. The default setting for useMoles is false.
 *
 * To run the simulation execute the following line in shell (works with the box and cell centered spatial discretization method):
 * <tt>./test_ccco2 </tt> or <tt>./test_boxco2 </tt>
 */
template <class TypeTag >
class HeterogeneousProblem : public ImplicitPorousMediaProblem<TypeTag>
{

    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        lPhaseIdx = Indices::wPhaseIdx,
        gPhaseIdx = Indices::nPhaseIdx
    };
    enum {
        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx
    };
    enum {
        BrineIdx = FluidSystem::BrineIdx,
        CO2Idx = FluidSystem::CO2Idx
    };
    enum {
        conti0EqIdx = Indices::conti0EqIdx,
        contiCO2EqIdx = conti0EqIdx + CO2Idx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef typename GET_PROP_TYPE(TypeTag, CO2Table) CO2Table;
    typedef Dumux::CO2<Scalar, CO2Table> CO2;
    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };
    enum { dofCodim = isBox ? dim : 0 };
    //! property that defines whether mole or mass fractions are used
    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    /*!
     * \brief The constructor
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    HeterogeneousProblem(TimeManager &timeManager,
                     const GridView &gridView)
        : ParentType(timeManager, GridCreator::grid().leafGridView()),
          //Boundary Id Setup:
          injectionTop_(1),
          injectionBottom_(2),
          dirichletBoundary_(3),
          noFlowBoundary_(4),
          intersectionToVertexBC_(*this)
    {
            nTemperature_       = GET_RUNTIME_PARAM(TypeTag, int, FluidSystem.NTemperature);
            nPressure_          = GET_RUNTIME_PARAM(TypeTag, int, FluidSystem.NPressure);
            pressureLow_        = GET_RUNTIME_PARAM(TypeTag, Scalar, FluidSystem.PressureLow);
            pressureHigh_       = GET_RUNTIME_PARAM(TypeTag, Scalar, FluidSystem.PressureHigh);
            temperatureLow_     = GET_RUNTIME_PARAM(TypeTag, Scalar, FluidSystem.TemperatureLow);
            temperatureHigh_    = GET_RUNTIME_PARAM(TypeTag, Scalar, FluidSystem.TemperatureHigh);
            depthBOR_           = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.DepthBOR);
            name_               = GET_RUNTIME_PARAM(TypeTag, std::string, Problem.Name);
            injectionRate_      = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.InjectionRate);

        /* Alternative syntax:
         * typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
         * const Dune::ParameterTree &tree = ParameterTree::tree();
         * nTemperature_       = tree.template get<int>("FluidSystem.nTemperature");
         *
         * + We see what we do
         * - Reporting whether it was used does not work
         * - Overwriting on command line not possible
        */

        // set the spatial parameters by reading the DGF grid file
        this->spatialParams().setParams();

        eps_ = 1e-6;

        // initialize the tables of the fluid system
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        //stating in the console whether mole or mass fractions are used
        if(useMoles)
        {
            std::cout<<"problem uses mole fractions"<<std::endl;
        }
        else
        {
            std::cout<<"problem uses mass fractions"<<std::endl;
        }
    }

    /*!
     * \brief User defined output after the time integration
     *
     * Will be called diretly after the time integration.
     */
    void postTimeStep()
    {
        // Calculate storage terms
        PrimaryVariables storageL, storageG;
        this->model().globalPhaseStorage(storageL, lPhaseIdx);
        this->model().globalPhaseStorage(storageG, gPhaseIdx);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout<<"Storage: liquid=[" << storageL << "]"
                     << " gas=[" << storageG << "]\n";
        }
    }

    /*!
     * \brief Append all quantities of interest which can be derived
     *        from the solution of the current time step to the VTK
     *        writer.
     */
    void addOutputVtkFields()
    {
        typedef Dune::BlockVector<Dune::FieldVector<double, 1> > ScalarField;

         // get the number of degrees of freedom
         unsigned numDofs = this->model().numDofs();
         unsigned numElements = this->gridView().size(0);

         //create required scalar fields
         ScalarField *Kxx = this->resultWriter().allocateManagedBuffer(numElements);
         ScalarField *cellPorosity = this->resultWriter().allocateManagedBuffer(numElements);
         ScalarField *boxVolume = this->resultWriter().allocateManagedBuffer(numDofs);
         (*boxVolume) = 0;

         //Fill the scalar fields with values
         ScalarField *rank = this->resultWriter().allocateManagedBuffer(numElements);

         FVElementGeometry fvGeometry;
         VolumeVariables volVars;

         for (const auto& element : elements(this->gridView()))
         {
             int eIdx = this->elementMapper().index(element);
             (*rank)[eIdx] = this->gridView().comm().rank();
             fvGeometry.update(this->gridView(), element);


             for (int scvIdx = 0; scvIdx < fvGeometry.numScv; ++scvIdx)
             {
                 int dofIdxGlobal = this->model().dofMapper().subIndex(element, scvIdx, dofCodim);
                 volVars.update(this->model().curSol()[dofIdxGlobal],
                                *this,
                                element,
                                fvGeometry,
                                scvIdx,
                                false);
                 (*boxVolume)[dofIdxGlobal] += fvGeometry.subContVol[scvIdx].volume;
             }
             (*Kxx)[eIdx] = this->spatialParams().intrinsicPermeability(element, fvGeometry, /*element data*/ 0);
             (*cellPorosity)[eIdx] = this->spatialParams().porosity(element, fvGeometry, /*element data*/ 0);
         }

         //pass the scalar fields to the vtkwriter
         this->resultWriter().attachDofData(*Kxx, "Kxx", false); //element data
         this->resultWriter().attachDofData(*cellPorosity, "cellwisePorosity", false); //element data
         this->resultWriter().attachDofData(*boxVolume, "boxVolume", isBox);
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
    const std::string name() const
    { return name_; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * \param globalPos The position
     *
     * This problem assumes a geothermal gradient with
     * a surface temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return temperature_(globalPos);
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
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param vertex The vertex for which the boundary type is set
     */
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        intersectionToVertexBC_.boundaryTypes(values, vertex);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param intersection specifies the intersection at which boundary
     *           condition is to set
     */
    void boundaryTypes(BoundaryTypes &values, const Intersection &intersection) const
    {
#if DUNE_VERSION_NEWER(DUNE_ALUGRID, 2, 6)
        int boundaryId = intersection.impl().boundaryId();
#else
        int boundaryId = intersection.boundaryId();
#endif
        if (boundaryId < 1 || boundaryId > 4)
        {
            std::cout<<"invalid boundaryId: "<<boundaryId<<std::endl;
            DUNE_THROW(Dune::InvalidStateException, "Invalid " << boundaryId);
        }
        if (boundaryId == dirichletBoundary_)
            values.setAllDirichlet();
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
        initial_(values, globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
      * \param values Stores the Neumann values for the conservation equations in
     *               \f$ [ \textnormal{unit of conserved quantity} / (m^(dim-1) \cdot s )] \f$
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local index of the sub-control volume
     * \param boundaryFaceIdx The index of the boundary face
     *
     * The \a values store the mass flux of each phase normal to the boundary.
     * Negative values indicate an inflow.
     *
     * Depending on whether useMoles is set on true or false, the flux has to be given either in
     * kg/(m^2*s) or mole/(m^2*s) in the input file!! Convertion with molar mass obtained from fluid system FluidSystem::molarMass(nCompIdx)
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &intersection,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
#if DUNE_VERSION_NEWER(DUNE_ALUGRID, 2, 6)
        int boundaryId = intersection.impl().boundaryId();
#else
        int boundaryId = intersection.boundaryId();
#endif
        values = 0;
        if (boundaryId == injectionBottom_)
        {
            values[contiCO2EqIdx] = -injectionRate_; //see above: either give in kg/(m^2*s) or mole/(m^2*s) depending on useMoles
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
    void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    /*!
     * \brief Returns the initial phase state for a control volume.
     *
     * \param vertex The vertex
     * \param vIdxGlobal The global index of the vertex
     * \param globalPos The global position
     */
    int initialPhasePresence(const Vertex &vertex,
                             int &vIdxGlobal,
                             const GlobalPosition &globalPos) const
    { return Indices::wPhaseOnly; }

    // \}

private:
    /*!
     * \brief Evaluates the initial values for a control volume
     *
     * The internal method for the initial condition
     *
     * \param values Stores the initial values for the conservation equations in
     *               \f$ [ \textnormal{unit of primary variables} ] \f$
     * \param globalPos The global position
     */
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        Scalar temp = temperature_(globalPos);
        Scalar densityW = FluidSystem::Brine::liquidDensity(temp, 1e7);

        Scalar pl =  1e5 - densityW*this->gravity()[dimWorld-1]*(depthBOR_ - globalPos[dimWorld-1]);
        Scalar moleFracLiquidCO2 = 0.00;
        Scalar moleFracLiquidBrine = 1.0 - moleFracLiquidCO2;

        Scalar meanM =
            FluidSystem::molarMass(BrineIdx)*moleFracLiquidBrine +
            FluidSystem::molarMass(CO2Idx)*moleFracLiquidCO2;
        if(useMoles) // mole-fraction formulation
        {
            values[Indices::switchIdx] = moleFracLiquidCO2;
        }
        else // mass-fraction formulation
        {
            Scalar massFracLiquidCO2 = moleFracLiquidCO2*FluidSystem::molarMass(CO2Idx)/meanM;
            values[Indices::switchIdx] = massFracLiquidCO2;
        }
        values[Indices::pressureIdx] = pl;
    }

    Scalar temperature_(const GlobalPosition globalPos) const
    {
        Scalar T = 283.0 + (depthBOR_ - globalPos[dimWorld-1])*0.03;
        return T;
    }

    Scalar depthBOR_;
    Scalar injectionRate_;
    Scalar eps_;

    int nTemperature_;
    int nPressure_;

    std::string name_ ;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;

    int injectionTop_;
    int injectionBottom_;
    int dirichletBoundary_;
    int noFlowBoundary_;

    const IntersectionToVertexBC<TypeTag> intersectionToVertexBC_;
};
} //end namespace

#endif
