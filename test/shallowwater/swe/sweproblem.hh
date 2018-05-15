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
 * \ingroup SWE Tests
 * \brief A  simple dam break test for the SWEs.
 */
#ifndef DUMUX_SWE_TEST_PROBLEM_HH
#define DUMUX_SWE_TEST_PROBLEM_HH

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/shallowwater/properties.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/shallowwater/swe/problem.hh>
#include <dumux/shallowwater/swe/model.hh>
#include "swetestspatialparams.hh"
#include <dumux/shallowwater/numericalfluxes/exactriemannsolver.hh>
#include <dumux/shallowwater/numericalfluxes/fluxrotation.hh>
#include <dumux/shallowwater/numericalfluxes/boundaryFluxes.hh>

namespace Dumux
{
/*!
 * \ingroup SweTests
 * \brief A simple dambreak test for the shallow water equations
 */
template <class TypeTag>
class SweTestProblem;


// Specify the properties for the problem
namespace Properties
{
    NEW_TYPE_TAG(SweTestTypeTag, INHERITS_FROM(CCTpfaModel, Swe, SweTestSpatialParams));

// Use 2d YaspGrid
SET_TYPE_PROP(SweTestTypeTag, Grid, Dune::YaspGrid<2>);

// Set the physical problem to be solved
SET_TYPE_PROP(SweTestTypeTag, Problem,SweTestProblem<TypeTag>);

SET_TYPE_PROP(SweTestTypeTag, SpatialParams, SweTestSpatialParams<TypeTag>);

} // end namespace Dumux



/*!
 * \ingroup Shallow water equations model
 * \ingroup ImplicitTestProblems
 *
 * \brief A simple dambreak test
 *
 * This problem uses the \ref SweModel
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_swe -parameterFile test_swes.input -TimeManager.TEnd 10</tt>
 *
 * where the initial time step is 0.01 seconds, and the end of the
 * simulation time is 10 seconds
 */
template <class TypeTag>
// evtl von SweProblem ableiten?
class SweTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;


    enum {
        // copy some indices for convenience
        massBalanceIdx = Indices::massBalanceIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx,
        // Grid and world dimension
        dimWorld = GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    struct boundaryInfo{
        std::vector<double> time;
        std::vector<double> value;
        std::vector<int> boundaryType;
    };

    std::vector<std::vector<double>> boundary_boxes;
    std::vector<int> bd_types;
    std::vector<boundaryInfo> boundaryValuesVector;
    std::vector<double> hA_boundarySum;
    std::map<int,double>hA_boundaryMap;

public:
    /*!
     * \brief Constructor
     *
     * \param timeManager The Dumux TimeManager for simulation management.
     * \param gridView The grid view on the spatial domain of the problem
     */
    SweTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \name set Boundary values
     */
    // \{

    /*!
     * \brief save the boundary value for a given time.
     *
     * \param globalPos The position for which the boundary type is set
     */

    int setBoundaryValues(std::vector<int> bd_types, int number_of_bd){

        std::ifstream infile;
        std::string line;
        std::string splittedLine;
        int read_bdId;
        double read_bdTime;
        double read_bdValue;
        this->boundaryValuesVector.resize(number_of_bd);
        this->hA_boundarySum.resize(number_of_bd);

        //read the boundary file
        infile.open("boundary.dat");
        if (!infile) {
            std::cerr << "Unable to open boundary file " << "bounday.dat" << std::endl;
        }
        while(std::getline(infile,line))
        {
            std::stringstream ssin(line);
            while(ssin.good()){
                while(ssin >> read_bdId >> read_bdTime >> read_bdValue){
                    this->boundaryValuesVector[read_bdId].time.push_back(read_bdTime);
                    this->boundaryValuesVector[read_bdId].value.push_back(read_bdValue);
                }
            }
        }
        this->number_of_bd = number_of_bd;
        this->bd_types = bd_types;
        return 0;
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
     * \param globalPos The position for which the boundary type is set
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;

        bcTypes.setAllNeumann();
        return bcTypes;
    }

    PrimaryVariables neumann(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);

        //we need the Riemann invariants to compute the values depending of the boundary type
        //since we use a weak imposition we do not have a dirichlet value. We impose fluxes
        //based on q,h, etc. computed with the Riemann invariants

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& nxy = scvf.unitOuterNormal();

        Scalar cellStatesLeft[4] = {0.0};
        Scalar cellStatesRight[4] = {0.0};

        cellStatesLeft[0]  = insideVolVars.getH();
        cellStatesLeft[1]  = insideVolVars.getU();
        cellStatesLeft[2]  = insideVolVars.getV();
        cellStatesLeft[3]  = insideVolVars.getBottom();


        //first case reflecting wall boundary
        cellStatesRight[0] = cellStatesLeft[0];
        cellStatesRight[1] = -cellStatesLeft[1];
        cellStatesRight[2] = -cellStatesLeft[2];
        cellStatesRight[3] = cellStatesLeft[3];

        int bdType = 0; //0 = Neumann, 1 = definded q, 2 = defined h
        double bdValue = 0.0;

        //call boundaryfluxes for computing the Riemann invariants
        boundaryFluxes(nxy,scvf.area(),minHBoundary_,
                       cellStatesLeft,cellStatesRight,bdType,bdValue);


        stateRotation(nxy,cellStatesLeft);
        stateRotation(nxy,cellStatesRight);

        Scalar riemannFlux[3] = {0.0};
        computeExactRiemann(riemannFlux,cellStatesLeft[0],cellStatesRight[0],
                            cellStatesLeft[1],cellStatesRight[1],
                            cellStatesLeft[2],cellStatesRight[2],insideVolVars.getGravity());


        rotateFluxBack(nxy,riemannFlux);

        values[massBalanceIdx] = riemannFlux[0];
        values[velocityXIdx]   = riemannFlux[1];
        values[velocityYIdx]   = riemannFlux[2];



        return values;
    }

    //! set neumann condition for phases (flux, [kg/(m^2 s)])
    /*
    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);

        //we need the Riemann invariants to compute the values depending of the boundary type
        //since we use a weak imposition we do not have a dirichlet value. We impose fluxes
        //based on q,h, etc. computed with the Riemann invariants
        values[massBalanceIdx] = 0.0;
        values[velocityXIdx] = 0.0;
        values[velocityYIdx] = 0.0;

        return values;
    }
    */

    //! do some preprocessing
    void preTimeStep(const SolutionVector& curSol,
                      const GridVariables& gridVariables,
                      const Scalar timeStepSize)
    {
        // compute the mass in the entire domain to make sure the tracer is conserved
        Scalar tracerMass = 0.0;

        // bulk elements
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            for (auto&& scv : scvs(fvGeometry))
            {
                //check if the element is a boundary element and store the boudary condition type
                //and the boundary condition value (sum up the overall values over all partions
                // for each boundary)
            }
        }
    }


    //! do some postprocessing
    void postTimeStep(const SolutionVector& curSol,
                      const GridVariables& gridVariables,
                      const Scalar timeStepSize)
    {
        // compute the mass in the entire domain to make sure the tracer is conserved
        Scalar tracerMass = 0.0;

        // bulk elements
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            for (auto&& scv : scvs(fvGeometry))
            {

            }
        }
    }

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial values for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     *
     * \param globalPos The position for which the boundary type is set
     */
    PrimaryVariables initial(const Element& element) const
    {
        PrimaryVariables values(0.0);

        values[massBalanceIdx] = 0.001;
        values[velocityXIdx] = 0.0;
        values[velocityYIdx] = 0.0;

        auto someInitSol = element.geometry().center();
        if (someInitSol[0] < 10.001)
        {
            values[massBalanceIdx] = 4.0;
        }

        return values;
    };

    // \}

private:


    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;
    static constexpr Scalar minHBoundary_ = 1.0E-6;
};

} //end namespace Dumux

#endif
