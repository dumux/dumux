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
#include <dune/grid/common/gridenums.hh> // for GhostEntity
#include <dumux/discretization/cellcentered/godunov/properties.hh>
#include <dumux/shallowwater/properties.hh>
#include <dumux/shallowwater/swe/problem.hh>
#include <dumux/shallowwater/swe/model.hh>
#include "swetestspatialparams.hh"
#include <dumux/shallowwater/numericalfluxes/exactriemannsolver.hh>
#include <dumux/shallowwater/numericalfluxes/fluxrotation.hh>
#include <dumux/shallowwater/numericalfluxes/boundaryFluxes.hh>

#include <dune/grid/uggrid.hh>

#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<map>

namespace Dumux
{
/*!
 * \ingroup SweTests
 * \brief Shallow water equations
 */
template <class TypeTag> class SweTestProblem;

// Specify the properties for the problem
namespace Properties
{

NEW_TYPE_TAG(SweTestTypeTag, INHERITS_FROM(GodunovModel, Swe));

SET_TYPE_PROP(SweTestTypeTag, Grid, Dune::UGGrid<2>);
//SET_TYPE_PROP(SweTestTypeTag, Grid, Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>);


// Set the physical problem to be solved
SET_TYPE_PROP(SweTestTypeTag, Problem,SweTestProblem<TypeTag>);

// Set the spatial parameters
SET_PROP(SweTestTypeTag, SpatialParams)
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = SweTestSpatialParams<FVGridGeometry, Scalar>;
};

} // end namespace Properties



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
// evtl von SweProblem ableiten? class SweTestProblem : public SweProblem<TypeTag>
class SweTestProblem : public SweProblem<TypeTag>
{
    //using ParentType = SweProblem<TypeTag>;
    using ParentType = SweProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using NeumannFluxes = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using Indices = typename ModelTraits::Indices;


    enum {
        // copy some indices for convenience
        massBalanceIdx = Indices::massBalanceIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx,
        // Grid and world dimension
        dimWorld = GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;


    struct BoundaryBox{
        int id;
        double x0,x1,y0,y1;
        std::string boundaryfilename;
    };

    struct BoundaryValues{
        int id;
        std::string type;
        double bdvalue_old;
        double relaxfactor = 0.05;
        std::vector<double> x;
        std::vector<double> y;
    };

public:
    /*!
     * \brief Constructor
     *
     * \param timeManager The Dumux TimeManager for simulation management.
     * \param gridView The grid view on the spatial domain of the problem
     */
    SweTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                   std::shared_ptr<typename ParentType::SpatialParams> spatialParams)
    : ParentType(fvGridGeometry, spatialParams)
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

    int setBoundaryValues(){
        this->readBoundaryBoxFile("boundary_boxes.dat");
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

    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);

        //we need the Riemann invariants to compute the values depending of the boundary type
        //since we use a weak imposition we do not have a dirichlet value. We impose fluxes
        //based on q,h, etc. computed with the Riemann invariants

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& nxy = scvf.unitOuterNormal();
        const auto& ip = scvf.ipGlobal();

        Scalar cellStatesLeft[4] = {0.0};
        Scalar cellStatesRight[4] = {0.0};

        cellStatesLeft[0]  = insideVolVars.getH();
        cellStatesLeft[1]  = insideVolVars.getU();
        cellStatesLeft[2]  = insideVolVars.getV();
        cellStatesLeft[3]  = insideVolVars.getBottom();

        //first case reflecting wall boundary (no-flow)
        cellStatesRight[0] = cellStatesLeft[0];
        cellStatesRight[1] = -cellStatesLeft[1];
        cellStatesRight[2] = -cellStatesLeft[2];
        cellStatesRight[3] = cellStatesLeft[3];

        int bdType = 0; //0 = Neumann, 1 = definded q, 2 = defined h
        double bdValue = 0.0;
        auto myBoundaryId = this->getBoundaryId(ip[0],ip[1]); //use the integration point

        if (myBoundaryId != 0){
            auto boundaryValues = this->boundaryValuesMap_.at(myBoundaryId);

            if (boundaryValues.type == "hq-curve")
            {
                bdValue = getWQBoundaryValue(time_,myBoundaryId);
                bdType = 2;
            }
            if (boundaryValues.type == "depth")
            {
                bdValue = getBoundaryValue(time_,myBoundaryId);
                bdType = 2;
            }
            if (boundaryValues.type == "discharge")
            {
                bdValue = getBoundaryValue(time_,myBoundaryId);
                using std::abs;
                if (abs(bdValue) < (1.0E-20))
                {
                    bdValue = 0.0;
                    bdType = 0;
                }else{
                    bdValue = (bdValue / this->hBoundarySum_[boundaryValues.id])
                               * hBoundarySegmentMap_.at(scvf.index());
                    bdType = 1;
                }
            }
        }else{
            bdType = 0;
            bdValue = 0.0;
        }


        //call boundaryfluxes for computing the Riemann invariants
        auto localBD = this->minHBoundary_;

        boundaryFluxes(nxy,scvf.area(),localBD,
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


    std::map<std::string,std::vector<double>> & xdmfGetVariable(const SolutionVector& curSol,
                      const GridVariables& gridVariables,
                      const Scalar time)
    {
        //get all cell data for XDMF
        auto fvGeometry = localView(this->fvGridGeometry());
        int size = this->fvGridGeometry().gridView().size(0);

        this->XDMFCellData_["h"] = std::vector<double>(size);
        this->XDMFCellData_["u"] = std::vector<double>(size);
        this->XDMFCellData_["v"] = std::vector<double>(size);
        this->XDMFCellData_["z"] = std::vector<double>(size);
        this->XDMFCellData_["theta"] = std::vector<double>(size);

        // bulk elements
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            //get local index of the element
            const auto eIdx = this->fvGridGeometry().elementMapper().index(element);

            //check if the actual element is a ghost
            if (element.partitionType() != Dune::GhostEntity)
            {
                for (auto&& scv : scvs(fvGeometry))
                {
                    auto h = elemVolVars[scv].getH();
                    auto u = elemVolVars[scv].getU();
                    auto v = elemVolVars[scv].getV();
                    auto z = elemVolVars[scv].getBottom();

                    //TODO here we assume that we have only one scv this is ugly
                    XDMFCellData_["h"][eIdx] = h;
                    XDMFCellData_["u"][eIdx] = u;
                    XDMFCellData_["v"][eIdx] = v;
                    XDMFCellData_["z"][eIdx] = z;
                }
            }
        }

        return XDMFCellData_;
    }


    //! do some preprocessing
    void preTimeStep(const SolutionVector& curSol,
                      const GridVariables& gridVariables,
                      const Scalar time,
                      const Scalar timeStepSize)
    {
        // compute the mass in the entire domain to make sure the tracer is conserved
        Scalar tracerMass = 0.0;
        time_ = time;
        timeStepSize_ = timeStepSize;
        double isNoGhost = true;

        this->hBoundarySegmentMap_.clear();

        //clear the existing map and create a new one
        std::fill(this->hBoundarySum_.begin(), this->hBoundarySum_.end(), 0);

        // bulk elements
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            //check if the actual element is a ghost
            if (element.partitionType() == Dune::GhostEntity)
            {
                isNoGhost = false;
            }else
            {
                isNoGhost = true;
            }

            for (auto&& scv : scvs(fvGeometry))
            {
                auto h = elemVolVars[scv].getH();
                auto u = elemVolVars[scv].getU();
                auto v = elemVolVars[scv].getV();

                for (auto&& scvf : scvfs(fvGeometry))
                {
                    if (scvf.boundary()){
                       if (this->hBoundarySegmentMap_.find(scvf.index()) == this->hBoundarySegmentMap_.end())
                       {
                            if (h > this->minHBoundary_){
                                this->hBoundarySegmentMap_[scvf.index()] = scvf.area() * h;

                                //save for sum
                                if (isNoGhost){
                                    auto ip = scvf.ipGlobal();
                                    auto myBoundaryId = this->getBoundaryId(ip[0],ip[1]); //use the integration point
                                    //if not no-flow
                                    if (myBoundaryId > 0)
                                    {
                                        auto boundaryValues = this->boundaryValuesMap_.at(myBoundaryId);
                                        this->hBoundarySum_[boundaryValues.id] += scvf.area() * h;
                                    }
                                }
                            }
                       }
                    }

                }
            }
        }
        for(std::vector<int>::size_type i = 0; i != this->hBoundarySum_.size(); i++){
            this->hBoundarySum_[i] = this->fvGridGeometry().gridView().comm().sum(this->hBoundarySum_[i]);
        }
    }


    //! do some postprocessing
    void postTimeStep(const SolutionVector& curSol,
                      const GridVariables& gridVariables,
                      const Scalar timeStepSize)
    {
        // compute the mass in the entire domain to make sure the tracer is conserved
        Scalar tracerMass = 0.0;
        bool isNoGhost = false;

        //clear the existing map and create a new one
        std::fill(this->hBoundaryFlux_.begin(), this->hBoundaryFlux_.end(), 0);

        // bulk elements
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            //check if the actual element is a ghost
            if (element.partitionType() == Dune::GhostEntity)
            {
                isNoGhost = false;
            }else
            {
                isNoGhost = true;
            }

            for (auto&& scv : scvs(fvGeometry))
            {
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    auto h = elemVolVars[scv].getH();
                    auto u = elemVolVars[scv].getU();
                    auto v = elemVolVars[scv].getV();

                    if (scvf.boundary()){
                        if (this->hBoundarySegmentMap_.find(scvf.index()) == this->hBoundarySegmentMap_.end())
                        {
                            if (isNoGhost){

                                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                                const auto& insideVolVars = elemVolVars[insideScv];
                                const auto& nxy = scvf.unitOuterNormal();
                                const auto& ip = scvf.ipGlobal();

                                Scalar cellStatesLeft[4] = {0.0};
                                Scalar cellStatesRight[4] = {0.0};

                                cellStatesLeft[0]  =  elemVolVars[scv].getH();
                                cellStatesLeft[1]  =  elemVolVars[scv].getU();
                                cellStatesLeft[2]  =  elemVolVars[scv].getV();
                                cellStatesLeft[3]  =  elemVolVars[scv].getBottom();

                                //first case reflecting wall boundary (no-flow)
                                cellStatesRight[0] = cellStatesLeft[0];
                                cellStatesRight[1] = -cellStatesLeft[1];
                                cellStatesRight[2] = -cellStatesLeft[2];
                                cellStatesRight[3] = cellStatesLeft[3];

                                int bdType = 0; //0 = Neumann, 1 = definded q, 2 = defined h
                                double bdValue = 0.0;

                                auto myBoundaryId = this->getBoundaryId(ip[0],ip[1]); //use the integration point
                                auto boundaryValues = this->boundaryValuesMap_.at(myBoundaryId);

                                if (myBoundaryId != 0){
                                    if (boundaryValues.type == "hq-curve")
                                    {
                                        bdValue = getWQBoundaryValue(time_,myBoundaryId);
                                        bdType = 2;
                                    }
                                    if (boundaryValues.type == "depth")
                                    {
                                        bdValue = getBoundaryValue(time_,myBoundaryId);
                                        bdType = 2;
                                    }
                                    if (boundaryValues.type == "discharge")
                                    {
                                        bdValue = getBoundaryValue(time_,myBoundaryId);
                                        using std::abs;
                                        if (abs(bdValue) < (1.0E-20))
                                        {
                                            bdValue = 0.0;
                                            bdType = 0;
                                        }else{
                                            bdValue = (bdValue / this->hBoundarySum_[boundaryValues.id])
                                                      * hBoundarySegmentMap_.at(scvf.index());
                                            bdType = 1;
                                        }
                                    }

                                    //call boundaryfluxes for computing the Riemann invariants
                                    boundaryFluxes(nxy,scvf.area(),0.2,
                                                   cellStatesLeft,cellStatesRight,bdType,bdValue);

                                    stateRotation(nxy,cellStatesLeft);
                                    stateRotation(nxy,cellStatesRight);

                                    Scalar riemannFlux[3] = {0.0};
                                    computeExactRiemann(riemannFlux,cellStatesLeft[0],cellStatesRight[0],
                                        cellStatesLeft[1],cellStatesRight[1],
                                        cellStatesLeft[2],cellStatesRight[2],insideVolVars.getGravity());

                                    rotateFluxBack(nxy,riemannFlux);

                                    //Store the computed flux of the Riemann solver
                                    this->hBoundaryFlux_[boundaryValues.id] += scvf.area() * riemannFlux[0];
                                }
                            }
                        }
                    }
                }
            }
        }
        for(std::vector<int>::size_type i = 0; i != this->hBoundaryFlux_.size(); i++){
            this->hBoundaryFlux_[i] = this->fvGridGeometry().gridView().comm().sum(this->hBoundaryFlux_[i]);
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
        auto elemId = this->fvGridGeometry().elementMapper().index(element);

        values[massBalanceIdx] = hInit_[elemId];
        values[velocityXIdx] = uInit_[elemId];
        values[velocityYIdx] = vInit_[elemId];

        return values;
    };

        /*!
     * \brief Evaluate the source term for all balance equations within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * For this method, the \a values parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
     NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);

        const auto& globalPos = scv.center();
        const auto& volVars = elemVolVars[scv];
        //const auto K = this->spatialParams().permeability(element, scv, EmptyElementSolution{});


        return source;
    }



    // \}
    void readBoundaryBoxFile(std::string filename)
    {
        double x0,x1,y0,y1;
        int id;
        std::string boundaryfilename;
        std::ifstream infile;
        std::string line;
        std::string comment ("#");
        std::cout << "Debug, start reading boundary file " << std::flush;;

        infile.open(filename);
        if (!infile) {
            std::cerr << "Unable to open boundary file " << filename << std::endl;
            //TODO exit the programm
        }
        while(std::getline(infile,line))
        {
            std::cout << "Debug, start reading boundary file " << std::endl;
            std::stringstream ssin(line);
            while(ssin.good()){
                //check if comment line
                std::size_t found = line.find(comment);
                if (found == std::string::npos){
                    while(ssin >> x0 >> x1 >> y0 >> y1 >> id >> boundaryfilename ){
                        BoundaryBox myBox;
                        myBox.x0 = x0;
                        myBox.x1 = x1;
                        myBox.y0 = y0;
                        myBox.y1 = y1;
                        myBox.id = id;
                        myBox.boundaryfilename = boundaryfilename;
                        boundaryBoxesVec_.push_back(myBox);
                        numberOfBoundaryfiles_ = std::max(id,numberOfBoundaryfiles_);
                    }
                }else{
                    ssin.ignore(256);
                }
            }
        }
        this->readBoundaryValuesFiles();
    }

    int getBoundaryId(double x, double y) const
    {
        int searchedId = 0;

        //loop over all boxes
        for(std::vector<int>::size_type i = 0; i != boundaryBoxesVec_.size(); i++) {
            auto x0 = std::min(boundaryBoxesVec_[i].x0,boundaryBoxesVec_[i].x1);
            auto x1 = std::max(boundaryBoxesVec_[i].x0,boundaryBoxesVec_[i].x1);
            auto y0 = std::min(boundaryBoxesVec_[i].y0,boundaryBoxesVec_[i].y1);
            auto y1 = std::max(boundaryBoxesVec_[i].y0,boundaryBoxesVec_[i].y1);
            int id = boundaryBoxesVec_[i].id;

            searchedId = 0;
            if (((x >= x0)&&(x <= x1)) && ((y >= y0)&&(y <= y1))){
                searchedId = id;
                return searchedId;
            }
        }
        return searchedId;

    }


    //get the boundary value for a given time and id
    double getBoundaryValue(auto m_time, int idPosition) const
    {
        auto boundaryValues = this->boundaryValuesMap_.at(idPosition);
        int steps = boundaryValues.x.size();
        double bd_value;

        if  (boundaryValues.type == "hq-curve")
        {
            std::cerr << "This function does not work for bdtype: hq-curve " << std::endl;
        }else{
            for (int j=0; j< steps-1 ;++j){
                double t1 = boundaryValues.x[j];
                double t2 = boundaryValues.x[j+1];
                double v1 = boundaryValues.y[j];
                double v2 = boundaryValues.y[j+1];

                if ((t1 <= m_time)&&(m_time <= t2)){
                    if (v1 == v2){
                        bd_value = v2;
                    }else{
                        bd_value = ((m_time-t1)/((t2-t1)/(v2-v1))) + v1;
                    }
                    break;
                }
            }
        }

        //some extra cases outside the time space we retrun first/last values
        if (m_time < boundaryValues.x[0]) bd_value = boundaryValues.y[0];
        if (m_time > boundaryValues.x.back()) bd_value = boundaryValues.y.back();

        return bd_value;
    }

    //get the boundary value for a given time and id for hq we also give the oldvalue
    double getWQBoundaryValue(double sim_q, int idPosition) const
    {
        auto boundaryValues = this->boundaryValuesMap_.at(idPosition);
        int steps = boundaryValues.x.size();
        double bd_value;

        if  (boundaryValues.type == "hq-curve")
        {
            for (int j=0; j< steps-1 ;++j){
                //time means here discharge,
                double q1 = boundaryValues.x[j];
                double q2 = boundaryValues.x[j+1];
                double w1 = boundaryValues.y[j];
                double w2 = boundaryValues.y[j+1];

                if ((q1 <= sim_q)&&(sim_q <= q2)){
                    bd_value = ((sim_q-q1)/((q2-q1)/(w2-w1))) + w1;
                    double h_new = bd_value;
                    bd_value = boundaryValues.bdvalue_old +
                               boundaryValues.relaxfactor *
                               (h_new - boundaryValues.bdvalue_old);
                    break;
                }
            }
        }else{
            std::cerr << "This function only works for bdtype: hq-curve " << std::endl;               //TODO exit the programm
        }

        //some extra cases outside the time space we retrun first/last values
        if (sim_q < boundaryValues.x[0]) bd_value = boundaryValues.y[0];
        if (sim_q > boundaryValues.x.back()) bd_value = boundaryValues.y.back();

        return bd_value;
    }

    //set the input data from map (we get this map from the XDMF/HDF5-Reader)
    void setInputData(auto elementdata)
    {
        //read the init data
        if ((elementdata.find("h") != elementdata.end())&&
            (elementdata.find("u") != elementdata.end())&&
            (elementdata.find("v") != elementdata.end()))
        {
            hInit_ = elementdata["h"];
            uInit_ = elementdata["u"];
            vInit_ = elementdata["v"];
        }else{
            std::cout << "Can not find initial data in elementdata"<< std::endl;
        }
        std::cout << "alive in inputdata with h[0] " << hInit_[0] << std::endl;
    }

private:

    std::map<int,double> waterdischarge_;
    std::map<int,double> waterdepth_;

    std::map<int,double> hBoundarySegmentMap_;
    std::vector<double> hBoundarySum_;
    std::vector<double> hBoundaryFlux_;

    std::map<std::string, std::vector<double>> XDMFCellData_;

    std::vector<BoundaryBox> boundaryBoxesVec_;
    std::map<int,BoundaryValues> boundaryValuesMap_;

    int numberOfBoundaryfiles_ = 0;
    double time_ = 0;
    double timeStepSize_ = 0;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;
    static constexpr Scalar minHBoundary_ = 1.0E-6;

    /*Initial data */
    std::vector<double> hInit_;
    std::vector<double> uInit_;
    std::vector<double> vInit_;


    /*!
     * \brief read boundary value files.
     *
     *
     * \param
     */
    //Read in boundary files automatically!
    void readBoundaryValuesFiles()
    {
        std::ifstream infile;
        std::vector<BoundaryValues> boundaryValues;
        std::string line;
        std::string comment ("#");
        std::string bdtypeString ("bdtype");
        std::string bdidString ("bdid");
        std::string initialvalueString ("initialvalue");

        std::string typeDepth ("depth");
        std::string typeDischarge ("discharge");
        std::string typeHQcurver ("hq-curve");
        std::string dummy;

        std::vector<double> x;
        std::vector<double> y;
        double inx,iny;
        std::string actualfile;
        int maxBoundaryId = 1;


        //open an read all files n = numberOfBoundaryfiles there might be more boxes as files!
        for (std::vector<int>::size_type i = 0; i != boundaryBoxesVec_.size(); i++){

            infile.open(this->boundaryBoxesVec_[i].boundaryfilename);
            if (!infile) {
                std::cerr << "\nError in readBoundaryValuesFiles() Unable to open boundary file "
                          << this->boundaryBoxesVec_[i].boundaryfilename << std::endl;
                //TODO exit the programm
            }
            BoundaryValues myValues;

            while(std::getline(infile,line))
            {

                std::stringstream ssin(line);
                while(ssin.good()){
                    //check if comment line
                    std::size_t found = line.find(comment);
                    if (found == std::string::npos){

                        //check for bdtype
                        found = line.find(bdtypeString);
                        if (found != std::string::npos){
                            ssin >> dummy >>  myValues.type;
                        }

                        //check for bdid
                        found = line.find(bdidString);
                        if (found != std::string::npos){
                            ssin >> dummy >>  myValues.id;
                            maxBoundaryId = std::max(myValues.id,maxBoundaryId);
                        }

                        //check for  initialvalue
                        found = line.find(initialvalueString);
                        if (found != std::string::npos){
                            ssin >> dummy >>  myValues.bdvalue_old >> myValues.relaxfactor;
                        }
                        //else read in variables
                        while(ssin >> inx >> iny){
                            myValues.x.push_back(inx);
                            myValues.y.push_back(iny);
                        }
                    }else{
                        ssin.ignore(256);
                    }
                }
            }
            infile.close();
            this->boundaryValuesMap_[myValues.id] = myValues;
        }

        //resize hBoundarySum_
        maxBoundaryId = this->fvGridGeometry().gridView().comm().max(maxBoundaryId);
        this->hBoundarySum_.resize(maxBoundaryId +1);
        this->hBoundaryFlux_.resize(maxBoundaryId +1);

    }
};

} //end namespace Dumux

#endif
