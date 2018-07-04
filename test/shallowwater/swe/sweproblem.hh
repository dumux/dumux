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

    //TODO map that saves the water flux for each boundary

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
        double x = 4.0;
        double y = 4.2;

        //TODO check the boundary type
        auto myBoundaryId = this->getBoundaryId(x,y); //TODO we need global x and y!
        auto boundaryValues = this->boundaryValuesMap_.at(myBoundaryId);

        if (boundaryValues.type == "hq-curve")
        {
            bdValue = getWQBoundaryValue(time_,myBoundaryId);
        }else{
            bdValue = getBoundaryValue(time_,myBoundaryId);
        }

        if (boundaryValues.type == "hq-curve") bdType = 2;
        if (boundaryValues.type == "depth") bdType = 2;
        if (boundaryValues.type == "discharge") bdType = 1;

        //TODO scale the boundary?
        //auto segIndex = is.boundarySegmentIndex();
        //bdValue = (bdValue / this->hA_boundarySum[bd_id])* this->hA_boundaryMap.at(segIndex);


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
                      const Scalar time,
                      const Scalar timeStepSize)
    {
        // compute the mass in the entire domain to make sure the tracer is conserved
        Scalar tracerMass = 0.0;
        time_ = time;
        timeStepSize_ = timeStepSize;

        // bulk elements
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVariables.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, curSol);

            //TODO loop over all boundary scvf
                //1) ha_boundarySum berechnen
                //2) parallel_sum von ha_boundarySum berechnen (overlapping?) ask if ghost element!
                //3) make a map ha_boundary[segIndex]

            for (auto&& scv : scvs(fvGeometry))
            {
                //check if the element is a boundary element and store the boudary condition type
                //and the boundary condition value (sum up the overall values over all partions
                // for each boundary)
            }

            for (auto&& scvf : scvfs(fvGeometry))
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
    void readBoundaryBoxFile(std::string filename)
    {
        double x0,x1,y0,y1;
        int id;
        std::string boundaryfilename;
        std::ifstream infile;
        std::string line;
        std::string comment ("#");

        infile.open(filename);
        if (!infile) {
            std::cerr << "Unable to open boundary file " << filename << std::endl;
            //TODO exit the programm
        }
        while(std::getline(infile,line))
        {
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

        //TODO read boundaryBoxesFiles
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


private:

    std::map<int,double> waterdischarge_;
    std::map<int,double> waterdepth_;

    std::vector<BoundaryBox> boundaryBoxesVec_;
    std::map<int,BoundaryValues> boundaryValuesMap_;

    int numberOfBoundaryfiles_ = 0;
    double time_ = 0;
    double timeStepSize_ = 0;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;
    static constexpr Scalar minHBoundary_ = 1.0E-6;

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
    }
};

} //end namespace Dumux

#endif
