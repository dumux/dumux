// $Id: outputToFile.hh 3736 2010-06-15 09:52:10Z pnuske $
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief A collection of functions for output purposes.
 *          These output files are meant for visualization with another program (matlab, gnuplot...)
 *
 */
#ifndef DUMUX_OUTPUTTOFILE_HH
#define DUMUX_OUTPUTTOFILE_HH

#include <dumux/common/valgrind.hh>
#include <dumux/common/propertysystem.hh>


#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
namespace Dumux
{
namespace Properties
{
    NEW_PROP_TAG(Scalar);
    NEW_PROP_TAG(Problem);
    NEW_PROP_TAG(GridView);
    NEW_PROP_TAG(DofMapper);
    NEW_PROP_TAG(ElementSolutionVector);
    NEW_PROP_TAG(SolutionVector);
    NEW_PROP_TAG(FVElementGeometry);
    NEW_PROP_TAG(TwoPIAIndices);
    NEW_PROP_TAG(NumEq);
    NEW_PROP_TAG(MaterialLaw);
    NEW_PROP_TAG(ElementVolumeVariables);

    NEW_PROP_TAG(AwnSurface);
    NEW_PROP_TAG(AwnSurfaceParams);


}

/*!
 * \brief A collection of functions for output purposes.
 *          These output files are meant for visualization with another program (matlab, gnuplot...)
 */
template<class TypeTag>
class OutputToFile
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(DofMapper)) DofMapper;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSolutionVector)) ElementSolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementVolumeVariables)) ElementVolumeVariables;



    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIAIndices)) Indices;
    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),

        // copy some indices for convenience
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        pcIdx = Indices::pcIdx, // in ia models cap pressure is a primary Variable

        pW = Indices::pW,
        sN = Indices::sN,

        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*
     * \brief A function that puts Sw, pc, awn data (for one or more points in the simulation domain) into columns and saves it to a text file.
     *          Each line contains the resulting values of a time step.
     *          These output files are meant for visualization with another program (matlab, gnuplot...)
     */
    void smallResultWriter(Problem & problem)
    {
        //BEWARE globalIdx is not the same as paraview index !!! therefore selcection of the node has to be done via coordinates :-(
        static bool FirstFlag = true;

        typedef typename GridView::template Codim<0>::Entity Element;
        typedef typename GridView::template Codim<0>::Iterator ElementIterator;
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(DofMapper)) DofMapper;
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(ElementSolutionVector)) ElementSolutionVector;
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(SolutionVector)) SolutionVector;

        ElementSolutionVector tmpSol;

        SolutionVector & globalSol = problem.model().curSol();
    //                SolutionVector & oldGlobalSol = problem.model().prevSol();

        const DofMapper &dofMapper = problem.model().dofMapper();
        Scalar satWIdx[10];
        Scalar pcIdx[10];
        Scalar awnIdx[10];
        Scalar pcBounding[10];
        Scalar dSw[10];
        Scalar dt[10];
        Scalar dSwdt[10];
    //                Scalar x1(3.), y1(3.);
        Scalar x1(-101.), y1(1.), z1(0.);
    //                Scalar x1(-99.8), y1(1.), z1(30.);
        int numKnots; // for output loop (more than one node can be put to the file)
        int scvIdx;

        FVElementGeometry fvElemGeom;
        ElementVolumeVariables elemVolVars;

        ElementIterator endit = problem.gridView().template end<0>();
        for (ElementIterator elementIt = problem.gridView().template begin<0>() ; elementIt != endit; ++elementIt)
        {

            fvElemGeom.update(problem.gridView(), *elementIt);
            elemVolVars.update(problem, *elementIt, fvElemGeom, false);
            Scalar SwOld = 0;

            int numVerts = elementIt->template count<dim>();
            for (int i = 0; i < numVerts; ++i)
            {
                scvIdx = i;

                const MaterialLawParams &pcParamsDrainageInit =
                                problem.spatialParameters().materialLawParamsDrainage(*elementIt, fvElemGeom, scvIdx);

                numKnots = 0;
                int globalIdx = dofMapper.map(*elementIt, i, dim);

                const GlobalPosition & globalPos = fvElemGeom.subContVol[i].global;

                //BEWARE globalIdx is not the same as paraview index !!!
    //                        if(globalPos[0]<x1+eps_ and globalPos[0]>x1-eps_ and globalPos[1]<y1+eps_ and globalPos[1]>y1-eps_) // 2D
                if(globalPos[0]<x1+eps_ and globalPos[0]>x1-eps_ and globalPos[1]<y1+eps_ and globalPos[1]>y1-eps_ and globalPos[2]<z1+eps_ and globalPos[2]>z1-eps_)// 3D
                {
                    satWIdx[numKnots]   = elemVolVars[i].saturation(wPhaseIdx);
                    pcIdx[numKnots]     = elemVolVars[i].capillaryPressure();
                    pcBounding[numKnots]= MaterialLaw::pC(pcParamsDrainageInit, satWIdx[numKnots]); // this is the pc on the bounding curve
                    awnIdx[numKnots]    = elemVolVars[i].awn();

                    dSw[numKnots]         = (elemVolVars[i].saturation(wPhaseIdx) - SwOld );
    //                            Scalar Sw = problem.model().localJacobian().curElemDat()[i].saturation(wPhaseIdx);
    //                            cout<<"Sw: " << Sw << "SwOld: " << SwOld << "dSw " << dSw[numKnots] <<"\n";
                    dt[numKnots]         = problem.timeManager().timeStepSize();
                    dSwdt[numKnots]        = dSw[numKnots] / dt[numKnots];
                }
                numKnots++;
    //                        if (globalIdx == 4) use COORDINATES !!!
    //                        {
    //
    ////                            problem.model().localJacobian().restrictToElement(tmpSol, oldGlobalSol);
    ////                            problem.model().localJacobian().setPreviousSolution(tmpSol);
    ////
    ////                            Scalar SwOld = problem.model().localJacobian().prevElemDat()[i].saturation(wPhaseIdx);
    ////
    ////                            problem.model().localJacobian().restrictToElement(tmpSol, globalSol);
    ////                            problem.model().localJacobian().setCurrentSolution(tmpSol);
    //
    //                            satWIdx[numKnots]      = problem.model().localJacobian().curElemDat()[i].saturation(wPhaseIdx);
    //                            pcIdx[numKnots]     = problem.model().localJacobian().curElemDat()[i].capillaryPressure();
    //                            pcBounding[numKnots]= MaterialLaw::pC(pcParamsDrainageInit, satWIdx[numKnots]); // this is where it should be
    //                            awnIdx[numKnots]    = problem.model().localJacobian().curElemDat()[i].awn();
    //
    //                            dSw[numKnots]         = (problem.model().localJacobian().curElemDat()[i].saturation(wPhaseIdx)
    //                                                        - SwOld );
    ////                            Scalar Sw = problem.model().localJacobian().curElemDat()[i].saturation(wPhaseIdx);
    ////                            cout<<"Sw: " << Sw << "SwOld: " << SwOld << "dSw " << dSw[numKnots] <<"\n";
    //                            dt[numKnots]         = problem.model().problem().timeManager().timeStepSize();
    //                            dSwdt[numKnots]        = dSw[numKnots] / dt[numKnots];
    //                        }
    //                        numKnots++;
    //                        if (globalIdx == 51) use COORDINATES !!!
    //                        {
    //
    ////                            problem.model().localJacobian().restrictToElement(tmpSol, oldGlobalSol);
    ////                            problem.model().localJacobian().setPreviousSolution(tmpSol);
    ////
    ////                            Scalar SwOld = problem.model().localJacobian().prevElemDat()[i].saturation(wPhaseIdx);
    ////
    ////                            problem.model().localJacobian().restrictToElement(tmpSol, globalSol);
    ////                            problem.model().localJacobian().setCurrentSolution(tmpSol);
    //
    //                            satWIdx[numKnots]      = problem.model().localJacobian().curElemDat()[i].saturation(wPhaseIdx);
    //                            pcIdx[numKnots]     = problem.model().localJacobian().curElemDat()[i].capillaryPressure();
    //                            pcBounding[numKnots]= MaterialLaw::pC(pcParamsDrainageInit, satWIdx[numKnots]); // this is where it should be
    //                            awnIdx[numKnots]    = problem.model().localJacobian().curElemDat()[i].awn();
    //
    //                            dSw[numKnots]         = (problem.model().localJacobian().curElemDat()[i].saturation(wPhaseIdx)
    //                                                        - SwOld );
    ////                            Scalar Sw = problem.model().localJacobian().curElemDat()[i].saturation(wPhaseIdx);
    ////                            cout<<"Sw: " << Sw << "SwOld: " << SwOld << "dSw " << dSw[numKnots] <<"\n";
    //                            dt[numKnots]         = problem.model().problem().timeManager().timeStepSize();
    //                            dSwdt[numKnots]        = dSw[numKnots] / dt[numKnots];
    //                        }
    //                        numKnots++;
            }
        }

    //                Scalar time = problem.timeManager().time();
    //                Scalar timeStep = problem.timeManager().timeStepSize();
    //                Scalar co2Mass = mass[0];
    //                Scalar brineMass = mass[1];
    //                Scalar co2Flux = flux[1];
        int timeStepIndex = problem.timeManager().timeStepIndex();
    //                Scalar CPUtime = problem.timeManager().timer.elapsed();
    //                Dune::Timer timer;
    //                Scalar CPUtime = timer.elapsed();

        std::string dataFileName = problem.name();
        std::string fileName = (boost::format("%s.dat")%dataFileName).str();


        std::string dataFileName2 = problem.name();
        std::string fileName2 = (boost::format("%sdSwdt.dat")%dataFileName).str();


        std::ofstream dataFile;
        std::ofstream dataFile2;

        if (FirstFlag == true)
        {
            dataFile.open(fileName.c_str());
            dataFile << "# This is an DuMuX output file for further processing with the preferred graphics program of your choice. \n";
            dataFile << "# These Numbers are printed in order to check whether the data points are one the bounding curve, as they should be.\n";

            dataFile << "# This output file was written from "<< __FILE__ << ", line " <<__LINE__ << "\n";
            dataFile << "# This output file was generated at " << __TIME__ <<", "<< __DATE__<< "\n";
            dataFile << "\n";
            dataFile << "# Header\n";
    //                    dataFile << "# \ttime \ttimestep \ttotal_CO2_mass \ttotal_brine_mass \tCO2_flux_across_z_=_80_m \tCPUtime" <<std::endl;
            dataFile << "#timestep \tSwIdx0 \t\tpcIdx0 \t\tpcBoundingIdx0 " << std::endl;
            dataFile.close();


            dataFile2.open(fileName2.c_str());
            dataFile2 << "# This is an DuMuX output file for further processing with the preferred graphics program of your choice. \n";
            dataFile2 << "# These Numbers are to be printed in order to be compared with the results of benchmark11 of Class et.al.(2009) A benchmark study on problems related to CO2 storage in geological formations.\n";

            dataFile2 << "# This output file was written from "<< __FILE__ << ", line " <<__LINE__ << "\n";
            dataFile2 << "# This output file was generated at " << __TIME__ <<", "<< __DATE__<< "\n";
            dataFile2 << "\n";
            dataFile2 << "# Header\n";
    //                    dataFile << "# \ttime \ttimestep \ttotal_CO2_mass \ttotal_brine_mass \tCO2_flux_across_z_=_80_m \tCPUtime" <<std::endl;
            dataFile2 << "#timestep \tdSwIdx0 \t\tdtIdx0 \t\tdSwdtIdx0 " << std::endl;
            dataFile2.close();

            FirstFlag = false ;
        }
        dataFile.open(fileName.c_str() , std::ios::app);
    //                dataFile <<timeStepIndex<<"\t" << time <<"\t"<< timeStep <<"\t\t"<< co2Mass <<"\t\t"<< brineMass <<"\t\t"<< co2Flux <<"\t\t\t"<<"\t"<<CPUtime<< "\n";
        dataFile << timeStepIndex ;
        for (int i=0; i<numKnots; i++)
        {
            dataFile <<"\t"<<satWIdx[i] <<"\t"<< pcIdx[i]  <<"\t"<< pcBounding[i] <<"\t"<< awnIdx[i] ;
        }
        dataFile <<"\n";
        dataFile.close();

        dataFile2.open(fileName2.c_str() , std::ios::app);
        dataFile2 << timeStepIndex ;
        for (int i=0; i<numKnots; i++)
        {
    //                    dataFile2 <<"\t"<< dSw[i] <<"\t"<< dt[i]  <<"\t"<< dSwdt[i] ;
        }
        dataFile2 <<"\n";
        dataFile2.close();
    return;
    }


    /*
     * \brief A function that puts Sw, pc or Sw, pc, awn into a file.
     *          Right now only the derivatis of the awn surface are evaluated.
     *          This is for checking whether the parameters (e.g. a1...6 of the awn surface) are put in correctly, The regularization is smooth etc...
     *          These output files are meant for visualization with another program (matlab, gnuplot...)
     */
    void smallMaterialLawWriter(Problem & problem)
    {
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(AwnSurface)) AwnSurface;
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(AwnSurfaceParams)) AwnSurfaceParams;

        int scvIdx;
        int numSteps = 21; // this must not be to big -> segfault (if it has to be big: may be make std::)vector
        Scalar Swmax = 1. , Swmin =0.;
        Scalar pcmax = 12000., pcmin =0.;
        Scalar deltaSw = (double) (Swmax - Swmin) / numSteps;
        Scalar deltapc = (double) (pcmax - pcmin) / numSteps;
        const int sizeArray = numSteps*numSteps ;
        double Sw[sizeArray] ;
        double pc[sizeArray] ;
        double dawndpc[sizeArray], dawndSw[sizeArray] ;
        double SwTemp(0.), pcTemp(0.), awnTemp(0.), dawndpcTemp(0.), dawndSwTemp(0.);

        AwnSurfaceParams awnSurfaceParams_;
        MaterialLawParams pcParamsDrainageInit;


        FVElementGeometry fvElemGeom;

        ElementIterator endit = problem.gridView().template end<0>();
        for (ElementIterator elementIt = problem.gridView().template begin<0>() ; elementIt != endit; ++elementIt)
        {

            fvElemGeom.update(problem.gridView(), *elementIt);

            int numVerts = elementIt->template count<dim>();


            for (int i = 0; i < numVerts; ++i)
            {
                scvIdx = i;
                awnSurfaceParams_ = this->spatialParameters().awnSurfaceParams(*elementIt, fvElemGeom, scvIdx);

                pcParamsDrainageInit =
                                problem.spatialParameters().materialLawParamsDrainage(*elementIt, fvElemGeom, scvIdx);
            }
        }



        static bool FirstFlag = true;

        std::string dataFileName = this->name();
        std::string fileName = (boost::format("%s_dawndSw_dawndpc.dat")%dataFileName).str();
        std::ofstream dataFile;

        for(int horse=0; horse <numSteps; horse++ )
        {
            for (int deer=0; deer < numSteps; deer++)
            {
                if (FirstFlag == true)
                {
                    dataFile.open(fileName.c_str());
                    dataFile << "# This is an DuMuX output file for further processing with the preferred graphics program of your choice. \n";
                    dataFile << "# These Numbers are for comparison whether the input material laws (coefficients) are 'correct', i.e. visualization with another program.\n";

                    dataFile << "# This output file was written from "<< __FILE__ << ", line " <<__LINE__ << "\n";
                    dataFile << "# This output file was generated at " << __TIME__ <<", "<< __DATE__<< "\n";
                    dataFile << "\n";
                    dataFile << "# Header\n";
                    dataFile << "# \tSw \t\tpc \t\t dawndSw \t\t dawndpc " <<std::endl;
                    dataFile.close();
                    FirstFlag = false ;
                }


                dawndSwTemp = AwnSurface::dawndSw(awnSurfaceParams_, SwTemp, pcTemp ) ;
                dawndpcTemp = AwnSurface::dawndpc(awnSurfaceParams_, SwTemp, pcTemp ) ;

                Sw[horse*numSteps + deer]         = SwTemp;
                pc[horse*numSteps + deer]         = pcTemp;

                dawndSw[horse*numSteps + deer]     = dawndSwTemp;
                dawndpc[horse*numSteps + deer]     = dawndpcTemp;

                SwTemp += deltaSw;

            }

            SwTemp=0.;
            pcTemp += deltapc;
        }
        dataFile.open(fileName.c_str() , std::ios::app);
        for (int mule =0 ; mule < sizeArray; mule++)
        {
            dataFile << Sw[mule] <<"\t" << pc[mule] << "\t" << dawndSw[mule] << "\t" << dawndpc[mule] << "\n";
        }
        dataFile.close();
        exit (1);
        return ;
    }

private:
        static const Scalar eps_=1e-6 ;

};
}
#endif
