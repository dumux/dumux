// $Id$
#ifndef DUNE_TIMELOOPSUBPROBS_HH
#define DUNE_TIMELOOPSUBPROBS_HH

#include "dumux/timedisc/timestep.hh"
#include "dumux/timedisc/rungekuttastep.hh"

namespace Dune
{
template<class G, class Model>
class TimeLoopSubProbs
{
typedef    typename G::ctype ct;
    enum
    {    n = G::dimension, dimworld = G::dimensionworld};

    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename G::HostGridType HG;
    typedef typename HG::template Codim<0>::EntityPointer HostEntityPointer;

public:

    void calculatemacrodispersion(Model&, const int, int&, FieldVector<ct,n>&, FieldVector<ct,n>&);
    void calculateconvectivecorrection(Model&, const int, int&, int&, FieldVector<ct,n>&, FieldVector<ct,n>&, FieldVector<ct,n>&);

    void executeD(Model& model, const int& subProblemNumber, int& coarseIndex, FieldVector<ct,n>& lowerLeft, FieldVector<ct,n>& upperRight, int& nDt = *(new int(1e2)), double& tEnd = *(new double(1e6)), bool firstRun=false)
    {
        // initialize solution with initial values
        model.initial();

        // now do the time steps
        double t = tStart_;
        int k = 0;

        if (firstRun)
        {
            while(t < tEnd_)
            {
                k++;

                timeStep_.execute(model, t, dt_, maxDt_, tEnd_, cFLFactor_);

                // generate output
                t += dt_;
                t = std::min(t, tEnd_);
                std::cout << ",subproblem timestep: " << k << "\t t=" << t
                << "\t dt_=" << dt_ << std::endl;

                if (t == tEnd_)
                {
                    for (int i = 0; i < model.variables().diffSize; i++)
                    {
                        meanSat_ += (*model)[i];
                    }
                    meanSat_ /= model.variables().diffSize;
                    std::cout<<"meansat = "<<meanSat_<<std::endl;
                    if (meanSat_ < minMaxMeanSat_)
                    {
                        model.initial();
                        k=0;
                        t=tStart_;
                        tEnd_+=2*dt_;
                        dt_=1e100;
                    }
                }

            }
        }
        else
        {
            bool reachedtEnd = false;
            while (k < nDt)
            {
                k++;

                timeStep_.execute(model, t, dt_, maxDt_, tEnd_, cFLFactor_);

                calculatemacrodispersion(model, subProblemNumber, coarseIndex, lowerLeft, upperRight);

                // generate output
                t += dt_;
                t = std::min(t, tEnd_);

                if (t==tEnd_)
                {
                    reachedtEnd = true;
                }

                if (reachedtEnd && meanSat_ <= 0.99)
                {
                    tEnd_+=dt_;
                }
                if (reachedtEnd && meanSat_> 0.99)
                {
                    break;
                }

                //                std::cout << ",subproblem timestep: " << k << "\t t=" << t
                //                << "\t dt_=" << dt_ << std::endl;

                //                std::cout<<"saturation ="<<model.variables().saturation<<"pressure = "<<model.variables().pressure<<std::endl;
                //                VTKWriter<typename G::LeafGridView> vtkwriter(model.variables().grid.leafView());
                //                char fname[128];
                //                sprintf(fname, "dispsubprob-%05d-%05d-%05d",coarseIndex, subProblemNumber, k);
                ////                vtkwriter.addCellData(model.variables().saturation, "saturation");
                //                vtkwriter.addCellData(model.variables().pressure, "total pressure p~");
                //                vtkwriter.write(fname, VTKOptions::ascii);
            }
            std::cout << ",subproblem timestep: " << k << "\t t=" << t
            << "\t dt_=" << dt_ << std::endl;
        }

        if (firstRun)
        {
            nDt=k;
            tEnd=tEnd_;
        }
        return;
    }

    void executeM(Model& model, const int& subProblemNumber, int& coarseIndex, int& numberInSelf, FieldVector<ct, n>& faceGlobalCoarse, FieldVector<ct,n>& lowerLeft, FieldVector<ct,n>& upperRight, int& nDt = *(new int(1e2)), double& tEnd = *(new double(1e6)), bool firstRun=false)
    {
        typedef typename G::LevelGridView GV;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename GV::IndexSet IS;

        typedef typename HG::LevelGridView HGV;

        typedef typename HG::template Codim<0>::EntityPointer HostEntityPointer;
        typedef typename HGV::IndexSet HIS;
        int fineLev = model.variables().diffLevel;
        int size = model.variables().diffSize;

        const GV& gridViewCoarse(model.variables().grid.levelView(0));
        const GV& gridViewFine(model.variables().grid.levelView(fineLev));
        const HGV& gridViewHost(gridViewCoarse.grid().getHostGrid().levelView(0));

        //    const IS& coarseIset = gridViewCoarse.indexSet();
        const HIS& coarseIsetHost = gridViewHost.indexSet();

        // initialize solution with initial values
        model.initial();

        // now do the time steps
        double t = tStart_;
        int k = 0;

        if (firstRun)
        {
            while(t < tEnd_)
            {
                k++;

                timeStep_.execute(model, t, dt_, maxDt_, tEnd_, cFLFactor_);

                // generate output
                t += dt_;
                t = std::min(t, tEnd_);
                std::cout << ",subproblem timestep: " << k << "\t t=" << t
                << "\t dt_=" << dt_ << std::endl;

                if (t == tEnd_)
                {
                    meanSat_=0;
                    int count = 0;
                    Iterator endCIt = gridViewCoarse.template end<0>();
                    for (Iterator cIt = gridViewCoarse.template begin<0>(); cIt != endCIt; ++cIt)
                    {
                        const HostEntityPointer& hostCoarseIt = model.variables().grid.template getHostEntity<0>(*cIt);

                        int coarseInd = coarseIsetHost.index(*hostCoarseIt);
                        Iterator eEndIt = gridViewFine.template end<0>();
                        Iterator eBeginIt = gridViewFine.template begin<0>();
                        for (Iterator it = eBeginIt; it != eEndIt; ++it)
                        {

                            int index = model.variables().diffMapper.map(*it);
                            if (coarseInd == coarseIndex)
                            {
                                meanSat_ += (*model)[index];
                                count++;
                            }
                            //                            std::cout<<"coarseInd = "<<coarseInd<<"coarseIndex = "<<coarseIndex<<std::endl;
                        }
                    }
                    meanSat_ /= count;
                    std::cout<<"meansat = "<<meanSat_<<std::endl;
                    if (meanSat_ < minMaxMeanSat_)
                    {
                        model.initial();
                        k=0;
                        t=tStart_;
                        tEnd_+=4*dt_;
                        dt_=1e100;
                    }
                }

            }
        }
        else
        {
            bool reachedtEnd = false;
            while (k < nDt)
            {
                k++;

                timeStep_.execute(model, t, dt_, maxDt_, tEnd_, cFLFactor_);

                calculateconvectivecorrection(model, subProblemNumber, coarseIndex, numberInSelf, faceGlobalCoarse, lowerLeft, upperRight);

                // generate output
                t += dt_;
                t = std::min(t, tEnd_);

                if (t==tEnd_)
                {
                    reachedtEnd = true;
                }

                if (reachedtEnd && meanSat_ <= 0.99)
                {
                    tEnd_+=dt_;
                }
                if (reachedtEnd && meanSat_> 0.99)
                {
                    break;
                }

                //                std::cout << ",subproblem timestep: " << k << "\t t=" << t
                //                << "\t dt_=" << dt_ << std::endl;

                //                                std::cout<<"saturation ="<<model.variables().saturation<<"pressure = "<<model.variables().pressure<<std::endl;
                //                                VTKWriter<typename G::LeafGridView> vtkwriter(model.variables().grid.leafView());
                //                                char fname[128];
                //                                sprintf(fname, "convsubprob-%05d-%05d-%05d",coarseIndex, numberInSelf, k);
                //                                vtkwriter.addCellData(model.variables().saturation, "saturation");
                //                                vtkwriter.addCellData(model.variables().pressure, "total pressure p~");
                //                                vtkwriter.write(fname, VTKOptions::ascii);
            }
            std::cout << ",subproblem timestep: " << k << "\t t=" << t
            << "\t dt_=" << dt_ << std::endl;
        }

        if (firstRun)
        {
            nDt=k;
            tEnd=tEnd_;
        }
        return;
    }

    TimeLoopSubProbs(double te = 1e6,const double mMMeanSat = 0.95,
            const double cfl = 0.95, TimeStep<G,
            Model>& tist = *(new RungeKuttaStep<G, Model> (1))) :
    tEnd_(te), minMaxMeanSat_(mMMeanSat), cFLFactor_(cfl), timeStep_(tist), dt_(1e100), tStart_(0), maxDt_(1e100), timeStepCounter_(0), meanSat_(0)
    {
    }

private:
    double tEnd_;
    double minMaxMeanSat_;
    const double cFLFactor_;
    TimeStep<G, Model>& timeStep_;
    double dt_;
    const double tStart_;
    const double maxDt_;
    int timeStepCounter_;
    double meanSat_;
};
template<class G, class Model>void TimeLoopSubProbs<G,Model>::calculatemacrodispersion(Model& model, const int subProblemNumber, int& coarseIndex, FieldVector<ct,n>& lowerLeft, FieldVector<ct,n>& upperRight)
{

    typedef typename G::LevelGridView GV;
    typedef typename GV::template Codim<0>::Iterator Iterator;
    typedef typename IntersectionIteratorGetter<G,LevelTag>::IntersectionIterator
    IntersectionIterator;
    typedef typename G::HostGridType HG;
    typedef typename HG::template Codim<0>::EntityPointer HostEntityPointer;

    typedef typename Entity::Geometry Geometry;

    int fineLev = model.variables().diffLevel;
    int size = model.variables().diffSize;

    const GV& gridView(model.variables().grid.levelView(fineLev));

    meanSat_=0;

    FieldVector<double, n> meanVelocity(0);
    for (int i = 0; i < size; i++)
    {
        meanSat_ += (*model)[i];
    }
    meanSat_ /= size;

    //        std::cout<<"meanSat = "<<meanSat_<<"meanvel = "<<meanVelocity<<std::endl;

    FieldVector<double, n> meanVf(0);//mean of finescale velocity v times fine scale fractional flow function f
    double fwSmean;// fractional flow function values of mean saturation

    double satIn=0;
    double satOut=0;
    int count = 0;

    Iterator eEndIt = gridView.template end<0>();
    Iterator eBeginIt = gridView.template begin<0>();
    for (Iterator it = eBeginIt; it != eEndIt; ++it)
    {
        // element geometry
        const Geometry& geometry = it->geometry();

        GeometryType gt = geometry.type();
        // cell center in reference element
        const FieldVector<ct,n>& local = ReferenceElements<ct,n>::general(gt).position(0, 0);

        // get global coordinate of cell center
        const FieldVector<ct,n> global = geometry.global(local);

        const HostEntityPointer& hostIt = model.variables().grid.template getHostEntity<0>(*it);

        int index = model.variables().diffMapper.map(*it);
        double sat = (*model)[index];

        if (index == 1)
        fwSmean = model.diffproblem.materialLaw.fractionalW(meanSat_, global, *hostIt, local);

        FieldVector<double,n*2> fluxVector(0);

        IntersectionIterator
        endit = IntersectionIteratorGetter<G, LevelTag>::end(*it);
        for (IntersectionIterator
                is = IntersectionIteratorGetter<G, LevelTag>::begin(*it); is
                !=endit; ++is)
        {
            // get geometry type of face
            Dune::GeometryType gtf = is->intersectionSelfLocal().type();

            // center in face's reference element
            const Dune::FieldVector<ct,n-1>&
            faceLocal = Dune::ReferenceElements<ct,n-1>::general(gtf).position(0,0);

            // center of face in global coordinates
            Dune::FieldVector<ct,dimworld> faceGlobal = is->intersectionGlobal().global(faceLocal);

            double faceVol = is->intersectionGlobal().volume();

            // get normal vector
            FieldVector<double,n> unitOuterNormal = is->unitOuterNormal(faceLocal);

            int faceNumber = is->numberInSelf();

            fluxVector[faceNumber] = model.variables().velocity[index][faceNumber] * unitOuterNormal;
            fluxVector[faceNumber] *= faceVol;

            if (is->boundary())
            {
                // center of face inside volume reference element
                const Dune::FieldVector<ct,n>&
                facelocalDim = Dune::ReferenceElements<ct,n>::general(gtf).position(faceNumber,1);

                //get boundary type
                BoundaryConditions::Flags bctype = model.diffproblem.bctypeSat(faceGlobal, *it, facelocalDim);

                if (bctype == BoundaryConditions::dirichlet)
                {
                    switch(subProblemNumber)
                    {
                        case 1:
                        if (faceGlobal[0]> lowerLeft[0])
                        {
                            satOut+=sat;
                            count++;
                        }
                        else
                        {
                            satIn+=sat;
                            count++;
                        }
                        break;
                        case 2:
                        if (faceGlobal[1]> lowerLeft[1])
                        {
                            satOut+=sat;
                            count++;
                        }
                        else
                        {
                            satIn+=sat;
                            count++;
                        }
                        break;
                        default:
                        break;
                    }

                }

            }
        }
        // calculate velocity on reference element
        FieldVector<ct,n> refVelocity;
        refVelocity[0] = 0.5*(fluxVector[1] - fluxVector[0]);
        refVelocity[1] = 0.5*(fluxVector[3] - fluxVector[2]);

        // get the transposed Jacobian of the element mapping
        const FieldMatrix<ct,n,n>& jacobianInv = geometry.jacobianInverseTransposed(local);
        FieldMatrix<ct,n,n> jacobianT(jacobianInv);
        jacobianT.invert();

        // calculate the element velocity by the Piola transformation
        FieldVector<ct,n> elementVelocity(0);
        jacobianT.umtv(refVelocity, elementVelocity);
        elementVelocity /= geometry.integrationElement(local);

        for (int i=0; i<n; i++)
        {
            meanVelocity[i]+=elementVelocity[i];
            //                    std::cout<<"velocity "<<model.variables().velocity[i][j]<<std::endl;

            meanVf[i] += (elementVelocity[i] * model.diffproblem.materialLaw.fractionalW(sat, global, *hostIt, local));
        }
    }
    double dist;
    switch (subProblemNumber)
    {
        case 1:
        dist=upperRight[0] - lowerLeft[0];
        break;
        case 2:
        dist=upperRight[1] - lowerLeft[1];
        break;
    }
    meanVelocity /= size;
    meanVf /= size;
    meanVelocity *= fwSmean;
    //        std::cout<<"meanveltimesfwmean ="<<meanvelocity<<std::endl;
    FieldVector<double, n> D = meanVelocity - meanVf;

    //divide by pressure gradient = 1/distance
    D *= dist;

    //calculate saturation gradient
    count /= 2;
    satIn /= count;
    satOut /= count;
    double satGrad = (satIn - satOut)/dist;
    //        std::cout<<"satin = "<<satin<<"satout = "<<satout<<"satgrad ="<<satgrad<<std::endl;

    //divide D by saturationgradient
    D /= satGrad;
    //                    std::cout<<"D ="<<D<<std::endl;

    //write D and satmean into soil
    switch (subProblemNumber)
    {
        case 1:
        model.diffproblem.soil.getDispersion()[coarseIndex][timeStepCounter_][0][0]= D[0];
        model.diffproblem.soil.getDispersion()[coarseIndex][timeStepCounter_][1][1]= D[1];
        model.diffproblem.soil.getDispersionSat()[coarseIndex][timeStepCounter_]= meanSat_;
        break;
        case 2:
        model.diffproblem.soil.getDispersion()[coarseIndex][timeStepCounter_][0][0]+= D[0];
        model.diffproblem.soil.getDispersion()[coarseIndex][timeStepCounter_][1][1]+= D[1];
        model.diffproblem.soil.getDispersionSat()[coarseIndex][timeStepCounter_]+= meanSat_;
        model.diffproblem.soil.getDispersionSat()[coarseIndex][timeStepCounter_]*=0.5;
        break;
        default:
        break;
    }
    timeStepCounter_++;
}

template<class G, class Model>void TimeLoopSubProbs<G,Model>::calculateconvectivecorrection(Model& model, const int subProblemNumber, int& coarseIndex, int& numberInSelfCoarse, FieldVector<ct, n>& faceGlobalCoarse,FieldVector<ct,n>& lowerLeft, FieldVector<ct,n>& upperRight)
{
    //std::cout<<"numberinself = "<<numberInSelfCoarse<<"coarseIndex = "<<coarseIndex<<std::endl;

    typedef typename Entity::Geometry Geometry;
    typedef typename G::LevelGridView GV;
    typedef typename GV::template Codim<0>::Iterator Iterator;
    typedef typename GV::IndexSet IS;
    typedef typename IntersectionIteratorGetter<G,LevelTag>::IntersectionIterator
    IntersectionIterator;
    typedef typename G::HostGridType HG;
    typedef typename HG::LevelGridView HGV;

    typedef typename G::template Codim<0>::EntityPointer EntityPointer;
    typedef typename HG::template Codim<0>::EntityPointer HostEntityPointer;
    typedef typename HG::Traits::template Codim<0>::Entity HostEntity;
    typedef typename HGV::IndexSet HIS;
    int fineLev = model.variables().diffLevel;
    int size = model.variables().diffSize;

    const GV& gridViewCoarse(model.variables().grid.levelView(0));
    const GV& gridViewFine(model.variables().grid.levelView(fineLev));
    const HGV& gridViewHost(gridViewCoarse.grid().getHostGrid().levelView(0));

    //    const IS& coarseIset = gridViewCoarse.indexSet();
    const HIS& coarseIsetHost = gridViewHost.indexSet();

    meanSat_=0;

    double meanVelocity=0;

    //    std::cout<<"meanSat = "<<meanSat_<<"meanvel = "<<meanVelocity<<std::endl;

    double meanVf=0;//mean of finescale velocity v times fine scale fractional flow function f
    double fwSmean;// fractional flow function values of mean saturation

    Dune::BlockVector<Dune::FieldMatrix<double,2,2> > D(2);
    Dune::FieldVector<double,2> meanSatHelp(0);

    double pressureIn=0;
    double pressureOut=0;
    double satIn=0;
    double satOut=0;
    int count = 0;
    int cellCounter = 0;
    int coarseCellNumber =0;

    Iterator eEndIt = gridViewFine.template end<0>();
    Iterator eBeginIt = gridViewFine.template begin<0>();
    for (Iterator it = eBeginIt; it != eEndIt; ++it)
    {
        // element geometry
        const Geometry& geometry = it->geometry();

        GeometryType gt = geometry.type();
        // cell center in reference element
        const FieldVector<ct,n>& local = ReferenceElements<ct,n>::general(gt).position(0, 0);

        // get global coordinate of cell center
        const FieldVector<ct,n> global = geometry.global(local);

        const HostEntityPointer& hostIt = model.variables().grid.template getHostEntity<0>(*it);

        int index = model.variables().diffMapper.map(*it);
        double sat = (*model)[index];
        double pressure = model.variables().pressure[index];

        HostEntity& hostEntity = *hostIt;
        HostEntityPointer fatherPointer = hostEntity.father();
        int fatherLevel = fatherPointer.level();
//        std::cout<<"fatherLevel = "<<fatherLevel<<std::endl;
        while (fatherLevel != 0)
        {
            HostEntity& fatherEntity = *fatherPointer;
            fatherPointer = fatherEntity.father();
            fatherLevel = fatherPointer.level();
//            std::cout<<"fatherLevel = "<<fatherLevel<<std::endl;
        }

        int coarseInd = coarseIsetHost.index(*fatherPointer);
//                std::cout<<"coarseInd = "<<coarseInd<<"coarseIndex = "<<coarseIndex<<std::endl;
        if (coarseInd == coarseIndex)
        {
            satOut+=sat;
            pressureOut+=pressure;
            meanSatHelp[0] += sat;
            //                std::cout<<"pressureOut = "<<pressureOut<<std::endl;
        }
        else
        {
            satIn+=sat;
            pressureIn+=pressure;
            meanSatHelp[1] += sat;
            meanSat_ += sat;
            //                std::cout<<"meanSat = "<<meanSat_<<std::endl;
            //                std::cout<<"pressureIn = "<<pressureIn<<std::endl;
            cellCounter++;
        }

        double meanSat = meanSat_/cellCounter;
        fwSmean = model.diffproblem.materialLaw.fractionalW(meanSat, global, *hostIt, local);

        IntersectionIterator
        endit = IntersectionIteratorGetter<G, LevelTag>::end(*it);
        for (IntersectionIterator
                is = IntersectionIteratorGetter<G, LevelTag>::begin(*it); is
                !=endit; ++is)
        {
            // get geometry type of face
            Dune::GeometryType gtf = is->intersectionSelfLocal().type();

            // center in face's reference element
            const Dune::FieldVector<ct,n-1>&
            faceLocal = Dune::ReferenceElements<ct,n-1>::general(gtf).position(0,0);

            // center of face in global coordinates
            Dune::FieldVector<ct,dimworld> faceGlobal = is->intersectionGlobal().global(faceLocal);

            double faceVol = is->intersectionGlobal().volume();

            // get normal vector
            FieldVector<double,n> unitOuterNormal = is->unitOuterNormal(faceLocal);

            int faceNumber = is->numberInSelf();

            if (coarseInd != coarseIndex)
            {
                switch(subProblemNumber)
                {
                    case 1:
                    if (faceGlobal[0] == faceGlobalCoarse[0])
                    {
                        meanVelocity += model.variables().velocity[index][faceNumber][0];
                        meanVf += (model.variables().velocity[index][faceNumber][0] * model.diffproblem.materialLaw.fractionalW(sat, global, *hostIt, local));
                        count++;
                    }
                    break;
                    case 2:
                    if (faceGlobal[0] == faceGlobalCoarse[0])
                    {
                        meanVelocity -= model.variables().velocity[index][faceNumber][0];
                        meanVf -= (model.variables().velocity[index][faceNumber][0] * model.diffproblem.materialLaw.fractionalW(sat, global, *hostIt, local));
                        count++;
                    }
                    break;
                    case 3:
                    if (faceGlobal[1] == faceGlobalCoarse[1])
                    {
                        meanVelocity += model.variables().velocity[index][faceNumber][1];
                        meanVf += (model.variables().velocity[index][faceNumber][1] * model.diffproblem.materialLaw.fractionalW(sat, global, *hostIt, local));
                        count++;
                    }
                    break;
                    case 4:
                    if (faceGlobal[1] == faceGlobalCoarse[1])
                    {
                        meanVelocity -= model.variables().velocity[index][faceNumber][1];
                        meanVf -= (model.variables().velocity[index][faceNumber][1] * model.diffproblem.materialLaw.fractionalW(sat, global, *hostIt, local));
                        count++;
                    }
                    break;
                }
            }
        }
    }
    meanSatHelp/=cellCounter;
    Iterator endCIt = gridViewCoarse.template end<0>();
    for (Iterator cIt = gridViewCoarse.template begin<0>(); cIt != endCIt; ++cIt)
    {
        const HostEntityPointer& hostCoarseIt = model.variables().grid.template getHostEntity<0>(*cIt);

        int coarseInd = coarseIsetHost.index(*hostCoarseIt);

        int timesize = model.diffproblem.soil.getDispersionSat().M();
        //subtract DgradS
        for (int i=0;i<timesize;i++)
        {
            if (model.diffproblem.soil.getDispersionSat()[coarseIndex][i]> meanSatHelp[coarseCellNumber])
            {
                double satdiff1 = model.diffproblem.soil.getDispersionSat()[coarseIndex][i] - meanSatHelp[coarseCellNumber];
                double satdiff2 = 1e100;
                if (i)
                {
                    satdiff2 = meanSatHelp[coarseCellNumber] - model.diffproblem.soil.getDispersionSat()[coarseIndex][i-1];
                }
                if (satdiff1 < satdiff2)
                {
                    D[coarseCellNumber]=model.diffproblem.soil.getDispersion()[coarseInd][i];
                    //                        std::cout<<"result1 = "<<problem.soil.getDispersion()[indexI][i]<<std::endl;
                }
                else
                {
                    D[coarseCellNumber]=model.diffproblem.soil.getDispersion()[coarseInd][i-1];
                    //                        std::cout<<"result2 = "<<problem.soil.getDispersion()[indexI][i-1]<<std::endl;
                    //                        std::cout<<"satGrad = "<<satGradient<<std::endl;
                }
                if (i==(timesize-1))
                {
                    D[coarseCellNumber]=model.diffproblem.soil.getDispersion()[coarseInd][i];
                    break;
                }
                break;
            }
        }

        coarseCellNumber++;
    }
    double dist;
    switch (subProblemNumber)
    {
        case 1:
        dist= upperRight[0] - faceGlobalCoarse[0];
        break;
        case 2:
        dist=faceGlobalCoarse[0] - lowerLeft[0];
        break;
        case 3:
        dist=upperRight[1] - faceGlobalCoarse[1];
        break;
        case 4:
        dist=faceGlobalCoarse[1] - lowerLeft[1];
        break;
    }
//    std::cout<<"cellCounter = "<<cellCounter<<"count ="<<count<<std::endl;
    //    std::cout<<"meanVelocity = "<<meanVelocity<<std::endl;
    meanVelocity /= count;
    meanSat_ /= cellCounter;
    meanVf /= count;
    //    std::cout<<"fwSmean ="<<fwSmean<<std::endl;
    meanVelocity *= fwSmean;
    //    std::cout<<"meanveltimesfwmean ="<<meanVelocity<<std::endl;
    double M = meanVelocity - meanVf;
//        std::cout<<"m = "<<M<<std::endl;

    //calculate saturation gradient
    satIn /= cellCounter;
    satOut /= cellCounter;
    pressureIn /= cellCounter;
    pressureOut /= cellCounter;

    double satGrad = (satIn - satOut)/dist;
    double pressureGrad = (pressureIn - pressureOut)/dist;
    //    std::cout<<"pressureIn = "<<pressureIn<<"pressureOut = "<<pressureOut<<std::endl;
//        std::cout<<"satgrad ="<<satGrad<<"pressureGrad = "<<pressureGrad<<std::endl;

    //scale by pressure gradient = 1/dist!!!
    if (pressureGrad != 0)
    {
        M /= pressureGrad;
    }
//    std::cout<<"D000 = "<<D[0][0][0]<<"D100 = "<<D[1][0][0]<<"D011 = "<<D[0][1][1]<<"D111 = "<<D[1][1][1]<<std::endl;

    //write M and satmean into soil
    switch (subProblemNumber)
    {
        case 1:
        model.diffproblem.soil.getM()[coarseIndex][timeStepCounter_][numberInSelfCoarse]= M-0.5*(D[0][0][0]+D[1][0][0])*satGrad;
        model.diffproblem.soil.getMSat()[coarseIndex][timeStepCounter_][numberInSelfCoarse]= meanSat_;
        break;
        case 2:
        model.diffproblem.soil.getM()[coarseIndex][timeStepCounter_][numberInSelfCoarse]= M-0.5*(D[0][0][0]+D[1][0][0])*satGrad;
        model.diffproblem.soil.getMSat()[coarseIndex][timeStepCounter_][numberInSelfCoarse]= meanSat_;
        break;
        case 3:
        model.diffproblem.soil.getM()[coarseIndex][timeStepCounter_][numberInSelfCoarse]= M-0.5*(D[0][1][1]+D[1][1][1])*satGrad;
        model.diffproblem.soil.getMSat()[coarseIndex][timeStepCounter_][numberInSelfCoarse]= meanSat_;
        break;
        case 4:
        model.diffproblem.soil.getM()[coarseIndex][timeStepCounter_][numberInSelfCoarse]= M-0.5*(D[0][1][1]+D[1][1][1])*satGrad;
        model.diffproblem.soil.getMSat()[coarseIndex][timeStepCounter_][numberInSelfCoarse]= meanSat_;
        break;
        default:
        break;
    }
    timeStepCounter_++;

}

}
#endif
