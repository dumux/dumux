// $Id$
#ifndef DUNE_TIMELOOPSUBPROBS_HH
#define DUNE_TIMELOOPSUBPROBS_HH

#include "dumux/timedisc/timestep.hh"
#include "dumux/timedisc/rungekuttastep.hh"

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Model, class CoarseScaleParameterType>
class TimeLoopSubProbs
{
typedef    typename Model::Scalar Scalar;
    enum
    {   dim = Grid::dimension, dimWorld = Grid::dimensionworld};

    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::HostGridType HostGrid;
    typedef typename HostGrid::template Codim<0>::EntityPointer HostElementPointer;

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:

    void calculateMacroDispersion(Model&, const int, const int&, const GlobalPosition&, const GlobalPosition&);

    void calculateFluxCorrection(Model&, const int,const int&, const int&, const GlobalPosition&, const GlobalPosition&, const GlobalPosition&);

    void execute(Model& model, const int& subProblemNumber,const int& globalIdxCoarseCurrent, const GlobalPosition& lowerLeft, const GlobalPosition& upperRight, bool correctionType, const int& faceNumberCurrent, const GlobalPosition& globalPosFaceCoarseCurrent, Scalar& tEnd = *(new Scalar(1e6)))
    {
        typedef typename Grid::LevelGridView GridView;
        typedef typename GridView::template Codim<0>::Iterator ElementIterator;
        typedef typename GridView::IndexSet IndexSet;

        typedef typename HostGrid::LevelGridView HostGridView;
        typedef typename HostGrid::Traits::template Codim<0>::Entity HostElement;
        typedef typename HostGrid::template Codim<0>::EntityPointer HostElementPointer;
        typedef typename HostGridView::IndexSet HostIndexSet;

        // initialize solution with initial values
        model.initial();

        // now do the time steps
        Scalar t = tStart_;
        int k = 0;

        bool reachedTEnd = false;
        bool stop = true;
        while (stop)
        {
            k++;

            if ((t+dt_)>tEnd_)
            {
                tEnd_+=dt_;
            }

            timeStep_.execute(model, t, dt_, maxDt_, tEnd_, cFLFactor_);

            if (correctionType)
            {
                calculateMacroDispersion(model, subProblemNumber, globalIdxCoarseCurrent, lowerLeft, upperRight);
            }
            else
            {
                calculateFluxCorrection(model, subProblemNumber, globalIdxCoarseCurrent, faceNumberCurrent, globalPosFaceCoarseCurrent, lowerLeft, upperRight);
            }

            // generate output
            t += dt_;
            t = std::min(t, tEnd_);

            if (t==tEnd_)
            {
                reachedTEnd = true;
            }

            if (reachedTEnd && meanSat_ <= minMaxMeanSat_)
            {
                tEnd_+=2*dt_;
                reachedTEnd = false;
            }
            if (meanSat_> minMaxMeanSat_)
            {
                stop = false;
            }
            //            if (correctionType)
            //            {
            //                //                VTKWriter<typename Grid::LeafGridView> vtkwriter(model.variables.grid.leafView());
            //                //                char fname[128];
            //                //                sprintf(fname, "dispsubprob-%02d-%02d-%03d",globalIdxCoarseCurrent, subProblemNumber, k);
            //                //                vtkwriter.addCellData((*model), "saturation");
            //                //                //                                              vtkwriter.addCellData(model.variables.pressure, "total pressure p~");
            //                //                vtkwriter.write(fname, VTKOptions::ascii);
            //            }
            //            else
            //            {
            //                VTKWriter<typename Grid::LeafGridView> vtkwriter(model.variables.grid.leafView());
            //                char fname[128];
            //                sprintf(fname, "convsubprob-%02d-%02d-%03d",globalIdxCoarseCurrent, faceNumberCurrent, k);
            //                vtkwriter.addCellData((*model), "saturation");
            //                //												vtkwriter.addCellData(model.variables.pressure, "total pressure p~");
            //                vtkwriter.write(fname, VTKOptions::ascii);
            //            }
        }
        //			std::cout<<"saturation ="<<model.variables.saturation<<"pressure = "<<model.variables.pressure<<std::endl;
        if (correctionType)
        {
            std::cout << ",dispersion subproblem timestep: " << k << "\t t=" << t
            << "\t dt_=" << dt_ << std::endl;
        }
        else
        {
            std::cout << ",convection subproblem timestep: " << k << "\t t=" << t
            << "\t dt_=" << dt_ << std::endl;
        }

        return;
    }

    TimeLoopSubProbs(CoarseScaleParameterType& coarseParams, const Scalar mMMeanSat = 0.99,
            const Scalar cfl = 0.99, TimeStep<Grid,
            Model>& tist = *(new RungeKuttaStep<Grid, Model> (1))) :
    coarseParameters_(coarseParams),tEnd_(1e5), minMaxMeanSat_(mMMeanSat), cFLFactor_(cfl), timeStep_(tist), dt_(1e100), tStart_(0), maxDt_(1e100), timeStepCounter_(0), meanSat_(0), eps_(1e-6)
    {
    }

private:
    CoarseScaleParameterType& coarseParameters_;
    bool isBoundary_;
    Scalar tEnd_;
    Scalar minMaxMeanSat_;
    const Scalar cFLFactor_;
    TimeStep<Grid, Model>& timeStep_;
    Scalar dt_;
    const Scalar tStart_;
    const Scalar maxDt_;
    int timeStepCounter_;
    Scalar meanSat_;
    Scalar eps_;
};

template<class Grid, class Model,class CoarseScaleParameterType>void TimeLoopSubProbs<Grid, Model, CoarseScaleParameterType>::calculateMacroDispersion(Model& model, const int subProblemNumber, const int& globalIdxCoarseCurrent, const GlobalPosition& lowerLeft, const GlobalPosition& upperRight)
{
    typedef typename Element::Geometry Geometry;
    typedef typename Grid::LevelGridView GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename Grid::HostGridType HostGrid;
    typedef typename HostGrid::LevelGridView HostGridView;

    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename HostGrid::template Codim<0>::EntityPointer HostElementPointer;
    typedef typename HostGrid::Traits::template Codim<0>::Entity HostElement;
    typedef typename HostGridView::IndexSet HostIndexSet;
    int fineLev = model.variables.diffLevel;

    const GridView& gridViewFine(model.variables.grid.levelView(fineLev));

    meanSat_=0;
    Scalar fWMean=0;// fractional flow function values of mean saturation
    FieldVector<Scalar, dim> velocityMean(0);
    FieldVector<Scalar, dim> vTimesFWMean(0);//mean of finescale velocity v times fine scale fractional flow function f

    Scalar elementVolumeCoarse = 1;
    Scalar elementVolumeFaceCoarse=0;

    Scalar satIn=0;
    Scalar satOut=0;

    ElementIterator eItEnd = gridViewFine.template end<0>();
    ElementIterator eItBegin = gridViewFine.template begin<0>();
    for (ElementIterator eIt = eItBegin; eIt != eItEnd; ++eIt)
    {
        // element geometry
        const Geometry& geometry = eIt->geometry();

        GeometryType gt = geometry.type();

        // cell center in reference element
        const LocalPosition& localPos = ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        const GlobalPosition& globalPos = geometry.global(localPos);

        const HostElementPointer& eItHost = model.variables.grid.template getHostEntity<0>(*eIt);

        int globalIdx = model.variables.diffMapper.map(*eIt);

        Scalar elementVolume = geometry.integrationElement(localPos)*Dune::ReferenceElements<Scalar,dim>::general(gt).volume();

        Scalar sat = (*model)[globalIdx];

        HostElement& hostElement = *eItHost;
        HostElementPointer fatherPointer = hostElement.father();
        int fatherLevel = fatherPointer.level();
        //      std::cout<<"fatherLevel = "<<fatherLevel<<std::endl;
        while (fatherLevel != 0)
        {
            HostElement& fatherElement = *fatherPointer;
            fatherPointer = fatherElement.father();
            fatherLevel = fatherPointer.level();
            //          std::cout<<"fatherLevel = "<<fatherLevel<<std::endl;
        }

        const LocalPosition &localPosCoarse = Dune::ReferenceElements<Scalar,dim>::general(fatherPointer->geometry().type()).position(0, 0);

        elementVolumeCoarse = fatherPointer->geometry().integrationElement(localPosCoarse)*Dune::ReferenceElements<Scalar,dim>::general(gt).volume();
        meanSat_ += sat*elementVolume;
        //        fWMean += model.diffProblem.materialLaw.fractionalW(sat, globalPos, *eItHost, localPos)*elementVolume;
        fWMean=model.diffProblem.materialLaw.fractionalW((meanSat_/elementVolumeCoarse), globalPos, *eItHost, localPos);

        FieldVector<Scalar,dim*2> fluxVector(0);

        IntersectionIterator
        isItEnd = gridViewFine.template iend(*eIt);
        for (IntersectionIterator
                isIt = gridViewFine.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {
            // get geometry type of face
            Dune::GeometryType faceGT = isIt->intersectionSelfLocal().type();

            // center in face's reference element
            const Dune::FieldVector<Scalar,dim-1>&
            faceLocal = Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

            // center of face in global coordinates
            const GlobalPosition& globalPosFace = isIt->intersectionGlobal().global(faceLocal);

            Scalar faceArea = isIt->intersectionGlobal().volume();

            int faceNumber = isIt->numberInSelf();

            // get normal vector
            FieldVector<Scalar,dim> unitOuterNormal = isIt->unitOuterNormal(faceLocal);

            fluxVector[faceNumber] = model.variables.velocity[globalIdx][faceNumber] * unitOuterNormal;
            fluxVector[faceNumber] *= faceArea;

            if (isIt->boundary())
            {
                // center of face inside volume reference element
                const LocalPosition&
                localPosFace = Dune::ReferenceElements<Scalar,dim>::general(faceGT).position(faceNumber,1);

                //get boundary type
                BoundaryConditions::Flags bctype = model.diffProblem.bctypeSat(globalPosFace, *eIt, localPosFace);

                if (bctype == BoundaryConditions::dirichlet)
                {
                    Scalar eps = 1e-6;
                    elementVolumeFaceCoarse += elementVolume;
                    switch(subProblemNumber)
                    {
                        case 1:
                        if (globalPosFace[0]> lowerLeft[0]+eps)
                        {
                            satOut+=sat*elementVolume;

                        }
                        else
                        {
                            satIn+=sat*elementVolume;
                        }
                        break;
                        case 3:
                        if (globalPosFace[1]> lowerLeft[1]-eps)
                        {
                            satOut+=sat*elementVolume;
                        }
                        else
                        {
                            satIn+=sat*elementVolume;
                        }
                        break;
                        default:
                        break;
                    }

                }

            }
        }
        // calculate velocity on reference element
        FieldVector<Scalar,dim> refVelocity;
        refVelocity[0] = 0.5*(fluxVector[1] - fluxVector[0]);
        refVelocity[1] = 0.5*(fluxVector[3] - fluxVector[2]);

        // get the transposed Jacobian of the element mapping
        const FieldMatrix<Scalar,dim,dim>& jacobianInv = geometry.jacobianInverseTransposed(localPos);
        FieldMatrix<Scalar,dim,dim> jacobianT(jacobianInv);
        jacobianT.invert();

        // calculate the element velocity by the Piola transformation
        FieldVector<Scalar,dim> elementVelocity(0);
        jacobianT.umtv(refVelocity, elementVelocity);
        elementVelocity /= geometry.integrationElement(localPos);

        for (int i=0; i<dim; i++)
        {
            velocityMean[i] += elementVelocity[i]*elementVolume;
            vTimesFWMean[i] += (elementVelocity[i] * model.diffProblem.materialLaw.fractionalW(sat, globalPos, *eItHost, localPos))*elementVolume;
        }
    }
    Scalar dist;
    switch (subProblemNumber)
    {
        case 1:
        dist=upperRight[0] - lowerLeft[0];
        break;
        case 2:
        dist=upperRight[1] - lowerLeft[1];
        break;
    }

    velocityMean /= elementVolumeCoarse;
    vTimesFWMean /= elementVolumeCoarse;
    meanSat_ /= elementVolumeCoarse;
    //            std::cout<<velocityMean<<std::endl;
    //            std::cout<<vTimesFWMean<<std::endl;
    //        std::cout<<"fWMean = "<<fWMean<<std::endl;

    FieldVector<Scalar, dim> vMeanTimesFWMean = velocityMean;
    vMeanTimesFWMean *= fWMean;

    //                std::cout<<"vMeanTimesFWMean ="<<vMeanTimesFWMean<<std::endl;
    FieldVector<Scalar, dim> D = vMeanTimesFWMean - vTimesFWMean;
    //    std::cout<<"D ="<<D<<std::endl;

    //divide by pressure gradient = 1/distance
    D *= dist;
    //    std::cout<<"D ="<<D<<std::endl;
    //calculate saturation gradient
    satIn /= (0.5*elementVolumeFaceCoarse);
    satOut /= (0.5*elementVolumeFaceCoarse);
    Scalar satGrad = (satIn - satOut)/dist;
    //        std::cout<<"satin = "<<satin<<"satout = "<<satout<<"satgrad ="<<satgrad<<std::endl;

    //divide D by saturationgradient
    D /= satGrad;
//                         std::cout<<"D ="<<D<<std::endl;

    switch (subProblemNumber)
    {
        case 1:
        coarseParameters_.getDispersion().addEntry(D, globalIdxCoarseCurrent,timeStepCounter_);
        coarseParameters_.getDispersionSat().addEntry(meanSat_, globalIdxCoarseCurrent, timeStepCounter_);
        break;
        case 3:
        coarseParameters_.getDispersion().addEntry(D, globalIdxCoarseCurrent,timeStepCounter_);
        coarseParameters_.getDispersionSat().addEntry(meanSat_, globalIdxCoarseCurrent, timeStepCounter_);
        coarseParameters_.getDispersionSat()[globalIdxCoarseCurrent][timeStepCounter_]*=0.5;
        break;
        default:
        break;
    }

    timeStepCounter_++;
}

template<class Grid, class Model,class CoarseScaleParameterType>void TimeLoopSubProbs<Grid, Model, CoarseScaleParameterType>::calculateFluxCorrection(Model& model, const int subProblemNumber,const int& globalIdxCoarseCurrent,const int& faceNumberCurrent, const GlobalPosition& globalPosFaceCoarseCurrent,const GlobalPosition& lowerLeft, const GlobalPosition& upperRight)
{
    //std::cout<<"numberinself = "<<faceNumberCurrent<<"globalIdxCoarseCurrent = "<<globalIdxCoarseCurrent<<std::endl;

    typedef typename Element::Geometry Geometry;
    typedef typename Grid::LevelGridView GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename Grid::HostGridType HostGrid;
    typedef typename HostGrid::LevelGridView HostGridView;

    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename HostGrid::template Codim<0>::EntityPointer HostElementPointer;
    typedef typename HostGrid::Traits::template Codim<0>::Entity HostElement;
    typedef typename HostGridView::IndexSet HostIndexSet;
    int fineLev = model.variables.diffLevel;

    const GridView& gridViewCoarse(model.variables.grid.levelView(0));
    const GridView& gridViewFine(model.variables.grid.levelView(fineLev));
    const HostGridView& gridViewHost(gridViewCoarse.grid().getHostGrid().levelView(0));

    const IndexSet& indexSetCoarse = gridViewCoarse.indexSet();
    const HostIndexSet& indexSetCoarseHost = gridViewHost.indexSet();

    meanSat_=0;

    FieldVector<Scalar, dim> velocityMean(0);
    Scalar fWMean=0;// fractional flow function values of mean saturation
    FieldVector<Scalar, dim> vTimesFWMean(0);//mean of finescale velocity v times fine scale fractional flow function f

    Scalar pressureIn=0;
    Scalar pressureOut=0;
    Scalar satIn = 0;
    Scalar satOut = 0;

    Scalar elementVolumeCoarse = 1;
    Scalar elementVolumeFaceCoarse=0;

    ElementIterator eItEnd = gridViewFine.template end<0>();
    ElementIterator eItBegin = gridViewFine.template begin<0>();
    for (ElementIterator eIt = eItBegin; eIt != eItEnd; ++eIt)
    {
        // element geometry
        const Geometry& geometry = eIt->geometry();

        GeometryType gt = geometry.type();

        // cell center in reference element
        const LocalPosition& localPos = ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        const GlobalPosition& globalPos = geometry.global(localPos);

        const HostElementPointer& eItHost = model.variables.grid.template getHostEntity<0>(*eIt);

        int globalIdx = model.variables.diffMapper.map(*eIt);

        Scalar elementVolume = geometry.integrationElement(localPos)*Dune::ReferenceElements<Scalar,dim>::general(gt).volume();

        Scalar sat = (*model)[globalIdx];
        Scalar pressure = model.variables.pressure[globalIdx];

        HostElement& hostElement = *eItHost;
        HostElementPointer fatherPointerHost = hostElement.father();
        Element& element = *eIt;
        ElementPointer fatherPointer = element.father();
        int fatherLevel = fatherPointerHost.level();
        //		std::cout<<"fatherLevel = "<<fatherLevel<<std::endl;
        while (fatherLevel != 0)
        {
            HostElement& fatherElementHost = *fatherPointerHost;
            Element& fatherElement = *fatherPointer;
            fatherPointerHost = fatherElementHost.father();
            fatherPointer = fatherElement.father();
            fatherLevel = fatherPointerHost.level();
            //			std::cout<<"fatherLevel = "<<fatherLevel<<std::endl;
        }

        int globalIdxCoarse = indexSetCoarse.index(*fatherPointer);
        int globalIdxCoarseHost = indexSetCoarseHost.index(*fatherPointerHost);
        //        std::cout<<"globalIdxCoarse = "<<globalIdxCoarse<<std::endl;
        //        std::cout<<"globalIdxCoarseHost = "<<globalIdxCoarseHost<<std::endl;

        const LocalPosition &localPosCoarse = Dune::ReferenceElements<Scalar,dim>::general(fatherPointerHost->geometry().type()).position(0, 0);

        //				std::cout<<"globalIdxCoarseHost = "<<globalIdxCoarseHost<<"globalIdxCoarseCurrent = "<<globalIdxCoarseCurrent<<std::endl;
        if (globalIdxCoarseHost != globalIdxCoarseCurrent)
        {
            satIn+=sat*elementVolume;
            pressureIn+=pressure*elementVolume;
        }
        if (globalIdxCoarseHost == globalIdxCoarseCurrent)
        {
            satOut+=sat*elementVolume;
            pressureOut+=pressure*elementVolume;
            elementVolumeCoarse = fatherPointerHost->geometry().integrationElement(localPosCoarse)*Dune::ReferenceElements<Scalar,dim>::general(gt).volume();
            meanSat_ += sat*elementVolume;
            fWMean=model.diffProblem.materialLaw.fractionalW((meanSat_/elementVolumeCoarse), globalPos, *eItHost, localPos);
        }

        FieldVector<Scalar, dim*2> fluxVector;

        IntersectionIterator
        isItEnd = gridViewFine.template iend(*eIt);
        for (IntersectionIterator
                isIt = gridViewFine.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {
            if (globalIdxCoarseHost == globalIdxCoarseCurrent)
            {
                // get geometry type of face
                Dune::GeometryType faceGT = isIt->intersectionSelfLocal().type();

                // center in face's reference element
                const Dune::FieldVector<Scalar,dim-1>&
                faceLocal = Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                // center of face in global coordinates
                const GlobalPosition& globalPosFace = isIt->intersectionGlobal().global(faceLocal);

                Scalar faceArea = isIt->intersectionGlobal().volume();

                int faceNumber = isIt->numberInSelf();

                // get normal vector
                FieldVector<Scalar,dim> unitOuterNormal = isIt->unitOuterNormal(faceLocal);

                fluxVector[faceNumber] = model.variables.velocity[globalIdx][faceNumber] * unitOuterNormal;
                fluxVector[faceNumber] *= faceArea;

                if (isIt->neighbor())
                {
                    // neighbor's properties
                    ElementPointer neighborPointer = isIt->outside();
                    ElementPointer neighborPointerCoarse = (*neighborPointer).father();
                    int neighborFatherLevel = neighborPointerCoarse.level();
                    while(neighborFatherLevel != 0)
                    {
                        Element& neighborFatherElement = *neighborPointerCoarse;
                        neighborPointerCoarse = neighborFatherElement.father();
                        neighborFatherLevel = neighborPointerCoarse.level();
                    }

                    int globalIdxNeighborCoarse = indexSetCoarse.index(*neighborPointerCoarse);

                    if(globalIdxCoarse != globalIdxNeighborCoarse && globalIdxCoarseHost == globalIdxCoarseCurrent)
                    {
                        satIn += sat*elementVolume;
                        pressureIn += pressure*elementVolume;
                        elementVolumeFaceCoarse += elementVolume;
                    }

                }
                if (isIt->boundary())
                {
                    // center of face inside volume reference element
                    const LocalPosition&
                    localPosFace = Dune::ReferenceElements<Scalar,dim>::general(faceGT).position(faceNumber,1);

                    //get boundary type
                    BoundaryConditions::Flags bctype = model.diffProblem.bctypeSat(globalPosFace, *eIt, localPosFace);

                    if (bctype == BoundaryConditions::dirichlet)
                    {
                        satOut += sat*elementVolume;
                        pressureOut += pressure*elementVolume;
                        elementVolumeFaceCoarse += elementVolume;
                    }
                }
            }
        }
        if (globalIdxCoarseHost == globalIdxCoarseCurrent)
        {
            // calculate velocity on reference element
            FieldVector<Scalar,dim> refVelocity;
            refVelocity[0] = 0.5*(fluxVector[1] - fluxVector[0]);
            refVelocity[1] = 0.5*(fluxVector[3] - fluxVector[2]);

            // get the transposed Jacobian of the element mapping
            const FieldMatrix<Scalar,dim,dim>& jacobianInv = geometry.jacobianInverseTransposed(localPos);
            FieldMatrix<Scalar,dim,dim> jacobianT(jacobianInv);
            jacobianT.invert();

            // calculate the element velocity by the Piola transformation
            FieldVector<Scalar,dim> elementVelocity(0);
            jacobianT.umtv(refVelocity, elementVelocity);
            elementVelocity /= geometry.integrationElement(localPos);

            for (int i=0; i<dim; i++)
            {
                velocityMean[i] += elementVelocity[i]*elementVolume;
                vTimesFWMean[i] += (elementVelocity[i] * model.diffProblem.materialLaw.fractionalW(sat, globalPos, *eItHost, localPos))*elementVolume;
            }
        }
    }

    Dune::FieldVector<Scalar,dim> DIn(0);
    Dune::FieldVector<Scalar,dim> DOut(0);

    meanSat_ /= elementVolumeCoarse;
    //    std::cout<<"meanSat = "<<meanSat_<<std::endl;
    //    meanSatHelp /= elementVolumeCoarse;

    satIn /= elementVolumeCoarse;
    satOut /= elementVolumeCoarse;

    //    satIn /= (0.5*elementVolumeFaceCoarse);
    //    satOut /= (0.5*elementVolumeFaceCoarse);

    FieldVector<Scalar, dim> satGrad(0);
    GlobalPosition distVec(0);
    switch (subProblemNumber)
    {
        case 1:
        distVec[0]= upperRight[0] - globalPosFaceCoarseCurrent[0];
        satGrad = distVec;
        satGrad *= (satIn-satOut);
        break;
        case 2:
        distVec[0]=globalPosFaceCoarseCurrent[0] - lowerLeft[0];
        satGrad = distVec;
        satGrad *= (satOut-satIn);
        break;
        case 3:
        distVec[1]=upperRight[1] - globalPosFaceCoarseCurrent[1];
        satGrad = distVec;
        satGrad *= (satIn-satOut);
        break;
        case 4:
        distVec[1]=globalPosFaceCoarseCurrent[1] - lowerLeft[1];
        satGrad = distVec;
        satGrad *= (satOut-satIn);
        break;
    }
    Scalar dist = distVec.two_norm();

    //        std::cout<<"distVec = "<<distVec<<"dist = "<<dist<<std::endl;

    satGrad /= (dist*dist);

    ElementIterator eItCoarseEnd = gridViewCoarse.template end<0>();
    for (ElementIterator eItCoarse = gridViewCoarse.template begin<0>(); eItCoarse != eItCoarseEnd; ++eItCoarse)
    {
        const HostElementPointer& eItCoarseHost = model.variables.grid.template getHostEntity<0>(*eItCoarse);

        int globalIdxCoarseHost = indexSetCoarseHost.index(*eItCoarseHost);

        if (globalIdxCoarseHost == globalIdxCoarseCurrent)
        {

            int colNum = coarseParameters_.getDispersion()[globalIdxCoarseHost].size();

            //subtract DgradS
            for (int i=0;i<colNum;i++)
            {
                if (coarseParameters_.getDispersionSat()[globalIdxCoarseHost][i]> satOut)
                {
                    double satdiff1 = coarseParameters_.getDispersionSat()[globalIdxCoarseHost][i] - satOut;
                    double satdiff2 = 1e100;
                    if (i)
                    {
                        satdiff2 = satOut - coarseParameters_.getDispersionSat()[globalIdxCoarseHost][i-1];
                    }
                    if (satdiff1 < satdiff2)
                    {
                        DOut=coarseParameters_.getDispersion()[globalIdxCoarseHost][i];
                        //                    std::cout<<"D = "<<D<<std::endl;
                    }
                    else
                    {
                        DOut=coarseParameters_.getDispersion()[globalIdxCoarseHost][i-1];
                        //                    std::cout<<"D = "<<D<<std::endl;
                        //                        std::cout<<"satGrad = "<<satGradient<<std::endl;
                    }
                    if (i==(colNum-1))
                    {
                        DOut=coarseParameters_.getDispersion()[globalIdxCoarseHost][i];
                        break;
                    }
                    break;
                }
            }
        }
        else
        {
            int colNum = coarseParameters_.getDispersion()[globalIdxCoarseHost].size();

            //subtract DgradS
            for (int i=0;i<colNum;i++)
            {
                if (coarseParameters_.getDispersionSat()[globalIdxCoarseHost][i]> satIn)
                {
                    double satdiff1 = coarseParameters_.getDispersionSat()[globalIdxCoarseHost][i] - satIn;
                    double satdiff2 = 1e100;
                    if (i)
                    {
                        satdiff2 = satIn - coarseParameters_.getDispersionSat()[globalIdxCoarseHost][i-1];
                    }
                    if (satdiff1 < satdiff2)
                    {
                        DIn=coarseParameters_.getDispersion()[globalIdxCoarseHost][i];
                        //                    std::cout<<"D = "<<D<<std::endl;
                    }
                    else
                    {
                        DIn=coarseParameters_.getDispersion()[globalIdxCoarseHost][i-1];
                        //                    std::cout<<"D = "<<D<<std::endl;
                        //                        std::cout<<"satGrad = "<<satGradient<<std::endl;
                    }
                    if (i==(colNum-1))
                    {
                        DIn=coarseParameters_.getDispersion()[globalIdxCoarseHost][i];
                        break;
                    }
                    break;
                }
            }
        }
    }

    velocityMean /= elementVolumeCoarse;
    vTimesFWMean /= elementVolumeCoarse;
    //        std::cout<<"vTimesFWMean ="<<vTimesFWMean<<std::endl;
    //        std::cout<<"velocityMean = "<<velocityMean<<std::endl;

    //    fWMean/=elementVolumeCoarse;
    FieldVector<Scalar, dim> vMeanTimesFWMean = velocityMean;
    vMeanTimesFWMean *= fWMean;
    //        std::cout<<"meanveltimesfwmean ="<<vMeanTimesFWMean<<std::endl;

    FieldVector<Scalar, dim> fluxCorrection = vMeanTimesFWMean - vTimesFWMean;
    //                std::cout<<"m1 = "<<fluxCorrection<<std::endl;

    pressureIn /= elementVolumeCoarse;
    pressureOut /= elementVolumeCoarse;

    //    pressureIn /= (0.5*elementVolumeFaceCoarse);
    //    pressureOut /= (0.5*elementVolumeFaceCoarse);

    Scalar pressureGrad = (pressureIn - pressureOut)/dist;
    //    			std::cout<<"pressureGrad = "<<pressureGrad<<std::endl;

    //scale by pressure gradient = 1/dist!!!
    if (pressureGrad != 0)
    {
        fluxCorrection /= pressureGrad;
    }
    //        		std::cout<<"m2 = "<<fluxCorrection<<std::endl;


    //    FieldVector<Scalar, dim> dispersiveFlux(0);

    Scalar dispersiveFlux = DIn*satGrad;
    dispersiveFlux += DOut*satGrad;
    dispersiveFlux *= 0.5;

    //    std::cout<<"satGrad = "<<satGrad<<std::endl;
    //    std::cout<<"dispersiveFlux = "<<D<<std::endl;
    //    for (int i = 0; i<dim; i++)
    //    {
    //        dispersiveFlux[i] = D[i]*satGrad[i];
    //    }

    //    fluxCorrection -= dispersiveFlux;
    //    std::cout<<"dispersiveFlux = "<<dispersiveFlux<<std::endl;

    FieldVector<Scalar, dim> normal(distVec);
    normal /= dist;
    //    std::cout<<"normal = "<<normal<<std::endl;

    //                std::cout<<"fluxCorrection = "<<fluxCorrection<<std::endl;

    FieldVector<Scalar, 2*dim> newEntry(0);
    newEntry[faceNumberCurrent] = (fluxCorrection*normal-dispersiveFlux);
    coarseParameters_.getFluxCorr().addEntry(newEntry, globalIdxCoarseCurrent, timeStepCounter_);

    newEntry = 0;
    newEntry[faceNumberCurrent] = meanSat_;
    coarseParameters_.getFluxCorrSat().addEntry(newEntry, globalIdxCoarseCurrent, timeStepCounter_);

    timeStepCounter_++;

}

}
#endif
