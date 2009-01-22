// $Id$
#ifndef DUNE_TIMELOOPSUBPROBS_HH
#define DUNE_TIMELOOPSUBPROBS_HH

#include "dumux/timedisc/timestep.hh"
#include "dumux/timedisc/rungekuttastep.hh"

namespace Dune
{
template<class Grid, class Model>
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

    void calculateBoundaryFluxCorrection(Model&, const int,const int&,const int&, const GlobalPosition&, const GlobalPosition&, const GlobalPosition&);

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
                if (!isBoundary_)
                {
                    calculateFluxCorrection(model, subProblemNumber, globalIdxCoarseCurrent, faceNumberCurrent, globalPosFaceCoarseCurrent, lowerLeft, upperRight);
                }
                else
                {
                    calculateBoundaryFluxCorrection(model, subProblemNumber, globalIdxCoarseCurrent, faceNumberCurrent, globalPosFaceCoarseCurrent, lowerLeft, upperRight);
                }
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
//                VTKWriter<typename Grid::LeafGridView> vtkwriter(model.variables().grid.leafView());
//                char fname[128];
//                sprintf(fname, "dispsubprob-%02d-%02d-%03d",globalIdxCoarseCurrent, subProblemNumber, k);
//                vtkwriter.addCellData((*model), "saturation");
//                //                                              vtkwriter.addCellData(model.variables().pressure, "total pressure p~");
//                vtkwriter.write(fname, VTKOptions::ascii);
//            }
//            else
//            {
//                VTKWriter<typename Grid::LeafGridView> vtkwriter(model.variables().grid.leafView());
//                char fname[128];
//                sprintf(fname, "convsubprob-%02d-%02d-%03d",globalIdxCoarseCurrent, faceNumberCurrent, k);
//                vtkwriter.addCellData((*model), "saturation");
//                //												vtkwriter.addCellData(model.variables().pressure, "total pressure p~");
//                vtkwriter.write(fname, VTKOptions::ascii);
//            }
        }
        //			std::cout<<"saturation ="<<model.variables().saturation<<"pressure = "<<model.variables().pressure<<std::endl;
        std::cout << ",subproblem timestep: " << k << "\t t=" << t
        << "\t dt_=" << dt_ << std::endl;

        return;
    }

    TimeLoopSubProbs(bool isBoundary = false,const Scalar mMMeanSat = 0.99,
            const Scalar cfl = 0.99, TimeStep<Grid,
            Model>& tist = *(new RungeKuttaStep<Grid, Model> (1))) :
    isBoundary_(isBoundary),tEnd_(1e5), minMaxMeanSat_(mMMeanSat), cFLFactor_(cfl), timeStep_(tist), dt_(1e100), tStart_(0), maxDt_(1e100), timeStepCounter_(0), meanSat_(0)
    {
    }

private:
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
};

template<class Grid, class Model>void TimeLoopSubProbs<Grid,Model>::calculateMacroDispersion(Model& model, const int subProblemNumber, const int& globalIdxCoarseCurrent, const GlobalPosition& lowerLeft, const GlobalPosition& upperRight)
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
    int fineLev = model.variables().diffLevel;

    const GridView& gridViewCoarse(model.variables().grid.levelView(0));
    const GridView& gridViewFine(model.variables().grid.levelView(fineLev));

    meanSat_=0;
    Scalar fwSMean=0;// fractional flow function values of mean saturation

    FieldVector<Scalar, dim> meanVelocity(0);
    FieldVector<Scalar, dim> meanVf(0);//mean of finescale velocity v times fine scale fractional flow function f

    Scalar elementVolumeCoarse = 1;
    Scalar faceAreaCoarse=0;
    Scalar sumFaceArea = 0;

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

        const HostElementPointer& eItHost = model.variables().grid.template getHostEntity<0>(*eIt);

        int globalIdx = model.variables().diffMapper.map(*eIt);

        Scalar elementVolume = geometry.integrationElement(localPos)*Dune::ReferenceElements<Scalar,dim>::general(gt).volume();
        Scalar elementFaceArea = 0;

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
        fwSMean += model.diffProblem.materialLaw.fractionalW(sat, globalPos, *eItHost, localPos)*elementVolume;

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
            sumFaceArea += 0.5*faceArea;
            elementFaceArea += faceArea;

            int faceNumber = isIt->numberInSelf();

            // get normal vector
            FieldVector<Scalar,dim> unitOuterNormal = isIt->unitOuterNormal(faceLocal);

            fluxVector[faceNumber] = model.variables().velocity[globalIdx][faceNumber] * unitOuterNormal;
            fluxVector[faceNumber] *= faceArea;

            if (isIt->boundary())
            {
                sumFaceArea += 0.5*faceArea;

                // center of face inside volume reference element
                const LocalPosition&
                localPosFace = Dune::ReferenceElements<Scalar,dim>::general(faceGT).position(faceNumber,1);

                //get boundary type
                BoundaryConditions::Flags bctype = model.diffProblem.bctypeSat(globalPosFace, *eIt, localPosFace);

                if (bctype == BoundaryConditions::dirichlet)
                {
                    Scalar eps = 1e-6;
                    switch(subProblemNumber)
                    {
                        case 1:
                        if (globalPosFace[0]> lowerLeft[0]+eps)
                        {
                            satOut+=sat*faceArea;
                            faceAreaCoarse += faceArea;
                        }
                        else
                        {
                            satIn+=sat*faceArea;
                            faceAreaCoarse += faceArea;
                        }
                        break;
                        case 3:
                        if (globalPosFace[1]> lowerLeft[1]-eps)
                        {
                            satOut+=sat*faceArea;
                            faceAreaCoarse += faceArea;
                        }
                        else
                        {
                            satIn+=sat*faceArea;
                            faceAreaCoarse += faceArea;
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
            meanVelocity[i]+=elementVelocity[i]*elementFaceArea;
            meanVf[i] += (elementVelocity[i] * model.diffProblem.materialLaw.fractionalW(sat, globalPos, *eItHost, localPos))*elementFaceArea;
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
    meanVelocity /= sumFaceArea;
    meanVf /= sumFaceArea;
    fwSMean /= elementVolumeCoarse;
    meanSat_ /= elementVolumeCoarse;

    FieldVector<Scalar, dim> meanVelocityMeanFwS = meanVelocity;
    meanVelocityMeanFwS *= fwSMean;

    //        std::cout<<"meanveltimesfwmean ="<<meanvelocity<<std::endl;
    FieldVector<Scalar, dim> D = meanVf-meanVelocityMeanFwS;

    //divide by pressure gradient = 1/distance
    D *= dist;

    //calculate saturation gradient
    satIn /=faceAreaCoarse;
    satOut /= faceAreaCoarse;
    Scalar satGrad = (satIn - satOut)/dist;
    //        std::cout<<"satin = "<<satin<<"satout = "<<satout<<"satgrad ="<<satgrad<<std::endl;

    //divide D by saturationgradient
    D /= satGrad;
    //                    std::cout<<"D ="<<D<<std::endl;


    //write Dispersion coefficient and mean saturation into soil
    int vectorsize = model.diffProblem.soil.getDispersion()[globalIdxCoarseCurrent].size();
    if (vectorsize < (timeStepCounter_+1))
    {
        Dune::FieldVector<Scalar, dim> addVec(0);
        Dune::FieldVector<Scalar, 1> addVecSat(0);
        model.diffProblem.soil.getDispersion()[globalIdxCoarseCurrent].push_back(addVec);
        model.diffProblem.soil.getDispersionSat()[globalIdxCoarseCurrent].push_back(addVecSat);
    }

    switch (subProblemNumber)
    {
        case 1:
        model.diffProblem.soil.getDispersion()[globalIdxCoarseCurrent][timeStepCounter_][0]= D[0];
        model.diffProblem.soil.getDispersion()[globalIdxCoarseCurrent][timeStepCounter_][1]= D[1];
        model.diffProblem.soil.getDispersionSat()[globalIdxCoarseCurrent][timeStepCounter_]= meanSat_;
        break;
        case 2:
        model.diffProblem.soil.getDispersion()[globalIdxCoarseCurrent][timeStepCounter_][0]+= D[0];
        model.diffProblem.soil.getDispersion()[globalIdxCoarseCurrent][timeStepCounter_][1]+= D[1];
        model.diffProblem.soil.getDispersionSat()[globalIdxCoarseCurrent][timeStepCounter_]+= meanSat_;
        model.diffProblem.soil.getDispersionSat()[globalIdxCoarseCurrent][timeStepCounter_]*=0.5;
        break;
        default:
        break;
    }
    timeStepCounter_++;
}

template<class Grid, class Model>void TimeLoopSubProbs<Grid,Model>::calculateFluxCorrection(Model& model, const int subProblemNumber,const int& globalIdxCoarseCurrent,const int& faceNumberCurrent, const GlobalPosition& globalPosFaceCoarseCurrent,const GlobalPosition& lowerLeft, const GlobalPosition& upperRight)
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
    int fineLev = model.variables().diffLevel;

    const GridView& gridViewCoarse(model.variables().grid.levelView(0));
    const GridView& gridViewFine(model.variables().grid.levelView(fineLev));
    const HostGridView& gridViewHost(gridViewCoarse.grid().getHostGrid().levelView(0));

    //	const IndexSet& coarseIset = gridViewCoarse.indexSet();
    const HostIndexSet& indexSetCoarseHost = gridViewHost.indexSet();

    meanSat_=0;

    Scalar meanVelocity=0;
    //	std::cout<<"meanSat = "<<meanSat_<<"meanvel = "<<meanVelocity<<std::endl;

    Scalar meanVf=0;//mean of finescale velocity v times fine scale fractional flow function f
    Scalar fwSMean=0;// fractional flow function values of mean saturation

    Scalar pressureIn=0;
    Scalar pressureOut=0;

    Scalar elementVolumeCoarse = 1;
    Scalar faceAreaCoarse=0;

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

        const HostElementPointer& eItHost = model.variables().grid.template getHostEntity<0>(*eIt);

        int globalIdx = model.variables().diffMapper.map(*eIt);

        Scalar elementVolume = geometry.integrationElement(localPos)*Dune::ReferenceElements<Scalar,dim>::general(gt).volume();

        Scalar sat = (*model)[globalIdx];
        Scalar pressure = model.variables().pressure[globalIdx];

        HostElement& hostElement = *eItHost;
        HostElementPointer fatherPointer = hostElement.father();
        int fatherLevel = fatherPointer.level();
        //		std::cout<<"fatherLevel = "<<fatherLevel<<std::endl;
        while (fatherLevel != 0)
        {
            HostElement& fatherElement = *fatherPointer;
            fatherPointer = fatherElement.father();
            fatherLevel = fatherPointer.level();
            //			std::cout<<"fatherLevel = "<<fatherLevel<<std::endl;
        }

        int globalIdxCoarse = indexSetCoarseHost.index(*fatherPointer);

        const LocalPosition &localPosCoarse = Dune::ReferenceElements<Scalar,dim>::general(fatherPointer->geometry().type()).position(0, 0);

        //				std::cout<<"globalIdxCoarse = "<<globalIdxCoarse<<"globalIdxCoarseCurrent = "<<globalIdxCoarseCurrent<<std::endl;
        if (globalIdxCoarse != globalIdxCoarseCurrent)
        {
            elementVolumeCoarse = fatherPointer->geometry().integrationElement(localPosCoarse)*Dune::ReferenceElements<Scalar,dim>::general(gt).volume();
            pressureIn+=pressure*elementVolume;
            meanSat_ += sat*elementVolume;
            fwSMean += model.diffProblem.materialLaw.fractionalW(sat, globalPos, *eItHost, localPos)*elementVolume;
            //				std::cout<<"meanSat = "<<meanSat_<<std::endl;
        }
        if (globalIdxCoarse == globalIdxCoarseCurrent)
        {
            pressureOut+=pressure*elementVolume;
        }

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

            if (globalIdxCoarse != globalIdxCoarseCurrent)
            {
                switch(subProblemNumber)
                {
                    case 1:
                    if (globalPosFace[0] == globalPosFaceCoarseCurrent[0])
                    {
                        meanVelocity += model.variables().velocity[globalIdx][faceNumber][0]*faceArea;
                        meanVf += (model.variables().velocity[globalIdx][faceNumber][0] * model.diffProblem.materialLaw.fractionalW(sat, globalPos, *eItHost, localPos))*faceArea;
                        faceAreaCoarse += faceArea;
                    }
                    break;
                    case 2:
                    if (globalPosFace[0] == globalPosFaceCoarseCurrent[0])
                    {
                        meanVelocity += model.variables().velocity[globalIdx][faceNumber][0]*faceArea;
                        meanVf += (model.variables().velocity[globalIdx][faceNumber][0] * model.diffProblem.materialLaw.fractionalW(sat, globalPos, *eItHost, localPos))*faceArea;
                        faceAreaCoarse += faceArea;
                    }
                    break;
                    case 3:
                    if (globalPosFace[1] == globalPosFaceCoarseCurrent[1])
                    {
                        meanVelocity += model.variables().velocity[globalIdx][faceNumber][1]*faceArea;
                        meanVf += (model.variables().velocity[globalIdx][faceNumber][1] * model.diffProblem.materialLaw.fractionalW(sat, globalPos, *eItHost, localPos))*faceArea;
                        faceAreaCoarse += faceArea;
                    }
                    break;
                    case 4:
                    if (globalPosFace[1] == globalPosFaceCoarseCurrent[1])
                    {
                        meanVelocity += model.variables().velocity[globalIdx][faceNumber][1]*faceArea;
                        meanVf += (model.variables().velocity[globalIdx][faceNumber][1] * model.diffProblem.materialLaw.fractionalW(sat, globalPos, *eItHost, localPos))*faceArea;
                        faceAreaCoarse += faceArea;
                    }
                    break;
                }
            }
        }
    }

    Scalar dist;
    switch (subProblemNumber)
    {
        case 1:
        dist= upperRight[0] - globalPosFaceCoarseCurrent[0];
        break;
        case 2:
        dist=globalPosFaceCoarseCurrent[0] - lowerLeft[0];
        break;
        case 3:
        dist=upperRight[1] - globalPosFaceCoarseCurrent[1];
        break;
        case 4:
        dist=globalPosFaceCoarseCurrent[1] - lowerLeft[1];
        break;
    }

    meanSat_ /= elementVolumeCoarse;
    //	std::cout<<"meanSat_ = "<<meanSat_<<std::endl;

    meanVelocity /= faceAreaCoarse;
    meanVf /= faceAreaCoarse;
    //		std::cout<<"fwSMean ="<<fwSMean<<std::endl;
    //			std::cout<<"meanVelocity = "<<meanVelocity<<std::endl;

    fwSMean/=elementVolumeCoarse;
    Scalar meanVelocityMeanFwS = meanVelocity * fwSMean;
    //				std::cout<<"meanveltimesfwmean ="<<meanVelocity<<std::endl;

    Scalar fluxCorrection = meanVf - meanVelocityMeanFwS;
    //std::cout<<"m1 = "<<fluxCorrection<<std::endl;

    pressureIn /= elementVolumeCoarse;
    pressureOut /= elementVolumeCoarse;

    Scalar pressureGrad = (pressureIn - pressureOut)/dist;
    //		std::cout<<"pressureIn = "<<pressureIn<<"pressureOut = "<<pressureOut<<std::endl;
    //			std::cout<<"pressureGrad = "<<pressureGrad<<std::endl;

    //scale by pressure gradient = 1/dist!!!
    if (pressureGrad != 0)
    {
        fluxCorrection /= pressureGrad;
    }
    //		std::cout<<"m2 = "<<fluxCorrection<<std::endl;

    //write fluxCorrection and satmean into soil
    int vectorsize = model.diffProblem.soil.getM()[globalIdxCoarseCurrent].size();
    if (vectorsize < (timeStepCounter_+1))
    {
        Dune::FieldVector<Scalar, dim*2> addVec(0);
        model.diffProblem.soil.getM()[globalIdxCoarseCurrent].push_back(addVec);
        model.diffProblem.soil.getMSat()[globalIdxCoarseCurrent].push_back(addVec);
    }

    model.diffProblem.soil.getM()[globalIdxCoarseCurrent][timeStepCounter_][faceNumberCurrent]= fluxCorrection;
    model.diffProblem.soil.getMSat()[globalIdxCoarseCurrent][timeStepCounter_][faceNumberCurrent]= meanSat_;

    timeStepCounter_++;

}

template<class Grid, class Model>void TimeLoopSubProbs<Grid,Model>::calculateBoundaryFluxCorrection(Model& model, const int subProblemNumber,const int& globalIdxCoarseCurrent,const int& faceNumberCurrent, const GlobalPosition& globalPosFaceCoarseCurrent,const GlobalPosition& lowerLeft, const GlobalPosition& upperRight)
{
    //std::cout<<"numberinself = "<<faceNumberCurrent<<"globalIdxCoarseCurrent = "<<globalIdxCoarseCurrent<<std::endl;

    typedef typename Element::Geometry Geometry;
    typedef typename Grid::LevelGridView GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename IntersectionIteratorGetter<Grid,LevelTag>::IntersectionIterator
    IntersectionIterator;
    typedef typename Grid::HostGridType HostGrid;
    typedef typename HostGrid::LevelGridView HostGridView;

    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename HostGrid::template Codim<0>::EntityPointer HostElementPointer;
    typedef typename HostGrid::Traits::template Codim<0>::Entity HostElement;
    typedef typename HostGridView::IndexSet HostIndexSet;
    int fineLev = model.variables().diffLevel;

    const GridView& gridViewFine(model.variables().grid.levelView(fineLev));

    meanSat_=0;

    Scalar meanVelocity=0;

    //  std::cout<<"meanSat = "<<meanSat_<<"meanvel = "<<meanVelocity<<std::endl;

    Scalar meanVf=0;//mean of finescale velocity v times fine scale fractional flow function f
    Scalar fwSMean=0;// fractional flow function values of mean saturation

    Scalar pressure=0;

    Scalar elementVolumeCoarse = 1;
    Scalar faceAreaCoarse=0;

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

        const HostElementPointer& eItHost = model.variables().grid.template getHostEntity<0>(*eIt);

        int globalIdx = model.variables().diffMapper.map(*eIt);

        Scalar elementVolume = geometry.integrationElement(localPos)*Dune::ReferenceElements<Scalar,dim>::general(gt).volume();

        Scalar sat = (*model)[globalIdx];
        pressure += model.variables().pressure[globalIdx];

        meanSat_ += sat*elementVolume;
        //              std::cout<<"meanSat = "<<meanSat_<<std::endl;

        fwSMean = model.diffProblem.materialLaw.fractionalW(1.0, globalPos, *eItHost, localPos);

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

            switch(subProblemNumber)
            {
                case 1:
                if (globalPosFace[0] == globalPosFaceCoarseCurrent[0])
                {
                    //                      std::cout<<"globalPosFace = "<<globalPosFace<<"globalPosFaceCoarseCurrent = "<<globalPosFaceCoarseCurrent<<std::endl;
                    meanVelocity += model.variables().velocity[globalIdx][faceNumber][0]*faceArea;
                    meanVf += (model.variables().velocity[globalIdx][faceNumber][0] * model.diffProblem.materialLaw.fractionalW(1.0, globalPos, *eItHost, localPos))*faceArea;
                    faceAreaCoarse += faceArea;
                }
                break;
                case 2:
                if (globalPosFace[0] == globalPosFaceCoarseCurrent[0])
                {
                    //                      std::cout<<"globalPosFace = "<<globalPosFace<<"globalPosFaceCoarseCurrent = "<<globalPosFaceCoarseCurrent<<std::endl;
                    meanVelocity += model.variables().velocity[globalIdx][faceNumber][0]*faceArea;
                    meanVf += (model.variables().velocity[globalIdx][faceNumber][0] * model.diffProblem.materialLaw.fractionalW(1.0, globalPos, *eItHost, localPos))*faceArea;
                    faceAreaCoarse += faceArea;
                }
                break;
                case 3:
                if (globalPosFace[1] == globalPosFaceCoarseCurrent[1])
                {
                    //                      std::cout<<"globalPosFace = "<<globalPosFace<<"globalPosFaceCoarseCurrent = "<<globalPosFaceCoarseCurrent<<std::endl;
                    meanVelocity += model.variables().velocity[globalIdx][faceNumber][1]*faceArea;
                    meanVf += (model.variables().velocity[globalIdx][faceNumber][1] * model.diffProblem.materialLaw.fractionalW(1.0, globalPos, *eItHost, localPos))*faceArea;
                    faceAreaCoarse += faceArea;
                }
                break;
                case 4:
                if (globalPosFace[1] == globalPosFaceCoarseCurrent[1])
                {
                    //                      std::cout<<"globalPosFace = "<<globalPosFace<<"globalPosFaceCoarseCurrent = "<<globalPosFaceCoarseCurrent<<std::endl;
                    meanVelocity += model.variables().velocity[globalIdx][faceNumber][1]*faceArea;
                    meanVf += (model.variables().velocity[globalIdx][faceNumber][1] * model.diffProblem.materialLaw.fractionalW(1.0, globalPos, *eItHost, localPos))*faceArea;
                    faceAreaCoarse += faceArea;
                }
                break;

            }
        }
    }

    Scalar dist;
    switch (subProblemNumber)
    {
        case 1:
        dist= upperRight[0] - lowerLeft[0];
        break;
        case 2:
        dist= upperRight[0] - lowerLeft[0];
        break;
        case 3:
        dist=upperRight[1] - lowerLeft[1];
        break;
        case 4:
        dist=upperRight[1] - lowerLeft[1];
        break;
    }

    meanSat_ /= elementVolumeCoarse;
    pressure/= elementVolumeCoarse;

    meanVf /= faceAreaCoarse;
    //      std::cout<<"fwSMean ="<<fwSMean<<std::endl;

    meanVelocity /=faceAreaCoarse;

    Scalar meanVelocityMeanFwS = meanVelocity * fwSMean;
    //              std::cout<<"meanveltimesfwmean ="<<meanVelocity<<std::endl;

    Scalar fluxCorrection = meanVf - meanVelocityMeanFwS;
    //          std::cout<<"m1 = "<<fluxCorrection<<std::endl;

    Scalar pressureGrad = (1-pressure)/(0.5*dist);
    //      std::cout<<"pressureIn = "<<pressureIn<<"pressureOut = "<<pressureOut<<std::endl;
    //          std::cout<<"pressureGrad = "<<pressureGrad<<std::endl;

    //scale by pressure gradient = 1/dist!!!
    if (pressureGrad != 0)
    {
        fluxCorrection /= pressureGrad;
    }
    //      std::cout<<"m2 = "<<fluxCorrection<<std::endl;

    //write fluxCorrection and satmean into soil
    int vectorsize = model.diffProblem.soil.getM()[globalIdxCoarseCurrent].size();
    if (vectorsize < (timeStepCounter_+1))
    {
        Dune::FieldVector<Scalar, dim*2> addVec(0);
        model.diffProblem.soil.getM()[globalIdxCoarseCurrent].push_back(addVec);
        model.diffProblem.soil.getMSat()[globalIdxCoarseCurrent].push_back(addVec);
    }

    model.diffProblem.soil.getM()[globalIdxCoarseCurrent][timeStepCounter_][faceNumberCurrent]= fluxCorrection;
    model.diffProblem.soil.getMSat()[globalIdxCoarseCurrent][timeStepCounter_][faceNumberCurrent]= meanSat_;

    timeStepCounter_++;

}

}
#endif
