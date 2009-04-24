// $Id$
#ifndef DUNE_TIMELOOPSUBPROBS_HH
#define DUNE_TIMELOOPSUBPROBS_HH

#include "dumux/timedisc/timestep.hh"
#include "dumux/timedisc/rungekuttastep.hh"

/**
 * @file
 * @brief  Timeloop class for solving a local fine-scale problem within a local-global upscaling approach.
 * @author Markus Wolff
 */

namespace Dune
{
//! \ingroup MultiMulti
/*!  Timeloop class for solving a local fine-scale problem within a local-global upscaling approach.
 * Local fine-scale problems are solved in order to obtain coarse-scale quantities.
 * Local problems consisting of one coarse cell and
 * local problems consisting of two coarse cells are solved.
 * The coarse scale quantities are calculated according to the following paper:
 *
 * Efendiev and Durlofsky, "A Generalized Convection-Diffusion Model for Subgrid Transport in Porous Media", SIAM Multiscale Modeling and Simulation, Vol. 3, 504 - 526, 2003.*/
/*! Template parameters are:
 *  - Grid                       a Dune::subgrid grid type
 *  - Model                      the model type (here of type IMPESSubProbs)
 *  - CoarseScaleParameterType   a class type defining the course scale quantities
 */

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

    //!calculates a saturation dependent macro dispersion coefficient from the solution of the local fine-scale problems
    void calculateMacroDispersion(Model&, const int, const int&, const GlobalPosition&, const GlobalPosition&);

    //!calculates a saturation dependent convective flux correction from the solution of the local fine-scale problems
    void calculateFluxCorrection(Model&, const int,const int&, const int&, const GlobalPosition&, const GlobalPosition&, const GlobalPosition&);

    //!start calculations for a subproblem
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

        //bool variables for timeloop control
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

            if (correctionType)// correctionType = true -> subproblem for macro dispersion
            {
                calculateMacroDispersion(model, subProblemNumber, globalIdxCoarseCurrent, lowerLeft, upperRight);
            }
            else// correctionType = false -> subproblem for convective flux correction
            {
                calculateFluxCorrection(model, subProblemNumber, globalIdxCoarseCurrent, faceNumberCurrent, globalPosFaceCoarseCurrent, lowerLeft, upperRight);
            }

            // determine time t
            t += dt_;
            t = std::min(t, tEnd_);

            if (t==tEnd_)
            {
                reachedTEnd = true;
            }

            //check whether the subdomain is sufficiently saturated with the infiltrating fluid, if tEnd is reached.
            //if that is not the case, enlarge tEnd
            if (reachedTEnd && meanSat_ <= minMaxMeanSat_)
            {
                tEnd_+=2*dt_;
                reachedTEnd = false;
            }
            //if that is the case, stop the time loop
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

        if (correctionType)// correctionType = true -> subproblem for macro dispersion
        {
            std::cout << ",dispersion subproblem timestep: " << k << "\t t=" << t
            << "\t dt_=" << dt_ << std::endl;
        }
        else// correctionType = false -> subproblem for convective flux correction
        {
            std::cout << ",convection subproblem timestep: " << k << "\t t=" << t
            << "\t dt_=" << dt_ << std::endl;
        }

        return;
    }

    /**
     * \param coarseParams object of type CoarseScaleParameterType containing the coarse-scale parameters to be computed
     * \param mMMeanSat minimum averaged saturation which has to be reached within a sub domain
     * \param cfl cfl-security-factor
     * TimeStep object defining the time-discretisation scheme
     */
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

//!calculates a saturation dependent macro dispersion coefficient from the solution of the local fine-scale problems
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

    //get the gridView
    const GridView& gridViewFine(model.variables.grid.levelView(fineLev));

    //some variable declarations
    meanSat_=0;
    Scalar fWMean=0;// fractional flow function values dependent on mean saturation
    FieldVector<Scalar, dim> velocityMean(0); //fine scale velocity averaged over current coarse cell
    FieldVector<Scalar, dim> vTimesFWMean(0);//mean of fine-scale velocity v times fine scale fractional flow function f

    Scalar elementVolumeCoarse = 1;
    Scalar elementVolumeFaceCoarse=0;//volume of the fine cells at a coarse face

    Scalar satIn=0;//saturation averaged over the coarse inlet face
    Scalar satOut=0;//saturation averaged over the coarse outlet face

    //start iteration over fine cells
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

        //get host element (global)
        const HostElementPointer& eItHost = model.variables.grid.template getHostEntity<0>(*eIt);


        int globalIdx = model.variables.diffMapper.map(*eIt);

        Scalar elementVolume = geometry.integrationElement(localPos)*Dune::ReferenceElements<Scalar,dim>::general(gt).volume();

        //get saturation
        Scalar sat = (*model)[globalIdx];

        //get the coarse scale entity the fine scale entity belongs to
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

        //weight saturation with fine cell volume
        meanSat_ += sat*elementVolume;

        //get fWMean as function of the average Saturation
        fWMean=model.diffProblem.materialLaw.fractionalW((meanSat_/elementVolumeCoarse), globalPos, *eItHost, localPos);

        //fluxVector for calculation of the element velocity
        FieldVector<Scalar,dim*2> fluxVector(0);

        //start iteration over faces of the current fine-scale element
        IntersectionIterator
        isItEnd = gridViewFine.template iend(*eIt);
        for (IntersectionIterator
                isIt = gridViewFine.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {
            // get geometry type of face
            Dune::GeometryType faceGT = isIt->geometryInInside().type();

            // center in face's reference element
            const Dune::FieldVector<Scalar,dim-1>&
            faceLocal = Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

            // center of face in global coordinates
            const GlobalPosition& globalPosFace = isIt->geometry().global(faceLocal);

            Scalar faceArea = isIt->geometry().volume();

            int faceNumber = isIt->indexInInside();

            // get normal vector
            FieldVector<Scalar,dim> unitOuterNormal = isIt->unitOuterNormal(faceLocal);

            //store flux for calculation of the element velocity
            fluxVector[faceNumber] = model.variables.velocity[globalIdx][faceNumber] * unitOuterNormal;
            fluxVector[faceNumber] *= faceArea;

            if (isIt->boundary())
            {
                // center of face inside volume reference element
                const LocalPosition&
                localPosFace = Dune::ReferenceElements<Scalar,dim>::general(faceGT).position(faceNumber,1);

                //get boundary type
                BoundaryConditions::Flags bctype = model.diffProblem.bctypeSat(globalPosFace, *eIt, localPosFace);

                //get saturation at inlet and outlet -> inlet and outlet correspond to dirichlet boundaries
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
        }//end of intersection iterator
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

        //weight by element volume
        for (int i=0; i<dim; i++)
        {
            velocityMean[i] += elementVelocity[i]*elementVolume;
            vTimesFWMean[i] += (elementVelocity[i] * model.diffProblem.materialLaw.fractionalW(sat, globalPos, *eItHost, localPos))*elementVolume;
        }
    }//end of element iterator

    //determine the distance between inlet and outlet
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

    //fine scale velocities were weighted with the element volume!
    velocityMean /= elementVolumeCoarse;
    vTimesFWMean /= elementVolumeCoarse;
    meanSat_ /= elementVolumeCoarse;
    //            std::cout<<velocityMean<<std::endl;
    //            std::cout<<vTimesFWMean<<std::endl;
    //        std::cout<<"fWMean = "<<fWMean<<std::endl;

    //calculate the term consisting of the product of the averaged quantities
    FieldVector<Scalar, dim> vMeanTimesFWMean = velocityMean;
    vMeanTimesFWMean *= fWMean;
    //                std::cout<<"vMeanTimesFWMean ="<<vMeanTimesFWMean<<std::endl;

    //get macro dispersion
    FieldVector<Scalar, dim> D = vMeanTimesFWMean - vTimesFWMean;
    //    std::cout<<"D ="<<D<<std::endl;

    //scale macrodispersion with pressure gradient = 1/distance
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

    //store D dependent on the averaged saturation
    switch (subProblemNumber)
    {
        case 1:
        coarseParameters_.getDispersion().addEntry(D, globalIdxCoarseCurrent,timeStepCounter_);
        coarseParameters_.getDispersionSat().addEntry(meanSat_, globalIdxCoarseCurrent, timeStepCounter_);
        break;
        case 3:
        coarseParameters_.getDispersion().addEntry(D, globalIdxCoarseCurrent,timeStepCounter_);
        coarseParameters_.getDispersionSat().addEntry(meanSat_, globalIdxCoarseCurrent, timeStepCounter_);
        coarseParameters_.getDispersionSat()[globalIdxCoarseCurrent][timeStepCounter_]*=0.5;//arithmetic mean of the averaged saturation of the subproblem in x- and y-direction
        break;
        default:
        break;
    }

    //update timeStepCounter
    timeStepCounter_++;
}//end of dispersion function

//!calculates a saturation dependent convective flux correction from the solution of the local fine-scale problems
template<class Grid, class Model,class CoarseScaleParameterType>void TimeLoopSubProbs<Grid, Model, CoarseScaleParameterType>::calculateFluxCorrection(Model& model, const int subProblemNumber,const int& globalIdxCoarseCurrent,const int& faceNumberCurrent, const GlobalPosition& globalPosFaceCoarseCurrent,const GlobalPosition& lowerLeft, const GlobalPosition& upperRight)
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

    //gridViews for coarse and fine subgrid and for coarse global grid
    const GridView& gridViewCoarse(model.variables.grid.levelView(0));
    const GridView& gridViewFine(model.variables.grid.levelView(fineLev));
    const HostGridView& gridViewHost(gridViewCoarse.grid().getHostGrid().levelView(0));

    //index sets for coarse subgrid and coarse global grid
    const IndexSet& indexSetCoarse = gridViewCoarse.indexSet();
    const HostIndexSet& indexSetCoarseHost = gridViewHost.indexSet();

    //some variable declarations
    meanSat_=0;

    FieldVector<Scalar, dim> velocityMean(0);//fine scale velocity averaged over current coarse cell
    Scalar fWMean=0;// fractional flow function values dependent on mean saturation
    FieldVector<Scalar, dim> vTimesFWMean(0);//mean of fine-scale velocity v times fine scale fractional flow function f

    Scalar pressureIn=0;//mean pressure of the coarse inlet cell
    Scalar pressureOut=0;//mean pressure of the current coarse cell
    Scalar satIn = 0;//mean saturation of the coarse inlet cell
    Scalar satOut = 0;//mean saturation of the current coarse cell

    Scalar elementVolumeCoarse = 1;//just initialise with value != 0
    Scalar elementVolumeFaceCoarse=0;//sum of volumes of the fine cells at the interface

    //start iteration over the fine grid elements
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

        //get saturation and pressure
        Scalar sat = (*model)[globalIdx];
        Scalar pressure = model.variables.pressure[globalIdx];

        //get the coarse scale entity the fine scale entity belongs to
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

        //decide on which coarse scale element the fine scale element is located
        if (globalIdxCoarseHost != globalIdxCoarseCurrent)
        {
            //weight saturation and pressure with fine cell volume
            satIn+=sat*elementVolume;
            pressureIn+=pressure*elementVolume;
        }
        if (globalIdxCoarseHost == globalIdxCoarseCurrent)
        {
            //weight saturation and pressure with fine cell volume
            satOut+=sat*elementVolume;
            pressureOut+=pressure*elementVolume;
            //get volume of coarse element
            elementVolumeCoarse = fatherPointerHost->geometry().integrationElement(localPosCoarse)*Dune::ReferenceElements<Scalar,dim>::general(gt).volume();
            meanSat_ += sat*elementVolume;
            //f as function of the mean saturation: is updated until meanSat_ is completely calculated
            fWMean=model.diffProblem.materialLaw.fractionalW((meanSat_/elementVolumeCoarse), globalPos, *eItHost, localPos);
        }

        //fluxVector for calculation of the element velocity
        FieldVector<Scalar, dim*2> fluxVector;

        //start iteration over faces of the current fine-scale element
        IntersectionIterator
        isItEnd = gridViewFine.template iend(*eIt);
        for (IntersectionIterator
                isIt = gridViewFine.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        {
            //make sure that the fine scale element is located at the current coarse scale element.
            if (globalIdxCoarseHost == globalIdxCoarseCurrent)
            {
                // get geometry type of face
                Dune::GeometryType faceGT = isIt->geometryInInside().type();

                // center in face's reference element
                const Dune::FieldVector<Scalar,dim-1>&
                faceLocal = Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

                // center of face in global coordinates
                const GlobalPosition& globalPosFace = isIt->geometry().global(faceLocal);

                Scalar faceArea = isIt->geometry().volume();

                int faceNumber = isIt->indexInInside();

                // get normal vector
                FieldVector<Scalar,dim> unitOuterNormal = isIt->unitOuterNormal(faceLocal);

                //store flux for calculation of the element velocity
                fluxVector[faceNumber] = model.variables.velocity[globalIdx][faceNumber] * unitOuterNormal;
                fluxVector[faceNumber] *= faceArea;
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

            //weight by element volume
            for (int i=0; i<dim; i++)
            {
                velocityMean[i] += elementVelocity[i]*elementVolume;
                vTimesFWMean[i] += (elementVelocity[i] * model.diffProblem.materialLaw.fractionalW(sat, globalPos, *eItHost, localPos))*elementVolume;
            }
        }
    }

    //dispersion coefficient of coarse inlet and outlet cell
    Dune::FieldVector<Scalar,dim> DIn(0);
    Dune::FieldVector<Scalar,dim> DOut(0);

    //divide values through total volume
    meanSat_ /= elementVolumeCoarse;
    satIn /= elementVolumeCoarse;
    satOut /= elementVolumeCoarse;

    //calculate the saturation gradient -> gradient between two coarse cells
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
    satGrad /= (dist*dist);

    //determine dispersion coefficients of the two neighbouring coarse cells
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
    }//end iteration over coarse cells -> dispersion coefficients

    //divide by total volume (one coarse cell) -> fine scale velocities were weighted with the element volume!
    velocityMean /= elementVolumeCoarse;
    vTimesFWMean /= elementVolumeCoarse;
    //        std::cout<<"vTimesFWMean ="<<vTimesFWMean<<std::endl;
    //        std::cout<<"velocityMean = "<<velocityMean<<std::endl;

    //calculate the term consisting of the product of the averaged quantities
    FieldVector<Scalar, dim> vMeanTimesFWMean = velocityMean;
    vMeanTimesFWMean *= fWMean;
    //        std::cout<<"meanveltimesfwmean ="<<vMeanTimesFWMean<<std::endl;

    //calculate fluctuation
    FieldVector<Scalar, dim> fluxCorrection = vMeanTimesFWMean - vTimesFWMean;
    //                std::cout<<"m1 = "<<fluxCorrection<<std::endl;

    //divide by total volume (one coarse cell) to get the average
    pressureIn /= elementVolumeCoarse;
    pressureOut /= elementVolumeCoarse;

    //calculate the pressure gradient between the two coarse cells
    Scalar pressureGrad = (pressureIn - pressureOut)/dist;
    //    			std::cout<<"pressureGrad = "<<pressureGrad<<std::endl;

    //scale by pressure gradient
    if (pressureGrad != 0)
    {
        fluxCorrection /= pressureGrad;
    }
    //        		std::cout<<"m2 = "<<fluxCorrection<<std::endl;

    //arithmetic average of the dispersion coefficients
    Scalar dispersiveFlux = DIn*satGrad;
    dispersiveFlux += DOut*satGrad;
    dispersiveFlux *= 0.5;

    //get absolute value of the normal to make the flux correction a scalar
    FieldVector<Scalar, dim> normal(distVec);
    normal /= dist;
    //    std::cout<<"normal = "<<normal<<std::endl;

    //                std::cout<<"fluxCorrection = "<<fluxCorrection<<std::endl;

    //store convective flux correction
    FieldVector<Scalar, 2*dim> newEntry(0);
    newEntry[faceNumberCurrent] = (fluxCorrection*normal-dispersiveFlux);
    coarseParameters_.getFluxCorr().addEntry(newEntry, globalIdxCoarseCurrent, timeStepCounter_);

    newEntry = 0;
    newEntry[faceNumberCurrent] = meanSat_;
    coarseParameters_.getFluxCorrSat().addEntry(newEntry, globalIdxCoarseCurrent, timeStepCounter_);

    timeStepCounter_++;
}//end convective fluxCorrection function

}
#endif
