// $Id$
#ifndef DUNE_UPSCALINGPREPROCESS_HH
#define DUNE_UPSCALINGPREPROCESS_HH

#include <dumux/material/property_baseclasses.hh>
#include "variableclasssubprobs.hh"

#include "dumux/upscaledsaturation/subproblems/dispersionsubproblem.hh"
#include "dumux/upscaledsaturation/subproblems/convectionsubproblem.hh"
#include "dumux/upscaledsaturation/preprocess/fvdiffusionsubprobs.hh"
#include "dumux/upscaledsaturation/preprocess/fvtransportsubprobs.hh"
#include "dumux/upscaledsaturation/preprocess/impessubprobs.hh"
#include "dumux/upscaledsaturation/preprocess/timeloopsubprobs.hh"
#include "dumux/upscaledsaturation/variablematrix.hh"
#include "dumux/io/exportcorrection.hh"

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar, class CoarseScaleParameterType> class UpscalingPreprocess
{
    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };

typedef    typename Grid::LevelGridView GridView;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename Grid::template Codim<0>::HierarchicIterator HierarchicElementIterator;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

    void calcDispersiveFluxCorrection(int,int);
    void calcConvectiveFluxCorrection(int,int);
    void exportData(int);

    void XProblem1(Dune::SubGrid<dim,Grid>&, const int&, const GlobalPosition&, const GlobalPosition&, const int&, bool, const int, const GlobalPosition&);
    void YProblem1(Dune::SubGrid<dim,Grid>&, const int&, const GlobalPosition&, const GlobalPosition&, const int&, bool, const int, const GlobalPosition&);
    void XProblem2(Dune::SubGrid<dim,Grid>&, const int&, const GlobalPosition&, const GlobalPosition&, const int&, bool, const int, const GlobalPosition&);
    void YProblem2(Dune::SubGrid<dim,Grid>&, const int&, const GlobalPosition&, const GlobalPosition&, const int&, bool, const int, const GlobalPosition&);

public:
    UpscalingPreprocess(Grid& grid,Fluid& wettingPhase , Fluid& nonWettingPhase, Matrix2p<Grid, Scalar>& soil, CoarseScaleParameterType& coarseParams, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid, Scalar>), Scalar tEnd = 1e6)
    :grid_(grid), wettingPhase_(wettingPhase), nonWettingPhase_(nonWettingPhase),soil_(soil),coarseParameters_(coarseParams),materialLaw_(materialLaw),eps_(1e-6)
    {}

#ifdef PARALLEL
    UpscalingPreprocess(Grid& grid, Dune::CollectiveCommunication<MPI_Comm>& comm,Fluid& wettingPhase , Fluid& nonWettingPhase, Matrix2p<Grid, Scalar>& soil, CoarseScaleParameterType& coarseParams, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid, Scalar>), Scalar tEnd = 1e6)
    :grid_(grid), communicator_(comm),wettingPhase_(wettingPhase), nonWettingPhase_(nonWettingPhase),soil_(soil),coarseParameters_(coarseParams),materialLaw_(materialLaw),eps_(1e-6)
    {}
#endif

    void preprocessexecute(int coarseLev = 0,int fineLev=-1)
    {
        calcDispersiveFluxCorrection(coarseLev,fineLev);
        calcConvectiveFluxCorrection(coarseLev,fineLev);
        exportData(coarseLev);

        return;
    }

private:
    Grid& grid_;
    Fluid& wettingPhase_;
    Fluid& nonWettingPhase_;
    Matrix2p<Grid, Scalar>& soil_;
    CoarseScaleParameterType& coarseParameters_;
    TwoPhaseRelations<Grid, Scalar>& materialLaw_;
    Scalar eps_;
#ifdef PARALLEL
    const Dune::CollectiveCommunication<MPI_Comm>& communicator_;
#endif
};

template<class Grid, class Scalar, class CoarseScaleParameterType>void UpscalingPreprocess<Grid,Scalar,CoarseScaleParameterType>::calcDispersiveFluxCorrection(int coarseLev,int fineLev)
{
    bool calcDispersion = true;

    if (fineLev<0)
    fineLev=grid_.maxLevel();

    const GridView& gridView = grid_.levelView(coarseLev);

    const IndexSet& indexSetCoarse = grid_.levelIndexSet(coarseLev);
    int size = indexSetCoarse.size(0);

    coarseParameters_.getDispersion().resizeRows(size);
    coarseParameters_.getDispersionSat().resizeRows(size);

#ifdef PARALLEL
    int numProc = communicator_.size();
#endif

    ElementIterator eItCoarseEnd = gridView.template end<0>();
    for (ElementIterator eItCoarse = gridView.template begin<0>(); eItCoarse != eItCoarseEnd; ++eItCoarse)
    {
        int globalIdxCoarse = indexSetCoarse.index(*eItCoarse);
        //            std::cout<<"coarse cell "<<globalIdxCoarse<<std::endl;

#ifdef PARALLEL
        int rank = communicator_.rank();
        bool takeEntity = false;

        for (int i = 0;i<numProc;i++)
        {
            if (rank == i)
            {
                if (globalIdxCoarse >= size/numProc*i && globalIdxCoarse < size/numProc*(i+1))
                {
                    takeEntity = true;
                }
            }
        }
        if (!takeEntity)
        {
            continue;
        }
#endif

        Dune::SubGrid<dim,Grid> subGrid(grid_);

        subGrid.createBegin();

        HierarchicElementIterator eItEnd = eItCoarse-> hend(fineLev);
        for (HierarchicElementIterator eIt = eItCoarse->hbegin(fineLev); eIt != eItEnd; ++eIt)
        {
            if (eIt->level() == fineLev)
            {
                subGrid.addPartial(eIt);
            }
        }
        subGrid.createEnd();

        GlobalPosition lowerLeft = eItCoarse->geometry().corner(0);
        GlobalPosition upperRight = eItCoarse->geometry().corner(3);

        XProblem1(subGrid, fineLev, lowerLeft, upperRight,globalIdxCoarse, calcDispersion);
        YProblem1(subGrid, fineLev, lowerLeft, upperRight, globalIdxCoarse, calcDispersion);
    }

#ifdef PARALLEL
    communicator_.barrier();
#endif
    return;

}

template<class Grid, class Scalar, class CoarseScaleParameterType>void UpscalingPreprocess<Grid, Scalar, CoarseScaleParameterType>::calcConvectiveFluxCorrection(int coarseLev,int fineLev)
{
    bool calcConvection = false;

    if (fineLev<0)
    fineLev=grid_.maxLevel();

    const GridView& gridView = grid_.levelView(coarseLev);

    const IndexSet& indexSetCoarse = grid_.levelIndexSet(coarseLev);
    int size = indexSetCoarse.size(0);

    coarseParameters_.getFluxCorr().resizeRows(size);
    coarseParameters_.getFluxCorrSat().resizeRows(size);

#ifdef PARALLEL
    int numProc = communicator_.size();
#endif

    ElementIterator eItCoarseEnd = gridView.template end<0>();
    for (ElementIterator eItCoarse = gridView.template begin<0>(); eItCoarse != eItCoarseEnd; ++eItCoarse)
    {
        int globalIdxCoarse = indexSetCoarse.index(*eItCoarse);

#ifdef PARALLEL
        int rank = communicator_.rank();
        bool takeEntity = false;

        for (int i = 0;i<numProc;i++)
        {
            if (rank == i)
            {
                if (globalIdxCoarse >= size/numProc*i && globalIdxCoarse < size/numProc*(i+1))
                {
                    takeEntity = true;
                }
            }
        }
        if (!takeEntity)
        {
            continue;
        }
#endif

        IntersectionIterator isItCoarseEnd = gridView.template iend(*eItCoarse);
        for (IntersectionIterator isItCoarse = gridView.ibegin(*eItCoarse); isItCoarse != isItCoarseEnd; ++isItCoarse)
        {
            if ( isItCoarse->neighbor() )
            {
                GeometryType gtCoarse = isItCoarse->intersectionSelfLocal().type();
                FieldVector<Scalar,dim-1> faceLocalCoarse = ReferenceElements<Scalar,dim-1>::general(gtCoarse).position(0,0);
                const GlobalPosition& globalPosFaceCoarse = isItCoarse->intersectionGlobal().global(faceLocalCoarse);

                // local number of facet
                int numberInSelf = isItCoarse->numberInSelf();
                ElementPointer neighborPointer = isItCoarse->outside();

                Dune::SubGrid<dim,Grid> subGrid(grid_);
                subGrid.createBegin();

                HierarchicElementIterator eItEnde = eItCoarse-> hend(fineLev);
                for (HierarchicElementIterator eIt = eItCoarse->hbegin(fineLev); eIt != eItEnde; ++eIt)
                {
                    if (eIt->level() == fineLev)
                    {
                        subGrid.addPartial(eIt);
                    }
                }

                HierarchicElementIterator eItEnd = neighborPointer-> hend(fineLev);
                for (HierarchicElementIterator eIt = neighborPointer->hbegin(fineLev); eIt != eItEnd; ++eIt)
                {
                    if (eIt->level() == fineLev)
                    {
                        subGrid.addPartial(eIt);
                    }
                }
                subGrid.createEnd();

                GlobalPosition lowerLeftHelp = eItCoarse->geometry().corner(0);
                GlobalPosition upperRightHelp = eItCoarse->geometry().corner(3);

                if (globalPosFaceCoarse[0] <= lowerLeftHelp[0] + eps_ && globalPosFaceCoarse[0] >= lowerLeftHelp[0] - eps_)
                {
                    GlobalPosition lowerLeft = neighborPointer->geometry().corner(0);
                    GlobalPosition upperRight = upperRightHelp;

                    XProblem1(subGrid, fineLev, lowerLeft, upperRight,globalIdxCoarse,calcConvection, numberInSelf, globalPosFaceCoarse);
                }
                if (globalPosFaceCoarse[0] <= upperRightHelp[0] + eps_ && globalPosFaceCoarse[0] >= upperRightHelp[0] - eps_)
                {
                    GlobalPosition lowerLeft = lowerLeftHelp;
                    GlobalPosition upperRight = neighborPointer->geometry().corner(3);

                    XProblem2(subGrid, fineLev, lowerLeft, upperRight,globalIdxCoarse,calcConvection, numberInSelf,globalPosFaceCoarse);
                }
                if (globalPosFaceCoarse[1] <= lowerLeftHelp[1] + eps_ && globalPosFaceCoarse[1] >= lowerLeftHelp[1] - eps_)
                {
                    GlobalPosition lowerLeft = neighborPointer->geometry().corner(0);
                    GlobalPosition upperRight = upperRightHelp;

                    YProblem1(subGrid, fineLev, lowerLeft, upperRight, globalIdxCoarse,calcConvection, numberInSelf,globalPosFaceCoarse);
                }
                if (globalPosFaceCoarse[1] <= upperRightHelp[1] + eps_ && globalPosFaceCoarse[1] >= upperRightHelp[1] - eps_)
                {
                    GlobalPosition lowerLeft = lowerLeftHelp;
                    GlobalPosition upperRight = neighborPointer->geometry().corner(3);

                    YProblem2(subGrid, fineLev, lowerLeft, upperRight, globalIdxCoarse,calcConvection, numberInSelf,globalPosFaceCoarse);
                }
            }
        }
    }
#ifdef PARALLEL
    communicator_.barrier();
#endif

    return;
}

template<class Grid, class Scalar, class CoarseScaleParameterType>void UpscalingPreprocess<Grid, Scalar, CoarseScaleParameterType>::XProblem1(Dune::SubGrid<dim,Grid>& subGrid,const int& fineLev,const GlobalPosition& lowerLeft,const GlobalPosition& upperRight,const int& globalIdxCoarse, bool correctionType, const int numberInSelf = 0, const GlobalPosition& globalPosFaceCoarse = *(new GlobalPosition(0)))
{
    typedef Dune::SubGrid<dim,Grid> SubGrid;
    typedef Dune::VariableClassSubProbs<SubGrid, Scalar> VC;
    typedef Dune::ConvSubProblemX1<SubGrid, Scalar, VC> ConvSubProblem;
    typedef Dune::FVDiffSubProbs<SubGrid, Scalar, VC> Diffusion;
    typedef Dune::FVTransSubProbs<SubGrid, Scalar,VC> Transport;
    typedef Dune::IMPESSubProbs<SubGrid, Diffusion, Transport,VC> IMPES;

    VC subProblemVariables(subGrid, fineLev);
    ConvSubProblem subProblem(subProblemVariables, wettingPhase_, nonWettingPhase_, soil_, materialLaw_, lowerLeft, upperRight);

    Diffusion subDiffusion(subGrid,subProblem);
    Transport subTransport(subGrid,subProblem);
    IMPES subImpes(subDiffusion,subTransport);
    TimeLoopSubProbs<SubGrid, IMPES,CoarseScaleParameterType> subTimeloop(coarseParameters_);
    subTimeloop.execute(subImpes, 1, globalIdxCoarse, lowerLeft, upperRight, correctionType, numberInSelf, globalPosFaceCoarse);

    std::cout<<
#ifdef PARALLEL
    "Process "<< communicator_.rank() <<" : "<<""<<
#endif
    "---------- end subproblem x1, coarse cell "<<(globalIdxCoarse+1)<<" interface "<< numberInSelf <<" ---------- "<<std::endl;
}

template<class Grid, class Scalar, class CoarseScaleParameterType>void UpscalingPreprocess<Grid, Scalar, CoarseScaleParameterType>::XProblem2(Dune::SubGrid<dim,Grid>& subGrid,const int& fineLev,const GlobalPosition& lowerLeft,const GlobalPosition& upperRight,const int& globalIdxCoarse, bool correctionType, const int numberInSelf = 0, const GlobalPosition& globalPosFaceCoarse= *(new GlobalPosition(0)))
{
    typedef Dune::SubGrid<dim,Grid> SubGrid;
    typedef Dune::VariableClassSubProbs<SubGrid, Scalar> VC;
    typedef Dune::ConvSubProblemX2<SubGrid, Scalar, VC> ConvSubProblem;
    typedef Dune::FVDiffSubProbs<SubGrid, Scalar, VC> Diffusion;
    typedef Dune::FVTransSubProbs<SubGrid, Scalar,VC> Transport;
    typedef Dune::IMPESSubProbs<SubGrid, Diffusion, Transport,VC> IMPES;

    VC subProblemVariables(subGrid, fineLev);
    ConvSubProblem subProblem(subProblemVariables, wettingPhase_, nonWettingPhase_, soil_, materialLaw_, lowerLeft, upperRight);

    Diffusion subDiffusion(subGrid,subProblem);
    Transport subTransport(subGrid,subProblem);
    IMPES subImpes(subDiffusion,subTransport);
    TimeLoopSubProbs<SubGrid, IMPES,CoarseScaleParameterType> subTimeloop(coarseParameters_);
    subTimeloop.execute(subImpes, 2, globalIdxCoarse,lowerLeft, upperRight, correctionType, numberInSelf, globalPosFaceCoarse);

    std::cout<<
#ifdef PARALLEL
    "Process "<< communicator_.rank() <<" : "<<""<<
#endif
    "---------- end subproblem x2, coarse cell "<<(globalIdxCoarse+1)<<" interface "<< numberInSelf <<" ---------- "<<std::endl;
}

template<class Grid, class Scalar, class CoarseScaleParameterType>void UpscalingPreprocess<Grid, Scalar, CoarseScaleParameterType>::YProblem1(Dune::SubGrid<dim,Grid>& subGrid,const int& fineLev,const GlobalPosition& lowerLeft,const GlobalPosition& upperRight,const int& globalIdxCoarse, bool correctionType, const int numberInSelf = 0, const GlobalPosition& globalPosFaceCoarse= *(new GlobalPosition(0)))
{
    typedef Dune::SubGrid<dim,Grid> SubGrid;
    typedef Dune::VariableClassSubProbs<SubGrid, Scalar> VC;
    typedef Dune::ConvSubProblemY1<SubGrid, Scalar, VC> ConvSubProblem;
    typedef Dune::FVDiffSubProbs<SubGrid, Scalar, VC> Diffusion;
    typedef Dune::FVTransSubProbs<SubGrid, Scalar,VC> Transport;
    typedef Dune::IMPESSubProbs<SubGrid, Diffusion, Transport,VC> IMPES;

    VC subProblemVariables(subGrid, fineLev);
    ConvSubProblem subProblem(subProblemVariables, wettingPhase_, nonWettingPhase_, soil_, materialLaw_, lowerLeft, upperRight);

    Diffusion subDiffusion(subGrid,subProblem);
    Transport subTransport(subGrid,subProblem);
    IMPES subImpes(subDiffusion,subTransport);
    TimeLoopSubProbs<SubGrid, IMPES,CoarseScaleParameterType> subTimeloop(coarseParameters_);

    subTimeloop.execute(subImpes, 3, globalIdxCoarse, lowerLeft, upperRight, correctionType, numberInSelf,globalPosFaceCoarse);

    std::cout<<
#ifdef PARALLEL
    "Process "<< communicator_.rank() <<" : "<<""<<
#endif
    "---------- end subproblem y1, coarse cell "<<(globalIdxCoarse+1)<<" interface "<< numberInSelf <<" ---------- "<<std::endl;
}
template<class Grid, class Scalar, class CoarseScaleParameterType>void UpscalingPreprocess<Grid, Scalar, CoarseScaleParameterType>::YProblem2(Dune::SubGrid<dim,Grid>& subGrid,const int& fineLev,const GlobalPosition& lowerLeft,const GlobalPosition& upperRight,const int& globalIdxCoarse, bool correctionType, const int numberInSelf = 0, const GlobalPosition& globalPosFaceCoarse= *(new GlobalPosition(0)))
{
    typedef Dune::SubGrid<dim,Grid> SubGrid;
    typedef Dune::VariableClassSubProbs<SubGrid, Scalar> VC;
    typedef Dune::ConvSubProblemY2<SubGrid, Scalar, VC> ConvSubProblem;
    typedef Dune::FVDiffSubProbs<SubGrid, Scalar, VC> Diffusion;
    typedef Dune::FVTransSubProbs<SubGrid, Scalar,VC> Transport;
    typedef Dune::IMPESSubProbs<SubGrid, Diffusion, Transport,VC> IMPES;

    VC subProblemVariables(subGrid, fineLev);
    ConvSubProblem subProblem(subProblemVariables, wettingPhase_, nonWettingPhase_, soil_, materialLaw_, lowerLeft, upperRight);

    Diffusion subDiffusion(subGrid,subProblem);
    Transport subTransport(subGrid,subProblem);
    IMPES subImpes(subDiffusion,subTransport);
    TimeLoopSubProbs<SubGrid, IMPES,CoarseScaleParameterType> subTimeloop(coarseParameters_);

    subTimeloop.execute(subImpes, 4, globalIdxCoarse,lowerLeft, upperRight, correctionType, numberInSelf,globalPosFaceCoarse);

    std::cout<<
#ifdef PARALLEL
    "Process "<< communicator_.rank() <<" : "<<""<<
#endif
    "---------- end subproblem y2, coarse cell "<<(globalIdxCoarse+1)<<" interface "<< numberInSelf <<" ---------- "<<std::endl;
}

template<class Grid, class Scalar, class CoarseScaleParameterType>void UpscalingPreprocess<Grid, Scalar, CoarseScaleParameterType>::exportData(int coarseLev)
{
#ifdef PARALLEL

    VariableMatrix<Scalar, Dune::FieldVector<Scalar, dim> >& dispersion = coarseParameters_.getDispersion();
    VariableMatrix<Scalar, Dune::FieldVector<Scalar, 1> >& dispersionSat = coarseParameters_.getDispersionSat();
    VariableMatrix<Scalar, Dune::FieldVector<Scalar, dim*2> >& fluxCorr = coarseParameters_.getFluxCorr();
    VariableMatrix<Scalar, Dune::FieldVector<Scalar, dim*2> >& fluxCorrSat = coarseParameters_.getFluxCorrSat();

    for (int i=0;i<dispersion.rowSize();i++)
    {
        int columnSize = dispersion[i].size();

        communicator_.sum(&columnSize,1);

        if (dispersion[i].size() == 0)
        {
            dispersion.resizeColumnX(columnSize,i);
            for (int j = 0; j < columnSize; j++)
            {
                dispersion[i][j]=0;
            }
        }
    }
    for (int i=0;i<dispersion.rowSize();i++)
    {
        for (unsigned int j=0;j<dispersion[i].size();j++)
        {
            FieldVector<Scalar,dim> entrySum = dispersion[i][j];
            communicator_.sum(&entrySum,1);
            dispersion[i][j]=entrySum;
        }
    }

    for (int i=0;i<dispersionSat.rowSize();i++)
    {
        int columnSize = dispersionSat[i].size();

        communicator_.sum(&columnSize,1);

        if (dispersionSat[i].size() == 0)
        {
            dispersionSat.resizeColumnX(columnSize,i);
            for (int j = 0; j < columnSize; j++)
            {
                dispersionSat[i][j]=0;
            }
        }
    }
    for (int i=0;i<dispersionSat.rowSize();i++)
    {
        for (unsigned int j=0;j<dispersionSat[i].size();j++)
        {
            FieldVector<Scalar,1> entrySum = dispersionSat[i][j];
            communicator_.sum(&entrySum,1);
            dispersionSat[i][j]=entrySum;
        }
    }

    for (int i=0;i<fluxCorr.rowSize();i++)
    {
        int columnSize = fluxCorr[i].size();

        communicator_.sum(&columnSize,1);

        if (fluxCorr[i].size() == 0)
        {
            fluxCorr.resizeColumnX(columnSize,i);
            for (int j = 0; j < columnSize; j++)
            {
                fluxCorr[i][j]=0;
            }
        }
    }
    for (int i=0;i<fluxCorr.rowSize();i++)
    {
        for (unsigned int j=0;j<fluxCorr[i].size();j++)
        {
            FieldVector<Scalar,dim*2> entrySum = fluxCorr[i][j];
            communicator_.sum(&entrySum,1);
            fluxCorr[i][j]=entrySum;
        }
    }

    for (int i=0;i<fluxCorrSat.rowSize();i++)
    {
        int columnSize = fluxCorrSat[i].size();

        communicator_.sum(&columnSize,1);

        if (fluxCorrSat[i].size() == 0)
        {
            fluxCorrSat.resizeColumnX(columnSize,i);
            for (int j = 0; j < columnSize; j++)
            {
                fluxCorrSat[i][j]=0;
            }
        }
    }
    for (int i=0;i<fluxCorrSat.rowSize();i++)
    {
        for (unsigned int j=0;j<fluxCorrSat[i].size();j++)
        {
            FieldVector<Scalar,dim*2> entrySum = fluxCorrSat[i][j];
            communicator_.sum(&entrySum,1);
            fluxCorrSat[i][j]=entrySum;
        }
    }

    if (communicator_.rank()==0)
    {
        exportToFile(dispersion.getMatrix(), "dataDispersion");
        exportToFile(dispersionSat.getMatrix(), "dataDispersionSat");
        exportToFile(fluxCorr.getMatrix(),"dataM");
        exportToFile(fluxCorrSat.getMatrix(),"dataMSat");
    }
#else
    std::vector<std::vector<Dune::FieldVector<Scalar, dim> > >& dispersionData = coarseParameters_.getDispersion().writeMatrix();
    std::vector<std::vector<Dune::FieldVector<Scalar, 1> > >& dispersionSatData = coarseParameters_.getDispersionSat().writeMatrix();

    std::vector<std::vector<Dune::FieldVector<Scalar, 2*dim> > >& mData = coarseParameters_.getFluxCorr().writeMatrix();
    std::vector<std::vector<Dune::FieldVector<Scalar, 2*dim> > >& mSatData = coarseParameters_.getFluxCorrSat().writeMatrix();

    exportToFile(dispersionData, "dataDispersion");
    exportToFile(dispersionSatData, "dataDispersionSat");
    exportToFile(mData,"dataM");
    exportToFile(mSatData,"dataMSat");
#endif

    return;
}
}
#endif
