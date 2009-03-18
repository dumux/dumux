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

/**
 * @file
 * @brief  Class for running a preprocessing for a local-global upscaling approach
 * @author Markus Wolff
 *
 * \defgroup MultiMulti
 */

namespace Dune
{
//! \ingroup MultiMulti
/*! Class for running a preprocessing for a local-global upscaling approach.
 * Local fine-scale problems are solved in order to obtain coarse-scale quantities.
 * Local problems consisting of one coarse cell and
 * local problems consisting of two coarse cells are solved.
 * The calculation of the coarse-scale quantities is defined within the time loop.*/
/*! Template parameters are:
 *
 * - Grid                       a DUNE grid type
 * - Scalar                     a scalar type (usually double)
 * - CoarseScaleParameterType   a class type defining the course scale quantities
 */

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

    //!solves the local fine scale problems to obtain coarse scale dispersion coefficients -> int, int are the different levels
    void calcDispersiveFluxCorrection(int,int);

    //!solves the local fine scale problems to obtain coarse scale flux correction -> int, int are the different levels
    void calcConvectiveFluxCorrection(int,int);

    //!export of the calculated data into *.dat files -> int is the coarse Level
    void exportData(int);

    void XProblem1(Dune::SubGrid<dim,Grid>&, const int&, const GlobalPosition&, const GlobalPosition&, const int&, bool, const int, const GlobalPosition&);
    void YProblem1(Dune::SubGrid<dim,Grid>&, const int&, const GlobalPosition&, const GlobalPosition&, const int&, bool, const int, const GlobalPosition&);
    void XProblem2(Dune::SubGrid<dim,Grid>&, const int&, const GlobalPosition&, const GlobalPosition&, const int&, bool, const int, const GlobalPosition&);
    void YProblem2(Dune::SubGrid<dim,Grid>&, const int&, const GlobalPosition&, const GlobalPosition&, const int&, bool, const int, const GlobalPosition&);

public:
    //! For sequential preprocessing.
    /**
     * \param grid grid object of type Grid
     * \param wettingPhase phase object of a type derived from class Fluid
     * \param nonWettingPhase phase object of a type derived from class Fluid
     * \param soil opject of a type derived from class Matrix2p. Contains the global soil parameters.
     * \param coarseParams object of type CoarseScaleParameterType containing the coarse-scale parameters to be computed
     * \param materialLaw object of type TwoPhaseRelations. Contains the global material laws.
     */
    UpscalingPreprocess(Grid& grid,Fluid& wettingPhase , Fluid& nonWettingPhase, Matrix2p<Grid, Scalar>& soil, CoarseScaleParameterType& coarseParams, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid, Scalar>))
    :grid_(grid), wettingPhase_(wettingPhase), nonWettingPhase_(nonWettingPhase),soil_(soil),coarseParameters_(coarseParams),materialLaw_(materialLaw)
    {}

#ifdef PARALLEL
    //! For parallel preprocessing.
    /**
     * \param grid grid object of type Grid
     * \param communicator object of type Dune::CollectiveCommunication. Default communicator type MPI_Comm.
     * \param wettingPhase phase object of a type derived from class Fluid
     * \param nonWettingPhase phase object of a type derived from class Fluid
     * \param soil opject of a type derived from class Matrix2p. Contains the global soil parameters.
     * \param coarseParams object of type CoarseScaleParameterType containing the coarse-scale parameters to be computed.
     * \param materialLaw object of type TwoPhaseRelations. Contains the global material laws.
     */
    UpscalingPreprocess(Grid& grid, Dune::CollectiveCommunication<MPI_Comm>& comm, Fluid& wettingPhase , Fluid& nonWettingPhase, Matrix2p<Grid, Scalar>& soil, CoarseScaleParameterType& coarseParams, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid, Scalar>))
    :grid_(grid), communicator_(comm), numProc_(comm.size()), wettingPhase_(wettingPhase), nonWettingPhase_(nonWettingPhase),soil_(soil),coarseParameters_(coarseParams),materialLaw_(materialLaw)
    {
        int rank = comm.rank();
        maxRank_ = comm.max(rank);
    }
#endif

    //! start a preprocessing procedure
    void preprocessexecute(int coarseLev = 0,int fineLev=-1)
    {
        calcDispersiveFluxCorrection(coarseLev,fineLev);
        calcConvectiveFluxCorrection(coarseLev,fineLev);
        exportData(coarseLev);

        return;
    }

private:
    Grid& grid_;//! the grid of the global domain
#ifdef PARALLEL
    const Dune::CollectiveCommunication<MPI_Comm>& communicator_;
    const int numProc_;
    int maxRank_;
#endif
    Fluid& wettingPhase_;//! the wetting phase fluid properties
    Fluid& nonWettingPhase_;//! the non-wetting phase fluid properties
    Matrix2p<Grid, Scalar>& soil_;//! the soil properties, defined for the global domain
    CoarseScaleParameterType& coarseParameters_;//! container of the coarse scale model parameter variables
    TwoPhaseRelations<Grid, Scalar>& materialLaw_;//! the material laws, defined for the global domain
};

//!solves the local fine scale problems to obtain coarse scale dispersion coefficients
template<class Grid, class Scalar, class CoarseScaleParameterType>void UpscalingPreprocess<Grid,Scalar,CoarseScaleParameterType>::calcDispersiveFluxCorrection(int coarseLev,int fineLev)
{
    //identification of the kind of preprocessing calculation in the timeloop: true = dispersion coefficient
    bool calcDispersion = true;

    //is a fine level given in the function call?
    if (fineLev<0)
    fineLev=grid_.maxLevel();

    //get some grid features
    const GridView& gridView = grid_.levelView(coarseLev);
    const IndexSet& indexSetCoarse = grid_.levelIndexSet(coarseLev);
    int size = indexSetCoarse.size(0);

    //set number of rows of the matrixes containing the coarse scale model parameter, which will be calculated in the preprocessing.
    coarseParameters_.getDispersion().resizeRows(size);
    coarseParameters_.getDispersionSat().resizeRows(size);

    //start iteration over the coarse scale elements
    ElementIterator eItCoarseEnd = gridView.template end<0>();
    for (ElementIterator eItCoarse = gridView.template begin<0>(); eItCoarse != eItCoarseEnd; ++eItCoarse)
    {
        int globalIdxCoarse = indexSetCoarse.index(*eItCoarse);

        //if PARALLEL is defined, the coarse cells will be distributed on the available processes
#ifdef PARALLEL
        int rank = communicator_.rank();
        bool takeEntity = false;

        for (int i = 0;i<numProc_;i++)
        {
            if (rank == i)
            {
                //all coarse cells with globalIdxCoarse < (size-1)
                if (globalIdxCoarse >= size/numProc_*i && globalIdxCoarse < size/numProc_*(i+1))
                {
                    takeEntity = true;
                }
            }
        }
        if (rank == maxRank_)
        {
            //coarse cell with globalIdxCoarse == (size-1)
            if (globalIdxCoarse == (size-1))
            {
                takeEntity = true;
            }
        }
        if (!takeEntity)
        {
            continue;
        }
#endif

        //create a subgrid which contains one coarse cell
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

        //get coordinates of the lower left and the upper right corner of the subdomain
        GlobalPosition lowerLeft = eItCoarse->geometry().corner(0);
        GlobalPosition upperRight = eItCoarse->geometry().corner(3);

        //call of the function, which solves the local problem in x-direction.
        XProblem1(subGrid, fineLev, lowerLeft, upperRight,globalIdxCoarse, calcDispersion);

        //call of the function, which solves the local problem in y-direction.
        YProblem1(subGrid, fineLev, lowerLeft, upperRight, globalIdxCoarse, calcDispersion);
    }

#ifdef PARALLEL
    communicator_.barrier();
#endif
    return;

}

//!solves the local fine scale problems to obtain coarse scale flux correction
template<class Grid, class Scalar, class CoarseScaleParameterType>void UpscalingPreprocess<Grid, Scalar, CoarseScaleParameterType>::calcConvectiveFluxCorrection(int coarseLev,int fineLev)
{
    //identification of the kind of preprocessing calculation in the timeloop: false = flux correction
    bool calcConvection = false;

    //is a fine level given in the function call?
    if (fineLev<0)
    fineLev=grid_.maxLevel();

    //get some grid features
    const GridView& gridView = grid_.levelView(coarseLev);
    const IndexSet& indexSetCoarse = grid_.levelIndexSet(coarseLev);
    int size = indexSetCoarse.size(0);

    //set number of rows of the matrixes containing the coarse scale model parameter, which will be calculated in the preprocessing.
    coarseParameters_.getFluxCorr().resizeRows(size);
    coarseParameters_.getFluxCorrSat().resizeRows(size);

    //start iteration over the coarse scale elements
    ElementIterator eItCoarseEnd = gridView.template end<0>();
    for (ElementIterator eItCoarse = gridView.template begin<0>(); eItCoarse != eItCoarseEnd; ++eItCoarse)
    {
        int globalIdxCoarse = indexSetCoarse.index(*eItCoarse);

        //if PARALLEL is defined, the coarse cells will be distributed on the available processes
#ifdef PARALLEL
        int rank = communicator_.rank();
        bool takeEntity = false;

        for (int i = 0;i<numProc_;i++)
        {
            if (rank == i)
            {
                //all coarse cells with globalIdxCoarse < (size-1)
                if (globalIdxCoarse >= size/numProc_*i && globalIdxCoarse < size/numProc_*(i+1))
                {
                    takeEntity = true;
                }
            }
        }
        if (rank == maxRank_)
        {
            //coarse cell with globalIdxCoarse == (size-1)
            if (globalIdxCoarse == (size-1))
            {
                takeEntity = true;
            }
        }

        if (!takeEntity)
        {
            continue;
        }
#endif

        //start iteration over the coarse cell interfaces -> local subdomains consist of two neighbouring coarse cells
        IntersectionIterator isItCoarseEnd = gridView.template iend(*eItCoarse);
        for (IntersectionIterator isItCoarse = gridView.ibegin(*eItCoarse); isItCoarse != isItCoarseEnd; ++isItCoarse)
        {
            if ( isItCoarse->neighbor() )
            {
                GeometryType gtCoarse = isItCoarse->intersectionSelfLocal().type();
                FieldVector<Scalar,dim-1> faceLocalCoarse = ReferenceElements<Scalar,dim-1>::general(gtCoarse).position(0,0);
                const GlobalPosition& globalPosFaceCoarse = isItCoarse->intersectionGlobal().global(faceLocalCoarse);

                // local number of face
                int faceNumberCoarse = isItCoarse->numberInSelf();
                ElementPointer neighborPointer = isItCoarse->outside();

                //create a subgrid which contains coarse cell and its neighbouring cell on the interface "faceNumberCoarse"
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

                GlobalPosition lowerLeft(0);
                GlobalPosition upperRight(0);

                //check which is the current interface
                switch (faceNumberCoarse)
                {
                    case 0://left face
                    //get coordinates of the lower left and the upper right corner of the domain
                    lowerLeft = neighborPointer->geometry().corner(0);
                    upperRight = eItCoarse->geometry().corner(3);

                    XProblem1(subGrid, fineLev, lowerLeft, upperRight,globalIdxCoarse,calcConvection, faceNumberCoarse, globalPosFaceCoarse);
                    break;
                    case 1://right face
                    //get coordinates of the lower left and the upper right corner of the domain
                    lowerLeft = eItCoarse->geometry().corner(0);
                    upperRight = neighborPointer->geometry().corner(3);

                    XProblem2(subGrid, fineLev, lowerLeft, upperRight,globalIdxCoarse,calcConvection, faceNumberCoarse,globalPosFaceCoarse);
                    break;
                    case 2://lower face
                    //get coordinates of the lower left and the upper right corner of the domain
                    lowerLeft = neighborPointer->geometry().corner(0);
                    upperRight = eItCoarse->geometry().corner(3);

                    YProblem1(subGrid, fineLev, lowerLeft, upperRight, globalIdxCoarse,calcConvection, faceNumberCoarse,globalPosFaceCoarse);
                    break;
                    case 3://upper face
                    //get coordinates of the lower left and the upper right corner of the domain
                    lowerLeft = eItCoarse->geometry().corner(0);;
                    upperRight = neighborPointer->geometry().corner(3);

                    YProblem2(subGrid, fineLev, lowerLeft, upperRight, globalIdxCoarse,calcConvection, faceNumberCoarse,globalPosFaceCoarse);
                    break;
                }
            }
        }
    }
#ifdef PARALLEL
    communicator_.barrier();
#endif

    return;
}

template<class Grid, class Scalar, class CoarseScaleParameterType>void UpscalingPreprocess<Grid, Scalar, CoarseScaleParameterType>::XProblem1(Dune::SubGrid<dim,Grid>& subGrid,const int& fineLev,const GlobalPosition& lowerLeft,const GlobalPosition& upperRight,const int& globalIdxCoarse, bool correctionType, const int faceNumberCoarse = 0, const GlobalPosition& globalPosFaceCoarse = *(new GlobalPosition(0)))
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
    subTimeloop.execute(subImpes, 1, globalIdxCoarse, lowerLeft, upperRight, correctionType, faceNumberCoarse, globalPosFaceCoarse);

    std::cout<<
#ifdef PARALLEL
    "Process "<< communicator_.rank() <<" : "<<""<<
#endif
    "---------- end subproblem x1, coarse cell "<<(globalIdxCoarse+1)<<" interface "<< faceNumberCoarse <<" ---------- "<<std::endl;
}

template<class Grid, class Scalar, class CoarseScaleParameterType>void UpscalingPreprocess<Grid, Scalar, CoarseScaleParameterType>::XProblem2(Dune::SubGrid<dim,Grid>& subGrid,const int& fineLev,const GlobalPosition& lowerLeft,const GlobalPosition& upperRight,const int& globalIdxCoarse, bool correctionType, const int faceNumberCoarse = 0, const GlobalPosition& globalPosFaceCoarse= *(new GlobalPosition(0)))
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
    subTimeloop.execute(subImpes, 2, globalIdxCoarse,lowerLeft, upperRight, correctionType, faceNumberCoarse, globalPosFaceCoarse);

    std::cout<<
#ifdef PARALLEL
    "Process "<< communicator_.rank() <<" : "<<""<<
#endif
    "---------- end subproblem x2, coarse cell "<<(globalIdxCoarse+1)<<" interface "<< faceNumberCoarse <<" ---------- "<<std::endl;
}

template<class Grid, class Scalar, class CoarseScaleParameterType>void UpscalingPreprocess<Grid, Scalar, CoarseScaleParameterType>::YProblem1(Dune::SubGrid<dim,Grid>& subGrid,const int& fineLev,const GlobalPosition& lowerLeft,const GlobalPosition& upperRight,const int& globalIdxCoarse, bool correctionType, const int faceNumberCoarse = 0, const GlobalPosition& globalPosFaceCoarse= *(new GlobalPosition(0)))
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

    subTimeloop.execute(subImpes, 3, globalIdxCoarse, lowerLeft, upperRight, correctionType, faceNumberCoarse,globalPosFaceCoarse);

    std::cout<<
#ifdef PARALLEL
    "Process "<< communicator_.rank() <<" : "<<""<<
#endif
    "---------- end subproblem y1, coarse cell "<<(globalIdxCoarse+1)<<" interface "<< faceNumberCoarse <<" ---------- "<<std::endl;
}
template<class Grid, class Scalar, class CoarseScaleParameterType>void UpscalingPreprocess<Grid, Scalar, CoarseScaleParameterType>::YProblem2(Dune::SubGrid<dim,Grid>& subGrid,const int& fineLev,const GlobalPosition& lowerLeft,const GlobalPosition& upperRight,const int& globalIdxCoarse, bool correctionType, const int faceNumberCoarse = 0, const GlobalPosition& globalPosFaceCoarse= *(new GlobalPosition(0)))
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

    subTimeloop.execute(subImpes, 4, globalIdxCoarse,lowerLeft, upperRight, correctionType, faceNumberCoarse,globalPosFaceCoarse);

    std::cout<<
#ifdef PARALLEL
    "Process "<< communicator_.rank() <<" : "<<""<<
#endif
    "---------- end subproblem y2, coarse cell "<<(globalIdxCoarse+1)<<" interface "<< faceNumberCoarse <<" ---------- "<<std::endl;
}

template<class Grid, class Scalar, class CoarseScaleParameterType>void UpscalingPreprocess<Grid, Scalar, CoarseScaleParameterType>::exportData(int coarseLev)
{
#ifdef PARALLEL

    VariableMatrix<Scalar, Dune::FieldVector<Scalar, dim> >& dispersion = coarseParameters_.getDispersion();
    VariableMatrix<Scalar, Dune::FieldVector<Scalar, 1> >& dispersionSat = coarseParameters_.getDispersionSat();
    VariableMatrix<Scalar, Dune::FieldVector<Scalar, dim*2> >& fluxCorr = coarseParameters_.getFluxCorr();
    VariableMatrix<Scalar, Dune::FieldVector<Scalar, dim*2> >& fluxCorrSat = coarseParameters_.getFluxCorrSat();

    //resize the matrixes of all processes to have the size of the assembled matrix
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
    //assemble the complete matrix of all matrices of all processes
    for (int i=0;i<dispersion.rowSize();i++)
    {
        for (unsigned int j=0;j<dispersion[i].size();j++)
        {
            FieldVector<Scalar,dim> entrySum = dispersion[i][j];
            communicator_.sum(&entrySum,1);
            dispersion[i][j]=entrySum;
        }
    }

    //resize the matrixes of all processes to have the size of the assembled matrix
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
    //assemble the complete matrix of all matrices of all processes
    for (int i=0;i<dispersionSat.rowSize();i++)
    {
        for (unsigned int j=0;j<dispersionSat[i].size();j++)
        {
            FieldVector<Scalar,1> entrySum = dispersionSat[i][j];
            communicator_.sum(&entrySum,1);
            dispersionSat[i][j]=entrySum;
        }
    }

    //resize the matrixes of all processes to have the size of the assembled matrix
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
    //assemble the complete matrix of all matrices of all processes
    for (int i=0;i<fluxCorr.rowSize();i++)
    {
        for (unsigned int j=0;j<fluxCorr[i].size();j++)
        {
            FieldVector<Scalar,dim*2> entrySum = fluxCorr[i][j];
            communicator_.sum(&entrySum,1);
            fluxCorr[i][j]=entrySum;
        }
    }

    //resize the matrixes of all processes to have the size of the assembled matrix
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
    //assemble the complete matrix of all matrices of all processes
    for (int i=0;i<fluxCorrSat.rowSize();i++)
    {
        for (unsigned int j=0;j<fluxCorrSat[i].size();j++)
        {
            FieldVector<Scalar,dim*2> entrySum = fluxCorrSat[i][j];
            communicator_.sum(&entrySum,1);
            fluxCorrSat[i][j]=entrySum;
        }
    }

    //call function which writes data into files -> *.dat
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

    //call function which writes data into files -> *.dat
    exportToFile(dispersionData, "dataDispersion");
    exportToFile(dispersionSatData, "dataDispersionSat");
    exportToFile(mData,"dataM");
    exportToFile(mSatData,"dataMSat");
#endif

    return;
}
}
#endif
