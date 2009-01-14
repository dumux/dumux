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
#include "dumux/io/exporttodgf.hh"
//#include <dune/istl/matrix.hh>
//#include <dune/common/fvector.hh>

namespace Dune
{

template<class G, class RT> class UpscalingPreprocess
{
    enum
    {
        dim = G::dimension, dimworld = G::dimensionworld
    };

typedef    typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename G::ctype ct;
    typedef typename G::LevelGridView GV;
    typedef typename GV::IndexSet IS;
    typedef typename GV::template Codim<0>::Iterator Iterator;
    typedef typename GV::IntersectionIterator IntersectionIterator;

    typedef typename G::template Codim<0>::HierarchicIterator HierarchicIterator;
    typedef typename G::template Codim<0>::EntityPointer EntityPointer;

    void calcdispersivecorrection(int,int);
    void calcconvectivecorrection(int,int);
    void exportdata(int);

    void TestProblem(Dune::SubGrid<dim,G>&, int&, FieldVector<ct,dim>&, FieldVector<ct,dim>&, int&, int&,FieldVector<ct, dim>& ,int&);
    void XProblem1(Dune::SubGrid<dim,G>&, int&, FieldVector<ct,dim>&, FieldVector<ct,dim>&, int&, int&, FieldVector<ct, dim>&, int&);
    void YProblem1(Dune::SubGrid<dim,G>&, int&, FieldVector<ct,dim>&, FieldVector<ct,dim>&, int&, int&, FieldVector<ct, dim>&, int&);
    void XProblem2(Dune::SubGrid<dim,G>&, int&, FieldVector<ct,dim>&, FieldVector<ct,dim>&, int&, int&, FieldVector<ct, dim>&, int&);
    void YProblem2(Dune::SubGrid<dim,G>&, int&, FieldVector<ct,dim>&, FieldVector<ct,dim>&, int&, int&, FieldVector<ct, dim>&, int&);

public:
    UpscalingPreprocess(G& g, Fluid& wP , Fluid& nwP, Matrix2p<G, RT>& s, TwoPhaseRelations<G, RT>& law = *(new TwoPhaseRelations<G,RT>), double tEnd = 1e6)
    :grid_(g), wettingPhase_(wP), nonWettingPhase_(nwP),soil_(s),materialLaw_(law), tEnd_(tEnd), nDt_(0), firstRun_(true)
    {}

    void preprocessexecute(int coarseLev = 0,int fineLev=-1)
    {
        calcdispersivecorrection(coarseLev,fineLev);
        calcconvectivecorrection(coarseLev,fineLev);
        exportdata(coarseLev);
    }

private:
    G& grid_;
    Fluid& wettingPhase_;
    Fluid& nonWettingPhase_;
    Matrix2p<G, RT>& soil_;
    TwoPhaseRelations<G, RT>& materialLaw_;
    double tEnd_;
    int nDt_;
    bool firstRun_;
};

template<class G, class RT>void UpscalingPreprocess<G,RT>::calcdispersivecorrection(int coarseLev,int fineLev)
{
    typedef Dune::SubGrid<dim,G> SG;
    typedef Dune::VariableClassSubProbs<SG, RT> VC;
    typedef Dune::DiffSubProblemX<SG, RT, VC> DiffSubProblemX;
    typedef Dune::DiffSubProblemY<SG, RT, VC> DiffSubProblemY;
    typedef Dune::FVDiffSubProbs<SG, RT, VC> Diffusion;
    typedef Dune::FVTransSubProbs<SG, RT,VC> Transport;
    typedef Dune::IMPESSubProbs<SG, Diffusion, Transport,VC> IMPES;

    if (fineLev<0)
    fineLev=grid_.maxLevel();

    const GV& gridView = grid_.levelView(coarseLev);

    const IS& coarseIset = grid_.levelIndexSet(coarseLev);
    int size = coarseIset.size(0);

    Iterator endCIt = gridView.template end<0>();
    for (Iterator cIt = gridView.template begin<0>(); cIt != endCIt; ++cIt)
    {
        int coarseIndex = coarseIset.index(*cIt);
        //            std::cout<<"coarse cell "<<coarseIndex<<std::endl;

        Dune::SubGrid<dim,G> subGrid(grid_);

        subGrid.createBegin();

        HierarchicIterator endIt = cIt-> hend(fineLev);
        for (HierarchicIterator it = cIt->hbegin(fineLev); it != endIt; ++it)
        {
            if (it->level() == fineLev)
            {
                subGrid.addPartial(it);
            }
        }
        subGrid.createEnd();

        FieldVector<ct,dim> lowerLeft = cIt->geometry().corner(0);
        FieldVector<ct,dim> upperRight = cIt->geometry().corner(3);

        //            std::cout<<"LowerLeft = "<<LowerLeft<<" UpperRight = "<< UpperRight<<std::endl;

        if(firstRun_)
        {
            VC subProblemVariablesTest(subGrid, fineLev);

            DiffSubProblemX subProblemTest(subProblemVariablesTest, wettingPhase_, nonWettingPhase_, soil_, materialLaw_,lowerLeft,upperRight);

            Diffusion subDiffusionTest(subGrid,subProblemTest,fineLev);
            Transport subTransportTest(subGrid,subProblemTest,fineLev);
            IMPES subImpesTest(subDiffusionTest,subTransportTest);
            TimeLoopSubProbs<SG, IMPES> subTimeloopTest(tEnd_);
            subTimeloopTest.executeD(subImpesTest, 1, coarseIndex, lowerLeft, upperRight, nDt_,tEnd_ , firstRun_);

            soil_.getDispersion().setSize(size,nDt_);
            soil_.getDispersionSat().setSize(size,nDt_);

            soil_.getDispersion()=0;
            soil_.getDispersionSat()=0;

            firstRun_=false;
            std::cout<<"-------------------- end dispersion testsubproblem -------------------- "<<std::endl;
        }

        VC subProblemVariables1(subGrid, fineLev);
        VC subProblemVariables2(subGrid, fineLev);

        DiffSubProblemX subProblem1(subProblemVariables1, wettingPhase_, nonWettingPhase_, soil_, materialLaw_, lowerLeft, upperRight);
        DiffSubProblemY subProblem2(subProblemVariables2, wettingPhase_, nonWettingPhase_, soil_, materialLaw_, lowerLeft, upperRight);

        Diffusion subDiffusion1(subGrid,subProblem1,fineLev);
        Transport subTransport1(subGrid,subProblem1,fineLev);
        IMPES subImpes1(subDiffusion1,subTransport1);
        TimeLoopSubProbs<SG, IMPES> subTimeloop1(tEnd_);
        subTimeloop1.executeD(subImpes1, 1, coarseIndex, lowerLeft, upperRight, nDt_);
        //        std::cout<<"saturation ="<<subProblemVariables1.saturation<<"pressure = "<<subProblemVariables1.pressure<<std::endl;

        std::cout<<"---------- end subproblem 1, coarse cell "<<(coarseIndex+1)<<" ---------- "<<std::endl;

        Diffusion subDiffusion2(subGrid,subProblem2,fineLev);
        Transport subTransport2(subGrid,subProblem2,fineLev);
        IMPES subImpes2(subDiffusion2,subTransport2);
        TimeLoopSubProbs<SG, IMPES> subTimeloop2(tEnd_);
        subTimeloop2.executeD(subImpes2, 2, coarseIndex, lowerLeft, upperRight, nDt_);
        //        std::cout<<"saturation ="<<subProblemVariables2.saturation<<"pressure = "<<subProblemVariables2.pressure<<std::endl;

        std::cout<<"---------- end subproblem 2, coarse cell "<<(coarseIndex+1)<<" ---------- "<<std::endl;
    }

    //        Dune::Matrix<Dune::FieldVector<double,1> > Dxx(size,nDt_);
    //        Dune::Matrix<Dune::FieldVector<double,1> > Dyy(size,nDt_);
    //        Dxx=0;
    //        Dyy=0;
    //    Dune::Matrix<Dune::FieldVector<double,1> > data(size,3*nDt_);
    //    data=0;
    //
    //    for (int i=0;i<size;i++)
    //    {
    //        for (int j=0;j<nDt_;j++)
    //        {
    //            //                Dxx[i][j] = soil_.getDispersion()[i][j][0][0];
    //            //                Dyy[i][j] = soil_.getDispersion()[i][j][1][1];
    //            data[i][j] = soil_.getDispersion()[i][j][0][0];
    //            data[i][j+nDt_] = soil_.getDispersion()[i][j][1][1];
    //            data[i][j+2*nDt_] = soil_.getDispersionSat()[i][j];
    //
    //            //            std::cout<<soil_.getDispersionSat()[i][j]<<std::endl;
    //            //            std::cout<<soil_.getDispersion()[i][j]<<std::endl;
    //        }
    //    }
    //    //        exportToDGF(grid_.levelView(coarseLev),soil_.getDispersionSat(),nDt_,"saturationdata");
    //    //        exportToDGF(grid_.levelView(coarseLev),Dxx,nDt_,"Dxxdata");
    //    //        exportToDGF(grid_.levelView(coarseLev),Dyy,nDt_,"Dyydata");
    //    exportToDGF(grid_.levelView(coarseLev),data,(3*nDt_),"data");

    return;
}

template<class G, class RT>void UpscalingPreprocess<G,RT>::calcconvectivecorrection(int coarseLev,int fineLev)
{
    if (fineLev<0)
    fineLev=grid_.maxLevel();

    firstRun_ = true;

    tEnd_ *=2;

    const GV& gridView = grid_.levelView(coarseLev);

    const IS& coarseIset = grid_.levelIndexSet(coarseLev);
    int size = coarseIset.size(0);

    Iterator endCIt = gridView.template end<0>();
    for (Iterator cIt = gridView.template begin<0>(); cIt != endCIt; ++cIt)
    {
        int coarseIndex = coarseIset.index(*cIt);
        //            std::cout<<"coarse cell "<<coarseIndex<<std::endl;

        IntersectionIterator endCIIt = gridView.template iend(*cIt);
        for (IntersectionIterator cIIt = gridView.ibegin(*cIt); cIIt != endCIIt; ++cIIt)
        {
            if ( cIIt->neighbor() )
            {
                GeometryType gtCoarse = cIIt->intersectionSelfLocal().type();
                FieldVector<ct,dim-1> faceLocalCoarse = ReferenceElements<ct,dim-1>::general(gtCoarse).position(0,0);
                FieldVector<ct,dim> faceGlobalCoarse = cIIt->intersectionGlobal().global(faceLocalCoarse);

                // local number of facet
                int numberInSelf = cIIt->numberInSelf();
                EntityPointer outside = cIIt->outside();

                Dune::SubGrid<dim,G> subGrid(grid_);
                subGrid.createBegin();

                HierarchicIterator endIt = cIt-> hend(fineLev);
                for (HierarchicIterator it = cIt->hbegin(fineLev); it != endIt; ++it)
                {
                    if (it->level() == fineLev)
                    {
                        subGrid.addPartial(it);
                    }
                }

                HierarchicIterator endItNB = outside-> hend(fineLev);
                for (HierarchicIterator itNB = outside->hbegin(fineLev); itNB != endItNB; ++itNB)
                {
                    if (itNB->level() == fineLev)
                    {
                        subGrid.addPartial(itNB);
                    }
                }
                subGrid.createEnd();

                FieldVector<ct,dim> lowerLeftHelp = cIt->geometry().corner(0);
                FieldVector<ct,dim> upperRightHelp = cIt->geometry().corner(3);

                if (faceGlobalCoarse[0] == lowerLeftHelp[0])
                {
                    //                    std::cout<<"faceGlobalCoarse[0] ="<<faceGlobalCoarse[0]<<"lowerLeftHelp[0] = "<<lowerLeftHelp[0]<<std::endl;
                    FieldVector<ct, dim> lowerLeft = outside->geometry().corner(0);
                    FieldVector<ct,dim> upperRight = cIt->geometry().corner(3);

//                    std::cout<<"lowerLeftx1 = "<<lowerLeft<<"upperRight = "<<upperRight<<std::endl;

                    XProblem1(subGrid, fineLev, lowerLeft, upperRight,coarseIndex, numberInSelf, faceGlobalCoarse,size);
                }
                if (faceGlobalCoarse[0] == upperRightHelp[0])
                {
                    //                    std::cout<<"faceGlobalCoarse[0] ="<<faceGlobalCoarse[0]<<"upperRightHelp[0] = "<<upperRightHelp[0]<<std::endl;
                    FieldVector<ct, dim> lowerLeft = lowerLeftHelp;
                    FieldVector<ct,dim> upperRight = outside->geometry().corner(3);
//                    std::cout<<"lowerLeftx2 = "<<lowerLeft<<"upperRight = "<<upperRight<<std::endl;

                    if(firstRun_)
                    {
                        TestProblem(subGrid, fineLev, lowerLeft, upperRight,coarseIndex, numberInSelf, faceGlobalCoarse, size);
                    }
                    XProblem2(subGrid, fineLev, lowerLeft, upperRight,coarseIndex, numberInSelf,faceGlobalCoarse, size);
                }
                if (faceGlobalCoarse[1] == lowerLeftHelp[1])
                {
                    //                    std::cout<<"faceGlobalCoarse[1] ="<<faceGlobalCoarse[1]<<"lowerLeftHelp[1] = "<<lowerLeftHelp[1]<<std::endl;
                    FieldVector<ct, dim> lowerLeft = outside->geometry().corner(0);
                    FieldVector<ct,dim> upperRight = upperRightHelp;
//                    std::cout<<"lowerLefty1 = "<<lowerLeft<<"upperRight = "<<upperRight<<std::endl;

                    YProblem1(subGrid, fineLev, lowerLeft, upperRight, coarseIndex, numberInSelf,faceGlobalCoarse, size);
                }
                if (faceGlobalCoarse[1] == upperRightHelp[1])
                {
                    //                    std::cout<<"faceGlobalCoarse[1] ="<<faceGlobalCoarse[1]<<"upperRightHelp[1] = "<<upperRightHelp[1]<<std::endl;
                    FieldVector<ct, dim> lowerLeft = lowerLeftHelp;
                    FieldVector<ct,dim> upperRight = outside->geometry().corner(3);
//                    std::cout<<"lowerLefty2 = "<<lowerLeft<<"upperRight = "<<upperRight<<std::endl;

                    YProblem2(subGrid, fineLev, lowerLeft, upperRight, coarseIndex, numberInSelf,faceGlobalCoarse, size);
                }
            }
        }
    }
    //    Dune::Matrix<Dune::FieldVector<double,1> > data(size,8*nDt_);
    //    data=0;
    //
    //    for (int i=0;i<size;i++)
    //    {
    //        for (int j=0;j<nDt_;j++)
    //        {
    //            //            data[i][j] = soil_.getDispersion()[i][j][0][0];
    //            //            data[i][j+nDt_] = soil_.getDispersion()[i][j][1][1];
    //            //            data[i][j+2*nDt_] = soil_.getDispersionSat()[i][j];
    //            data[i][j] = soil_.getM()[i][j][0];
    //            data[i][j+nDt_] = soil_.getM()[i][j][1];
    //            data[i][j+2*nDt_] = soil_.getM()[i][j][2];
    //            data[i][j+3*nDt_] = soil_.getM()[i][j][3];
    //            data[i][j+4*nDt_] = soil_.getMSat()[i][j][0];
    //            data[i][j+5*nDt_] = soil_.getMSat()[i][j][1];
    //            data[i][j+6*nDt_] = soil_.getMSat()[i][j][2];
    //            data[i][j+7*nDt_] = soil_.getMSat()[i][j][3];
    //            std::cout<<soil_.getMSat()[i][j]<<std::endl;
    //            std::cout<<soil_.getM()[i][j]<<std::endl;
    //        }
    //    }
    //
    //    exportToDGF(gridView,data,8*nDt_,"convectiondata");

    return;
}

template<class G, class RT>void UpscalingPreprocess<G, RT>::TestProblem(Dune::SubGrid<dim,G>& subGrid,int& fineLev,FieldVector<ct,dim>& lowerLeft,FieldVector<ct,dim>& upperRight,int& coarseIndex, int& numberInSelf, FieldVector<ct, dim>& faceGlobalCoarse, int& size)
{
    typedef Dune::SubGrid<dim,G> SG;
    typedef Dune::VariableClassSubProbs<SG, RT> VC;
    typedef Dune::ConvSubProblemX2<SG, RT, VC> ConvSubProblem;
    typedef Dune::FVDiffSubProbs<SG, RT, VC> Diffusion;
    typedef Dune::FVTransSubProbs<SG, RT,VC> Transport;
    typedef Dune::IMPESSubProbs<SG, Diffusion, Transport,VC> IMPES;

    VC subProblemVariablesTest(subGrid, fineLev);

    ConvSubProblem subProblemTest(subProblemVariablesTest, wettingPhase_, nonWettingPhase_, soil_, materialLaw_,lowerLeft,upperRight);

    Diffusion subDiffusionTest(subGrid,subProblemTest,fineLev);
    Transport subTransportTest(subGrid,subProblemTest,fineLev);
    IMPES subImpesTest(subDiffusionTest,subTransportTest);
    TimeLoopSubProbs<SG, IMPES> subTimeloopTest(tEnd_);
    subTimeloopTest.executeM(subImpesTest, 1, coarseIndex, numberInSelf, faceGlobalCoarse, lowerLeft, upperRight, nDt_, tEnd_ , firstRun_);

    soil_.getM().setSize(size,nDt_);
    soil_.getMSat().setSize(size,nDt_);

    soil_.getM()=0;
    soil_.getMSat()=0;

    firstRun_=false;
    std::cout<<"-------------------- end testsubproblem -------------------- "<<std::endl;
}

template<class G, class RT>void UpscalingPreprocess<G, RT>::XProblem1(Dune::SubGrid<dim,G>& subGrid,int& fineLev,FieldVector<ct,dim>& lowerLeft,FieldVector<ct,dim>& upperRight,int& coarseIndex, int& numberInSelf, FieldVector<ct,dim>& faceGlobalCoarse, int& size)
{
    typedef Dune::SubGrid<dim,G> SG;
    typedef Dune::VariableClassSubProbs<SG, RT> VC;
    typedef Dune::ConvSubProblemX1<SG, RT, VC> ConvSubProblem;
    typedef Dune::FVDiffSubProbs<SG, RT, VC> Diffusion;
    typedef Dune::FVTransSubProbs<SG, RT,VC> Transport;
    typedef Dune::IMPESSubProbs<SG, Diffusion, Transport,VC> IMPES;

    VC subProblemVariables(subGrid, fineLev);
    ConvSubProblem subProblem(subProblemVariables, wettingPhase_, nonWettingPhase_, soil_, materialLaw_, lowerLeft, upperRight);

    Diffusion subDiffusion(subGrid,subProblem,fineLev);
    Transport subTransport(subGrid,subProblem,fineLev);
    IMPES subImpes(subDiffusion,subTransport);
    TimeLoopSubProbs<SG, IMPES> subTimeloop(tEnd_);
    subTimeloop.executeM(subImpes, 1, coarseIndex, numberInSelf, faceGlobalCoarse, lowerLeft, upperRight, nDt_);

    std::cout<<"---------- end subproblem x1, coarse cell "<<(coarseIndex+1)<<" interface "<< numberInSelf <<" ---------- "<<std::endl;
}

template<class G, class RT>void UpscalingPreprocess<G, RT>::XProblem2(Dune::SubGrid<dim,G>& subGrid,int& fineLev,FieldVector<ct,dim>& lowerLeft,FieldVector<ct,dim>& upperRight,int& coarseIndex, int& numberInSelf, FieldVector<ct, dim>& faceGlobalCoarse, int& size)
{
    typedef Dune::SubGrid<dim,G> SG;
    typedef Dune::VariableClassSubProbs<SG, RT> VC;
    typedef Dune::ConvSubProblemX2<SG, RT, VC> ConvSubProblem;
    typedef Dune::FVDiffSubProbs<SG, RT, VC> Diffusion;
    typedef Dune::FVTransSubProbs<SG, RT,VC> Transport;
    typedef Dune::IMPESSubProbs<SG, Diffusion, Transport,VC> IMPES;

    VC subProblemVariables(subGrid, fineLev);
    ConvSubProblem subProblem(subProblemVariables, wettingPhase_, nonWettingPhase_, soil_, materialLaw_, lowerLeft, upperRight);

    Diffusion subDiffusion(subGrid,subProblem,fineLev);
    Transport subTransport(subGrid,subProblem,fineLev);
    IMPES subImpes(subDiffusion,subTransport);
    TimeLoopSubProbs<SG, IMPES> subTimeloop(tEnd_);
    subTimeloop.executeM(subImpes, 2, coarseIndex, numberInSelf, faceGlobalCoarse,lowerLeft, upperRight, nDt_);

    std::cout<<"---------- end subproblem x2, coarse cell "<<(coarseIndex+1)<<" interface "<< numberInSelf <<" ---------- "<<std::endl;
}

template<class G, class RT>void UpscalingPreprocess<G, RT>::YProblem1(Dune::SubGrid<dim,G>& subGrid,int& fineLev,FieldVector<ct,dim>& lowerLeft,FieldVector<ct,dim>& upperRight,int& coarseIndex, int& numberInSelf, FieldVector<ct, dim>& faceGlobalCoarse, int& size)
{
    typedef Dune::SubGrid<dim,G> SG;
    typedef Dune::VariableClassSubProbs<SG, RT> VC;
    typedef Dune::ConvSubProblemY1<SG, RT, VC> ConvSubProblem;
    typedef Dune::FVDiffSubProbs<SG, RT, VC> Diffusion;
    typedef Dune::FVTransSubProbs<SG, RT,VC> Transport;
    typedef Dune::IMPESSubProbs<SG, Diffusion, Transport,VC> IMPES;

    VC subProblemVariables(subGrid, fineLev);
    ConvSubProblem subProblem(subProblemVariables, wettingPhase_, nonWettingPhase_, soil_, materialLaw_, lowerLeft, upperRight);

    Diffusion subDiffusion(subGrid,subProblem,fineLev);
    Transport subTransport(subGrid,subProblem,fineLev);
    IMPES subImpes(subDiffusion,subTransport);
    TimeLoopSubProbs<SG, IMPES> subTimeloop(tEnd_);

    subTimeloop.executeM(subImpes, 3, coarseIndex, numberInSelf,faceGlobalCoarse, lowerLeft, upperRight, nDt_);

    std::cout<<"---------- end subproblem y1, coarse cell "<<(coarseIndex+1)<<" interface "<< numberInSelf <<" ---------- "<<std::endl;
}
template<class G, class RT>void UpscalingPreprocess<G, RT>::YProblem2(Dune::SubGrid<dim,G>& subGrid,int& fineLev,FieldVector<ct,dim>& lowerLeft,FieldVector<ct,dim>& upperRight,int& coarseIndex, int& numberInSelf, FieldVector<ct, dim>& faceGlobalCoarse, int& size)
{
    typedef Dune::SubGrid<dim,G> SG;
    typedef Dune::VariableClassSubProbs<SG, RT> VC;
    typedef Dune::ConvSubProblemY2<SG, RT, VC> ConvSubProblem;
    typedef Dune::FVDiffSubProbs<SG, RT, VC> Diffusion;
    typedef Dune::FVTransSubProbs<SG, RT,VC> Transport;
    typedef Dune::IMPESSubProbs<SG, Diffusion, Transport,VC> IMPES;

    VC subProblemVariables(subGrid, fineLev);
    ConvSubProblem subProblem(subProblemVariables, wettingPhase_, nonWettingPhase_, soil_, materialLaw_, lowerLeft, upperRight);

    Diffusion subDiffusion(subGrid,subProblem,fineLev);
    Transport subTransport(subGrid,subProblem,fineLev);
    IMPES subImpes(subDiffusion,subTransport);
    TimeLoopSubProbs<SG, IMPES> subTimeloop(tEnd_);

    subTimeloop.executeM(subImpes, 4, coarseIndex, numberInSelf,faceGlobalCoarse,lowerLeft, upperRight, nDt_);

    std::cout<<"---------- end subproblem y2, coarse cell "<<(coarseIndex+1)<<" interface "<< numberInSelf <<" ---------- "<<std::endl;
}

template<class G, class RT>void UpscalingPreprocess<G,RT>::exportdata(int coarseLev)
{
    int size = soil_.getDispersion().N();
    int dispersionSize = soil_.getDispersion().M();

    Dune::Matrix<Dune::FieldVector<double,1> > data(size,11*nDt_);
    data=0;

    for (int i=0;i<size;i++)
    {
        for (int j=0;j<nDt_;j++)
        {
            if (j<dispersionSize)
            {
                data[i][j] = soil_.getDispersion()[i][j][0][0];
                data[i][j+nDt_] = soil_.getDispersion()[i][j][1][1];
                data[i][j+2*nDt_] = soil_.getDispersionSat()[i][j];
            }
            data[i][j+3*nDt_] = soil_.getM()[i][j][0];
            data[i][j+4*nDt_] = soil_.getM()[i][j][1];
            data[i][j+5*nDt_] = soil_.getM()[i][j][2];
            data[i][j+6*nDt_] = soil_.getM()[i][j][3];
            data[i][j+7*nDt_] = soil_.getMSat()[i][j][0];
            data[i][j+8*nDt_] = soil_.getMSat()[i][j][1];
            data[i][j+9*nDt_] = soil_.getMSat()[i][j][2];
            data[i][j+10*nDt_] = soil_.getMSat()[i][j][3];
            //            std::cout<<soil_.getMSat()[i][j]<<std::endl;
            //            std::cout<<soil_.getM()[i][j]<<std::endl;
        }
    }

    exportToDGF(grid_.levelView(coarseLev),data,11*nDt_,"data");
}
}
#endif
