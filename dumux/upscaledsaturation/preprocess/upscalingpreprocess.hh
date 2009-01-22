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
#include "dumux/io/exportcorrection.hh"
//#include <dune/istl/matrix.hh>
//#include <dune/common/fvector.hh>

namespace Dune {

template<class Grid, class Scalar> class UpscalingPreprocess {
	enum {
		dim = Grid::dimension, dimWorld = Grid::dimensionworld
	};

typedef	typename Grid::LevelGridView GridView;
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

	void XProblem1(Dune::SubGrid<dim,Grid>&, const int&, const GlobalPosition&, const GlobalPosition&, const int&, bool, const int, const GlobalPosition&, bool);
	void YProblem1(Dune::SubGrid<dim,Grid>&, const int&, const GlobalPosition&, const GlobalPosition&, const int&, bool, const int, const GlobalPosition&, bool);
	void XProblem2(Dune::SubGrid<dim,Grid>&, const int&, const GlobalPosition&, const GlobalPosition&, const int&, bool, const int, const GlobalPosition&, bool);
	void YProblem2(Dune::SubGrid<dim,Grid>&, const int&, const GlobalPosition&, const GlobalPosition&, const int&, bool, const int, const GlobalPosition&, bool);

public:
	UpscalingPreprocess(Grid& grid, Fluid& wettingPhase , Fluid& nonWettingPhase, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>), Scalar tEnd = 1e6)
	:grid_(grid), wettingPhase_(wettingPhase), nonWettingPhase_(nonWettingPhase),soil_(soil),materialLaw_(materialLaw)
	{}

	void preprocessexecute(int coarseLev = 0,int fineLev=-1)
	{
		calcDispersiveFluxCorrection(coarseLev,fineLev);
		calcConvectiveFluxCorrection(coarseLev,fineLev);
		exportData(coarseLev);
	}

private:
	Grid& grid_;
	Fluid& wettingPhase_;
	Fluid& nonWettingPhase_;
	Matrix2p<Grid, Scalar>& soil_;
	TwoPhaseRelations<Grid, Scalar>& materialLaw_;
	bool firstRun_;
};

template<class Grid, class Scalar>void UpscalingPreprocess<Grid,Scalar>::calcDispersiveFluxCorrection(int coarseLev,int fineLev)
{
	bool calcDispersion = true;

	if (fineLev<0)
	fineLev=grid_.maxLevel();

	const GridView& gridView = grid_.levelView(coarseLev);

	const IndexSet& indexSetCoarse = grid_.levelIndexSet(coarseLev);
	int size = indexSetCoarse.size(0);

	soil_.getDispersion().resize(size);
	soil_.getDispersionSat().resize(size);

	ElementIterator eItCoarseEnd = gridView.template end<0>();
	for (ElementIterator eItCoarse = gridView.template begin<0>(); eItCoarse != eItCoarseEnd; ++eItCoarse)
	{
		int globalIdxCoarse = indexSetCoarse.index(*eItCoarse);
		//            std::cout<<"coarse cell "<<globalIdxCoarse<<std::endl;

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
	return;
}

template<class Grid, class Scalar>void UpscalingPreprocess<Grid,Scalar>::calcConvectiveFluxCorrection(int coarseLev,int fineLev)
{
	bool calcConvection = false;

	if (fineLev<0)
	fineLev=grid_.maxLevel();

	const GridView& gridView = grid_.levelView(coarseLev);

	const IndexSet& indexSetCoarse = grid_.levelIndexSet(coarseLev);
	int size = indexSetCoarse.size(0);

	soil_.getM().resize(size);
	soil_.getMSat().resize(size);

	ElementIterator eItCoarseEnd = gridView.template end<0>();
	for (ElementIterator eItCoarse = gridView.template begin<0>(); eItCoarse != eItCoarseEnd; ++eItCoarse)
	{
		int globalIdxCoarse = indexSetCoarse.index(*eItCoarse);
		//			std::cout<<"coarse cell "<<globalIdxCoarse<<std::endl;

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

				if (globalPosFaceCoarse[0] == lowerLeftHelp[0])
				{
					//					std::cout<<"globalPosFaceCoarse[0] ="<<globalPosFaceCoarse[0]<<"lowerLeftHelp[0] = "<<lowerLeftHelp[0]<<std::endl;
					GlobalPosition lowerLeft = neighborPointer->geometry().corner(0);
					GlobalPosition upperRight = upperRightHelp;

					//					std::cout<<"lowerLeftx1 = "<<lowerLeft<<"upperRight = "<<upperRight<<std::endl;

					XProblem1(subGrid, fineLev, lowerLeft, upperRight,globalIdxCoarse,calcConvection, numberInSelf, globalPosFaceCoarse);
				}
				if (globalPosFaceCoarse[0] == upperRightHelp[0])
				{
					//					std::cout<<"globalPosFaceCoarse[0] ="<<globalPosFaceCoarse[0]<<"upperRightHelp[0] = "<<upperRightHelp[0]<<std::endl;
					GlobalPosition lowerLeft = lowerLeftHelp;
					GlobalPosition upperRight = neighborPointer->geometry().corner(3);
					//					std::cout<<"lowerLeftx2 = "<<lowerLeft<<"upperRight = "<<upperRight<<std::endl;

					XProblem2(subGrid, fineLev, lowerLeft, upperRight,globalIdxCoarse,calcConvection, numberInSelf,globalPosFaceCoarse);
				}
				if (globalPosFaceCoarse[1] == lowerLeftHelp[1])
				{
					//					std::cout<<"globalPosFaceCoarse[1] ="<<globalPosFaceCoarse[1]<<"lowerLeftHelp[1] = "<<lowerLeftHelp[1]<<std::endl;
					GlobalPosition lowerLeft = neighborPointer->geometry().corner(0);
					GlobalPosition upperRight = upperRightHelp;
					//					std::cout<<"lowerLefty1 = "<<lowerLeft<<"upperRight = "<<upperRight<<std::endl;

					YProblem1(subGrid, fineLev, lowerLeft, upperRight, globalIdxCoarse,calcConvection, numberInSelf,globalPosFaceCoarse);
				}
				if (globalPosFaceCoarse[1] == upperRightHelp[1])
				{
					//					std::cout<<"globalPosFaceCoarse[1] ="<<globalPosFaceCoarse[1]<<"upperRightHelp[1] = "<<upperRightHelp[1]<<std::endl;
					GlobalPosition lowerLeft = lowerLeftHelp;
					GlobalPosition upperRight = neighborPointer->geometry().corner(3);
					//					std::cout<<"lowerLefty2 = "<<lowerLeft<<"upperRight = "<<upperRight<<std::endl;

					YProblem2(subGrid, fineLev, lowerLeft, upperRight, globalIdxCoarse,calcConvection, numberInSelf,globalPosFaceCoarse);
				}
			}
			else
			{
				if (isItCoarse->boundary())
				{
					bool isBoundary = true;

					GeometryType gtCoarse = isItCoarse->intersectionSelfLocal().type();
					FieldVector<Scalar,dim-1> faceLocalCoarse = ReferenceElements<Scalar,dim-1>::general(gtCoarse).position(0,0);
					const GlobalPosition& globalPosFaceCoarse = isItCoarse->intersectionGlobal().global(faceLocalCoarse);

					// local number of facet
					int numberInSelf = isItCoarse->numberInSelf();

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
					subGrid.createEnd();

					GlobalPosition lowerLeft = eItCoarse->geometry().corner(0);
					GlobalPosition upperRight = eItCoarse->geometry().corner(3);

					if (globalPosFaceCoarse[0] == lowerLeft[0])
					{
						XProblem1(subGrid, fineLev, lowerLeft, upperRight,globalIdxCoarse,calcConvection, numberInSelf, globalPosFaceCoarse,isBoundary);
					}
					if (globalPosFaceCoarse[0] == upperRight[0])
					{
						XProblem2(subGrid, fineLev, lowerLeft, upperRight,globalIdxCoarse,calcConvection, numberInSelf,globalPosFaceCoarse, isBoundary);
					}
					if (globalPosFaceCoarse[1] == lowerLeft[1])
					{
						YProblem1(subGrid, fineLev, lowerLeft, upperRight, globalIdxCoarse,calcConvection, numberInSelf,globalPosFaceCoarse,isBoundary);
					}
					if (globalPosFaceCoarse[1] == upperRight[1])
					{
						YProblem2(subGrid, fineLev, lowerLeft, upperRight, globalIdxCoarse,calcConvection, numberInSelf,globalPosFaceCoarse,isBoundary);
					}
				}
			}
		}
	}

	return;
}

template<class Grid, class Scalar>void UpscalingPreprocess<Grid, Scalar>::XProblem1(Dune::SubGrid<dim,Grid>& subGrid,const int& fineLev,const GlobalPosition& lowerLeft,const GlobalPosition& upperRight,const int& globalIdxCoarse, bool correctionType, const int numberInSelf = 0, const GlobalPosition& globalPosFaceCoarse = *(new GlobalPosition(0)), bool isBoundadry = false)
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
	TimeLoopSubProbs<SubGrid, IMPES> subTimeloop(isBoundadry);
	subTimeloop.execute(subImpes, 1, globalIdxCoarse, lowerLeft, upperRight, correctionType, numberInSelf, globalPosFaceCoarse);

	std::cout<<"---------- end subproblem x1, coarse cell "<<(globalIdxCoarse+1)<<" interface "<< numberInSelf <<" ---------- "<<std::endl;
}

template<class Grid, class Scalar>void UpscalingPreprocess<Grid, Scalar>::XProblem2(Dune::SubGrid<dim,Grid>& subGrid,const int& fineLev,const GlobalPosition& lowerLeft,const GlobalPosition& upperRight,const int& globalIdxCoarse, bool correctionType, const int numberInSelf = 0, const GlobalPosition& globalPosFaceCoarse= *(new GlobalPosition(0)), bool isBoundadry = false)
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
	TimeLoopSubProbs<SubGrid, IMPES> subTimeloop(isBoundadry);
	subTimeloop.execute(subImpes, 2, globalIdxCoarse,lowerLeft, upperRight, correctionType, numberInSelf, globalPosFaceCoarse);

	std::cout<<"---------- end subproblem x2, coarse cell "<<(globalIdxCoarse+1)<<" interface "<< numberInSelf <<" ---------- "<<std::endl;
}

template<class Grid, class Scalar>void UpscalingPreprocess<Grid, Scalar>::YProblem1(Dune::SubGrid<dim,Grid>& subGrid,const int& fineLev,const GlobalPosition& lowerLeft,const GlobalPosition& upperRight,const int& globalIdxCoarse, bool correctionType, const int numberInSelf = 0, const GlobalPosition& globalPosFaceCoarse= *(new GlobalPosition(0)), bool isBoundadry = false)
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
	TimeLoopSubProbs<SubGrid, IMPES> subTimeloop(isBoundadry);

	subTimeloop.execute(subImpes, 3, globalIdxCoarse, lowerLeft, upperRight, correctionType, numberInSelf,globalPosFaceCoarse);

	std::cout<<"---------- end subproblem y1, coarse cell "<<(globalIdxCoarse+1)<<" interface "<< numberInSelf <<" ---------- "<<std::endl;
}
template<class Grid, class Scalar>void UpscalingPreprocess<Grid, Scalar>::YProblem2(Dune::SubGrid<dim,Grid>& subGrid,const int& fineLev,const GlobalPosition& lowerLeft,const GlobalPosition& upperRight,const int& globalIdxCoarse, bool correctionType, const int numberInSelf = 0, const GlobalPosition& globalPosFaceCoarse= *(new GlobalPosition(0)), bool isBoundadry = false)
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
	TimeLoopSubProbs<SubGrid, IMPES> subTimeloop(isBoundadry);

	subTimeloop.execute(subImpes, 4, globalIdxCoarse,lowerLeft, upperRight, correctionType, numberInSelf,globalPosFaceCoarse);

	std::cout<<"---------- end subproblem y2, coarse cell "<<(globalIdxCoarse+1)<<" interface "<< numberInSelf <<" ---------- "<<std::endl;
}

template<class Grid, class Scalar>void UpscalingPreprocess<Grid,Scalar>::exportData(int coarseLev)
{
	std::vector<std::vector<Dune::FieldVector<Scalar, dim> > > dispersionData = soil_.getDispersion();
	std::vector<std::vector<Dune::FieldVector<Scalar, 1> > > dispersionSatData = soil_.getDispersionSat();

	std::vector<std::vector<Dune::FieldVector<Scalar, 2*dim> > > mData = soil_.getM();
	std::vector<std::vector<Dune::FieldVector<Scalar, 2*dim> > > mSatData = soil_.getMSat();

	exportToFile(dispersionData, "dataDispersion");
	exportToFile(dispersionSatData, "dataDispersionSat");
	exportToFile(mData,"dataM");
	exportToFile(mSatData,"dataMSat");
}
}
#endif
