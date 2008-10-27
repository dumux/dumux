#ifndef EXPORTTODGF_HH
#define EXPORTTODGF_HH

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <boost/format.hpp>

namespace Dune
{
//template<int dim> struct P1
//{
//	bool contains(Dune::GeometryType gt) {
//		return gt.dim() == 0;
//	}
//};

template<class GridView, class Data>
void exportToDGF(const GridView& gridView, const Data& data, int paramnumber = 1,
		std::string dataFileName = "data", bool cellWise = true)
{
	typedef typename GridView::Grid Grid;
	typedef typename Grid::ctype DT;

	enum
	{	n = GridView::dimension};


//	const Grid& gridHelp(gridView.grid());
	// mapper: one data element per vertex

	typedef typename GridView::template Codim< 0>::Iterator ElementIterator;
	typedef typename GridView::template Codim< n>::Iterator VertexIterator;
	typedef typename GridView::IntersectionIterator IntersectionIterator;

	typedef typename GridView::IndexSet IndexSet;
//	typedef MultipleCodimMultipleGeomTypeMapper<Grid,IndexSet,P1> VertexMapper;
//	VertexMapper vertexmapper(gridHelp,	gridHelp.leafIndexSet());

	const IndexSet& indexSet = gridView.indexSet();
	int size = indexSet.size(n);

	std::string dgfFileName = (boost::format("%s.dgf")%dataFileName).str();
	std::ofstream dataFile;

	dataFile.open(dgfFileName.c_str());
	dataFile <<
	"DGF\n"
	"VERTEX\n";

	if (!cellWise)
	{
		dataFile <<
		"PARAMETERS "<< paramnumber <<"\n";
	}

	Dune::BlockVector<Dune::FieldVector<double, n> > coordinates(size);

	VertexIterator vEndIt = gridView.template end<n>();
	for (VertexIterator vIt = gridView.template begin<n>(); vIt != vEndIt; ++vIt)
	{
//		int vertexIndex = vertexmapper.map(*vIt);
		int vertexIndex = indexSet.index(*vIt);
		coordinates[vertexIndex] = vIt->geometry()[0];
	}

	for (int vertexNum = 0; vertexNum < size; vertexNum++)
	{
		dataFile <<
		coordinates[vertexNum] << " ";
		if (!cellWise)
		{
			for (int i = 0; i<paramnumber;i++)
			{
				dataFile <<
				data[vertexNum][i] <<
				" ";
			}
		}
		dataFile <<
		"\n";
	}
	dataFile <<
	"#\n"
	"CUBE\n";

	if (cellWise)
	{
		dataFile <<
		"PARAMETERS "<< paramnumber <<"\n";
	}

	ElementIterator eEndIt = gridView.template end<0>();
	for (ElementIterator eIt = gridView.template begin<0>(); eIt != eEndIt; ++eIt)
	{
//		int elementIndex = vertexmapper.map(*eIt);
		int elementIndex = indexSet.index(*eIt);

		int vertexOnElement = eIt->geometry().corners();

		for (int i=0; i<vertexOnElement;i++)
		{
//			int vertexIndex = vertexmapper.map(*((*eIt).template entity<n>(i) ));
			int vertexIndex = indexSet.index(*((*eIt).template entity<n>(i) ));
			dataFile <<
			vertexIndex << " ";

		}
		if (cellWise)
		{
			for (int i = 0; i<paramnumber;i++)
			{
				dataFile <<
				data[elementIndex][i] <<
				" ";
			}
		}
		dataFile <<
		"\n";
	}
//	dataFile <<
//	"#\n"
//	"BOUNDARYSEGMENTS\n";
//
//	for (ElementIterator eIt = gridView.template begin<0>(); eIt != eEndIt; ++eIt)
//	{
//		IntersectionIterator endIIt = gridView.iend(*eIt);
//		for (IntersectionIterator iIt = gridView.ibegin(*eIt); iIt != endIIt; ++iIt)
//		{
//			if ((*iIt).boundary() && !(*iIt).neighbor())
//			{
//
//				dataFile <<
//				"2 ";
//
//				int vertexOnElement = (*eIt).geometry().corners();
//				int vertexOnIntersection = (*iIt).intersectionGlobal().corners();
//
//				for (int i=0; i<vertexOnElement;i++)
//				{
//					Dune::FieldVector<DT,n> globalE = eIt->geometry()[i];
//					for(int j = 0;j<vertexOnIntersection;j++)
//					{
//						Dune::FieldVector<DT,n> globalV = iIt->intersectionGlobal()[j];
//						if (globalE == globalV)
//						{
//							int vertexIndex = vertexmapper.index(*((*eIt).template entity<n>(i) ));
//
//							dataFile <<
//							vertexIndex << " ";
//						}
//					}
//				}
//				dataFile <<
//				"\n";
//			}
//		}
//	}

	dataFile <<
	"#\n"
	"BOUNDARYDOMAIN\n"
	"default 1\n"
	"#\n"
	"# " << dgfFileName;
	dataFile.close();
	return;
}
}
#endif
