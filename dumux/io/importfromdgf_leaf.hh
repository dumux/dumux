// $Id$
#ifndef IMPORTPARAMETERSFROMDGF_HH
#define IMPORTPARAMETERSFROMDGF_HH

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <boost/format.hpp>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>

namespace Dune
{

/** \todo Please doc me! */

template<class GV, class Data>
void importFromDGF(Data& data, std::string dataFileName, bool cellWise = true)
{
#ifndef DUMUX_NO_RESTART
    enum{n = GV::dimension};

//    typedef Dune::UGGrid<n> GInternal;
    typedef Dune::UGGrid<n> GInternal;
    typedef typename GInternal::LeafGridView GridView;
    typedef    typename GInternal::ctype DT;


    typedef typename GridView::template Codim< 0>::Iterator ElementIterator;
    typedef typename GridView::template Codim< n>::Iterator VertexIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef typename GridView::IndexSet IndexSet;

    std::string dgfFileName = (boost::format("%s.dgf")%dataFileName).str();

    Dune::GridPtr<GInternal> getData(dgfFileName);

    GInternal& dataGrid = *getData;

    const GridView& gridView(dataGrid.leafView());


    const IndexSet& indexSet = gridView.indexSet();

    int paramNum = 1;
    int elemNum = 1;
    if (cellWise)
    {
        ElementIterator eEndIt = gridView.template end<0>();
        for (ElementIterator eIt = gridView.template begin<0>(); eIt != eEndIt; ++eIt)
        {
            int elementIndex = indexSet.index(*eIt);
            std::vector<double> parameters = getData.parameters(*eIt);

            if (eIt == gridView.template begin<0>())
            {
                elemNum = indexSet.size(0);
                paramNum = parameters.size();
            }

            for (int i = 0; i<paramNum;i++)
            {
                data[elementIndex][i]=parameters[i];
            }
        }
    }
    if (!cellWise)
    {
        VertexIterator vEndIt = gridView.template end<n>();
        for (VertexIterator vIt = gridView.template begin<n>(); vIt != vEndIt; ++vIt)
        {
            std::vector<double> parameters = getData.parameters(*vIt);
            int vertexIndex = indexSet.index(*vIt);
            if (vIt == gridView.template begin<n>())
            {
                elemNum = indexSet.size(n);
                paramNum = parameters.size();
            }

            for (int i = 0; i<paramNum;i++)
            {
                data[vertexIndex][i]=parameters[i];
            }
        }
    }
    return;
#endif // DUMUX_NO_RESTART
}
}
#endif
