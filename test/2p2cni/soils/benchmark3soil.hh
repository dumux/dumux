#ifndef BENCHMARK3_SOIL_H
#define BENCHMARK3_SOIL_H

#include <dumux/material/matrixproperties.hh>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>

namespace Dune
{


template<class Grid, class Scalar>
class Benchmark3Soil: public HeterogeneousSoil<Grid, Scalar>
{
public:

	typedef typename Grid::ctype CoordScalar;
    enum {dim=Grid::dimension, dimWorld=Grid::dimensionworld};
    typedef	typename Grid::Traits::template Codim<0>::Entity Element;
    typedef	typename Grid::Traits::template Codim<dim>::Entity Vertex;
    typedef Dune::FieldVector<CoordScalar,dim>      LocalPosition;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
	typedef typename Grid::LeafGridView GV;
	typedef typename GV::template Codim<dim>::Iterator VertexIterator;
	typedef typename GV::IndexSet IndexSet;
    typedef typename GV::IndexSet IS;


	FieldMatrix<CoordScalar,dim,dim> &K(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const int idx=0)
	{
		const GV& gridView(_grid.leafView());
		const IndexSet& indexSet = gridView.indexSet();
		int index = indexSet.index(*((e).template entity<dim>(idx)));

			permloc_ = perm_[index]*1e-15;

			return permloc_;
	}

	double porosity(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const int idx=0)	const
	{
		double phi;
		const GV& gridView(_grid.leafView());
		const IndexSet& indexSet = gridView.indexSet();
		int index = indexSet.index(*((e).template entity<dim>(idx)));

		phi = poro_[index];
		return phi;
	}

	double Sr_w(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T = 283.15) const
	{
		return 0.2;
	}

	double Sr_n(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T = 283.15) const
	{
		return 0.05;
	}

	virtual double heatCap(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const int idx = 0) const
	{
		return 	(800 /* spec. heat cap. of sediment */
						* 2650 /* density of sediment */
						* (1-porosity(x, e, xi, idx)));
	}

	virtual double heatCond(const GlobalPosition &x, const Vertex& v, const LocalPosition &xi, const double sat, const int idx = 0) const
	{
		static const double ldry = 0.32;
		static const double lsat = 2.7;
		double l_wurz, Sw;
		Sw = sat;
		if(Sw<0.0) Sw = 0.0; /* ACHTUNG Regularisierung */
		if(Sw>1.) Sw = 1.; /* ACHTUNG Regularisierung */

		l_wurz = ldry + sqrt(Sw)*(lsat-ldry);

		if(isnan(l_wurz)) {
			std::cout <<"isnan heatcondwurzel \n"<<std::endl;
			l_wurz = 0.0;
		}
		return(l_wurz);
	}

	std::vector<double> paramRelPerm(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T = 283.15) const
	{
		std::vector<double> param(2);

		//Brooks-Corey parameters
		param[0] = 2.; // lambda
		param[1] = 10000.; // entry-pressure

		return param;
	}

	typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
	{
		return Matrix2p<Grid,Scalar>::brooks_corey;
	}

    void readSoilProperties()
	{
		const GV& gridView(_grid.leafView());
		const IndexSet& indexSet = gridView.indexSet();
		int numNodes=indexSet.size(dim);

		soilProps_.resize(numNodes);
		std::ifstream is(properties_file_);
		std::string str;

		for(int j=0; j<numNodes; j++){
		std::getline(is, str);
		std::stringstream strs(str);

		strs >> soilProps_[j][0] >> soilProps_[j][1] >> soilProps_[j][2] >>
		soilProps_[j][3] >> soilProps_[j][4];
		}
		return;
	}

    void setSoilProperties()
	{
		const GV& gridView(_grid.leafView());
		const IndexSet& indexSet = gridView.indexSet();

		int numNodes=indexSet.size(dim);
		perm_.resize(numNodes);
		poro_.resize(numNodes);

		VertexIterator vEndIt = gridView.template end<dim>();
		for (VertexIterator vIt = gridView.template begin<dim>(); vIt != vEndIt; ++vIt)
		{
			int vertexIndex = indexSet.index(*vIt);
			for (int j=0; j<numNodes; j++)
			{
				if (vIt->geometry().corner(0)[0]-0.05 <= soilProps_[j][0] && vIt->geometry().corner(0)[0]+0.05 >= soilProps_[j][0]
				    && vIt->geometry().corner(0)[1]-0.05 <= soilProps_[j][1] && vIt->geometry().corner(0)[1]+0.05 >= soilProps_[j][1]
				    && vIt->geometry().corner(0)[2]-3238.01 <= soilProps_[j][2] && vIt->geometry().corner(0)[2]-3237.914 >= soilProps_[j][2])
				{
					poro_[vertexIndex] = soilProps_[j][3];
					perm_[vertexIndex] = soilProps_[j][4];
				}
			}
		}
	}

	virtual bool readPropertiesFlag()
	{
		return read_;
	}
	Benchmark3Soil(const Grid& g, const bool read, const char* filename)
	:HeterogeneousSoil<Grid,Scalar>(g, read, filename), _grid(g), properties_file_(filename), read_(read)
	{
	  permTest_ = 0;
	  permloc_ = 0;
	  permlocWell_ = 0;
	  for (int i = 0; i < dim; i++)
		permTest_[i][i] = 1.0e-15;
	  for (int i = 0; i < dim; i++)
		permloc_[i][i] = 2.0e-14;
      for (int i = 0; i < dim; i++)
        permlocWell_[i][i] = 1.0e-12;
	}


	private:
	Dune::FieldMatrix<Scalar,dim,dim> permloc_, permlocWell_, permTest_;
	std::vector<double> perm_;
	std::vector<double> poro_;
	Dune::BlockVector<FieldVector<double, 5> > soilProps_;
	const Grid& _grid;
	int numNodes_;
	const char* properties_file_;
	bool read_;
};
} // end namespace
#endif
