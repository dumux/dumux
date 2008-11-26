#ifndef DUNE_VARIABLECLASS_HH
#define DUNE_VARIABLECLASS_HH

/**
 * @file
 * @brief  class including the variables
 * @author Markus Wolff
 */

namespace Dune
{

template<class G, class RT> class VariableClass
{
	template<int dim> struct ElementLayout
	{
		bool contains(Dune::GeometryType gt)
		{
			return gt.dim() == dim;
		}
	};

	enum{n=G::dimension};
	typedef typename G::ctype DT;
	typedef Dune::BlockVector< Dune::FieldVector<RT,1> > ScalarType;
	typedef Dune::BlockVector< Dune::FieldVector<RT,n> > SlopeType;
	typedef Dune::BlockVector< Dune::FieldVector<Dune::FieldVector<double, n>, 2*n> >
			VelType;

	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename G::Traits::LevelIndexSet IS;
	typedef Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout>
			ElementMapper;

public:
	G& grid;
	ElementMapper diffmapper;
	ElementMapper transmapper;
	int diffsize;
	int transsize;
	int difflevel;
	int translevel;


	ScalarType saturation;
	ScalarType pressure;
	VelType velocity;
	SlopeType slope;

	VariableClass(G& g, RT& initialsat = *(new RT(0)), RT& initalpress = *(new RT(0)), Dune::FieldVector<RT, n>& initialvel = *(new Dune::FieldVector<RT, n> (0)), int translev = -1, int difflev = -1)
	: grid(g), diffmapper(g, (difflev >= 0) ? g.levelIndexSet(difflev) : g.levelIndexSet(g.maxLevel())), transmapper(g, (translev >= 0) ? g.levelIndexSet(translev) : g.levelIndexSet(g.maxLevel())),
	  diffsize(diffmapper.size()),transsize(transmapper.size()),
	  difflevel((difflev >= 0) ? difflev : g.maxLevel()),translevel((translev >= 0) ? translev : g.maxLevel())
	{
		initsat(initialsat, transsize);
		initpress(initalpress, diffsize);
		initvel(initialvel, transsize);
		initslopes(0.0, transsize);
	}

	void initsat(RT& initialsat, int size)
	{
		saturation.resize(size);
		saturation=initialsat;
		return;
	}

	void initslopes(RT initialslope, int size)
		{
			slope.resize(size);
			slope=initialslope;
			return;
		}

	void initpress(RT& initialpress, int size)
	{
		pressure.resize(size);
		pressure=initialpress;
		return;
	}
	void initvel(Dune::FieldVector<RT, n>& initialvel, int size)
	{
		velocity.resize(size);
		velocity=initialvel;
		return;
	}

	ScalarType& sat() const
	{
		return saturation;
	}
	ScalarType& press() const
	{
		return pressure;
	}
	VelType& vel() const
	{
		return velocity;
	}

	const Dune::FieldVector<RT,1>& sat(const Dune::FieldVector<DT,n>& x,
			const Entity& e, const Dune::FieldVector<DT,n>& xi) const
	{
		return saturation[transmapper.map(e)];;
	}

	const Dune::FieldVector<RT,1>& press(const Dune::FieldVector<DT,n>& x,
			const Entity& e, const Dune::FieldVector<DT,n>& xi) const
	{
		return pressure[diffmapper.map(e)];
	}

	const Dune::FieldVector<DT,n>& vTotal(const Entity& e,
			const int numberInSelf) const
	{
		int elemId = transmapper.map(e);

		return (velocity[elemId][numberInSelf]);
	}

	/*! @brief prints the saturation to a VTK file
	 *
	 *  The file name is "<name>-<k>.vtu" where k is an integer number.
	 *  @param name specifies the name of the VTK file
	 *  @param k specifies a number
	 */

	void vtkout (const char* name, int k)
	{
			vtkoutdifflevel(name,k,difflevel,translevel);
	}

	void vtkoutdifflevel(const char* name, int k, int pressurelevel = 0, int satlevel = 0) const
	{
		if (pressurelevel == satlevel)
		{
			VTKWriter<typename G::LeafGridView> vtkwriter(grid.leafView());
			char fname[128];
			sprintf(fname, "%s-%05d", name, k);
			vtkwriter.addCellData(saturation, "saturation");
			vtkwriter.addCellData(pressure, "total pressure p~");
			vtkwriter.write(fname, VTKOptions::ascii);
		}
		else
		{
			Dune::VTKWriter<typename G::LevelGridView>
					vtkwriterpressure(grid.levelView(pressurelevel));
			char fname[128];
			sprintf(fname, "%s-press%05d", name, k);
			vtkwriterpressure.addCellData(pressure, "total pressure p~");
			vtkwriterpressure.write(fname, Dune::VTKOptions::ascii);

			Dune::VTKWriter<typename G::LevelGridView>
					vtkwritersaturation(grid.levelView(satlevel));
			sprintf(fname, "%s-%05d", name, k);
			vtkwritersaturation.addCellData(saturation, "saturation");
			vtkwritersaturation.write(fname, VTKOptions::ascii);
		}
		return;
	}
//	void vtkoutpressure(const char* name, int k) const
//	{
//		LeafVTKWriter<typename G::LeafGridView> vtkwriter(grid);
//		char fname[128];
//		sprintf(fname, "%s-press%05d", name, k);
//		vtkwriter.addCellData(pressure, "total pressure p~");
//		vtkwriter.write(fname, VTKOptions::ascii);
//	}
//	void vtkoutsamelevel(const char* name, int k, BlockVector<FieldVector<RT, 2> > uEx) const
//	{
//		ScalarType saturationEx(transsize);
//		ScalarType error(transsize);
//
//		for (int i=0; i<transsize; i++)
//		{
//			saturationEx[i] = uEx[i][0];
//			error[i]=uEx[i][1];
//		}
//		VTKWriter<typename G::LevelGridView> vtkwriter(grid);
//		char fname[128];
//		sprintf(fname, "%s-%05d", name, k);
//		vtkwriter.addCellData(saturation, "saturation");
//		vtkwriter.addCellData(pressure, "total pressure p~");
//		vtkwriter.addCellData(saturationEx, "saturation (exact solution)");
//		vtkwriter.addCellData(error, "error");
//		vtkwriter.write(fname, VTKOptions::ascii);
//	}
};
}
#endif
