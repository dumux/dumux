#ifndef DUNE_VARIABLECLASS_HH
#define DUNE_VARIABLECLASS_HH

/**
 * @file
 * @brief  class including the variables
 * @author Markus Wolff
 */

namespace Dune {

template<class G, class RT> class VariableClass {

	enum {n=G::dimension};
	typedef typename G::ctype DT;
	typedef Dune::BlockVector< Dune::FieldVector<RT,1> > ScalarType;
	typedef Dune::BlockVector< Dune::FieldVector<Dune::FieldVector<double, n>, 2*n> >
			VelType;
	typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
	ScalarType saturation;
	ScalarType pressure;
	VelType velocity;

	template<int dim> struct ElementLayout {
		bool contains(Dune::GeometryType gt) {
			return gt.dim() == dim;
		}
	};

	typedef typename G::Traits::LevelIndexSet IS;
	typedef Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout>
			ElementMapper;
	ElementMapper mapper;

	VariableClass(G& g, RT& initialsat = *(new RT(0)), RT& initalpress = *(new RT(0)), Dune::FieldVector<RT, n>& initialvel = *(new Dune::FieldVector<RT, n> (0)), int lev = 0) :
		grid(g), mapper(g, g.levelIndexSet(lev)), size(mapper.size()) {
		initsat(initialsat, size);
		initpress(initalpress, size);
		initvel(initialvel, size);
	}

	void initsat(RT& initialsat, int size) {
		saturation.resize(size);
		saturation=initialsat;
		return;
	}
	void initpress(RT& initialpress, int size) {
		pressure.resize(size);
		pressure=initialpress;
		return;
	}
	void initvel(Dune::FieldVector<RT, n>& initialvel, int size) {
		velocity.resize(size);
		velocity=initialvel;
		return;
	}

	ScalarType& sat() const {
		return saturation;
	}
	ScalarType& press() const {
		return pressure;
	}
	VelType& vel() const {
		return velocity;
	}

	RT& sat(const Dune::FieldVector<DT,n>& x, const Entity& e,
			const Dune::FieldVector<DT,n>& xi) const {
		return saturation[mapper.map(e)];;
	}

	RT& press(const Dune::FieldVector<DT,n>& x, const Entity& e,
			const Dune::FieldVector<DT,n>& xi) const {
		return pressure[mapper.map(e)];
	}

	const Dune::FieldVector<DT,n>& vTotal(const Entity& e,
			const int numberInSelf) const {
		int elemId = mapper.map(e);

		return (velocity[elemId][numberInSelf]);
	}

	/*! @brief prints the saturation to a VTK file
	 * 
	 *  The file name is "<name>-<k>.vtu" where k is an integer number.
	 *  @param name specifies the name of the VTK file
	 *  @param k specifies a number
	 */
	void vtkout(const char* name, int k, int pressurelevel = 0, int satlevel = 0) const {
		if (pressurelevel == satlevel) {
			VTKWriter<G> vtkwriter(grid);
			char fname[128];
			sprintf(fname, "%s-%05d", name, k);
			vtkwriter.addCellData(saturation, "saturation");
			vtkwriter.addCellData(pressure, "total pressure p~");
			vtkwriter.write(fname, VTKOptions::ascii);
		} else {
			Dune::VTKWriter<G, typename G::template Codim<0>::LevelIndexSet>
					vtkwriterpressure(grid, grid.levelIndexSet(pressurelevel));
			char fname[128];
			sprintf(fname, "%s-%05d", name, k);
			vtkwriterpressure.addCellData(pressure, "total pressure p~");
			vtkwriterpressure.write(fname, Dune::VTKOptions::ascii);

			VTKWriter<G, typename G::template Codim<0>::LevelIndexSet>
					vtkwritersaturation(grid, grid.levelIndexSet(satlevel));
			sprintf(fname, "%s-press%05d", name, k);
			vtkwritersaturation.addCellData(saturation, "saturation");
			vtkwritersaturation.write(fname, VTKOptions::ascii);
		}
		return;
	}

	G& grid;
private:
	int size;
};
}
#endif
