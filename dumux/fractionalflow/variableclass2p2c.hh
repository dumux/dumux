#ifndef DUNE_VARIABLECLASS2P2C_HH
#define DUNE_VARIABLECLASS2P2C_HH

/**
 * @file
 * @brief  class including the variables needed for decoupled 2p2c computations.
 * @author Jochen Fritz
 */

namespace Dune {

template<class G, class RT> class VariableClass2p2c {

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
	ScalarType totalConcentration;
	ScalarType wet_c1, nonwet_c1;
	ScalarType wet_c2, nonwet_c2;
	ScalarType volErr;

	template<int dim> struct ElementLayout 
	{
		bool contains(Dune::GeometryType gt) 
		{
			return gt.dim() == dim;
		}
	};

	typedef typename G::LevelGridView::IndexSet IS;
	typedef Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout>	ElementMapper;
	ElementMapper mapper;

	VariableClass2p2c(G& g, RT& initialsat = *(new RT(0)), RT& initalpress = *(new RT(0)), Dune::FieldVector<RT, n>& initialvel = *(new Dune::FieldVector<RT, n> (0)), int lev = 0) :
		grid(g), mapper(g, g.levelIndexSet(lev)), size(mapper.size()) 
	{
		initsat(initialsat, size);
		initpress(initalpress, size);
		initvel(initialvel, size);
		initC(size);
	}

	void initsat(RT& initialsat, int size) 
	{
		saturation.resize(size);
		saturation=initialsat;
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
	void initC(int size)
	{
		totalConcentration.resize(2*size);
		wet_c1.resize(size);
		wet_c2.resize(size);
		nonwet_c1.resize(size);
		nonwet_c2.resize(size);
		volErr.resize(size);
		
		totalConcentration = 0;
		wet_c1 = 0;
		wet_c2 = 0; 
		nonwet_c1 = 0;
		nonwet_c2 = 0;
		volErr = 0;
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

	const Dune::FieldVector<RT,1>& sat(const Dune::FieldVector<DT,n>& x, const Entity& e,
			const Dune::FieldVector<DT,n>& xi) const {
		return saturation[mapper.map(e)];;
	}

	const Dune::FieldVector<RT,1>& press(const Dune::FieldVector<DT,n>& x, const Entity& e,
			const Dune::FieldVector<DT,n>& xi) const {
		return pressure[mapper.map(e)];
	}

	const Dune::FieldVector<DT,n>& vTotal(const Entity& e,
			const int numberInSelf) const {
		int elemId = mapper.map(e);

		return (velocity[elemId][numberInSelf]);
	}

	/*! @brief writes all variables to a VTK File
	 * 
	 *  The file name is "<name>-<k>.vtu" where k is an integer number.
	 *  @param name specifies the name of the VTK file
	 *  @param k specifies a number
	 */
	void vtkout(const char* name, int k, int pressurelevel = 0, int satlevel = 0) const 
	{
		if (pressurelevel == satlevel) 
		{
			ScalarType C1, C2;
			C1.resize(size); C2.resize(size);
			for (int i = 0; i < size; i++)
			{
				C1[i] = totalConcentration[i];
				C2[i] = totalConcentration[i + size];
			}
			VTKWriter<G> vtkwriter(grid);
			char fname[128];
			sprintf(fname, "%s-%05d", name, k);
			vtkwriter.addCellData(saturation, "saturation");
			vtkwriter.addCellData(pressure, "total pressure p~");
			vtkwriter.addCellData(C1, "total concentration 1 [kg/m^3]");
			vtkwriter.addCellData(C2, "total concentration 2 [kg/m^3]");		
			vtkwriter.addCellData(volErr, "volume error [%]");
			vtkwriter.addCellData(wet_c1, "concentration 1 in wetting phase [kg/m^3]");
			vtkwriter.addCellData(nonwet_c1, "concentration 1 in non-wetting phase [kg/m^3]");
			vtkwriter.addCellData(wet_c2, "concentration 2 in wetting phase [kg/m^3]");
			vtkwriter.addCellData(nonwet_c2, "concentration 2 in non-wetting phase [kg/m^3]");
			vtkwriter.write(fname, VTKOptions::ascii);
		} 
		else 
		{
			Dune::VTKWriter<G, typename G::LevelGridView>
					vtkwriterpressure(grid.levelView(pressurelevel));
			char fname[128];
			sprintf(fname, "%s-%05d", name, k);
			vtkwriterpressure.addCellData(pressure, "total pressure p~");
			vtkwriterpressure.write(fname, Dune::VTKOptions::ascii);

			VTKWriter<G, typename G::LevelGridView>
					vtkwritersaturation(grid.levelView(satlevel));
			sprintf(fname, "%s-press%05d", name, k);
			vtkwritersaturation.addCellData(saturation, "saturation");
			vtkwritersaturation.write(fname, VTKOptions::ascii);
		}
		return;
	}
	void vtkoutpressure(const char* name, int k) const {
		VTKWriter<G> vtkwriter(grid);
		char fname[128];
		sprintf(fname, "%s-press%05d", name, k);
		vtkwriter.addCellData(pressure, "total pressure p~");
		vtkwriter.write(fname, VTKOptions::ascii);
	}	

	G& grid;
private:
	int size;
};
}
#endif
