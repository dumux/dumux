// $Id: variableclass.hh 670 2008-10-09 07:49:23Z markus $

#ifndef DUNE_SHALLOWVARIABLECLASS_HH
#define DUNE_SHALLOWVARIABLECLASS_HH

/**
 * @file
 * @brief  class including the variables
 * @author Markus Wolff
 */

namespace Dune
{

template<class G, class DT> class ShallowVariableClass
{
	enum
	{	dim=G::dimension};

	typedef Dune::BlockVector<Dune::FieldVector<DT,1> > ScalarType;
	typedef Dune::BlockVector<Dune::FieldVector<DT,dim> > VelType;
	typedef Dune::BlockVector<Dune::FieldVector<DT,dim+1> > SolutionType;
	typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
	template<int dim> struct ElementLayout
	{
		bool contains(Dune::GeometryType gt)
		{
			return gt.dim() == dim;
		}
	};

	typedef typename G::Traits::LeafIndexSet IS;

	typedef Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout>
			ElementMapper;
	ElementMapper elementmapper;

	ScalarType wDepth;
	VelType velocity;
	SolutionType globalSolution;

	G& grid;

	//Constructor: defines the initital values for the variables and the methods includes in the class

	ShallowVariableClass(G& g) :
		elementmapper(g, g.leafIndexSet()), grid(g), size(elementmapper.size())
	{
		sizeInitWDepth(size);
		sizeInitVel(size);
		sizeInitGlobalSolution(size);
	}

	//Definition of the methods to resize the vectors and fill with initial data. !!!Initial data is not the boundary condition, it´s just to initialize sth

	void sizeInitWDepth(int size)
	{
		wDepth.resize(size);
		wDepth=0;
		return;
	}

	void sizeInitVel(int size)
	{
		velocity.resize(size);
		velocity=0;
		return;
	}

	void sizeInitGlobalSolution(int size)
	{
		globalSolution.resize(size);
		globalSolution=0;
		return;
	}

	//Declaration of methods to return variable values

	ScalarType& returnWDepth()
	{
		return wDepth;
	}

	VelType& returnVel()
	{
		return velocity;
	}

	SolutionType& returnGlobalSol()
	{
		return globalSolution;
	}

	//Definition of the return methods

	const DT& returnWDepth(const Dune::FieldVector<DT,dim>& x, const Entity& e,
			const Dune::FieldVector<DT,dim>& xi)
	{
		return wDepth[elementmapper.map(e)];
	}

	const VelType& returnVel(const Dune::FieldVector<DT,dim>& x,
			const Entity& e, const Dune::FieldVector<DT,dim>& xi)
	{
		return velocity[elementmapper.map(e)];
	}

	const SolutionType& returnGlobalSol(const Dune::FieldVector<DT,dim>& x,
			const Entity& e, const Dune::FieldVector<DT,dim>& xi)
	{
		return globalSolution[elementmapper.map(e)];
	}

	void calcPrimVar()
	{
		for (i=0; i<size; i++)
		{

			if (globalSolution[i][0]> 0)
			{

				wDepth[i]=globalSolution[i][0];
				velocity[i][0]=globalSolution[i][1]/wDepth[i];
				velocity[i][1]=globalSolution[i][2]/wDepth[i];
			}
			else
			{
				wDepth[i]=0;
				velocity[i][0]=0;
				velocity[i][1]=0;
			}

			//std::cout<<"global Solution of cell"<<i<<"= "<<globalSolution[i]<<std::endl;

			//std::cout<<"wDepth = "<<wDepth[i]<<std::endl;
			//std::cout<<"velX = "<<velocity[i][0]<<std::endl;
			//std::cout<<"velY = "<<velocity[i][1]<<std::endl;

		}
		return;
	}

	void vtkout(const char* name, int k)
	{
		calcPrimVar();
		VTKWriter<typename G::LeafGridView> vtkwriter(grid.leafView());
		Dune::BlockVector<Dune::FieldVector<DT, 1> > vX(size);
		Dune::BlockVector<Dune::FieldVector<DT, 1> > vY(size);
		vX = 0;
		vY = 0;

		for (int i = 0; i < size; i++)
		{
			std::cout<<"wDepth = "<<wDepth[i]<<std::endl;
			std::cout<<"velX = "<<velocity[i][0]<<std::endl;
			std::cout<<"velY = "<<velocity[i][1]<<std::endl;
			vX[i] = velocity[i][0];
			vY[i] = velocity[i][1];
		}
		char fname[128];
		sprintf(fname, "%s-%05d", name, k);
		vtkwriter.addCellData(vX, "Velocity_X");
		vtkwriter.addCellData(vY, "Velocity_Y");
		vtkwriter.addCellData(wDepth, "wDepth");
		vtkwriter.write(fname, VTKOptions::ascii);

		return;
	}
private:
	int size;
	int i;
	double eps;
};
}
#endif
