// $Id$

#ifndef DUNE_IMPLICITFVTRANSPORT_HH
#define DUNE_IMPLICITFVTRANSPORT_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include "dumux/transport/transport.hh"
#include "dumux/transport/fv/computestorage.hh"
#include "dumux/transport/fv/computenumflux.hh"
#include "dumux/transport/fv/diffusivepart.hh"
#include "dumux/transport/transportproblem.hh"
#include "dumux/nonlinear/newtonmethodmatrix.hh"
#include "dumux/pardiso/pardiso.hh"

/**
 * @file
 * @brief  Finite Volume Transport Model
 * @author Yufei Cao, Bernd Flemisch
 */


namespace Dune
{
//! \ingroup transport
//! The finite volume model for the solution of the transport equation
template<class Grid, class Scalar, class VC, class Problem = TransportProblem<
        Grid, Scalar, VC> >
class ImplicitFVTransport: public Transport<Grid, Scalar, VC, Problem>
{
public:
    enum
    {
        dim = Grid::dimension
    };
    enum
    {
        dimWorld = Grid::dimensionworld
    };

    typedef typename VC::ScalarVectorType RepresentationType;
    typedef typename VC::VelType VelType;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LevelGridView GridView;
    typedef typename GridView::IndexSet IndexSet;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef BlockVector< Dune::FieldVector<Scalar,dim> > SlopeType;

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

    typedef Dune::FieldMatrix<Scalar, 1, 1> MB;
    typedef BCRSMatrix<MB> MatrixType;
    typedef BlockVector<Dune::FieldVector<Scalar, 1> > VectorType;

    typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
    typedef ImplicitFVTransport<Grid, Scalar, VC, Problem> ThisType;

    MatrixType A;
    VectorType f;
    VectorType u;
    VectorType uOldTimeStep;

    /*!
     *  \param t time
     *  \param dt time step size to estimate
     *  \param update vector to be filled with the update
     *
     *  This method calculates the update vector, i.e., the FV discretization
     *  of \f$\text{div}\, (f_\text{w} \boldsymbol{v}_t)\f$. The total velocity
     *  \f$\boldsymbol{v}_t)\f$ has to be given by the block vector \a velocity,
     *  containing values at each cell face. The fractional flow function \f$f_\text{w}\f$
     *  is realized by the numerical flux function \a numericalFlux such that the \f$i\f$-th entry
     *  of the vector update is obtained by calculating
     *  \f[ \sum_{j \in \mathcal{N}(i)} v_{ij}g(S_i, S_j) - v_{ji}g(S_j, S_i), \f]
     *  where \f$\mathcal{N}(i)\f$ denotes the set of neighbors of the cell \f$i\f$ and
     *  \f$g\f$ stands for \a numericalFlux. The normal velocities \f$v_{ij}\f$ and \f$v_{ji}\f$
     *  are given by
     *  \f{align*} v_{ij} = \int_{\Gamma_{ij}} \max(\boldsymbol{v}_\text{t}{\cdot}\boldsymbol{n}_{ij}, \, 0), \qquad
     *  v_{ji} = \int_{\Gamma_{ij}} \min(\boldsymbol{v}_\text{t}{\cdot}\boldsymbol{n}_{ij}, \, 0), \f}
     *  where \f$\boldsymbol{n}_{ij}\f$ denotes the unit normal vector from cell \f$i\f$ towards cell \f$j\f$.
     *
     *  Additionally to the \a update vector, the recommended time step size \a dt is calculated
     *  employing the usual CFL condition.
     */
    int update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar& cFLFac);

    void initializeMatrix();

    void assemble();

    void initialTransport();

    void solve();

    void update();

    const Scalar getDt()
    {
	return dt_;
    }

    void setDt(const Scalar dt)
    {
	dt_ = dt;
    }

    MatrixType& matrix()
    {
        return A;
    }

    VectorType& rhs()
    {
        return f;
    }

    VectorType& sol()
    {
        return u;
    }

    /*! @brief constructor
     *
     * @param grid a DUNE grid object
     * @param problem an object of class TransportProblem or derived
     * @param time time step
     * @param diffPart an object of class DiffusivePart or derived. This determines the diffusive flux incorporated in the transport.
     * @param numFl an object of class Numerical Flux or derived
     * @param storage an object of class ComputeStorage or derived
     */
    ImplicitFVTransport(Grid& grid, Problem& problem, const Scalar time,
            DiffusivePart<Grid,Scalar>& diffPart = *(new DiffusivePart<Grid, Scalar>), 
            ComputeNumFlux<Grid,Scalar>& numFl = *(new ComputeNumFlux<Grid,Scalar>),
            ComputeStorage<Scalar>& storage = *(new ComputeStorage<Scalar>)):
    Transport<Grid, Scalar, VC, Problem>(grid, problem), 
    A(grid.size(this->level(), 0), grid.size(this->level(), 0), (2*dim+1)*grid.size(this->level(), 0), BCRSMatrix<MB>::random),
    f(grid.size(this->level(), 0)), u(grid.size(this->level(), 0)), uOldTimeStep(grid.size(this->level(), 0)),
    grid_(grid), dt_(time), diffusivePart_(diffPart), numFlux_(numFl), storage_(storage)
    {
	initializeMatrix();
    }

private:
    const Grid& grid_;
    Scalar dt_;
    const DiffusivePart<Grid, Scalar>& diffusivePart_;
    const ComputeNumFlux<Grid, Scalar>& numFlux_;
    const ComputeStorage<Scalar>& storage_;
};

template<class Grid, class Scalar, class VC, class Problem>
int ImplicitFVTransport<Grid, Scalar, VC, Problem>::update(const Scalar t, Scalar& dt,
        RepresentationType& updateVec, Scalar& cFLFac = 1)
{
	update(); 

	updateVec = u; 
	updateVec -= uOldTimeStep;

	// divide by the fixed timestep 
	updateVec /= dt_; 

	dt = dt_;

	return 0;
}

template<class Grid, class Scalar, class VC, class Problem>
void ImplicitFVTransport<Grid, Scalar, VC, Problem>::initializeMatrix()
{

    const GridView& gridView = this->grid_.levelView(this->level());

    // determine matrix row sizes
    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = this->transProblem.variables.diffMapper.map(*eIt);

        // initialize row size
        int rowSize = 1;

        // run through all intersections with neighbors
        IntersectionIterator
        isItEnd = gridView.template iend(*eIt);
        for (IntersectionIterator
                isIt = gridView.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        if (isIt->neighbor())
        rowSize++;
        A.setrowsize(globalIdxI, rowSize);
    }
    A.endrowsizes();

    // determine position of matrix entries
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxI = this->transProblem.variables.diffMapper.map(*eIt);

        // add diagonal index
        A.addindex(globalIdxI, globalIdxI);

        // run through all intersections with neighbors
        IntersectionIterator
        isItEnd = gridView.template iend(*eIt);
        for (IntersectionIterator
                isIt = gridView.template ibegin(*eIt); isIt
                !=isItEnd; ++isIt)
        if (isIt->neighbor())
        {
            // access neighbor
            ElementPointer outside = isIt->outside();
            int globalIdxJ = this->transProblem.variables.diffMapper.map(*outside);

            // add off diagonal index
            A.addindex(globalIdxI, globalIdxJ);
        }
    }
    A.endindices();

    return;
}


template<class Grid, class Scalar, class VC, class Problem>
void ImplicitFVTransport<Grid, Scalar, VC, Problem>::assemble()
{
    // initialization: set matrix A to zero
    A = 0;

    const GridView& gridView = this->grid_.levelView(this->level());

    // loop over elements 
    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // cell geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // cell center in reference element
        const LocalPosition
        &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // cell center in global coordinates
        GlobalPosition globalPos = eIt->geometry().global(localPos);

        // cell volume, assume linear map here
        Scalar volume = eIt->geometry().integrationElement(localPos)
        *Dune::ReferenceElements<Scalar,dim>::general(gt).volume();

        // cell index
        int globalIdxI = this->transProblem.variables.transMapper.map(*eIt);

	// get saturation value at cell center
        Scalar satI = u[globalIdxI];
        Scalar satOldI = uOldTimeStep[globalIdxI];

        // get residual saturation 
        Scalar residualSatW = this->transProblem.soil().Sr_w(globalPos, *eIt, localPos);
        Scalar residualSatNW = this->transProblem.soil().Sr_n(globalPos, *eIt, localPos);

        // storage term
        Scalar storageU = storage_(satI) - storage_(satOldI);

        // set f_i
        f[globalIdxI] = storageU * this->transProblem.soil().porosity(globalPos, *eIt,localPos) * volume / dt_;

	// calculate epsilon 
	Scalar epsilon = 1e-5;

	// add epsilon*e_i to the solution 
        Scalar satIPlus = satI + epsilon;

	// storage term 
	Scalar storageUPlusEps = storage_(satIPlus)- storage_(satOldI);

	// subtract epsilon*e_i from the solution 
        Scalar satIMinus = satI - epsilon;

	// storage term 
	Scalar storageUMinusEps = storage_(satIMinus)- storage_(satOldI);

	// set A_ii
 	A[globalIdxI][globalIdxI] += 0.5/epsilon*(storageUPlusEps - storageUMinusEps);

        A[globalIdxI][globalIdxI] *= this->transProblem.soil().porosity(globalPos, *eIt,localPos) * volume;
        A[globalIdxI][globalIdxI] /= dt_;

        // run through all intersections with neighbors and boundary
        IntersectionIterator isItEnd = gridView.template iend(*eIt);
        for (IntersectionIterator isIt = gridView.template ibegin(*eIt); isIt != isItEnd; ++isIt)
        {
            // local number of facet
            int indexInInside = isIt->indexInInside();

            // get geometry type of face
            Dune::GeometryType faceGT = isIt->geometryInInside().type();

            // center in face's reference element
            const Dune::FieldVector<Scalar,dim-1>&
            faceLocal = Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

            // center of face inside volume reference element
            const LocalPosition&
            localPosFace = Dune::ReferenceElements<Scalar,dim>::general(faceGT).position(indexInInside,1);

            // get normal vector 
            Dune::FieldVector<Scalar,dimWorld> integrationOuterNormal = isIt->integrationOuterNormal(faceLocal);
            integrationOuterNormal *= Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).volume();

            // get face volume
            Scalar faceVol = isIt->geometry().volume();

            // in fact, satAvg is not used in diffusivePart_()!
            Scalar satAvg = 0;

            // handle interior face
            if (isIt->neighbor())
            {
                // access neighbor
                ElementPointer neighborPointer = isIt->outside();
                int globalIdxJ = this->transProblem.variables.transMapper.map(*neighborPointer);

                // compute factor in neighbor
                Dune::GeometryType neighborGT = neighborPointer->geometry().type();
                const LocalPosition&
                localPosNeighbor = Dune::ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);
               
                // cell center in global coordinates
                const GlobalPosition& globalPos = eIt->geometry().global(localPos);

                // neighbor cell center in global coordinates
                const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

                // distance vector between barycenters
                Dune::FieldVector<Scalar,dimWorld> distVec = globalPos - globalPosNeighbor;

                // compute distance between cell centers
                Scalar dist = distVec.two_norm();

	        // get saturation value at cell center
                Scalar satJ = u[globalIdxJ];

                // calculate the saturation gradient
                Dune::FieldVector<Scalar,dim> satGradient = distVec;
                satGradient *= (satJ - satI)/(dist*dist);

		// flux term 
		Scalar fluxU = numFlux_(*eIt, indexInInside, satI, satJ) -
                               diffusivePart_(*eIt, indexInInside, satAvg, satGradient, dt_, satI, satJ) 
                               * integrationOuterNormal;

                // add to f_i
                f[globalIdxI] += fluxU;

		// flux term 
		Scalar fluxUIPlusEps;
		if (satIPlus > 1.0-residualSatNW)
		{
			fluxUIPlusEps = fluxU; 
			epsilon *= 0.5;
		}
		else 
                {
                        // calculate the saturation gradient
                        Dune::FieldVector<Scalar,dim> satGradientIPlus = distVec;
                        satGradientIPlus *= (satJ - satIPlus)/(dist*dist);

 			fluxUIPlusEps = numFlux_(*eIt, indexInInside, satIPlus, satJ) -
                                       diffusivePart_(*eIt, indexInInside, satAvg, satGradientIPlus, dt_, satIPlus, satJ) 
                                       * integrationOuterNormal;
                }

		// flux term 
		Scalar fluxUIMinusEps;
		if (satIMinus < residualSatW)
		{
			fluxUIMinusEps = fluxU; 
			epsilon *= 0.5;
		}
		else 
                {
                        // calculate the saturation gradient
                        Dune::FieldVector<Scalar,dim> satGradientIMinus = distVec;
                        satGradientIMinus *= (satJ - satIMinus)/(dist*dist);

 			fluxUIMinusEps = numFlux_(*eIt, indexInInside, satIMinus, satJ) -
                                        diffusivePart_(*eIt, indexInInside, satAvg, satGradientIMinus, dt_, satIMinus, satJ) 
                                        * integrationOuterNormal;
                }

	 	A[globalIdxI][globalIdxI] += 0.5/epsilon*(fluxUIPlusEps - fluxUIMinusEps);

		// restore epsilon 
		if (satIPlus > 1.0-residualSatNW || satIMinus < residualSatW)
			epsilon *= 2.0;

		// add epsilon*e_j to the solution 
                Scalar satJPlus = satJ + epsilon;

		// flux term 
		Scalar fluxUJPlusEps;
		if (satJPlus > 1.0-residualSatNW)
		{
			fluxUJPlusEps = fluxU;
			epsilon *= 0.5;
		}
		else
                {
                        // calculate the saturation gradient
                        Dune::FieldVector<Scalar,dim> satGradientJPlus = distVec;
                	satGradientJPlus *= (satJPlus - satI)/(dist*dist);
 
 			fluxUJPlusEps = numFlux_(*eIt, indexInInside, satI, satJPlus) -
                                       diffusivePart_(*eIt, indexInInside, satAvg, satGradientJPlus, dt_, satI, satJPlus)
                                       * integrationOuterNormal;
                }

		// subtract epsilon*e_j from the solution 
                Scalar satJMinus = satJ - epsilon;

		// flux term 
		Scalar fluxUJMinusEps;
		if (satJMinus < residualSatW)
		{
			fluxUJMinusEps = fluxU;
			epsilon *= 0.5;
		}
		else 
                {
                	// calculate the saturation gradient
                	Dune::FieldVector<Scalar,dim> satGradientJMinus = distVec;
                	satGradientJMinus *= (satJMinus - satI)/(dist*dist);

 			fluxUJMinusEps  = numFlux_(*eIt, indexInInside, satI, satJMinus) -
                                        diffusivePart_(*eIt, indexInInside, satAvg, satGradientJMinus, dt_, satI, satJMinus)
                                        * integrationOuterNormal;
		}

		// add to A_ij
	 	A[globalIdxI][globalIdxJ] += 0.5/epsilon*(fluxUJPlusEps - fluxUJMinusEps);

		// restore epsilon 
		if (satJPlus > 1.0-residualSatNW || satJMinus < residualSatW)
			epsilon *= 2.0;
            }

            // handle boundary face
            if (isIt->boundary())
            {
                // center of face in global coordinates
                GlobalPosition globalPosFace = isIt->geometry().global(faceLocal);

                //get boundary type
                BoundaryConditions::Flags bctype = this->transProblem.bctypeSat(globalPosFace, *eIt, localPosFace);

                // Dirichlet boundary
                if (bctype == BoundaryConditions::dirichlet)
                {
                    Scalar satBound = this->transProblem.dirichletSat(globalPosFace, *eIt, localPosFace);

                    // distance vector between barycenters
                    Dune::FieldVector<Scalar,dimWorld> distVec = globalPos - globalPosFace;

                    // compute distance between cell centers
                    Scalar dist = distVec.two_norm();

                    // calculate the saturation gradient
                    Dune::FieldVector<Scalar,dim> satGradient = distVec;
                    satGradient *= (satBound - satI)/(dist*dist);      

                    // flux term 
	            Scalar fluxU = numFlux_(*eIt, indexInInside, satI, satBound) -
                                   diffusivePart_(*eIt, indexInInside, satAvg, satGradient, dt_, satI, satBound) 
                                   * integrationOuterNormal;

                    // add to f_i
                    f[globalIdxI] += fluxU;

                    // flux term 
	            Scalar fluxUIPlusEps;
		    if (satIPlus > 1.0-residualSatNW)
		    {
			fluxUIPlusEps = fluxU; 
			epsilon *= 0.5;
		    }
		    else 
                    {
                        // calculate the saturation gradient
                        Dune::FieldVector<Scalar,dim> satGradientIPlus = distVec;
                        satGradientIPlus *= (satBound - satIPlus)/(dist*dist);

 			fluxUIPlusEps = numFlux_(*eIt, indexInInside, satIPlus, satBound) -
                                        diffusivePart_(*eIt, indexInInside, satAvg, satGradientIPlus, dt_, satIPlus, satBound) 
                                        * integrationOuterNormal;
	            }

		    // flux term 
		    Scalar fluxUIMinusEps;
		    if (satIMinus < residualSatW)
		    {
	                fluxUIMinusEps = fluxU; 
			epsilon *= 0.5;
		    }
		    else 
                    {
                        // calculate the saturation gradient
                        Dune::FieldVector<Scalar,dim> satGradientIMinus = distVec;
                        satGradientIMinus *= (satBound - satIMinus)/(dist*dist);

 			fluxUIMinusEps = numFlux_(*eIt, indexInInside, satIMinus, satBound) -
                                         diffusivePart_(*eIt, indexInInside, satAvg, satGradientIMinus, dt_, satIMinus, satBound) 
                                         * integrationOuterNormal;
                    }                    

                    // get the flux through the boundary
                    Scalar velocity = this->transProblem.variables.vTotal(*eIt, indexInInside)*integrationOuterNormal;

		    // add to A_ii
                    // outflow boundary
                    if (velocity >= 0.0)
                        A[globalIdxI][globalIdxI] += 0.0;
                    // inflow boundary
                    else
	 	        A[globalIdxI][globalIdxI] += 0.5/epsilon*(fluxUIPlusEps - fluxUIMinusEps);         

		    // restore epsilon 
		    if (satIPlus > 1.0-residualSatNW || satIMinus < residualSatW)
			epsilon *= 2.0;
		}
                // Neumann boundary
                else
                {
                    // handling 1: assume the domain is big enough
                    // outflow boundary is constant Neumann boundary
		    // add to A_ii
	 	    A[globalIdxI][globalIdxI] += 0.0;     

                    // in fact, for this case, helpFactor is not used !
                    Scalar helpFactor = 0;      

                    // get the Neumann boundary value               
                    Scalar fluxU = this->transProblem.neumannSat(globalPosFace, *eIt, localPosFace, helpFactor) * faceVol;

                    // add to f_i
                    f[globalIdxI] += fluxU;
                }
            }
        }
    } // end grid traversal

    // for test
    //printmatrix(std::cout, A, "global stiffness matrix", "row", 11, 3);
    //printvector(std::cout, f, "right hand side", "row", 200, 1, 3);
}

template<class Grid, class Scalar, class VC, class Problem>
void ImplicitFVTransport<Grid, Scalar, VC, Problem>::initialTransport()
{
    const GridView& gridView = this->grid_.levelView(this->level());

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = gridView.template end<0>();
    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eItEnd; ++eIt)
    {
        // get geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // get cell center in reference element
        const LocalPosition
        &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        GlobalPosition globalPos = eIt->geometry().global(localPos);

        // initialize cell concentration
        this->transProblem.variables.saturation[this->transProblem.variables.transMapper.map(*eIt)] = this->transProblem.initSat(globalPos, *eIt, localPos);
    }

    u = this->transProblem.variables.saturation;

    return;
}


template<class Grid, class Scalar, class VC, class Problem>
void ImplicitFVTransport<Grid, Scalar, VC, Problem>:: solve()
    {
        Operator op(A);  // make operator out of matrix
        double red=1E-12;

#ifdef HAVE_PARDISO
        SeqPardiso<MatrixType,VectorType,VectorType> pardiso(A);
        LoopSolver<VectorType> solver(op, pardiso, red, 10, 2);
#else
        SeqILU0<MatrixType,VectorType,VectorType> ilu0(A,1.0);// a precondtioner
        BiCGSTABSolver<VectorType> solver(op,ilu0,red,10000,1);         // an inverse operator
#endif
        InverseOperatorResult r;
        solver.apply(u, f, r);

        return;
    }


template<class Grid, class Scalar, class VC, class Problem>
void ImplicitFVTransport<Grid, Scalar, VC, Problem>:: update()
    {
        uOldTimeStep = u;
        //NewtonMethodMatrix<Grid, ThisType> newtonMethod(this->grid_, *this);

        //In order to not enlarge the timestep in the Newton step
        double dtol = 1e-7;
        double rtol = 1e7;
        int maxIt = 20;
        double mindt = 1;
        int goodIt = 1;
        NewtonMethodMatrix<Grid, ThisType> newtonMethod(this->grid_, *this, dtol, rtol, maxIt, mindt, goodIt);

        newtonMethod.execute();

        return;
    }
}
#endif
