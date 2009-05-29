// $Id$

#ifndef DUNE_MPFAODIFFUSION_HH
#define DUNE_MPFAODIFFUSION_HH

#include <dune/common/helpertemplates.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dumux/pardiso/pardiso.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/common/typetraits.hh>
#include "dumux/diffusion/diffusion.hh"
#include "dumux/diffusion/diffusionproblem.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical diffusion model
 * @brief  MPFA O-method 
 * @brief  Remark1: only for 2-D quadrilateral grid.
 * @brief  Remark2: can use UGGrid or SGrid (YaspGrid); variable 'ch' is chosen to decide which grid will be used.
 * @brief  Remark3: without capillary pressure and gravity!
 * @author Yufei Cao
 */

namespace Dune
{
  //! \ingroup diffusion
  //! Base class for defining an instance of a numerical diffusion model.
  /*! An interface for defining a numerical diffusion model for the 
   *  solution of equations of the form 
   * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = 0, \f$, 
   * \f$p = g\f$ on \f$\Gamma_1\f$, and 
   * \f$-\lambda K \text{grad}\, p \cdot \mathbf{n} = J\f$ 
   * on \f$\Gamma_2\f$. Here, 
   * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability, 
   * and \f$\lambda\f$ the total mobility, possibly depending on the 
   * saturation.
	Template parameters are:

	- GridView      a DUNE gridView type
	- Scalar        type used for return values 
   */
  template<class GridView, class Scalar, class VC, class Problem = DiffusionProblem<GridView, Scalar, VC> >
  class MPFAODiffusion 
  : public Diffusion< GridView, Scalar, VC, Problem> 
  {	  
	  enum{dim = GridView::dimension};	
	  enum{dimWorld = GridView::dimensionworld};
	  typedef typename GridView::Traits::template Codim<0>::Entity Element;
	  typedef typename GridView::IndexSet IndexSet;
	  typedef typename GridView::template Codim<0>::Iterator ElementIterator;
	  typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
          typedef typename GridView::IntersectionIterator IntersectionIterator;

	  typedef Dune::FieldVector<Scalar,dim> LocalPosition;
	  typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
	  typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

	  typedef Dune::FieldMatrix<double,1,1> MB;
	  typedef Dune::BCRSMatrix<MB> MatrixType;
	  typedef BlockVector<FieldVector<Scalar, 1> > Vector;
	  
  public:
	typedef BlockVector< FieldVector<FieldVector<Scalar, dim>, 2*dim> > VelType;
	typedef BlockVector< Dune::FieldVector<Scalar,1> > RepresentationType;

        void assemble(const Scalar t); 

	void solve(); 

	void pressure(const Scalar t=0)
	{
                Dune::Timer timer;

                timer.reset();
		assemble(t);
                std::cout << "assembling MPFA O-matrix on level" << this->gridView.grid().maxLevel() << " took " << timer.elapsed() << " seconds" << std::endl;

                timer.reset();
		solve();
                std::cout << "solving MPFA O-matrix on level" << this->gridView.grid().maxLevel() << " took " << timer.elapsed() << " seconds" << std::endl;

		return;
	}
	
	void initializeMatrix();

	MPFAODiffusion(GridView& gridView, Problem& problem)
	  : Diffusion<GridView, Scalar, VC, Problem>(gridView, problem), 
	  M(problem.variables().gridSizeDiffusion(), problem.variables().gridSizeDiffusion(), (4*dim+(dim-1))*problem.variables().gridSizeDiffusion(), Dune::BCRSMatrix<MB>::random), 
	  f(problem.variables().gridSizeDiffusion()),
          solverName_("BiCGSTAB"),
	  preconditionerName_("SeqPardiso")
	{
		initializeMatrix();
	}

	MPFAODiffusion(GridView& gridView, Problem& problem, std::string solverName, std::string preconditionerName)
	  : Diffusion<GridView, Scalar, VC, Problem>(gridView, problem), 
	    M(problem.variables().gridSizeDiffusion(), problem.variables().gridSizeDiffusion(), (4*dim+(dim-1))*problem.variables().gridSizeDiffusion(), Dune::BCRSMatrix<MB>::random), 
	    f(problem.variables().gridSizeDiffusion()),
            solverName_(solverName),
	    preconditionerName_(preconditionerName) 
	{
		initializeMatrix();
	}

private:

	MatrixType M;
	RepresentationType f;
	std::string solverName_;
	std::string preconditionerName_;
  };


  template<class GridView, class Scalar, class VC, class Problem>
  void MPFAODiffusion<GridView, Scalar, VC, Problem>::initializeMatrix()
  {
            // input character to decide which grid to use; 's' for SGrid, 'u' for UGGrid
            char ch = 'u';

	    // determine matrix row sizes 
	    ElementIterator eItBegin = this->gridView.template begin<0>();
	    ElementIterator eItEnd = this->gridView.template end<0>();
	    for (ElementIterator eIt = eItBegin; eIt != eItEnd; ++eIt)
	      {
			// cell index
			int globalIdxI = this->diffProblem.variables().indexDiffusion(*eIt);
	
			// initialize row size
			int rowSize = 1;
			
			// run through all intersections with neighbors 
                        IntersectionIterator isItBegin = this->gridView.template ibegin(*eIt);
			IntersectionIterator isItEnd = this->gridView.template iend(*eIt);
			for (IntersectionIterator isIt = isItBegin; isIt!=isItEnd; ++isIt)
			  {
                            IntersectionIterator tempisIt = isIt;
			    IntersectionIterator tempisItBegin = isItBegin;

                            // 'nextisIt' iterates over next codimension 1 intersection neighboring with 'isIt'
                            IntersectionIterator nextisIt = ++tempisIt;

			    // get 'nextisIt'
			    switch (ch)
			      {
			        // for SGrid
			        case 's':
			        {
                                  if (nextisIt == isItEnd)
                                    nextisIt = isItBegin; 
			          else
			            {
			     	      nextisIt = ++tempisIt;
			       
				      if (nextisIt == isItEnd)
				        {
				          nextisIt = ++tempisItBegin;
				        }
			            }

				  break;
			        }
			        // for UGGrid
			        case 'u':
			        {
			  	  if (nextisIt == isItEnd)
				    nextisIt = isItBegin;

				  break;
			        }
			      }                        

			    if (isIt->neighbor()) 
			      rowSize++;			     

                            if (isIt->neighbor() && nextisIt->neighbor())
                              rowSize++;
			  } // end of 'for' IntersectionIterator			  

			// set number of indices in row globalIdxI to rowSize 			
			M.setrowsize(globalIdxI, rowSize);

	      } // end of 'for' ElementIterator

	    // indicate that size of all rows is defined 
	    M.endrowsizes();

	    // determine position of matrix entries 
	    for (ElementIterator eIt = eItBegin; eIt != eItEnd; ++eIt)
	      {
			// cell index
			int globalIdxI = this->diffProblem.variables().indexDiffusion(*eIt);
	
			// add diagonal index
			M.addindex(globalIdxI, globalIdxI);
		
			// run through all intersections with neighbors 
                        IntersectionIterator isItBegin = this->gridView.template ibegin(*eIt);
			IntersectionIterator isItEnd = this->gridView.template iend(*eIt);
			for (IntersectionIterator isIt = isItBegin; isIt!=isItEnd; ++isIt)
			  {
                            IntersectionIterator tempisIt = isIt;
			    IntersectionIterator tempisItBegin = isItBegin;

                            // 'nextisIt' iterates over next codimension 1 intersection neighboring with 'isIt'
			    // sequence of next is anticlockwise of 'isIt'
                            IntersectionIterator nextisIt = ++tempisIt;

 			    // get 'nextisIt'
			    switch (ch)
			      {
			        // for SGrid
			        case 's':
			        {
			          if (nextisIt == isItEnd)
			            {
			 	      nextisIt = isItBegin;
			            }
			          else
			            {
			  	      nextisIt = ++tempisIt;

				      if (nextisIt == isItEnd)
				        {
				          nextisIt = ++tempisItBegin;
				        }
			            }

				  break;
			        }
			        // for UGGrid
			        case 'u':
			        {
				  if (nextisIt == isItEnd)
				    nextisIt = isItBegin;

				  break;
			        }
			      }

			    if (isIt->neighbor()) 
			      {
					// access neighbor
					ElementPointer outside = isIt->outside();
					int globalIdxJ = this->diffProblem.variables().indexDiffusion(*outside);
		
					// add off diagonal index
					// add index (row,col) to the matrix 
					M.addindex(globalIdxI, globalIdxJ);
			      }

			    if (isIt->neighbor() && nextisIt->neighbor())
			      {
				// access the common neighbor of isIt's and nextisIt's outside
                                ElementPointer outside = isIt->outside();
				ElementPointer nextisItoutside = nextisIt->outside();
                                    
                                IntersectionIterator innerisItEnd = this->gridView.template iend(*outside);
                                IntersectionIterator innernextisItEnd = this->gridView.template iend(*nextisItoutside);

				for (IntersectionIterator innerisIt = this->gridView.template ibegin(*outside);
				     innerisIt!=innerisItEnd; ++innerisIt )
				  for (IntersectionIterator innernextisIt = this->gridView.template ibegin(*nextisItoutside);
				       innernextisIt!=innernextisItEnd; ++innernextisIt)
				    {
				      if (innerisIt->neighbor() && innernextisIt->neighbor())
					{
					  ElementPointer innerisItoutside = innerisIt->outside();
					  ElementPointer innernextisItoutside = innernextisIt->outside();

					  if (innerisItoutside == innernextisItoutside && innerisItoutside != isIt->inside())
					    {
					      int globalIdxJ = this->diffProblem.variables().indexDiffusion(*innerisItoutside);
					      
					      M.addindex(globalIdxI, globalIdxJ);
                                            }
					}
				    }
			      }			    
			  } // end of 'for' IntersectionIterator
	      } // end of 'for' ElementIterator

	    // indicate that all indices are defined, check consistency 
	    M.endindices();		

	    return;
  }

  
  // only for 2-D general quadrilateral
  template<class GridView, class Scalar, class VC, class Problem>
  void MPFAODiffusion<GridView, Scalar, VC, Problem>::assemble(const Scalar t=0)
	{
          // input character to decide which grid to use; 's' for SGrid, 'u' for UGGrid
          char ch = 'u';

	  // initialization: set global matrix M to zero
	  M = 0;

	  // introduce matrix R for vector rotation and R is initialized as zero matrix        
          FieldMatrix R(0);

	  // evaluate matrix R
	  if (dim==2)
	    for (int i=0; i<dim; ++i)
	      {
		R[0][1] = 1;
		R[1][0] = -1;
	      }
	      
	  // run through all elements
	  ElementIterator eItEnd = this->gridView.template end<0>();
	  for (ElementIterator eIt = this->gridView.template begin<0>(); eIt != eItEnd; ++eIt)
	    {
	                // get common geometry information for the following computation

			// cell 1 geometry type
			Dune::GeometryType gt1 = eIt->geometry().type();
			
			// cell 1 center in reference element
			const LocalPosition& 
                          localPos1 = Dune::ReferenceElements<Scalar,dim>::general(gt1).position(0,0);
			
			// get global coordinate of cell 1 center
			GlobalPosition globalPos1 = eIt->geometry().global(localPos1);

			// cell 1 volume
			double volume1 = eIt->geometry().integrationElement(localPos1)
			  *Dune::ReferenceElements<Scalar,dim>::general(gt1).volume();			  
			  			
			// cell 1 index
			int globalIdx1 = this->diffProblem.variables().indexDiffusion(*eIt);

			// evaluate right hand side
			f[globalIdx1] = this->diffProblem.sourcePress(globalPos1, *eIt, localPos1)*volume1;
	
			// get absolute permeability of cell 1
			FieldMatrix K1(this->diffProblem.soil().K(globalPos1,*eIt,localPos1));
	
			//compute total mobility of cell 1
			double lambda1;
                        double sat1 = this->diffProblem.variables().saturation()[globalIdx1];
	                lambda1 = this->diffProblem.materialLaw().mobTotal(sat1, globalPos1,*eIt,localPos1); 

			// if K1 is zero, no flux through cell1
			// for 2-D
			if (K1[0][0] == 0 && K1[0][1] == 0 && K1[1][0] == 0 && K1[1][1] == 0)
                          {
                            M[globalIdx1][globalIdx1] += 1.0;
                            continue;
                          }

			IntersectionIterator isItBegin = this->gridView.template ibegin(*eIt);
			IntersectionIterator isItEnd = this->gridView.template iend(*eIt);
			for (IntersectionIterator isIt = isItBegin; isIt!=isItEnd; ++isIt)
			  {
			    // intersection iterator 'nextisIt' is used to get geometry information
			    IntersectionIterator tempisIt = isIt;
			    IntersectionIterator tempisItBegin = isItBegin;

			    IntersectionIterator nextisIt = ++tempisIt;

			    // get 'nextisIt'
			    switch (ch)
			      {
			        // for SGrid
			        case 's':
			          {
			            if (nextisIt == isItEnd)
			              nextisIt = isItBegin;
			            else
			              {
			   	        nextisIt = ++tempisIt;

				        if (nextisIt == isItEnd)
				          {
				            nextisIt = ++tempisItBegin;
				          }
			              }

				    break;
			          }
			        // for UGGrid
			        case 'u':
			          {
				    if (nextisIt == isItEnd)
				      nextisIt = isItBegin;

				    break;
			          }
			      }			    
				    
			    // get geometry type of face 'isIt', i.e., the face between cell1 and cell2 (locally numbered)
			    Dune::GeometryType gtf12 = isIt->geometryInInside().type();
			    
			    // center in face's reference element
			    const Dune::FieldVector<Scalar,dim-1>& 
			      faceLocal12 = Dune::ReferenceElements<Scalar,dim-1>::general(gtf12).position(0,0);

			    // center of face inside volume reference element
			    const LocalPosition&
			      localPosFace12 = Dune::ReferenceElements<Scalar,dim>::general(gtf12).position(isIt->indexInInside(),1);
			      
			    // center of face in global coordinates, i.e., the midpoint of edge 'isIt'
			    GlobalPosition 
			      globalPosFace12 = isIt->geometry().global(faceLocal12);

			    // get face volume
			    double face12vol = isIt->geometry().volume();
			    
			    // get outer normal vector scaled with half volume of face 'isIt'
			    Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln1 
			      = isIt->integrationOuterNormal(faceLocal12);
			    integrationOuterNormaln1
			      *= Dune::ReferenceElements<Scalar,dim-1>::general(gtf12).volume()/2.0;


		            // get geometry type of 'nextisIt', i.e., face between cell1 and cell3 (locally numbered)
			    Dune::GeometryType gtf13 = nextisIt->geometryInInside().type();

			    // center in face's reference element
			    const Dune::FieldVector<Scalar,dim-1>& faceLocal13
                              = Dune::ReferenceElements<Scalar,dim-1>::general(gtf13).position(0,0);

			    // center of face inside volume reference element
			    const LocalPosition&
			      localPosFace13 = Dune::ReferenceElements<Scalar,dim>::general(gtf13).position(nextisIt->indexInInside(),1);

			    // center of face in global coordinates, i.e., the midpoint of edge 'nextisIt'
		            GlobalPosition globalPosFace13
                              = nextisIt->geometry().global(faceLocal13);

			    // get outer normal vector scaled with half volume of face 'nextisIt'		           
			    Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln3
			      = nextisIt->integrationOuterNormal(faceLocal13);
                            integrationOuterNormaln3
                              *= Dune::ReferenceElements<Scalar,dim-1>::general(gtf13).volume()/2.0;

			    // get the intersection node /bar^{x_3} between 'isIt' and 'nextisIt', denoted as 'corner1234'
			    // initialization of corner1234
			    GlobalPosition corner1234(0);	
				  				
			    // get the global coordinate of corner1234
                            for (int i=0; i<isIt->geometry().corners(); ++i)
			      {
			        GlobalPosition isItcorner = isIt->geometry().corner(i);

				for (int j=0; j<nextisIt->geometry().corners(); ++j)
				  {
				    GlobalPosition nextisItcorner = nextisIt->geometry().corner(j);

				    if (nextisItcorner == isItcorner)
				      {
					corner1234 = isItcorner;
					continue;
				      }
				  }
			      }

			    				    			  				
			    // handle interior face
			    if (isIt->neighbor())
			    {
			      // access neighbor cell 2 of 'isIt'
			      ElementPointer outside = isIt->outside();
                              int globalIdx2 = this->diffProblem.variables().indexDiffusion(*outside);
                  
			      // neighbor cell 2 geometry type
                              Dune::GeometryType gt2 = outside->geometry().type();

			      // neighbor cell 2 center in reference element
                              const LocalPosition&
				localPos2 = Dune::ReferenceElements<Scalar,dim>::general(gt2).position(0,0);

			      // get global coordinate of neighbor cell 2 center
                              GlobalPosition 
				globalPos2 = outside->geometry().global(localPos2);

			      // get absolute permeability of neighbor cell 2
                              FieldMatrix K2(this->diffProblem.soil().K(globalPos2, *outside, localPos2));

			      // get total mobility of neighbor cell 2
                              double lambda2;
                              double sat2 = this->diffProblem.variables().saturation()[globalIdx2];
	                      lambda2 = this->diffProblem.materialLaw().mobTotal(sat2, globalPos2, *outside, localPos2);		        

			      // 'nextisIt' is an interior face
			      if (nextisIt->neighbor())
				{
				  // get basic information of cell 1,2's neighbor cell 3,4
				  // neighbor cell 3
				  // access neighbor cell 3
				  ElementPointer nextisItoutside = nextisIt->outside();
				  int globalIdx3 = this->diffProblem.variables().indexDiffusion(*nextisItoutside);

				  // neighbor cell 3 geometry type
				  Dune::GeometryType gt3 = nextisItoutside->geometry().type();

				  // neighbor cell 3 center in reference element
				  const LocalPosition&
				    localPos3 = Dune::ReferenceElements<Scalar,dim>::general(gt3).position(0,0);

				  // get global coordinate of neighbor cell 3 center
				  GlobalPosition
				    globalPos3 = nextisItoutside->geometry().global(localPos3);

				  // get absolute permeability of neighbor cell 3
				  FieldMatrix K3(this->diffProblem.soil().K(globalPos3, *nextisItoutside, localPos3));

				  // get total mobility of neighbor cell 3
				  double lambda3;
                                  double sat3 = this->diffProblem.variables().saturation()[globalIdx3];
	                          lambda3 = this->diffProblem.materialLaw().mobTotal(sat3, globalPos3, *nextisItoutside, localPos3);
				  
				  // neighbor cell 4
				  GlobalPosition globalPos4(0);
				  FieldMatrix K4(0);
				  double lambda4;
				  int globalIdx4;

                                  IntersectionIterator innerisItEnd = this->gridView.template iend(*outside);
                                  IntersectionIterator innernextisItEnd = this->gridView.template iend(*nextisItoutside);
                                  for (IntersectionIterator innerisIt = this->gridView.template ibegin(*outside);
                                       innerisIt!=innerisItEnd; ++innerisIt )
                                    for (IntersectionIterator innernextisIt = this->gridView.template ibegin(*nextisItoutside); 
                                         innernextisIt!=innernextisItEnd; ++innernextisIt)
				      { 
                                        if (innerisIt->neighbor() && innernextisIt->neighbor())
                                          {
					    ElementPointer innerisItoutside = innerisIt->outside();
					    ElementPointer innernextisItoutside = innernextisIt->outside();
                                        
					    // find the common neighbor cell between cell 2 and cell 3, except cell 1
					    if (innerisItoutside == innernextisItoutside && innerisItoutside != isIt->inside())
					      {
					        // access neighbor cell 4
					        globalIdx4 = this->diffProblem.variables().indexDiffusion(*innerisItoutside);

					        // neighbor cell 4 geometry type
					        Dune::GeometryType gt4 = innerisItoutside->geometry().type();

					        // neighbor cell 4 center in reference element
					        const LocalPosition&
					        localPos4 = Dune::ReferenceElements<Scalar,dim>::general(gt4).position(0,0);

					        // get global coordinate of neighbor cell 4 center
					        globalPos4 = innerisItoutside->geometry().global(localPos4);

					        // get absolute permeability of neighbor cell 4
					        K4 += this->diffProblem.soil().K(globalPos4, *innerisItoutside, localPos4);

					        // get total mobility of neighbor cell 4
                                                double sat4 = this->diffProblem.variables().saturation()[globalIdx4];
	                                        lambda4 = this->diffProblem.materialLaw().mobTotal(sat4, globalPos4, *innerisItoutside, localPos4);	
                                              }
					  }
				      }

				  // computation of flux through the first half edge of 'isIt' and the flux
                                  // through the second half edge of 'nextisIt'
				    						
				  // get the information of the face 'isIt24' between cell2 and cell4 (locally numbered)
				  IntersectionIterator isIt24 = this->gridView.template ibegin(*outside);

				  for (IntersectionIterator innerisIt = this->gridView.template ibegin(*outside);
                                       innerisIt != innerisItEnd; ++innerisIt)
				    {
				      if (innerisIt->neighbor())
					{
					  if (innerisIt->outside() != isIt->inside())
					    {
					      for (int i=0; i<innerisIt->geometry().corners(); ++i)
						{
						  GlobalPosition innerisItcorner = innerisIt->geometry().corner(i);

						  if (innerisItcorner == corner1234)
						    {
						      isIt24 = innerisIt;
						      continue;
						    }
						}
					    }
					}				      
				    }

				  // get geometry type of face 'isIt24'
				  Dune::GeometryType gtf24 = isIt24->geometryInInside().type();

				  // center in face's reference element
				  const Dune::FieldVector<Scalar,dim-1>&
				    faceLocal24 = Dune::ReferenceElements<Scalar,dim-1>::general(gtf24).position(0,0);

				  // center of face in global coordinates, i.e., the midpoint of edge 'isIt24'
				  GlobalPosition
				    globalPosFace24 = isIt24->geometry().global(faceLocal24);

				  // get outer normal vector scaled with half volume of face 'isIt24'
				  Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln4
				    = isIt24->integrationOuterNormal(faceLocal24);
                                  integrationOuterNormaln4
                                    *= Dune::ReferenceElements<Scalar,dim-1>::general(gtf24).volume()/2.0;

				  // get the information of the face 'isIt34' between cell3 and cell4 (locally numbered)
				  IntersectionIterator isIt34 = this->gridView.template ibegin(*nextisItoutside);

				  for (IntersectionIterator innerisIt = this->gridView.template ibegin(*nextisItoutside);
                                       innerisIt != innernextisItEnd; ++innerisIt)
				    {
				      if (innerisIt->neighbor())
					{
					  if (innerisIt->outside() != isIt->inside())
					    {
					      for (int i=0; i<innerisIt->geometry().corners(); ++i)
						{
						  GlobalPosition innerisItcorner = innerisIt->geometry().corner(i);

						  if (innerisItcorner == corner1234)
						    {
						      isIt34 = innerisIt;
						      continue;
						    }
						}
					    }
					}				      
				    }

				  // get geometry type of face 'isIt34'
				  Dune::GeometryType gtf34 = isIt34->geometryInInside().type();

				  // center in face's reference element
				  const Dune::FieldVector<Scalar,dim-1>&
				    faceLocal34 = Dune::ReferenceElements<Scalar,dim-1>::general(gtf34).position(0,0);

				  // center of face in global coordinates, i.e., the midpoint of edge 'isIt34'
				  GlobalPosition
				    globalPosFace34 = isIt34->geometry().global(faceLocal34);

				  // get outer normal vector scaled with half volume of face 'isIt34'
				  Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln2
				    = isIt34->integrationOuterNormal(faceLocal34);
                                  integrationOuterNormaln2
                                    *= Dune::ReferenceElements<Scalar,dim-1>::general(gtf34).volume()/2.0;

				  // compute normal vectors nu11,nu21; nu12, nu22; nu13, nu23; nu14, nu24;
                                  FieldVector<Scalar,dim> nu11(0);
				  R.umv(globalPosFace13-globalPos1 ,nu11);

				  FieldVector<Scalar,dim> nu21(0);
				  R.umv(globalPos1-globalPosFace12, nu21);

				  FieldVector<Scalar,dim> nu12(0);
				  R.umv(globalPosFace24-globalPos2, nu12);

				  FieldVector<Scalar,dim> nu22(0);
				  R.umv(globalPosFace12-globalPos2, nu22);

				  FieldVector<Scalar,dim> nu13(0);
				  R.umv(globalPos3-globalPosFace13, nu13);

				  FieldVector<Scalar,dim> nu23(0);
				  R.umv(globalPos3-globalPosFace34, nu23);

				  FieldVector<Scalar,dim> nu14(0);
				  R.umv(globalPos4-globalPosFace24, nu14);

				  FieldVector<Scalar,dim> nu24(0);
				  R.umv(globalPosFace34-globalPos4, nu24);

				  // compute dF1, dF2, dF3, dF4 i.e., the area of quadrilateral made by normal vectors 'nu'
				  FieldVector<Scalar,dim> Rnu21(0);
				  R.umv(nu21, Rnu21);
				  double dF1 = fabs(nu11 * Rnu21);
				  
				  FieldVector<Scalar,dim> Rnu22(0);
				  R.umv(nu22, Rnu22);				  
				  double dF2 = fabs(nu12 * Rnu22);
				  
				  FieldVector<Scalar,dim> Rnu23(0);
				  R.umv(nu23, Rnu23);				  
				  double dF3 = fabs(nu13 * Rnu23);
				  			  
				  FieldVector<Scalar,dim> Rnu24(0);
				  R.umv(nu24, Rnu24);				  
				  double dF4 = fabs(nu14 * Rnu24);

				  // compute components needed for flux calculation, denoted as 'g'
				  FieldVector<Scalar,dim> K1nu11(0);
				  K1.umv(nu11, K1nu11);
				  FieldVector<Scalar,dim> K1nu21(0);
				  K1.umv(nu21, K1nu21);
				  FieldVector<Scalar,dim> K2nu12(0);
				  K2.umv(nu12, K2nu12);
				  FieldVector<Scalar,dim> K2nu22(0);
				  K2.umv(nu22, K2nu22);
				  FieldVector<Scalar,dim> K3nu13(0);
				  K3.umv(nu13, K3nu13);
				  FieldVector<Scalar,dim> K3nu23(0);
				  K3.umv(nu23, K3nu23);
				  FieldVector<Scalar,dim> K4nu14(0);
				  K4.umv(nu14, K4nu14);
				  FieldVector<Scalar,dim> K4nu24(0);
				  K4.umv(nu24, K4nu24);				 
                                  double g111 = lambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
				  double g121 = lambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
				  double g211 = lambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
				  double g221 = lambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                                  double g112 = lambda2 * (integrationOuterNormaln1 * K2nu12)/dF2;
				  double g122 = lambda2 * (integrationOuterNormaln1 * K2nu22)/dF2;
				  double g212 = lambda2 * (integrationOuterNormaln4 * K2nu12)/dF2;
				  double g222 = lambda2 * (integrationOuterNormaln4 * K2nu22)/dF2;
                                  double g113 = lambda3 * (integrationOuterNormaln2 * K3nu13)/dF3;
				  double g123 = lambda3 * (integrationOuterNormaln2 * K3nu23)/dF3;
				  double g213 = lambda3 * (integrationOuterNormaln3 * K3nu13)/dF3;
				  double g223 = lambda3 * (integrationOuterNormaln3 * K3nu23)/dF3;
                                  double g114 = lambda4 * (integrationOuterNormaln2 * K4nu14)/dF4;
				  double g124 = lambda4 * (integrationOuterNormaln2 * K4nu24)/dF4;
				  double g214 = lambda4 * (integrationOuterNormaln4 * K4nu14)/dF4;
				  double g224 = lambda4 * (integrationOuterNormaln4 * K4nu24)/dF4;
				  

				  // compute transmissibility matrix T = CA^{-1}B+F
				  Dune::FieldMatrix<Scalar,2*dim,2*dim> C(0), F(0), A(0), B(0);

				  // evaluate matrix C, F, A, B
				  C[0][0] = -g111;
				  C[0][2] = -g121;
				  C[1][1] = g114;
				  C[1][3] = g124;
				  C[2][1] = -g213;
				  C[2][2] = g223;
				  C[3][0] = g212;
				  C[3][3] = -g222;

				  F[0][0] = g111 + g121;
				  F[1][3] = -g114 - g124;
				  F[2][2] = g213 - g223;
				  F[3][1] = -g212 + g222;

				  A[0][0] = g111 + g112;
				  A[0][2] = g121;
				  A[0][3] = -g122;
				  A[1][1] = g114 + g113;
				  A[1][2] = -g123;
				  A[1][3] = g124;
				  A[2][0] = g211;
				  A[2][1] = -g213;
				  A[2][2] = g223 + g221;
				  A[3][0] = -g212;
				  A[3][1] = g214;
				  A[3][3] = g222 + g224;


				  B[0][0] = g111 + g121;
				  B[0][1] = g112 - g122;
				  B[1][2] = g113 - g123;
				  B[1][3] = g114 + g124;
				  B[2][0] = g211 + g221;
				  B[2][2] = -g213 + g223;
				  B[3][1] = -g212 + g222;
				  B[3][3] = g214 + g224;

				  // compute T
				  A.invert();
                                  F += B.leftmultiply(C.rightmultiply(A));
                                  Dune::FieldMatrix<Scalar,2*dim,2*dim> T(F);
	
				  // assemble the global matrix M and right hand side f			 
				  M[globalIdx1][globalIdx1] += T[0][0] + T[2][0];
				  M[globalIdx1][globalIdx2] += T[0][1] + T[2][1];
				  M[globalIdx1][globalIdx3] += T[0][2] + T[2][2];
				  M[globalIdx1][globalIdx4] += T[0][3] + T[2][3];
				 
				}
			      // 'nextisIt' is on the boundary			      
			      else
			      	{  
				  // computation of flux through the first half edge of 'isIt' and the flux
                                  // through the second half edge of 'nextisIt'
				   
                                  // get common geometry information for the following computation 				
				  // get the information of the face 'isIt24' between cell2 and cell4 (locally numbered)
				  IntersectionIterator isIt24 = this->gridView.template ibegin(*outside);
                                  IntersectionIterator innerisItEnd = this->gridView.template iend(*outside);
				  for (IntersectionIterator innerisIt = this->gridView.template ibegin(*outside);
                                       innerisIt != innerisItEnd; ++innerisIt)
				    {
				      if (innerisIt->boundary())
				    	{
					  for (int i=0; i<innerisIt->geometry().corners(); ++i)
					    {
					      GlobalPosition innerisItcorner = innerisIt->geometry().corner(i);

					      if (innerisItcorner == corner1234)
					        {
						  isIt24 = innerisIt;
						  continue;
					        }
					    }					     
					}				      
				    }

                                  // get geometry type of face 'isIt24'
				  Dune::GeometryType gtf24 = isIt24->geometryInInside().type();

				  // center in face's reference element
				  const Dune::FieldVector<Scalar,dim-1>&
				    faceLocal24 = Dune::ReferenceElements<Scalar,dim-1>::general(gtf24).position(0,0);

				  // center of face inside volume reference element
				  const LocalPosition&
				    localPosFace24 = Dune::ReferenceElements<Scalar,dim>::general(gtf24).position(isIt24->indexInInside(),1);

				  // center of face in global coordinates, i.e., the midpoint of edge 'isIt24'
				  GlobalPosition
				    globalPosFace24 = isIt24->geometry().global(faceLocal24);

				  // get outer normal vector scaled with half volume of face 'isIt24'
				  Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln4
				    = isIt24->integrationOuterNormal(faceLocal24);
                                  integrationOuterNormaln4
                                    *= Dune::ReferenceElements<Scalar,dim-1>::general(gtf24).volume()/2.0;

				  // get boundary condition for boundary face (nextisIt) center
				  BoundaryConditions::Flags nextisItbctype = this->diffProblem.bctypePress(globalPosFace13, *eIt, localPosFace13);

				  // 'nextisIt': Neumann boundary
				  if (nextisItbctype == BoundaryConditions::neumann)
				    {
				      // get Neumann boundary value of 'nextisIt'
				      double J3 = this->diffProblem.neumannPress(globalPosFace13, *eIt, localPosFace13);

				      // get boundary condition for boundary face (isIt24) center
				      BoundaryConditions::Flags isIt24bctype = 
                                        this->diffProblem.bctypePress(globalPosFace24, *outside, localPosFace24);  

                                      // 'isIt24': Neumann boundary
                                      if (isIt24bctype == BoundaryConditions::neumann)
                                        {
                                          // get neumann boundary value of 'isIt24'
				          double J4 = this->diffProblem.neumannPress(globalPosFace24, *outside, localPosFace24);

				          // compute normal vectors nu11,nu21; nu12, nu22;
                                          FieldVector<Scalar,dim> nu11(0);
				          R.umv(globalPosFace13-globalPos1 ,nu11);

				          FieldVector<Scalar,dim> nu21(0);
				          R.umv(globalPos1-globalPosFace12, nu21);

				          FieldVector<Scalar,dim> nu12(0);
				          R.umv(globalPosFace24-globalPos2, nu12);

				          FieldVector<Scalar,dim> nu22(0);
				          R.umv(globalPosFace12-globalPos2, nu22);
				  
				          // compute dF1, dF2 i.e., the area of quadrilateral made by normal vectors 'nu'
				          FieldVector<Scalar,dim> Rnu21(0);
				          R.umv(nu21, Rnu21);
				          double dF1 = fabs(nu11 * Rnu21);
				  
				          FieldVector<Scalar,dim> Rnu22(0);
				          R.umv(nu22, Rnu22);				  
				          double dF2 = fabs(nu12 * Rnu22);
				  
				          // compute components needed for flux calculation, denoted as 'g'
				          FieldVector<Scalar,dim> K1nu11(0);
				          K1.umv(nu11, K1nu11);
				          FieldVector<Scalar,dim> K1nu21(0);
				          K1.umv(nu21, K1nu21);
				          FieldVector<Scalar,dim> K2nu12(0);
				          K2.umv(nu12, K2nu12);
				          FieldVector<Scalar,dim> K2nu22(0);
				          K2.umv(nu22, K2nu22);	 
                                          double g111 = lambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
				          double g121 = lambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
				          double g211 = lambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
				          double g221 = lambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                                          double g112 = lambda2 * (integrationOuterNormaln1 * K2nu12)/dF2;
				          double g122 = lambda2 * (integrationOuterNormaln1 * K2nu22)/dF2;
				          double g212 = lambda2 * (integrationOuterNormaln4 * K2nu12)/dF2;
				          double g222 = lambda2 * (integrationOuterNormaln4 * K2nu22)/dF2;
                                  
				          // compute the matrix T & vector r in v = A^{-1}(Bu + r1) = Tu + r
				          Dune::FieldMatrix<Scalar,2*dim-1,2*dim-1> A(0);
                                          Dune::FieldMatrix<Scalar,2*dim-1,dim> B(0);
                                          Dune::FieldVector<Scalar,2*dim-1> r1(0), r(0);

				          // evaluate matrix A, B
				          A[0][0] = g111 + g112;
				          A[0][1] = g121;
				          A[0][2] = -g122;
				          A[1][0] = g211;
				          A[1][1] = g221;
				          A[2][0] = -g212;
				          A[2][2] = g222;
				 
				          B[0][0] = g111 + g121;
				          B[0][1] = g112 - g122;
				          B[1][0] = g211 + g221;
				          B[2][1] = g222 - g212;
 
                                          // evaluate vector r1
                                          r1[1] = -J3 * nextisIt->geometry().volume()/2.0;
                                          r1[2] = -J4 * isIt24->geometry().volume()/2.0;

                                          // compute T and r
				          A.invert();
                                          B.leftmultiply(A);      
                                          Dune::FieldMatrix<Scalar,2*dim-1,dim> T(B);
					  A.umv(r1, r);                                         
                                          			 			  
				          // assemble the global matrix M and right hand side f			 
				          M[globalIdx1][globalIdx1] += g111 + g121 - g111 * T[0][0] - g121 * T[1][0];
				          M[globalIdx1][globalIdx2] += -g111 * T[0][1] - g121 * T[1][1];
                                          f[globalIdx1] += g111 * r[0] + g121 * r[1];
                                       
                                        }
                                      // 'isIt24': Dirichlet boundary
                                      else
                                        {
                                          // get Dirichlet boundary value on 'isIt24'
				          double g4 = this->diffProblem.dirichletPress(globalPosFace24, *outside, localPosFace24);

				          // compute total mobility for Dirichlet boundary 'isIt24'
				          double alambda2;
					  alambda2 = lambda2;

				          // compute normal vectors nu11,nu21; nu12, nu22;
                                          FieldVector<Scalar,dim> nu11(0);
				          R.umv(globalPosFace13-globalPos1 ,nu11);

				          FieldVector<Scalar,dim> nu21(0);
				          R.umv(globalPos1-globalPosFace12, nu21);

				          FieldVector<Scalar,dim> nu12(0);
				          R.umv(globalPosFace24-globalPos2, nu12);

				          FieldVector<Scalar,dim> nu22(0);
				          R.umv(globalPosFace12-globalPos2, nu22);

				          // compute dF1, dF2 i.e., the area of quadrilateral made by normal vectors 'nu'
				          FieldVector<Scalar,dim> Rnu21(0);
				          R.umv(nu21, Rnu21);
				          double dF1 = fabs(nu11 * Rnu21);
				  
				          FieldVector<Scalar,dim> Rnu22(0);
				          R.umv(nu22, Rnu22);				  
				          double dF2 = fabs(nu12 * Rnu22);
				  
				          // compute components needed for flux calculation, denoted as 'g'
				          FieldVector<Scalar,dim> K1nu11(0);
				          K1.umv(nu11, K1nu11);
				          FieldVector<Scalar,dim> K1nu21(0);
				          K1.umv(nu21, K1nu21);
				          FieldVector<Scalar,dim> K2nu12(0);
				          K2.umv(nu12, K2nu12);
				          FieldVector<Scalar,dim> K2nu22(0);
				          K2.umv(nu22, K2nu22);	 
                                          double g111 = lambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
				          double g121 = lambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
				          double g211 = lambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
				          double g221 = lambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                                          double g112 = alambda2 * (integrationOuterNormaln1 * K2nu12)/dF2;
				          double g122 = alambda2 * (integrationOuterNormaln1 * K2nu22)/dF2;
                                  
				          // compute the matrix T & vector r in v = A^{-1}(Bu + r1) = Tu + r
				          FieldMatrix A(0), B(0); 
                                          Dune::FieldVector<Scalar,dim> r1(0), r(0);

				          // evaluate matrix A, B
				          A[0][0] = g111 + g112;
				          A[0][1] = g121;		
				          A[1][0] = g211;
				          A[1][1] = g221;
				         
				          B[0][0] = g111 + g121;
				          B[0][1] = g112 - g122;
				          B[1][0] = g211 + g221;
				         
                                          // evaluate vector r1
                                          r1[0] = g122 * g4;
                                          r1[1] = -J3 * nextisIt->geometry().volume()/2.0;
                                          
                                          // compute T and r
				          A.invert();
                                          B.leftmultiply(A);
                                          FieldMatrix T(B);
					  A.umv(r1, r);                                         
                                          			 			  
				          // assemble the global matrix M and right hand side f			 
				          M[globalIdx1][globalIdx1] += g111 + g121 - g111 * T[0][0] - g121 * T[1][0];
				          M[globalIdx1][globalIdx2] += -g111 * T[0][1] - g121 * T[1][1];
                                          f[globalIdx1] += g111 * r[0] + g121 * r[1];
                                                                                                                          
                                        }                                       				     
				    }
				  // 'nextisIt': Dirichlet boundary 
				  else
				    {
				      // get Dirichlet boundary value of 'nextisIt'
				      double g3 = this->diffProblem.dirichletPress(globalPosFace13, *eIt, localPosFace13);

				      // compute total mobility for Dirichlet boundary 'nextisIt'
				      double alambda1;
                                      alambda1 = lambda1;
				    
				      // get boundary condition for boundary face (isIt24) center
				      BoundaryConditions::Flags isIt24bctype = 
                                        this->diffProblem.bctypePress(globalPosFace24, *outside, localPosFace24);
                       
                                      // 'isIt24': Neumann boundary
                                      if (isIt24bctype == BoundaryConditions::neumann)
                                        {
                                          // get Neumann boundary value of 'isIt24'
				          double J4 = this->diffProblem.neumannPress(globalPosFace24, *outside, localPosFace24);

				          // compute normal vectors nu11,nu21; nu12, nu22;
                                          FieldVector<Scalar,dim> nu11(0);
				          R.umv(globalPosFace13-globalPos1 ,nu11);

				          FieldVector<Scalar,dim> nu21(0);
				          R.umv(globalPos1-globalPosFace12, nu21);

				          FieldVector<Scalar,dim> nu12(0);
				          R.umv(globalPosFace24-globalPos2, nu12);

				          FieldVector<Scalar,dim> nu22(0);
				          R.umv(globalPosFace12-globalPos2, nu22);

				          // compute dF1, dF2 i.e., the area of quadrilateral made by normal vectors 'nu'
				          FieldVector<Scalar,dim> Rnu21(0);
				          R.umv(nu21, Rnu21);
				          double dF1 = fabs(nu11 * Rnu21);
				  
				          FieldVector<Scalar,dim> Rnu22(0);
				          R.umv(nu22, Rnu22);				  
				          double dF2 = fabs(nu12 * Rnu22);
				  
				          // compute components needed for flux calculation, denoted as 'g'
				          FieldVector<Scalar,dim> K1nu11(0);
				          K1.umv(nu11, K1nu11);
				          FieldVector<Scalar,dim> K1nu21(0);
				          K1.umv(nu21, K1nu21);
				          FieldVector<Scalar,dim> K2nu12(0);
				          K2.umv(nu12, K2nu12);
				          FieldVector<Scalar,dim> K2nu22(0);
				          K2.umv(nu22, K2nu22);	 
                                          double g111 = alambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
				          double g121 = alambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
				          double g211 = alambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
				          double g221 = alambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                                          double g112 = lambda2 * (integrationOuterNormaln1 * K2nu12)/dF2;
				          double g122 = lambda2 * (integrationOuterNormaln1 * K2nu22)/dF2;
				          double g212 = lambda2 * (integrationOuterNormaln4 * K2nu12)/dF2;
				          double g222 = lambda2 * (integrationOuterNormaln4 * K2nu22)/dF2;
                                  
				          // compute the matrix T & vector r in v = A^{-1}(Bu + r1) = Tu + r
				          FieldMatrix A(0), B(0); 
                                          Dune::FieldVector<Scalar,dim> r1(0), r(0);

				          // evaluate matrix A, B
				          A[0][0] = g111 + g112;
				          A[0][1] = -g122;		
				          A[1][0] = -g212;
				          A[1][1] = g222;
				         
				          B[0][0] = g111 + g121;
				          B[0][1] = g112 - g122;
				          B[1][1] = g222 - g212;
				         
                                          // evaluate vector r1
                                          r1[0] = -g121 * g3;
                                          r1[1] = -J4 * isIt24->geometry().volume()/2.0;
                                          
                                          // compute T and r
				          A.invert();
                                          B.leftmultiply(A);
                                          FieldMatrix T(B);
					  A.umv(r1, r);                                         
                                          			 			  
				          // assemble the global matrix M and right hand side f			 
				          M[globalIdx1][globalIdx1] += (g111 + g121 - g111 * T[0][0]) + (g211 + g221 - g211 * T[0][0]);
				          M[globalIdx1][globalIdx2] += -g111 * T[0][1] - g211 * T[0][1];
                                          f[globalIdx1] += (g121 + g221) * g3 + (g111 + g211) * r[0];

                                        }
                                      // 'isIt24': Dirichlet boundary
                                      else
                                        {
                                          // get Dirichlet boundary value on 'isIt24'
				          double g4 = this->diffProblem.dirichletPress(globalPosFace24, *outside, localPosFace24);

				          // compute total mobility for Dirichlet boundary 'isIt24'
				          double alambda2;
                                          alambda2 = lambda2;

				          // compute normal vectors nu11,nu21; nu12, nu22;
                                          FieldVector<Scalar,dim> nu11(0);
				          R.umv(globalPosFace13-globalPos1 ,nu11);

				          FieldVector<Scalar,dim> nu21(0);
				          R.umv(globalPos1-globalPosFace12, nu21);

				          FieldVector<Scalar,dim> nu12(0);
				          R.umv(globalPosFace24-globalPos2, nu12);

				          FieldVector<Scalar,dim> nu22(0);
				          R.umv(globalPosFace12-globalPos2, nu22);

				          // compute dF1, dF2 i.e., the area of quadrilateral made by normal vectors 'nu'
				          FieldVector<Scalar,dim> Rnu21(0);
				          R.umv(nu21, Rnu21);
				          double dF1 = fabs(nu11 * Rnu21);
				  
				          FieldVector<Scalar,dim> Rnu22(0);
				          R.umv(nu22, Rnu22);				  
				          double dF2 = fabs(nu12 * Rnu22);
				  
				          // compute components needed for flux calculation, denoted as 'g'
				          FieldVector<Scalar,dim> K1nu11(0);
				          K1.umv(nu11, K1nu11);
				          FieldVector<Scalar,dim> K1nu21(0);
				          K1.umv(nu21, K1nu21);
				          FieldVector<Scalar,dim> K2nu12(0);
				          K2.umv(nu12, K2nu12);
				          FieldVector<Scalar,dim> K2nu22(0);
				          K2.umv(nu22, K2nu22);	 
                                          double g111 = alambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
				          double g121 = alambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
				          double g211 = alambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
				          double g221 = alambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                                          double g112 = alambda2 * (integrationOuterNormaln1 * K2nu12)/dF2;
				          double g122 = alambda2 * (integrationOuterNormaln1 * K2nu22)/dF2;
				          
				          // compute the matrix T & vector r 
				          FieldMatrix T(0); 
                                          Dune::FieldVector<Scalar,dim> r(0);
                                          
                                          double coe = g111 + g112; 

				          // evaluate matrix T
				          T[0][0] = g112 * (g111 + g121)/coe;
				          T[0][1] = -g111 * (g112 - g122)/coe;		
				          T[1][0] = g221 + g211 * (g112 - g121)/coe;
				          T[1][1] = -g211 * (g112 - g122)/coe;
				         				         			         
                                          // evaluate vector r
                                          r[0] = -(g4 * g122 * g111 + g3 * g112 * g121)/coe;
                                          r[1] = -g221 * g3 + (g3 * g211 * g121 - g4 * g211 * g122)/coe;                              			  
				          // assemble the global matrix M and right hand side f			 
				          M[globalIdx1][globalIdx1] += T[0][0] + T[1][0];
				          M[globalIdx1][globalIdx2] += T[0][1] + T[1][1];
                                          f[globalIdx1] -= r[0] + r[1];
              
                                        }	
                                    }
                                }			    
                            }
			    // handle boundary face 'isIt'
			    else
			    {
			      // get boundary condition for boundary face center of 'isIt'
			      BoundaryConditions::Flags isItbctype = this->diffProblem.bctypePress(globalPosFace12, *eIt, localPosFace12);
                              
			      // 'isIt' is on Neumann boundary
			      if (isItbctype == BoundaryConditions::neumann)
				{
				  // get Neumann boundary value
				  double J1 = this->diffProblem.neumannPress(globalPosFace12, *eIt, localPosFace12);
				  
                                  // evaluate right hand side
				  f[globalIdx1] -= face12vol*J1;

                                  // 'nextisIt' is on boundary
				  if (nextisIt->boundary())
				    {
				      // get boundary condition for boundary face center of 'nextisIt'
                                      BoundaryConditions::Flags nextisItbctype =
					this->diffProblem.bctypePress(globalPosFace13, *eIt, localPosFace13);

				      if (nextisItbctype == BoundaryConditions::dirichlet)
					{
				          // compute total mobility for Dirichlet boundary
				          double alambda1;
                                          alambda1 = lambda1;

					  // get Dirichlet boundary value
					  double g3 = this->diffProblem.dirichletPress(globalPosFace13, *eIt, localPosFace13);

				          // compute normal vectors nu11,nu21;
                                          FieldVector<Scalar,dim> nu11(0);
				          R.umv(globalPosFace13-globalPos1 ,nu11);

				          FieldVector<Scalar,dim> nu21(0);
				          R.umv(globalPos1-globalPosFace12, nu21);
	
				          // compute dF1, dF2 i.e., the area of quadrilateral made by normal vectors 'nu'
				          FieldVector<Scalar,dim> Rnu21(0);
				          R.umv(nu21, Rnu21);
				          double dF1 = fabs(nu11 * Rnu21);
				  				          				  
				          // compute components needed for flux calculation, denoted as 'g'
				          FieldVector<Scalar,dim> K1nu11(0);
				          K1.umv(nu11, K1nu11);
				          FieldVector<Scalar,dim> K1nu21(0);
				          K1.umv(nu21, K1nu21);				          
                                          double g111 = alambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
				          double g121 = alambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
				          double g211 = alambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
				          double g221 = alambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;
                                          				         
				          // assemble the global matrix M and right hand side f			 
				          M[globalIdx1][globalIdx1] += g221 - g211 * g121/g111;			         
                                          f[globalIdx1] -= (g211 * g121/g111 - g221) * g3 - (g211 * (-J1) * face12vol)/(2.0 * g111);

					}
				    }
                                  // 'nextisIt' is inside
                                  else
                                    {
				      // neighbor cell 3
				      // access neighbor cell 3
				      ElementPointer nextisItoutside = nextisIt->outside();
				      int globalIdx3 = this->diffProblem.variables().indexDiffusion(*nextisItoutside);

				      // neighbor cell 3 geometry type
				      Dune::GeometryType gt3 = nextisItoutside->geometry().type();

				      // neighbor cell 3 center in reference element
				      const LocalPosition&
				        localPos3 = Dune::ReferenceElements<Scalar,dim>::general(gt3).position(0,0);

				      // get global coordinate of neighbor cell 3 center
				      GlobalPosition
				        globalPos3 = nextisItoutside->geometry().global(localPos3);

				      // get absolute permeability of neighbor cell 3
				      FieldMatrix K3(this->diffProblem.soil().K(globalPos3, *nextisItoutside, localPos3));

				      // get total mobility of neighbor cell 3
				      double lambda3;
                                      double sat3 = this->diffProblem.variables().saturation()[globalIdx3];     
	                              lambda3 = this->diffProblem.materialLaw().mobTotal(sat3, globalPos3, *nextisItoutside, localPos3);
				      

				      // get the information of the face 'isIt34' between cell3 and cell4 (locally numbered)
				      IntersectionIterator isIt34 = this->gridView.template ibegin(*nextisItoutside);
                                      IntersectionIterator innernextisItEnd = this->gridView.template iend(*nextisItoutside);
				      for (IntersectionIterator innerisIt = this->gridView.template ibegin(*nextisItoutside);
                                           innerisIt != innernextisItEnd; ++innerisIt)
				        {
				          if (innerisIt->boundary())
				    	    {
					      for (int i=0; i<innerisIt->geometry().corners(); ++i)
					        {
					          GlobalPosition innerisItcorner = innerisIt->geometry().corner(i);

					          if (innerisItcorner == corner1234)
					            {
						      isIt34 = innerisIt;
						      continue;
					            }
					        }					     
					    }				      
				        }

				      // get geometry type of face 'isIt34'
				      Dune::GeometryType gtf34 = isIt34->geometryInInside().type();

				      // center in face's reference element
				      const Dune::FieldVector<Scalar,dim-1>&
				        faceLocal34 = Dune::ReferenceElements<Scalar,dim-1>::general(gtf34).position(0,0);

			              // center of face inside volume reference element
			              const LocalPosition&
			                localPosFace34 = Dune::ReferenceElements<Scalar,dim>::general(gtf34).position(isIt34->indexInInside(),1);

				      // center of face in global coordinates, i.e., the midpoint of edge 'isIt34'
				      GlobalPosition
				        globalPosFace34 = isIt34->geometry().global(faceLocal34);

				      // get outer normal vector scaled with half volume of face 'isIt34'
				      Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln2
				        = isIt34->integrationOuterNormal(faceLocal34);
                                      integrationOuterNormaln2
                                        *= Dune::ReferenceElements<Scalar,dim-1>::general(gtf34).volume()/2.0;

				      // get boundary condition for boundary face center of 'isIt34'
                                      BoundaryConditions::Flags isIt34bctype =
					this->diffProblem.bctypePress(globalPosFace34, *nextisItoutside, localPosFace34);

                                      // 'isIt34': Neumann boundary
                                      if (isIt34bctype == BoundaryConditions::neumann)
                                        {
				          // get Neumann boundary value
				          double J2 = this->diffProblem.neumannPress(globalPosFace34, *nextisItoutside, localPosFace34);    

				          // compute normal vectors nu11,nu21; nu13, nu23;
                                          FieldVector<Scalar,dim> nu11(0);
				          R.umv(globalPosFace13-globalPos1 ,nu11);

				          FieldVector<Scalar,dim> nu21(0);
				          R.umv(globalPos1-globalPosFace12, nu21);

				          FieldVector<Scalar,dim> nu13(0);
				          R.umv(globalPos3-globalPosFace13, nu13);

				          FieldVector<Scalar,dim> nu23(0);
				          R.umv(globalPos3-globalPosFace34, nu23);
	
				          // compute dF1, dF3 i.e., the area of quadrilateral made by normal vectors 'nu'
				          FieldVector<Scalar,dim> Rnu21(0);
				          R.umv(nu21, Rnu21);
				          double dF1 = fabs(nu11 * Rnu21);
				  				  
				          FieldVector<Scalar,dim> Rnu23(0);
				          R.umv(nu23, Rnu23);				  
				          double dF3 = fabs(nu13 * Rnu23);
				  
				          // compute components needed for flux calculation, denoted as 'g'
				          FieldVector<Scalar,dim> K1nu11(0);
				          K1.umv(nu11, K1nu11);
				          FieldVector<Scalar,dim> K1nu21(0);
				          K1.umv(nu21, K1nu21);				  
				          FieldVector<Scalar,dim> K3nu13(0);
				          K3.umv(nu13, K3nu13);
				          FieldVector<Scalar,dim> K3nu23(0);
				          K3.umv(nu23, K3nu23);				  				 
                                          double g111 = lambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
				          double g121 = lambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
				          double g211 = lambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
				          double g221 = lambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;                               
                                          double g113 = lambda3 * (integrationOuterNormaln2 * K3nu13)/dF3;
				          double g123 = lambda3 * (integrationOuterNormaln2 * K3nu23)/dF3;
				          double g213 = lambda3 * (integrationOuterNormaln3 * K3nu13)/dF3;
				          double g223 = lambda3 * (integrationOuterNormaln3 * K3nu23)/dF3;
                                  
				          // compute transmissibility matrix T = CA^{-1}B+F
				          Dune::FieldMatrix<Scalar,2*dim-1,2*dim-1> C(0), A(0);
                                          Dune::FieldMatrix<Scalar,2*dim-1,dim> F(0), B(0);

				          // evaluate matrix C, F, A, B
				          C[0][0] = -g111;
				          C[0][2] = -g121;
				          C[1][1] = -g113;
				          C[1][2] = g123;
				          C[2][1] = -g213;
				          C[2][2] = g223;

				          F[0][0] = g111 + g121;
				          F[1][1] = g113 - g123;
				          F[2][1] = g213 - g223;

				          A[0][0] = g111;
				          A[0][2] = g121;				  
				          A[1][1] = g113;
				          A[1][2] = -g123;				 
				          A[2][0] = g211;
				          A[2][1] = -g213;
				          A[2][2] = g223 + g221;
				  
				          B[0][0] = g111 + g121;
				          B[1][1] = g113 - g123;
				          B[2][0] = g211 + g221;
				          B[2][1] = g223 - g213;
				 
				          // compute T
				          A.invert();
                                          Dune::FieldMatrix<Scalar,2*dim-1,2*dim-1> CAinv(C.rightmultiply(A));
                                          F += B.leftmultiply(CAinv);
                                          Dune::FieldMatrix<Scalar,2*dim-1,dim> T(F);

                                          // compute vector r
			                  // evaluate r1
			                  Dune::FieldVector<Scalar,2*dim-1> r1(0);
			                  r1[0] = -J1 * face12vol/2.0;
                                          r1[1] = -J2 * isIt34->geometry().volume()/2.0;

				          // compute  r = CA^{-1}r1
				          Dune::FieldVector<Scalar,2*dim-1> r(0);
				          CAinv.umv(r1, r);
					 			 			    
				          // assemble the global matrix M and right hand side f			 
				          M[globalIdx1][globalIdx1] += T[2][0];
				          M[globalIdx1][globalIdx3] += T[2][1];
                                          f[globalIdx1] -= r[2];
 
                                        }
                                      // 'isIt34': Dirichlet boundary
                                      else
                                        {
					  // get Dirichlet boundary value
					  double g2 = this->diffProblem.dirichletPress(globalPosFace34, *nextisItoutside, localPosFace34);

				          // compute total mobility for Dirichlet boundary
				          double alambda3;
                                          alambda3 = lambda3;

				          // compute normal vectors nu11,nu21; nu13, nu23;
                                          FieldVector<Scalar,dim> nu11(0);
				          R.umv(globalPosFace13-globalPos1 ,nu11);

				          FieldVector<Scalar,dim> nu21(0);
				          R.umv(globalPos1-globalPosFace12, nu21);

				          FieldVector<Scalar,dim> nu13(0);
				          R.umv(globalPos3-globalPosFace13, nu13);

				          FieldVector<Scalar,dim> nu23(0);
				          R.umv(globalPos3-globalPosFace34, nu23);
				 
				          // compute dF1, dF3 i.e., the area of quadrilateral made by normal vectors 'nu'
				          FieldVector<Scalar,dim> Rnu21(0);
				          R.umv(nu21, Rnu21);
				          double dF1 = fabs(nu11 * Rnu21);
				  				  
				          FieldVector<Scalar,dim> Rnu23(0);
				          R.umv(nu23, Rnu23);				  
				          double dF3 = fabs(nu13 * Rnu23);
				  
				          // compute components needed for flux calculation, denoted as 'g'
				          FieldVector<Scalar,dim> K1nu11(0);
				          K1.umv(nu11, K1nu11);
				          FieldVector<Scalar,dim> K1nu21(0);
				          K1.umv(nu21, K1nu21);				  
				          FieldVector<Scalar,dim> K3nu13(0);
				          K3.umv(nu13, K3nu13);
				          FieldVector<Scalar,dim> K3nu23(0);
				          K3.umv(nu23, K3nu23);				  				 
                                          double g111 = lambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
				          double g121 = lambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
				          double g211 = lambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
				          double g221 = lambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;                               
				          double g213 = alambda3 * (integrationOuterNormaln3 * K3nu13)/dF3;
				          double g223 = alambda3 * (integrationOuterNormaln3 * K3nu23)/dF3;
                                  
				          // compute transmissibility matrix T = CA^{-1}B+F
				          FieldMatrix C(0), A(0), F(0), B(0); 

				          // evaluate matrix C, F, A, B
				          C[0][0] = -g111;
				          C[0][1] = -g121;
				          C[1][1] = g223;				          

				          F[0][0] = g111 + g121;
				          F[1][1] = g213 - g223;			         

				          A[0][0] = g111;
				          A[0][1] = g121;				  
				          A[1][0] = g211;
				          A[1][1] = g223 + g221;				 
				          
				          B[0][0] = g111 + g121;
				          B[1][0] = g211 + g221;
				          B[1][1] = g223 - g213;
				          
				          // compute T
				          A.invert();
                                          FieldMatrix CAinv(C.rightmultiply(A));
                                          F += B.leftmultiply(CAinv);
                                          FieldMatrix T(F);

                                          // compute vector r
			                  // evaluate r1, r2
			                  Dune::FieldVector<Scalar,dim> r1(0), r2(0);
                                          r1[1] = -g213 * g2;
                                          r2[0] = -J1 * face12vol/2.0;
                                          r2[1] = g213 * g2;

				          // compute  r = CA^{-1}r1
				          Dune::FieldVector<Scalar,dim> r(0);
				          CAinv.umv(r2, r);
                                          r += r1;
					 			 			    
				          // assemble the global matrix M and right hand side f			 
				          M[globalIdx1][globalIdx1] += T[1][0];
				          M[globalIdx1][globalIdx3] += T[1][1];
                                          f[globalIdx1] -= r[1];

                                        }                                      
                                    }
				}      				    							   
			      // 'isIt' is on Dirichlet boundary 
			      else
				{				  
                                  // get Dirichlet boundary value
				  double g1 = this->diffProblem.dirichletPress(globalPosFace12, *eIt, localPosFace12);

				  // compute total mobility for Dirichlet boundary
				  double alambda1;
                                  alambda1 = lambda1;

                                  // 'nextisIt' is on boundary
                                  if (nextisIt->boundary())
                                    {
				      // get boundary condition for boundary face (nextisIt) center
				      BoundaryConditions::Flags nextisItbctype 
                                        = this->diffProblem.bctypePress(globalPosFace13, *eIt, localPosFace13);

				      // 'nextisIt': Dirichlet boundary
				      if (nextisItbctype == BoundaryConditions::dirichlet)
				        {    
                                           // get Dirichlet boundary value of 'nextisIt'
				           double g3 = this->diffProblem.dirichletPress(globalPosFace13, *eIt,localPosFace13);

				          // compute total mobility for Dirichlet boundary
				          double alambda1;
                                          alambda1 = lambda1;
					  
				          // compute normal vectors nu11,nu21;
                                          FieldVector<Scalar,dim> nu11(0);
				          R.umv(globalPosFace13-globalPos1 ,nu11);

				          FieldVector<Scalar,dim> nu21(0);
				          R.umv(globalPos1-globalPosFace12, nu21);

				          // compute dF1 i.e., the area of quadrilateral made by normal vectors 'nu'
				          FieldVector<Scalar,dim> Rnu21(0);
				          R.umv(nu21, Rnu21);
				          double dF1 = fabs(nu11 * Rnu21);
				  				  				       
				          // compute components needed for flux calculation, denoted as 'g'
				          FieldVector<Scalar,dim> K1nu11(0);
				          K1.umv(nu11, K1nu11);
				          FieldVector<Scalar,dim> K1nu21(0);
				          K1.umv(nu21, K1nu21);						  				 
                                          double g111 = alambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
				          double g121 = alambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
				          double g211 = alambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
				          double g221 = alambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;                             
                                          
                                           // evaluate T1, T3, r1, r3
                                           double T1 = g111 + g121;
                                           double T3 = g211 + g221;
                                           double r1 = g111 * g1 + g121 * g3;
                                           double r3 = g211 * g1 + g221 * g3;

                                           // assemble matrix M and right hand side f
                                           M[globalIdx1][globalIdx1] += T1 + T3;
                                           f[globalIdx1] += r1 + r3;
                                        }
                                      // 'nextisIt': Neumann boundary
                                      else
                                        {
                                           // get Neumann boundary value of 'nextisIt'
				           double J3 = this->diffProblem.neumannPress(globalPosFace13, *eIt,localPosFace13);

				          // compute normal vectors nu11,nu21;
                                          FieldVector<Scalar,dim> nu11(0);
				          R.umv(globalPosFace13-globalPos1 ,nu11);

				          FieldVector<Scalar,dim> nu21(0);
				          R.umv(globalPos1-globalPosFace12, nu21);
	
				          // compute dF1 i.e., the area of quadrilateral made by normal vectors 'nu'
				          FieldVector<Scalar,dim> Rnu21(0);
				          R.umv(nu21, Rnu21);
				          double dF1 = fabs(nu11 * Rnu21);
				  				  				       
				          // compute components needed for flux calculation, denoted as 'g'
				          FieldVector<Scalar,dim> K1nu11(0);
				          K1.umv(nu11, K1nu11);
				          FieldVector<Scalar,dim> K1nu21(0);
				          K1.umv(nu21, K1nu21);						  				 
                                          double g111 = alambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
				          double g121 = alambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
				          double g211 = alambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
				          double g221 = alambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;   

                                           // evaluate T, r
                                           double T = g111 - g211 * g121/g221;
                                           double r = -T * g1 - g121 * (-J3) * nextisIt->geometry().volume()/ (2.0 * g221);

                                           // assemble matrix M and right hand side f
                                           M[globalIdx1][globalIdx1] += T;
                                           f[globalIdx1] -= r;
                                        }
                                    }
                                  // 'nextisIt' is inside
                                  else
                                    {
				      // neighbor cell 3
				      // access neighbor cell 3
				      ElementPointer nextisItoutside = nextisIt->outside();
				      int globalIdx3 = this->diffProblem.variables().indexDiffusion(*nextisItoutside);

				      // neighbor cell 3 geometry type
				      Dune::GeometryType gt3 = nextisItoutside->geometry().type();

				      // neighbor cell 3 center in reference element
				      const LocalPosition&
				        localPos3 = Dune::ReferenceElements<Scalar,dim>::general(gt3).position(0,0);

				      // get global coordinate of neighbor cell 3 center
				      GlobalPosition
				        globalPos3 = nextisItoutside->geometry().global(localPos3);

				      // get absolute permeability of neighbor cell 3
				      FieldMatrix K3(this->diffProblem.soil().K(globalPos3, *nextisItoutside, localPos3));

				      // get total mobility of neighbor cell 3
				      double lambda3;
                                      double sat3 = this->diffProblem.variables().saturation()[globalIdx3];
	                              lambda3 = this->diffProblem.materialLaw().mobTotal(sat3, globalPos3, *nextisItoutside, localPos3);
				      
				      // get the information of the face 'isIt34' between cell3 and cell4 (locally numbered)
				      IntersectionIterator isIt34 = this->gridView.template ibegin(*nextisItoutside);
                                      IntersectionIterator innernextisItEnd = this->gridView.template iend(*nextisItoutside);
				      for (IntersectionIterator innerisIt = this->gridView.template ibegin(*nextisItoutside);
                                           innerisIt != innernextisItEnd; ++innerisIt)
				        {
				          if (innerisIt->boundary())
				    	    {
					      for (int i=0; i<innerisIt->geometry().corners(); ++i)
					        {
					          GlobalPosition innerisItcorner = innerisIt->geometry().corner(i);

					          if (innerisItcorner == corner1234)
					            {
						      isIt34 = innerisIt;
						      continue;
					            }
					        }					     
					    }				      
				        }

				      // get geometry type of face 'isIt34'
				      Dune::GeometryType gtf34 = isIt34->geometryInInside().type();

				      // center in face's reference element
				      const Dune::FieldVector<Scalar,dim-1>&
				        faceLocal34 = Dune::ReferenceElements<Scalar,dim-1>::general(gtf34).position(0,0);

			              // center of face inside volume reference element
			              const LocalPosition&
			                localPosFace34 = Dune::ReferenceElements<Scalar,dim>::general(gtf34).position(isIt34->indexInInside(),1);

				      // center of face in global coordinates, i.e., the midpoint of edge 'isIt34'
				      GlobalPosition
				        globalPosFace34 = isIt34->geometry().global(faceLocal34);

				      // get outer normal vector scaled with half volume of face 'isIt34'
				      Dune::FieldVector<Scalar,dimWorld> integrationOuterNormaln2
				        = isIt34->integrationOuterNormal(faceLocal34);
                                      integrationOuterNormaln2
                                        *= Dune::ReferenceElements<Scalar,dim-1>::general(gtf34).volume()/2.0;

				      // get boundary condition for boundary face (isIt34) center
				      BoundaryConditions::Flags isIt34bctype 
                                        = this->diffProblem.bctypePress(globalPosFace34, *nextisItoutside, localPosFace34);

				      // 'isIt34': Dirichlet boundary
				      if (isIt34bctype == BoundaryConditions::dirichlet)
				        { 
                                           // get Dirichlet boundary value of 'isIt34'
				           double g2 = this->diffProblem.dirichletPress(globalPosFace34, *nextisItoutside,localPosFace34);

				          // compute total mobility for Dirichlet boundary
				          double alambda3;
                                          alambda3 = lambda3;

				          // compute normal vectors nu11,nu21; nu13, nu23;
                                          FieldVector<Scalar,dim> nu11(0);
				          R.umv(globalPosFace13-globalPos1 ,nu11);

				          FieldVector<Scalar,dim> nu21(0);
				          R.umv(globalPos1-globalPosFace12, nu21);

				          FieldVector<Scalar,dim> nu13(0);
				          R.umv(globalPos3-globalPosFace13, nu13);

				          FieldVector<Scalar,dim> nu23(0);
				          R.umv(globalPos3-globalPosFace34, nu23);

				          // compute dF1, dF3 i.e., the area of quadrilateral made by normal vectors 'nu'
				          FieldVector<Scalar,dim> Rnu21(0);
				          R.umv(nu21, Rnu21);
				          double dF1 = fabs(nu11 * Rnu21);
				  				  
				          FieldVector<Scalar,dim> Rnu23(0);
				          R.umv(nu23, Rnu23);				  
				          double dF3 = fabs(nu13 * Rnu23);
				  
				          // compute components needed for flux calculation, denoted as 'g'
				          FieldVector<Scalar,dim> K1nu11(0);
				          K1.umv(nu11, K1nu11);
				          FieldVector<Scalar,dim> K1nu21(0);
				          K1.umv(nu21, K1nu21);				  
				          FieldVector<Scalar,dim> K3nu13(0);
				          K3.umv(nu13, K3nu13);
				          FieldVector<Scalar,dim> K3nu23(0);
				          K3.umv(nu23, K3nu23);				  				 
                                          double g111 = alambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
				          double g121 = alambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
				          double g211 = alambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
				          double g221 = alambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;                              
				          double g213 = alambda3 * (integrationOuterNormaln3 * K3nu13)/dF3;
				          double g223 = alambda3 * (integrationOuterNormaln3 * K3nu23)/dF3;

				          // compute the matrix T & vector r 
				          FieldMatrix T(0); 
                                          Dune::FieldVector<Scalar,dim> r(0);
                                          
                                          double coe = g221 + g223; 

				          // evaluate matrix T
				          T[0][0] = g111 + g121 * (g223 - g211)/coe;
				          T[0][1] = -g121 * (g223 - g213)/coe;		
				          T[1][0] = g223 * (g211 + g221)/coe;
				          T[1][1] = -g221 * (g223 - g213)/coe;
				         				         			         
                                          // evaluate vector r                                       
                                          r[0] = -g111 * g1 + (g1 * g121 * g211 - g2 * g213 * g121)/coe;    
                                          r[1] = -(g1 * g211 * g223 + g2 * g221 * g213)/coe;                          			  
				          // assemble the global matrix M and right hand side f			 
				          M[globalIdx1][globalIdx1] += T[0][0] + T[1][0];
				          M[globalIdx1][globalIdx3] += T[0][1] + T[1][1];
                                          f[globalIdx1] -= r[0] + r[1];
   
                                        }
                                      // 'isIt34': Neumann boundary
                                      else
                                        {
                                           // get Neumann boundary value of 'isIt34'
				           double J2 = this->diffProblem.neumannPress(globalPosFace34, *nextisItoutside,localPosFace34);
				       
				          // compute normal vectors nu11,nu21; nu13, nu23;
                                          FieldVector<Scalar,dim> nu11(0);
				          R.umv(globalPosFace13-globalPos1 ,nu11);

				          FieldVector<Scalar,dim> nu21(0);
				          R.umv(globalPos1-globalPosFace12, nu21);

				          FieldVector<Scalar,dim> nu13(0);
				          R.umv(globalPos3-globalPosFace13, nu13);

				          FieldVector<Scalar,dim> nu23(0);
				          R.umv(globalPos3-globalPosFace34, nu23);
				 
				          // compute dF1, dF3 i.e., the area of quadrilateral made by normal vectors 'nu'
				          FieldVector<Scalar,dim> Rnu21(0);
				          R.umv(nu21, Rnu21);
				          double dF1 = fabs(nu11 * Rnu21);
				  				  
				          FieldVector<Scalar,dim> Rnu23(0);
				          R.umv(nu23, Rnu23);				  
				          double dF3 = fabs(nu13 * Rnu23);
				  
				          // compute components needed for flux calculation, denoted as 'g'
				          FieldVector<Scalar,dim> K1nu11(0);
				          K1.umv(nu11, K1nu11);
				          FieldVector<Scalar,dim> K1nu21(0);
				          K1.umv(nu21, K1nu21);				  
				          FieldVector<Scalar,dim> K3nu13(0);
				          K3.umv(nu13, K3nu13);
				          FieldVector<Scalar,dim> K3nu23(0);
				          K3.umv(nu23, K3nu23);				  				 
                                          double g111 = alambda1 * (integrationOuterNormaln1 * K1nu11)/dF1;
				          double g121 = alambda1 * (integrationOuterNormaln1 * K1nu21)/dF1;
				          double g211 = alambda1 * (integrationOuterNormaln3 * K1nu11)/dF1;
				          double g221 = alambda1 * (integrationOuterNormaln3 * K1nu21)/dF1;                              
                                          double g113 = lambda3 * (integrationOuterNormaln2 * K3nu13)/dF3;
				          double g123 = lambda3 * (integrationOuterNormaln2 * K3nu23)/dF3;
				          double g213 = lambda3 * (integrationOuterNormaln3 * K3nu13)/dF3;
				          double g223 = lambda3 * (integrationOuterNormaln3 * K3nu23)/dF3;
                                  
				          // compute the matrix T & vector r in v = A^{-1}(Bu + r1) = Tu + r
				          FieldMatrix A(0), B(0); 
                                          Dune::FieldVector<Scalar,dim> r1(0), r(0);

				          // evaluate matrix A, B
				          A[0][0] = g113;
				          A[0][1] = -g123;		
				          A[1][0] = -g213;
				          A[1][1] = g221 + g223;
				         
				          B[0][1] = g113 - g123;
				          B[1][0] = g211 + g221;
				          B[1][1] = g223 - g213;
				         
                                          // evaluate vector r1
                                          r1[0] = -J2 * isIt34->geometry().volume()/2.0;
                                          r1[1] = -g211 * g1;
                                          
                                          // compute T and r
				          A.invert();
                                          B.leftmultiply(A);
                                          FieldMatrix T(B);
					  A.umv(r1, r);                                         
                                          			 			  
				          // assemble the global matrix M and right hand side f			 
				          M[globalIdx1][globalIdx1] += (g111 + g121 - g121 * T[1][0]) + (g211 + g221 - g221 * T[1][0]);
				          M[globalIdx1][globalIdx3] += -g121 * T[1][1] - g221 * T[1][1];
                                          f[globalIdx1] += (g111 + g211) * g1 + (g121 + g221) * r[1];

                                        } 
                                     } 		   
			         } 
                             }
			  
			  } // end all intersections 	
		    
	    } // end grid traversal

	  // get the number of nonzero terms in the matrix
	  double num_nonzero = 0;
	  
	  // determine position of matrix entries 
	  for (ElementIterator eIt = this->gridView.template begin<0>(); eIt != this->gridView.template end<0>(); ++eIt)
	    {
	              // cell index
		      int globalIdxI = this->diffProblem.variables().indexDiffusion(*eIt);

		      if (M[globalIdxI][globalIdxI] != 0)
		        ++num_nonzero;
				
		      // run through all intersections with neighbors 
                      IntersectionIterator isItBegin = this->gridView.template ibegin(*eIt);
	              IntersectionIterator isItEnd = this->gridView.template iend(*eIt);
		      for (IntersectionIterator isIt = isItBegin; isIt!=isItEnd; ++isIt)
			{
                          IntersectionIterator tempisIt = isIt;
			  IntersectionIterator tempisItBegin = isItBegin;

                          // 'nextisIt' iterates over next codimension 1 intersection neighboring with 'isIt'
			  // sequence of next is anticlockwise of 'isIt'
                          IntersectionIterator nextisIt = ++tempisIt;

 			  // get 'nextisIt'
			  switch (ch)
			    {
			      // for SGrid
			      case 's':
			      {
			        if (nextisIt == isItEnd)
			          {
			 	    nextisIt = isItBegin;
			          }
			        else
			          {
			  	    nextisIt = ++tempisIt;

				    if (nextisIt == isItEnd)
				      {
				        nextisIt = ++tempisItBegin;
				      }
			          }

				break;
			      }
			      // for UGGrid
			      case 'u':
			      {
				if (nextisIt == isItEnd)
				  nextisIt = isItBegin;

				break;
			      }
			    }

			  if (isIt->neighbor()) 
			    {
			              // access neighbor
			              ElementPointer outside = isIt->outside();
				      int globalIdxJ = this->diffProblem.variables().indexDiffusion(*outside);
		
		                      if (M[globalIdxI][globalIdxJ] != 0)
				        ++num_nonzero;
			    }

			  if (isIt->neighbor() && nextisIt->neighbor())
			    {
		              // access the common neighbor of isIt's and nextisIt's outside
                              ElementPointer outside = isIt->outside();
		              ElementPointer nextisItoutside = nextisIt->outside();
                                    
                              IntersectionIterator innerisItEnd = this->gridView.template iend(*outside);
                              IntersectionIterator innernextisItEnd = this->gridView.template iend(*nextisItoutside);

			      for (IntersectionIterator innerisIt = this->gridView.template ibegin(*outside);
				   innerisIt!=innerisItEnd; ++innerisIt )
				for (IntersectionIterator innernextisIt = this->gridView.template ibegin(*nextisItoutside);
				     innernextisIt!=innernextisItEnd; ++innernextisIt)
				  {
				    if (innerisIt->neighbor() && innernextisIt->neighbor())
			              {
					ElementPointer innerisItoutside = innerisIt->outside();
					ElementPointer innernextisItoutside = innernextisIt->outside();

					if (innerisItoutside == innernextisItoutside && innerisItoutside != isIt->inside())
					  {
					    int globalIdxJ = this->diffProblem.variables().indexDiffusion(*innerisItoutside);
					      
					    if (M[globalIdxI][globalIdxJ] != 0)
				              ++num_nonzero; 
                                          }
			              }
				  }
			    }			    
			  } // end of 'for' IntersectionIterator
	      } // end of 'for' ElementIterator
	  
	  std::cout << "number of nonzero terms in the MPFA O-matrix on level " << this->gridView.grid().maxLevel() <<" nnmat: " << num_nonzero << std::endl;
	  
	  return;	    
	}
	

  template<class GridView, class Scalar, class VC, class Problem>
  void MPFAODiffusion<GridView, Scalar, VC, Problem>::solve()
  {
	Dune::MatrixAdapter<MatrixType,Vector,Vector> op(M);        // make linear operator from M
        Dune::InverseOperatorResult r;

	if (preconditionerName_ == "SeqILU0") {
                // preconditioner object
		SeqILU0<MatrixType,Vector,Vector> preconditioner(M, 1.0);
		if (solverName_ == "CG") {
                        // an inverse operator 
			CGSolver<Vector> solver(op, preconditioner, 1E-14, 1000, 1);
			solver.apply(this->diffProblem.variables().pressure(), f, r);
		} else if (solverName_ == "BiCGSTAB") {
			BiCGSTABSolver<Vector> solver(op, preconditioner, 1E-14, 1000, 1);
			solver.apply(this->diffProblem.variables().pressure(), f, r);
		} else
			DUNE_THROW(NotImplemented, "MPFAODiffusion :: solve : combination "
					<< preconditionerName_<< " and "<< solverName_ << ".");
	} else if (preconditionerName_ == "SeqPardiso") {
		SeqPardiso<MatrixType,Vector,Vector> preconditioner(M);
		if (solverName_ == "Loop") {
			LoopSolver<Vector> solver(op, preconditioner, 1E-14, 1000, 1);
			solver.apply(this->diffProblem.variables().pressure(), f, r);
		}else if (solverName_ == "BiCGSTAB") {
			BiCGSTABSolver<Vector> solver(op, preconditioner, 1E-14, 1000, 1);
			solver.apply(this->diffProblem.variables().pressure(), f, r);
		} else
			DUNE_THROW(NotImplemented, "MPFAODiffusion :: solve : combination "
					<< preconditionerName_<< " and "<< solverName_ << ".");
	} else
		DUNE_THROW(NotImplemented, "MPFAODiffusion :: solve : preconditioner "
				<< preconditionerName_ << ".");

	return;
  }
} // end of Dune namespace
#endif
