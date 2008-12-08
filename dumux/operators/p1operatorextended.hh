// $Id$ 

#ifndef DUNE_P1OPERATOR_HH
#define DUNE_P1OPERATOR_HH

#include<iostream>
#include<vector>
#include<set>
#include<map>
#include<stdio.h>
#include<stdlib.h>

#include<dune/common/timer.hh>
#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/mcmgmapper.hh>
#include<dune/grid/utility/intersectiongetter.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include"localstiffnessextended.hh"
#include"../functions/p1functionextended.hh" // for parallel extender class

/**
 * @file
 * @brief  defines a class for piecewise linear finite element functions
 * @author Peter Bastian
 */

/*! @defgroup DISC_Operators Operators
  @ingroup DISC
  @brief
  
  @section D1 Introduction
  <!--=================-->
  
  To be written
*/

namespace Dune
{
  /** @addtogroup DISC_Operators
   *
   * @{
   */
  /**
   * @brief defines a class for piecewise linear finite element functions
   *
   */

  // make a type to sort matrix entries
  typedef std::pair<int,int> P1OperatorLink;

  // template meta program for inserting indices
  template<int n, int c>
  struct P1Operator_meta {
	template<class Entity, class VMapper, class AMapper, class Refelem, class Matrix>
	static void addrowscube (const Entity& e, const VMapper& vertexmapper, const AMapper& allmapper, 
							 const Refelem& refelem, Matrix& A, std::vector<bool>& visited, 
							 int hangingnodes, std::set<P1OperatorLink>& links, std::map<int,int>& doubled2Original)
	{
	  if (refelem.type(0,0).isCube())
		{
		  for (int i=0; i<refelem.size(c); i++) // loop over subentities of codim c of e
			{
			  int index = allmapper.template map<c>(e,i);
			  if (!visited[index]) 
				{
				  int corners = refelem.size(i,c,n);
				  for (int j=0; j<corners/2; j++) // uses fact that diagonals are (0,corners-1), (1,corners-2) ...
					{
					  int alpha = vertexmapper.template map<n>(e,refelem.subEntity(i,c,j,n));
					  int beta = vertexmapper.template map<n>(e,refelem.subEntity(i,c,corners-1-j,n));
					  if (doubled2Original.find(alpha) != doubled2Original.end())
						  alpha = doubled2Original[alpha];
					  if (doubled2Original.find(beta) != doubled2Original.end())
						  beta = doubled2Original[beta];
					  A.incrementrowsize(alpha);
					  A.incrementrowsize(beta);
					  if (hangingnodes>0) // delete standard links
						{
						  links.erase(P1OperatorLink(alpha,beta));
						  links.erase(P1OperatorLink(beta,alpha));
						}
					}
				  visited[index] = true;
				}
			}
		}
	  if (refelem.type(0,0).isPyramid() && c==1)
		{
		  int index = allmapper.template map<c>(e,0);
		  if (!visited[index]) 
			{
			  int alpha = vertexmapper.template map<n>(e,0);
			  int beta = vertexmapper.template map<n>(e,2);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.incrementrowsize(alpha);
			  A.incrementrowsize(beta);
			  if (hangingnodes>0) // delete standard links
				{
				  links.erase(P1OperatorLink(alpha,beta));
				  links.erase(P1OperatorLink(beta,alpha));
				}
			  alpha = vertexmapper.template map<n>(e,1);
			  beta = vertexmapper.template map<n>(e,3);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.incrementrowsize(alpha);
			  A.incrementrowsize(beta);
			  if (hangingnodes>0) // delete standard links
				{
				  links.erase(P1OperatorLink(alpha,beta));
				  links.erase(P1OperatorLink(beta,alpha));
				}
			  visited[index] = true;
			}
		}	  
	  if (refelem.type(0,0).isPrism() && c==1)
		{
		  int index = allmapper.template map<c>(e,1);
		  if (!visited[index]) 
			{
			  int alpha = vertexmapper.template map<n>(e,0);
			  int beta = vertexmapper.template map<n>(e,4);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.incrementrowsize(alpha);
			  A.incrementrowsize(beta);
			  if (hangingnodes>0) // delete standard links
				{
				  links.erase(P1OperatorLink(alpha,beta));
				  links.erase(P1OperatorLink(beta,alpha));
				}
			  alpha = vertexmapper.template map<n>(e,1);
			  beta = vertexmapper.template map<n>(e,3);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.incrementrowsize(alpha);
			  A.incrementrowsize(beta);
			  if (hangingnodes>0) // delete standard links
				{
				  links.erase(P1OperatorLink(alpha,beta));
				  links.erase(P1OperatorLink(beta,alpha));
				}
			  visited[index] = true;
			}
		  index = allmapper.template map<c>(e,2);
		  if (!visited[index]) 
			{
			  int alpha = vertexmapper.template map<n>(e,1);
			  int beta = vertexmapper.template map<n>(e,5);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.incrementrowsize(alpha);
			  A.incrementrowsize(beta);
			  if (hangingnodes>0) // delete standard links
				{
				  links.erase(P1OperatorLink(alpha,beta));
				  links.erase(P1OperatorLink(beta,alpha));
				}
			  alpha = vertexmapper.template map<n>(e,2);
			  beta = vertexmapper.template map<n>(e,4);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.incrementrowsize(alpha);
			  A.incrementrowsize(beta);
			  if (hangingnodes>0) // delete standard links
				{
				  links.erase(P1OperatorLink(alpha,beta));
				  links.erase(P1OperatorLink(beta,alpha));
				}
			  visited[index] = true;
			}
		  index = allmapper.template map<c>(e,3);
		  if (!visited[index]) 
			{
			  int alpha = vertexmapper.template map<n>(e,0);
			  int beta = vertexmapper.template map<n>(e,5);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.incrementrowsize(alpha);
			  A.incrementrowsize(beta);
			  if (hangingnodes>0) // delete standard links
				{
				  links.erase(P1OperatorLink(alpha,beta));
				  links.erase(P1OperatorLink(beta,alpha));
				}
			  alpha = vertexmapper.template map<n>(e,2);
			  beta = vertexmapper.template map<n>(e,3);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.incrementrowsize(alpha);
			  A.incrementrowsize(beta);
			  if (hangingnodes>0) // delete standard links
				{
				  links.erase(P1OperatorLink(alpha,beta));
				  links.erase(P1OperatorLink(beta,alpha));
				}
			  visited[index] = true;
			}
		}	  
	  P1Operator_meta<n,c-1>::addrowscube(e,vertexmapper,allmapper,refelem,A,visited,hangingnodes,links,doubled2Original);
	  return;
	}
	template<class Entity, class VMapper, class AMapper, class Refelem, class Matrix>
	static void addindicescube (const Entity& e, const VMapper& vertexmapper, const AMapper& allmapper, 
				   const Refelem& refelem, Matrix& A, std::vector<bool>& visited, std::map<int,int>& doubled2Original)
	{
	  if (refelem.type(0,0).isCube())
		{
		  for (int i=0; i<refelem.size(c); i++)
			{
			  int index = allmapper.template map<c>(e,i);
			  if (!visited[index]) 
				{
				  int corners = refelem.size(i,c,n);
				  for (int j=0; j<corners/2; j++) // uses fact that diagonals are (0,corners-1), (1,corners-2) ...
					{
					  int alpha = vertexmapper.template map<n>(e,refelem.subEntity(i,c,j,n));
					  int beta = vertexmapper.template map<n>(e,refelem.subEntity(i,c,corners-1-j,n));
					  if (doubled2Original.find(alpha) != doubled2Original.end())
						  alpha = doubled2Original[alpha];
					  if (doubled2Original.find(beta) != doubled2Original.end())
						  beta = doubled2Original[beta];
					  A.addindex(alpha,beta);
					  A.addindex(beta,alpha);
					}
				  visited[index] = true;
				}
			}
		}
	  if (refelem.type(0,0).isPyramid() && c==1)
		{
		  int index = allmapper.template map<c>(e,0);
		  if (!visited[index]) 
			{
			  int alpha = vertexmapper.template map<n>(e,0);
			  int beta = vertexmapper.template map<n>(e,2);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.addindex(alpha,beta);
			  A.addindex(beta,alpha);
			  alpha = vertexmapper.template map<n>(e,1);
			  beta = vertexmapper.template map<n>(e,3);
			  A.addindex(alpha,beta);
			  A.addindex(beta,alpha);
			  visited[index] = true;
			}
		}	  
	  if (refelem.type(0,0).isPrism() && c==1)
		{
		  int index = allmapper.template map<c>(e,1);
		  if (!visited[index]) 
			{
			  int alpha = vertexmapper.template map<n>(e,0);
			  int beta = vertexmapper.template map<n>(e,4);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.addindex(alpha,beta);
			  A.addindex(beta,alpha);
			  alpha = vertexmapper.template map<n>(e,1);
			  beta = vertexmapper.template map<n>(e,3);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.addindex(alpha,beta);
			  A.addindex(beta,alpha);
			  visited[index] = true;
			}
		  index = allmapper.template map<c>(e,2);
		  if (!visited[index]) 
			{
			  int alpha = vertexmapper.template map<n>(e,1);
			  int beta = vertexmapper.template map<n>(e,5);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.addindex(alpha,beta);
			  A.addindex(beta,alpha);
			  alpha = vertexmapper.template map<n>(e,2);
			  beta = vertexmapper.template map<n>(e,4);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.addindex(alpha,beta);
			  A.addindex(beta,alpha);
			  visited[index] = true;
			}
		  index = allmapper.template map<c>(e,3);
		  if (!visited[index]) 
			{
			  int alpha = vertexmapper.template map<n>(e,0);
			  int beta = vertexmapper.template map<n>(e,5);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.addindex(alpha,beta);
			  A.addindex(beta,alpha);
			  alpha = vertexmapper.template map<n>(e,2);
			  beta = vertexmapper.template map<n>(e,3);
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  A.addindex(alpha,beta);
			  A.addindex(beta,alpha);
			  visited[index] = true;
			}
		}	  
	  P1Operator_meta<n,c-1>::addindicescube(e,vertexmapper,allmapper,refelem,A,visited,doubled2Original);
	  return;
	}
  };
  template<int n>
  struct P1Operator_meta<n,0> {
	template<class Entity, class VMapper, class AMapper, class Refelem, class Matrix>
	static void addrowscube (const Entity& e, const VMapper& vertexmapper, const AMapper& allmapper, 
							 const Refelem& refelem, Matrix& A, std::vector<bool>& visited,
							 int hangingnodes, std::set<P1OperatorLink>& links, std::map<int,int>& doubled2Original)
	{
	  if (!refelem.type(0,0).isCube()) return;
	  int corners = refelem.size(n);
	  for (int j=0; j<corners/2; j++) // uses fact that diagonals are (0,corners-1), (1,corners-2) ...
		{
		  int alpha = vertexmapper.template map<n>(e,refelem.subEntity(0,0,j,n));
		  int beta = vertexmapper.template map<n>(e,refelem.subEntity(0,0,corners-1-j,n));
		  if (doubled2Original.find(alpha) != doubled2Original.end())
			  alpha = doubled2Original[alpha];
		  if (doubled2Original.find(beta) != doubled2Original.end())
			  beta = doubled2Original[beta];
		  A.incrementrowsize(alpha);
		  A.incrementrowsize(beta);
		  if (hangingnodes>0) // delete standard links
			{
			  links.erase(P1OperatorLink(alpha,beta));
			  links.erase(P1OperatorLink(beta,alpha));
			}
		}
	  return;
	}
	template<class Entity, class VMapper, class AMapper, class Refelem, class Matrix>
	static void addindicescube (const Entity& e, const VMapper& vertexmapper, const AMapper& allmapper, 
				  const Refelem& refelem, Matrix& A, std::vector<bool>& visited, std::map<int,int>& doubled2Original)
	{
	  if (!refelem.type(0,0).isCube()) return;
	  int corners = refelem.size(n);
	  for (int j=0; j<corners/2; j++) // uses fact that diagonals are (0,corners-1), (1,corners-2) ...
		{
		  int alpha = vertexmapper.template map<n>(e,refelem.subEntity(0,0,j,n));
		  int beta = vertexmapper.template map<n>(e,refelem.subEntity(0,0,corners-1-j,n));
		  if (doubled2Original.find(alpha) != doubled2Original.end())
			  alpha = doubled2Original[alpha];
		  if (doubled2Original.find(beta) != doubled2Original.end())
			  beta = doubled2Original[beta];
		  A.addindex(alpha,beta);
		  A.addindex(beta,alpha);
		}
	  return;
	}
  };

  // handles case n=1, c=-1
  template<int n>
  struct P1Operator_meta<n,-1> {
	template<class Entity, class VMapper, class AMapper, class Refelem, class Matrix>
	static void addrowscube (const Entity& e, const VMapper& vertexmapper, const AMapper& allmapper, 
							 const Refelem& refelem, Matrix& A, std::vector<bool>& visited,
							 int hangingnodes, std::set<P1OperatorLink>& links, std::map<int,int>& doubled2Original)
	{
	  return;
	}
	template<class Entity, class VMapper, class AMapper, class Refelem, class Matrix>
	static void addindicescube (const Entity& e, const VMapper& vertexmapper, const AMapper& allmapper, 
				  const Refelem& refelem, Matrix& A, std::vector<bool>& visited, std::map<int,int>& doubled2Original)
	{
	  return;
	}
  };

    template<typename Grid, typename RT, typename Tag, bool hasHangingNodes>
    struct HangingNodesIdentifier
    {};
    
    template<typename Grid,typename RT,typename Tag>
    struct HangingNodesIdentifier<Grid,RT,Tag,false>
    {
      template<typename IndexSet, typename VM>
      inline void process(std::vector<unsigned char> S, std::vector<bool>& hanging, std::set<P1OperatorLink>& links, IndexSet& indexset, VM& vertexmapper) const
      {}
      
      inline std::size_t count() const
      {
	return 0;
      }
    };

    template<typename Grid,typename RT,typename Tag>
    class HangingNodesIdentifier<Grid,RT,Tag,true>
      : public IntersectionIteratorGetter<Grid,Tag>
    {
    public:
      typedef typename Grid::ctype DT;
      enum {n=Grid::dimension};
      typedef typename IntersectionIteratorGetter<Grid,Tag>::IntersectionIterator IntersectionIterator;
      typedef typename Grid::template Codim<0>::EntityPointer EEntityPointer;

      template<typename GridView, typename VM>
      void process(std::vector<unsigned char> S, std::vector<bool>& hanging, std::set<P1OperatorLink>& links, GridView& gridView, VM& vertexmapper)
      {
	
	typedef typename GridView::template Codim<0>::Iterator Iterator;
	// LOOP 1 : Prepare hanging node detection
	Iterator eendit = gridView.template end<0>();
	for (Iterator it = gridView.template begin<0>(); it!=eendit; ++it)
	  {
	    Dune::GeometryType gt = it->type();
	    const typename Dune::ReferenceElementContainer<DT,n>::value_type& 
	      refelem = ReferenceElements<DT,n>::general(gt);
	    
	    // compute S value in vertex
	    for (int i=0; i<refelem.size(n); i++)
	      {
		int alpha = vertexmapper.template map<n>(*it,i);
		if (S[alpha]>it->level()) S[alpha] = it->level(); // compute minimum
	      }
	  }

	// LOOP 2 : second stage of detecting hanging nodes
	for (Iterator it = gridView.template begin<0>(); it!=eendit; ++it)
	  {
	    Dune::GeometryType gt = it->type();
	    const typename Dune::ReferenceElementContainer<DT,n>::value_type& 
	      refelem = ReferenceElements<DT,n>::general(gt);
	    
	    // detect hanging nodes
	    IntersectionIterator endiit = end(*it);
	    for (IntersectionIterator iit = begin(*it); iit!=endiit; ++iit)
	      if (iit->neighbor())
		{
		  // check if neighbor is on lower level
		  const EEntityPointer outside = iit->outside();
		  if (it->level()<=outside->level()) continue;
		  
		  // loop over all vertices of this face
		  for (int j=0; j<refelem.size(iit->numberInSelf(),1,n); j++)
		    {
		      int alpha = vertexmapper.template map<n>(*it,refelem.subEntity(iit->numberInSelf(),1,j,n));
		      if (S[alpha]==it->level()) 
			hanging[alpha] = true;
		    }
		}
	  }
	
	// local to global maps
	int l2g[Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize];
	int fl2g[Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize];
	
	// LOOP 3 : determine additional links due to hanging nodes
	for (Iterator it = gridView.template begin<0>(); it!=eendit; ++it)
	  {
	    Dune::GeometryType gt = it->type();
	    const typename Dune::ReferenceElementContainer<DT,n>::value_type& 
	      refelem = ReferenceElements<DT,n>::general(gt);
	    
	    // build local to global map
	    bool hasHangingNodes = false; // flag set to true if this element has hanging nodes
	    for (int i=0; i<refelem.size(n); i++)
	      {
		l2g[i] = vertexmapper.template map<n>(*it,i);
		if (hanging[l2g[i]]) hasHangingNodes=true;
	      }
	    if (!hasHangingNodes) continue;
	    
	    // handle father element if hanging nodes were detected
	    // get father element
	    const EEntityPointer father = it->father();
	    
	    // build local to global map for father
	    for (int i=0; i<refelem.size(n); i++)
	      fl2g[i] = vertexmapper.template map<n>(*father,i);
	    
	    // a map that inverts l2g
	    std::map<int,int> g2l;
	    for (int i=0; i<refelem.size(n); i++)
	      g2l[l2g[i]] = i;
	    
	    // connect all fine nodes to all coarse nodes
	    for (int i=0; i<refelem.size(n); i++) // nodes in *it
	      for (int j=0; j<refelem.size(n); j++) // nodes in *father
		if (g2l.find(fl2g[j])==g2l.end())
		  {
		    links.insert(P1OperatorLink(l2g[i],fl2g[j]));
		    links.insert(P1OperatorLink(fl2g[j],l2g[i]));
		  }
	  }
	// count hanging nodes
	hangingnodes = 0;
	for (int i=0; i<vertexmapper.size(); i++) 
	  if (hanging[i]) hangingnodes++;
	
      }
      inline std::size_t count() const
      {
	return hangingnodes;
      }
    private:
      std::size_t hangingnodes;
    };
    

  
  /*! @brief A class for mapping a P1 function to a P1 function

  This class sets up a compressed row storage matrix with connectivity for P1 elements.
  It includes hanging nodes and is able to extend the matrix pattern arising
  from non-overlapping grids to minimum overlap.

  This class does not fill any entries into the matrix.

  The template parameter TypeTag describes what kind of Assembler we are. There two choices:
  <dt>LevelTag</dt><dd>We assemble on a grid level. </dd>
  <dt>LeafTag</dt><dd>We assemble on the leaf entities of the grid</dd>

  \tparam GV The grid view
   */
  template<typename TypeTag, class G, class RT, class GV, class LC, int m=1>
  class P1OperatorBase
  {
  public:    
	// export type used to store the matrix
	typedef FieldMatrix<RT,m,m> BlockType; 
	typedef BCRSMatrix<BlockType> RepresentationType;

	// mapper: one data element per vertex
	template<int dim>
	struct P1Layout
	{
	  bool contains (Dune::GeometryType gt)
	  {
              return gt.dim() == 0;
	  }
	}; 

	// mapper: one data element in every entity
	template<int dim>
	struct AllLayout
	{
	  bool contains (Dune::GeometryType gt)
	  {
		return true;
	  }
	}; 

	typedef typename G::ctype DT;
	enum {n=G::dimension};
      typedef typename GV::IndexSet IS;
	typedef typename G::template Codim<0>::Entity Entity;
	typedef typename GV::template Codim<0>::Iterator Iterator;
	typedef typename GV::template Codim<n>::Iterator VIterator;
	typedef typename G::template Codim<0>::EntityPointer EEntityPointer;
	typedef typename G::Traits::GlobalIdSet IDS;
	typedef typename IDS::IdType IdType;
	typedef std::set<IdType> GIDSet;
	typedef MultipleCodimMultipleGeomTypeMapper<G,IS,P1Layout> VM;
	typedef MultipleCodimMultipleGeomTypeMapper<G,IS,AllLayout> AM;

  private:

	// a function to compute the number of nonzeros
	// does not work for prisms and pyramids yet ?!
	int nnz (const IS& is)
	{
	  int s = 0;
	  s += is.size(n);   // vertices
	  s += 2*is.size(n-1); // edges

	  for (int c=0; c<n-1; c++)
		{
		  s += 2*is.size(GeometryType(GeometryType::cube,G::dimension-c))*(1<<(n-c-1));
		}

	  // hanging node correction
	  s += links.size();

	  return s;
	}
      
    
	// extra initialization function
	// 0) This method is executed before matrix is allocated
	// 1) determine hanging nodes as described in the paper
	// 2) generate a set with additional links
	//    The standard links are deleted later on
	bool init (const G& g, const GV& gridView, LC lc, bool extendoverlap)
	{
	  // parallel stuff we need to know
	  if (extendoverlap && g.overlapSize(0)>0)
		DUNE_THROW(GridError,"P1OperatorBase: extending overlap requires nonoverlapping grid");
	  extendOverlap = extendoverlap;
	  extraDOFs = 0;

	  // resize the S vector needed for detecting hanging nodes
	  watch.reset();
	  // resize hanging and initialize with false
	  hanging.resize(vertexmapper.size(), false);
	  // the number of levels never exceeds 100 ...
 	  std::vector<unsigned char> S(vertexmapper.size(), 100);

	  // Detect the hanging nodes
	  typedef HangingNodesIdentifier<G,RT,TypeTag,Dune::Capabilities::template hasHangingNodes<G>::v>
	    HangingNodesIdentifier;
	  HangingNodesIdentifier hangingNodes;
	  hangingNodes.process(S,hanging,links,gridView,vertexmapper);
	  
	  // count hanging nodes
	  hangingnodes = hangingNodes.count();

// 	  std::cout << "=== P1OperatorBase hanging node detection + add links " <<  watch.elapsed() << std::endl;

	  // compute additional links due to extended overlap
	  watch.reset();
	  if (extendOverlap)
		{
		  // set of neighbors in global ids for border vertices
		  std::map<int,GIDSet> borderlinks;

		  // compute extension
		  P1ExtendOverlap<G,GV,VM,LC> extender(lc);
#if HAVE_MPI
		  extender.extend(g,gv,vertexmapper,borderlinks,extraDOFs,gid2index,index2gid);
#else
		  extender.extend(g,gv,vertexmapper,borderlinks,extraDOFs,gid2index);
#endif
		  // put in extra links due to overlap 
		  // loop over all neighbors of border vertices
		  for (typename std::map<int,GIDSet>::iterator i=borderlinks.begin(); i!=borderlinks.end(); ++i)
			for (typename GIDSet::iterator j=(i->second).begin(); j!=(i->second).end(); ++j)
			  links.insert(P1OperatorLink(i->first,gid2index[*j]));
		  
		  // insert diagonal links for extra DOFs
		  for (int i=0; i<extraDOFs; i++)
			links.insert(P1OperatorLink(vertexmapper.size()+i,vertexmapper.size()+i));		  
		}

	  // Note: links contains now also connections that are standard.
	  // So below we have throw out these connections again!
// 	  std::cout << "=== P1OperatorBase parallel extend overlap " <<  watch.elapsed() << std::endl;

	  partitionType.resize(vertexmapper.size()+extraDOFs);
	  VIterator vendit = gridView.template end<n>();
	  for (VIterator it = gridView.template begin<n>(); it!=vendit; ++it)
	  {
		  partitionType[vertexmapper.map(*it)] = it->partitionType();
	  }
	  
	  VIterator ghostEndIt = gridView.template end<n>();
	  for (VIterator ghostIt = gridView.template begin<n>(); ghostIt!=ghostEndIt; ++ghostIt)
	  {
		  if (ghostIt->partitionType() != GhostEntity)
			  continue;
		  
		  int ghostIndex = vertexmapper.map(*ghostIt);
		  
		  const IdType& ghostGlobalId = index2gid[ghostIndex];
		  
		  VIterator realEndIt = gridView.template end<n>();
		  for (VIterator realIt = gridView.template begin<n>(); realIt!=realEndIt; ++realIt)
		  {
			  if (realIt->partitionType() == InteriorEntity)
				  continue;
			  
			  int realIndex = vertexmapper.map(*realIt);
			  
			  if (realIndex == ghostIndex)
				  continue;
			  
			  const IdType& realGlobalId = index2gid[realIndex];
			  
			  if (realGlobalId == ghostGlobalId) {
				  doubled2Original[std::max(realIndex, ghostIndex)] = std::min(realIndex, ghostIndex);
//				  if (realIt->partitionType() != GhostEntity)
					  break;
			  }
		  }
	  }

//	  for (typename std::map<int,int>::const_iterator it=doubled2Original.begin(); it!=doubled2Original.end(); ++it) 
//		  std::cout << grid.comm().rank() << ": doubled = " << it->first << ", original = " << it->second << std::endl;

	  return true;
	}


	// return number of rows/columns
	int size () const
	{
	  return vertexmapper.size()+extraDOFs;
	}

	struct MatEntry
	{
	  IdType first;
	  BlockType second;
	  MatEntry (const IdType& f, const BlockType& s) : first(f),second(s) {}
	  MatEntry () {}
	};

	// A DataHandle class to exchange matrix entries
	class MatEntryExchange : 
	  public CommDataHandleIF<MatEntryExchange,MatEntry> {
	  typedef typename RepresentationType::RowIterator rowiterator;
	  typedef typename RepresentationType::ColIterator coliterator;
	public:
	  //! export type of data for message buffer
	  typedef MatEntry DataType;

	  //! returns true if data for this codim should be communicated
	  bool contains (int dim, int codim) const
	  {
		return (codim==dim);
	  }

	  //! returns true if size per entity of given dim and codim is a constant
	  bool fixedsize (int dim, int codim) const
	  {
		return false;
	  }

	  /*! how many objects of type DataType have to be sent for a given entity

	  Note: Only the sender side needs to know this size. 
	  */
	  template<class EntityType>
	  size_t size (EntityType& e) const
	  {
		int i=vertexmapper.map(e);
		int n=0;
		for (coliterator j=A[i].begin(); j!=A[i].end(); ++j)
		  n++;
		return n;
	  }

	  //! pack data from user to message buffer
	  template<class MessageBuffer, class EntityType>
	  void gather (MessageBuffer& buff, const EntityType& e) const
	  {
		int i=vertexmapper.map(e);
		for (coliterator j=A[i].begin(); j!=A[i].end(); ++j)
		  {
			typename std::map<int,IdType>::const_iterator it=index2gid.find(j.index());
			if (it==index2gid.end())
			  DUNE_THROW(GridError,"MatEntryExchange::gather(): index not in map");
			buff.write(MatEntry(it->second,*j));
		  }
	  }

	  /*! unpack data from message buffer to user

	  n is the number of objects sent by the sender
	  */
	  template<class MessageBuffer, class EntityType>
	  void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
	  {
		int i=vertexmapper.map(e);
		
		for (size_t k=0; k<n; k++)
		  {
			MatEntry m;
			buff.read(m);
			typename std::map<IdType,int>::const_iterator it=gid2indexNew.find(m.first);
			if (it==gid2indexNew.end())
			  DUNE_THROW(GridError,"MatEntryExchange::scatter(): gid not in map");

			A[i][it->second] += m.second;
		  }
	  }

	  //! constructor
	  MatEntryExchange (const G& g, const std::map<IdType,int>& g2i, 
						const std::map<int,IdType>& i2g,
						const VM& vm,
						RepresentationType& a, std::vector<PartitionType>& pType, std::map<int,int>& d2O)
		: grid(g), gid2index(g2i), index2gid(i2g), vertexmapper(vm), A(a), partitionType(pType), doubled2Original(d2O)
	  {
		  for (typename std::map<IdType,int>::const_iterator it=gid2index.begin(); 
		  		it!=gid2index.end(); ++it)
			  gid2indexNew[it->first] = it->second; 
			  
		  for (typename std::map<int,IdType>::const_iterator it=index2gid.begin(); it!=index2gid.end(); ++it) 
			  if (doubled2Original.find(it->first) != doubled2Original.end())
				  gid2indexNew[it->second] = doubled2Original[it->first];
	  }
 
	private:
	  const G& grid;
	  std::map<IdType,int> gid2indexNew;
	  const std::map<IdType,int>& gid2index;
	  const std::map<int,IdType>& index2gid;
	  const VM& vertexmapper;
	  RepresentationType& A;
	  std::vector<PartitionType>& partitionType;
	  std::map<int,int>& doubled2Original;
	};

  public:

	P1OperatorBase (const G& g, const GV& gridView, LC lcomm, bool extendoverlap=false) 
            : grid(g),gv(gridView),lc(lcomm),vertexmapper(g,gridView.indexSet()),allmapper(g,gridView.indexSet()),links(),
              initialized(init(g,gridView,lc,extendoverlap)),A(size(),size(),nnz(gridView.indexSet()),RepresentationType::random)
		
	{
	  // Check for the TypeTag
	  dune_static_assert((is_same<TypeTag,LeafTag>::value 
              || is_same<TypeTag,LevelTag>::value), "TypeTag must be LeafTag or LevelTag");

	  // be verbose
//  	  std::cout << g.comm().rank() << ": " << "vector size = " << vertexmapper.size() << " + " << extraDOFs << std::endl;
 	  std::cout << g.comm().rank() << ": " << "making " << size() << "x" 
                    << size() << " matrix with " << nnz(gv.indexSet()) << " nonzeros" << std::endl;
//  	  std::cout << g.comm().rank() << ": " << "allmapper has size " << allmapper.size() << std::endl;
//  	  std::cout << g.comm().rank() << ": " << "vertexmapper has size " << vertexmapper.size() << std::endl;
//  	  std::cout << g.comm().rank() << ": " << "hanging nodes=" << hangingnodes << " links=" << links.size() << std::endl;
   
	  // set size of all rows to zero
	  for (size_t i=0; (int) i < (int) gv.indexSet().size(n); i++)
		A.setrowsize(i,0); 

	  // build needs a flag for all entities of all codims
	  std::vector<bool> visited(allmapper.size());
	  for (int i=0; i<allmapper.size(); i++) visited[i] = false;

	  // LOOP 4 : Compute row sizes
	  watch.reset();
	  Iterator eendit = gv.template end<0>();
	  for (Iterator it = gv.template begin<0>(); it!=eendit; ++it)
		{
		  Dune::GeometryType gt = it->type();
		  const typename Dune::ReferenceElementContainer<DT,n>::value_type& 
			refelem = ReferenceElements<DT,n>::general(gt);

		  // vertices, c=n
		  for (int i=0; i<refelem.size(n); i++)
			{
			  int index = allmapper.template map<n>(*it,i);
			  int alpha = vertexmapper.template map<n>(*it,i);
			  //			  std::cout << "index=" << index << " alpha=" << alpha << std::endl;
			  if (!visited[index]) 
				{
				  A.incrementrowsize(alpha);
				  visited[index] = true;
// 				  printf("increment row %04d\n",alpha);
				}
			}

		  // edges for all element types, c=n-1
		  for (int i=0; i<refelem.size(n-1); i++)
			{
			  int index = allmapper.template map<n-1>(*it,i);
			  int alphaOld, betaOld;
			  int alpha = alphaOld = vertexmapper.template map<n>(*it,refelem.subEntity(i,n-1,0,n));
			  int beta = betaOld = vertexmapper.template map<n>(*it,refelem.subEntity(i,n-1,1,n));
			  if (doubled2Original.find(alpha) != doubled2Original.end())
				  alpha = doubled2Original[alpha];
			  if (doubled2Original.find(beta) != doubled2Original.end())
				  beta = doubled2Original[beta];
			  if (!visited[index]) 
				{
				  A.incrementrowsize(alpha);
				  A.incrementrowsize(beta);
				  visited[index] = true;
				  if (hangingnodes>0 || extendOverlap) // delete standard links
					{
					  links.erase(P1OperatorLink(alpha,beta));
					  links.erase(P1OperatorLink(beta,alpha));
					}
				}
			}

		  // for codim n-2 to 0 we need a template metaprogram
		  if (!gt.isSimplex())
			P1Operator_meta<n,n-2>::addrowscube(*it,vertexmapper,allmapper,refelem,A,visited,hangingnodes+(extendOverlap),links,doubled2Original);
		}

	  // additional links due to hanging nodes
	  for (typename std::set<P1OperatorLink>::iterator i=links.begin(); i!=links.end(); ++i)
		A.incrementrowsize(i->first);

	  // now the row sizes have been set
	  A.endrowsizes();
// 	  std::cout << "=== P1OperatorBase compute row sizes " <<  watch.elapsed() << std::endl;

	  // clear the flags for the next round, actually that is not necessary because addindex takes care of this
	  for (int i=0; i<allmapper.size(); i++) visited[i] = false;

	  // LOOP 5 : insert the nonzeros
	  watch.reset();
	  for (Iterator it = gv.template begin<0>(); it!=eendit; ++it)
		{
		  Dune::GeometryType gt = it->type();
		  const typename Dune::ReferenceElementContainer<DT,n>::value_type&
			refelem = ReferenceElements<DT,n>::general(gt);
		  // 		  std::cout << "ELEM " << GeometryName(gt) << std::endl;

		  // vertices, c=n
		  for (int i=0; i<refelem.size(n); i++)
			{
			  int index = allmapper.template map<n>(*it,i);
			  int alpha = vertexmapper.template map<n>(*it,i);
			  // 			  std::cout << "vertex allindex " << index << std::endl; 
			  if (!visited[index]) 
				{
				  A.addindex(alpha,alpha);
				  visited[index] = true;
				}
			}

		  // edges for all element types, c=n-1
		  for (int i=0; i<refelem.size(n-1); i++)
			{
			  int index = allmapper.template map<n-1>(*it,i);
			  // 			  std::cout << "edge allindex " << index << std::endl; 
			  if (!visited[index]) 
				{
				  int alpha = vertexmapper.template map<n>(*it,refelem.subEntity(i,n-1,0,n));
				  int beta = vertexmapper.template map<n>(*it,refelem.subEntity(i,n-1,1,n));
				  if (doubled2Original.find(alpha) != doubled2Original.end())
					  alpha = doubled2Original[alpha];
				  if (doubled2Original.find(beta) != doubled2Original.end())
					  beta = doubled2Original[beta];
				  A.addindex(alpha,beta);
				  A.addindex(beta,alpha);
				  visited[index] = true;
// 				  printf("adding (%04d,%04d) index=%04d\n",alpha,beta,index);
// 				  printf("adding (%04d,%04d) index=%04d\n",beta,alpha,index);
				}
			}

		  // for codim n-2 to 0 we need a template metaprogram
		  if (!gt.isSimplex())
			P1Operator_meta<n,n-2>::addindicescube(*it,vertexmapper,allmapper,refelem,A,visited,doubled2Original);
		}

	  // additional links due to hanging nodes
	  for (typename std::set<P1OperatorLink>::iterator i=links.begin(); i!=links.end(); ++i) {
		  int indexI = i->first;
		  int indexJ = i->second;
		  if (doubled2Original.find(indexI) != doubled2Original.end())
			  indexI = doubled2Original[indexI];
		  if (doubled2Original.find(indexJ) != doubled2Original.end())
			  indexJ = doubled2Original[indexJ];
		A.addindex(indexI,indexJ);
	  }
	  
	  // now the matrix is ready for use
	  A.endindices();
// 	  std::cout << "=== P1OperatorBase index insertion " <<  watch.elapsed() << std::endl;

	  // delete additional links
	  links.clear();

// 	  std::cout << grid.comm().rank() << ": " << "matrix initialized" << std::endl;
	}

	//! return const reference to operator matrix
	const RepresentationType& operator* () const
	{
	  return A;
	}

	//! return reference to operator matrix
	RepresentationType& operator* ()
	{
	  return A;
	}

	//! makes matrix consistent in nonoverlapping case
	void sumEntries ()
	{
	  if (!extendOverlap) return;

	  // build forward map
	  //std::map<int,IdType> index2gid;
	  for (typename std::map<IdType,int>::iterator i=gid2index.begin(); i!=gid2index.end(); ++i)
		index2gid[i->second] = i->first;

	  // communicate matrix entries
	  MatEntryExchange datahandle(grid,gid2index,index2gid,vertexmapper,A,partitionType,doubled2Original);
	  lc.template communicate<MatEntryExchange>(datahandle, InteriorBorder_InteriorBorder_Interface, ForwardCommunication);
	}

  protected:
	Timer watch;
	const G& grid;	
	const GV gv;
	LC lc;
	VM vertexmapper;
	AM allmapper;
	std::vector<bool> hanging;
	std::set<P1OperatorLink> links; 
	int hangingnodes;
	bool extendOverlap;
	int extraDOFs;
	std::map<IdType,int> gid2index;
	std::map<int,IdType> index2gid;
	std::vector<PartitionType> partitionType;
	std::map<int,int> doubled2Original;
	bool initialized;
	RepresentationType A;
  };




  /*! @brief Extends P1OperatorBase by a generic methods to assemble global stiffness matrix from local stiffness matrices
   *
   *
   * The template parameter TypeTag describes what kind of Assembler we are. There two choices:
   * <dt>LevelTag</dt><dd>We assemble on a grid level. </dd>
   * <dt>LeafTag</dt><dd>We assemble on the leaf entities of the grid</dd>
   */
  template<typename TypeTag, class G, class RT, class GV, class LC, int m>
  class P1OperatorAssembler : public P1OperatorBase<TypeTag,G,RT,GV,LC,m>
  {
	typedef typename G::ctype DT;
	enum {n=G::dimension};
      typedef typename GV::IndexSet IS;
	typedef typename G::template Codim<0>::Entity Entity;
	typedef typename GV::template Codim<0>::Iterator Iterator;
	typedef typename GV::template Codim<n>::Iterator VIterator;
	typedef typename G::template Codim<0>::HierarchicIterator HierarchicIterator;
	typedef typename G::template Codim<0>::EntityPointer EEntityPointer;
	typedef typename P1FunctionExtended<GV,RT,LC,m>::RepresentationType VectorType;
	typedef typename VectorType::block_type VBlockType;
	typedef typename P1OperatorBase<TypeTag,G,RT,GV,LC,m>::RepresentationType MatrixType;
        typedef typename MatrixType::field_type MFieldType;
	typedef typename MatrixType::block_type MBlockType;
	typedef typename MatrixType::RowIterator rowiterator;
	typedef typename MatrixType::ColIterator coliterator;
	typedef typename P1OperatorBase<TypeTag,G,RT,GV,LC,m>::VM VM;
    typedef array<BoundaryConditions::Flags,m> BCBlockType;     // componentwise boundary conditions


	// A DataHandle class to exchange border rows
	class AccumulateBCFlags : 
      public CommDataHandleIF<AccumulateBCFlags,std::pair<BCBlockType,VBlockType> > {
	public:
	  //! export type of data for message buffer
	  typedef std::pair<BCBlockType,VBlockType> DataType; // BC flags are stored in a char per vertex

	  //! returns true if data for this codim should be communicated
	  bool contains (int dim, int codim) const
	  {
		return (codim==dim);
	  }

	  //! returns true if size per entity of given dim and codim is a constant
	  bool fixedsize (int dim, int codim) const
	  {
		return true;
	  }

	  /*! how many objects of type DataType have to be sent for a given entity

	  Note: Only the sender side needs to know this size. 
	  */
	  template<class EntityType>
	  size_t size (EntityType& e) const
	  {
		return 1;
	  }

	  //! pack data from user to message buffer
	  template<class MessageBuffer, class EntityType>
	  void gather (MessageBuffer& buff, const EntityType& e) const
	  {
		int alpha = vertexmapper.map(e);
		buff.write(DataType(essential[alpha],f[alpha])); 
	  }

	  /*! unpack data from message buffer to user

	  n is the number of objects sent by the sender
	  */
	  template<class MessageBuffer, class EntityType>
	  void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
	  {
		DataType x;
		buff.read(x);
		int alpha = vertexmapper.map(e);
		for (int i=0; i<m; i++)
		  if (x.first[i]>essential[alpha][i])
			{
			  essential[alpha][i] = x.first[i];
			  f[alpha][i] = x.second[i];
			}
	  }

	  //! constructor
	  AccumulateBCFlags (const G& g, const VM& vm, std::vector<BCBlockType>& e, VectorType& _f) 
		: grid(g), essential(e), vertexmapper(vm), f(_f)
	  {}
 
	private:
	  const G& grid;
	  std::vector<BCBlockType>& essential;
	  const VM& vertexmapper;
	  VectorType& f;
	};



  public:
	P1OperatorAssembler (const G& g, const GV& gridView, LC lcomm, bool extendoverlap=false)
	  : P1OperatorBase<TypeTag,G,RT,GV,LC,m>(g,gridView,lcomm,extendoverlap)
	{	}


	/*! @brief Assemble global stiffness matrix 

	This method takes an object that can compute local stiffness matrices and
    assembles the global linear system Au=f.
      
	  @param[in] loc    the local assembler providing element stiffness and boundary conditions for all elements  
	  @param[in,out] u  solution, contains initial values on input, Dirichlet values are set. The
                        type of boundary condition for a node is inferred from the values returned
                        by the local assembler. A node is of Neumann type if all elements referring
                        to that node report a Neumann boundary condition, it is set to Dirichlet
                        if a least one element reports a process or Dirichlet boundary condition. The difference
						between process and Dirichlet is that process always denotes a homogeneous Dirichlet
                        value. 
	  @param[in] f      right hand side is filled by this method

    Note that the rows corresponding to nodes at the Dirichlet boundary are filled
    with trivial equations of the form \f[1\cdot u_i = f_i \f] where \f$u_i\f$ and \f$f_i\f$ are both set to the 
    Dirichlet value at the \f$i\f$th node.

	 */
        template<class I>
	void assemble (LocalStiffness<I,G,RT,m>& loc, P1FunctionExtended<GV,RT,LC,m>& u, P1FunctionExtended<GV,RT,LC,m>& f)
	{
	  // check size
 	  if ((*u).N()!=this->A.M() || (*f).N()!=this->A.N())
 		DUNE_THROW(MathError,"P1OperatorAssembler::assemble(): size mismatch");		

	  // clear global stiffness matrix and right hand side
	  this->watch.reset();
	  this->A = 0;
	  *f = 0;
// 	  std::cout << "=== P1OperatorBase clear matrix " <<  this->watch.elapsed() << std::endl;

	  // allocate flag vector to hold flags for essential boundary conditions
	  std::vector<BCBlockType> essential(this->vertexmapper.size());
	  std::vector<VBlockType> dirichletIndexGlobal(this->vertexmapper.size());
	  for (typename std::vector<BCBlockType>::size_type i=0; i<essential.size(); i++) {
		essential[i].assign(BoundaryConditions::neumann);
		dirichletIndexGlobal[i] = 0;
	  }

	  // allocate flag vector to note hanging nodes whose row has been assembled
	  std::vector<unsigned char> treated(this->vertexmapper.size());
	  for (std::vector<unsigned char>::size_type
             i=0; i<treated.size(); i++) treated[i] = false;

	  // hanging node stuff
	  RT alpha[Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize][Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize];
	  
	  // local to global id mapping (do not ask vertex mapper repeatedly
	  int l2g[Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize];
	  int fl2g[Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::maxsize];

	  // run over all leaf elements
	  Iterator eendit = this->gv.template end<0>();
	  for (Iterator it = this->gv.template begin<0>(); it!=eendit; ++it)
		{
		  // parallelization
		  if (it->partitionType()==GhostEntity)
			continue;

		  // get access to shape functions for P1 elements
		  Dune::GeometryType gt = it->type();
		  const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type& 
			sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1);

		  // get local to global id map
		  for (int k=0; k<sfs.size(); k++)
			{
			  if (sfs[k].codim()!=n) DUNE_THROW(MathError,"expected codim==dim");
			  int alpha = this->vertexmapper.template map<n>(*it,sfs[k].entity());
			  if (this->doubled2Original.find(alpha) != this->doubled2Original.end())
				  alpha = this->doubled2Original[alpha];
			  l2g[k] = alpha;
			}

		  // build local stiffness matrix for P1 elements
		  // inludes rhs and boundary condition information
		  loc.template assemble<TypeTag>(*it); // assemble local stiffness matrix

		  // assemble constraints in hanging nodes
		  // as a by-product determine the interpolation factors WITH RESPECT TO NODES OF FATHER
		  bool hasHangingNodes = false;
		  for (int k=0; k<sfs.size(); k++) // loop over rows, i.e. test functions
			{
			  // process only hanging nodes
			  if (!this->hanging[l2g[k]]) continue;
			  hasHangingNodes = true;

			  // determine position of hanging node in father
			  EEntityPointer father=it->father(); // the father element
			  GeometryType gtf = father->type(); // fathers type
			  assert(gtf==gt); // in hanging node refinement the element type is preserved
			  const FieldVector<DT,n>& cpos=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1)[k].position();
			  FieldVector<DT,n> pos = it->geometryInFather().global(cpos); // map corner to father element

			  // evaluate interpolation factors and local to global mapping for father
			  for (int i=0; i<Dune::LagrangeShapeFunctions<DT,RT,n>::general(gtf,1).size(); ++i)
				{
				  alpha[i][k] = Dune::LagrangeShapeFunctions<DT,RT,n>::general(gtf,1)[i].evaluateFunction(0,pos);
				  fl2g[i] = this->vertexmapper.template map<n>(*father,i);
				  if (this->doubled2Original.find(fl2g[i]) != this->doubled2Original.end())
					  fl2g[i] = this->doubled2Original[fl2g[i]];
				}						  

			  // assemble the constraint row once
			  if (!treated[l2g[k]])
				{
				  bool throwflag=false;
				  int cnt=0;
				  for (int i=0; i<Dune::LagrangeShapeFunctions<DT,RT,n>::general(gtf,1).size(); ++i)
					if (std::abs(alpha[i][k])>1E-4)
					  {
						if (this->hanging[fl2g[i]])
						  {
							std::cout << "X i=" << father->geometry()[i]
									  << " k=" << it->geometry()[k]
									  << " alpha=" << alpha[i][k]
									  << " cnt=" << ++cnt
									  << std::endl;
							throwflag = true;
						  }
						this->A[l2g[k]][fl2g[i]] = MFieldType(); // Note: Interpolation is done a posteriori
					  }
				  if (throwflag)
					DUNE_THROW(GridError,"hanging node interpolated from hanging node");
				  for (int comp=0; comp<m; comp++)
					this->A[l2g[k]][l2g[k]][comp][comp] = static_cast<MFieldType>(1);
				  (*f)[l2g[k]] = 0;
				  (*u)[l2g[k]] = 0;
				  treated[l2g[k]] = true;
				}
			}

		  // accumulate local matrix into global matrix for non-hanging nodes
		  for (int i=0; i<sfs.size(); i++) // loop over rows, i.e. test functions
			{
			  // process only non-hanging nodes
			  if (this->hanging[l2g[i]]) continue;

			  // accumulate matrix
			  for (int j=0; j<sfs.size(); j++)
				{
				  // process only non-hanging nodes
				  if (this->hanging[l2g[j]]) continue;

				  // the standard entry
				  this->A[l2g[i]][l2g[j]] += loc.mat(i,j);
				}
			  
			  // essential boundary condition and rhs
			  for (int comp=0; comp<m; comp++)
				{
				  if (loc.bc(i)[comp]>essential[l2g[i]][comp])
					{
					  essential[l2g[i]][comp] = loc.bc(i)[comp];
					  dirichletIndexGlobal[l2g[i]][comp] = loc.dirichletIdx(i)[comp];
					  (*f)[l2g[i]][comp] = loc.rhs(i)[comp];
					}
				  if (essential[l2g[i]][comp]==BoundaryConditions::neumann)
					(*f)[l2g[i]][comp] += loc.rhs(i)[comp];
				}
			}

		  // add corrections for hanging nodes
		  if (hasHangingNodes)
			{
			  // a map that inverts l2g
			  std::map<int,int> g2l;
			  for (int i=0; i<sfs.size(); i++)
				g2l[l2g[i]] = i;

			  // a map that inverts fl2g
			  std::map<int,int> fg2l;
			  for (int i=0; i<sfs.size(); i++)
				fg2l[fl2g[i]] = i;

			  EEntityPointer father=it->father(); // DEBUG

			  // loop over nodes (rows) in father, assume father has same geometry type
			  for (int i=0; i<sfs.size(); ++i)
				{
				  // corrections to matrix
				  // loop over all nodes (columns) in father, assume father has same geometry type
				  for (int j=0; j<sfs.size(); ++j)
					{
					  // loop over hanging nodes
					  for (int k=0; k<sfs.size(); k++)
						{
						  // process only hanging nodes
						  if (!this->hanging[l2g[k]]) continue;

						  // first term, with i a coarse vertex
						  if ( std::abs(alpha[j][k])>1E-4 && g2l.find(fl2g[i])!=g2l.end() )
							{
							  accumulate(this->A[fl2g[i]][fl2g[j]],alpha[j][k],loc.mat(g2l[fl2g[i]],k));
							}

						  // first term, with i a fine non-hanging vertex
						  if ( std::abs(alpha[j][k])>1E-4 && fg2l.find(l2g[i])==fg2l.end() && !this->hanging[l2g[i]])
							{
							  accumulate(this->A[l2g[i]][fl2g[j]],alpha[j][k],loc.mat(i,k));
							}

						  // second term, with j a coarse vertex
						  if ( std::abs(alpha[i][k])>1E-4 && g2l.find(fl2g[j])!=g2l.end() )
							{
							  accumulate(this->A[fl2g[i]][fl2g[j]],alpha[i][k],loc.mat(k,g2l[fl2g[j]]));
							}

						  // second term, with j a fine non-hanging vertex
						  if ( std::abs(alpha[i][k])>1E-4 && fg2l.find(l2g[j])==fg2l.end() && !this->hanging[l2g[j]])
							{
							  accumulate(this->A[fl2g[i]][l2g[j]],alpha[i][k],loc.mat(k,j));
							}

						  // third term, loop over hanging nodes
						  for (int l=0; l<sfs.size(); l++)
							{
							  // process only hanging nodes
							  if (!this->hanging[l2g[l]]) continue;

							  if ( std::abs(alpha[i][k])>1E-4 && std::abs(alpha[j][l])>1E-4 )
								{
								  accumulate(this->A[fl2g[i]][fl2g[j]],alpha[i][k]*alpha[j][l],loc.mat(k,l));
								}
							}
						}
					}

				  // corrections to rhs
				  for (int comp=0; comp<m; comp++)
					if (essential[fl2g[i]][comp]==BoundaryConditions::neumann)
					  for (int k=0; k<sfs.size(); k++)
						{
						  // process only hanging nodes
						  if (!this->hanging[l2g[k]]) continue;
						  
						  if ( std::abs(alpha[i][k])>1E-4 && loc.bc(k)[comp]==BoundaryConditions::neumann )
							(*f)[fl2g[i]][comp] += alpha[i][k]*loc.rhs(k)[comp];
						}
				}
			}
		}

	  // accumulate matrix entries, should do also for systems ...
 	  if (this->extendOverlap)
 		this->sumEntries();

	  // send around boundary conditions
 	  if (this->extendOverlap)
		{
		  AccumulateBCFlags datahandle(this->grid,this->vertexmapper,essential,*f);
		  this->lc.template communicate<AccumulateBCFlags>(datahandle,InteriorBorder_InteriorBorder_Interface,ForwardCommunication);
		  // on a nonvoerlapping grid we have the correct BC at all interior and border vertices
		}
	  // what about the overlapping case ?

	  // muck up ghost vertices
	  VIterator vendit = this->gv.template end<n>();
	  for (VIterator it = this->gv.template begin<n>(); it!=vendit; ++it)
		if (it->partitionType()==GhostEntity)
		  {
			// index of this vertex
			int i=this->vertexmapper.map(*it);

			// muck up row
			coliterator endj=this->A[i].end();
			for (coliterator j=this->A[i].begin(); j!=endj; ++j)
			  {
				(*j) = MFieldType();
				if ((int) j.index()==i)
				  for (int comp=0; comp<m; comp++)
					(*j)[comp][comp] = static_cast<MFieldType>(1);
			  }
			(*f)[i] = 0;
		  }

	  // put in essential boundary conditions
	  rowiterator endi=this->A.end();
	  for (rowiterator i=this->A.begin(); i!=endi; ++i)
		{
		  // muck up extra rows
            if ((int) i.index()>=this->vertexmapper.size())
			{
			  coliterator endj=(*i).end();
			  for (coliterator j=(*i).begin(); j!=endj; ++j)
				{
				  (*j) = MFieldType();
				  if (j.index()==i.index())
					for (int comp=0; comp<m; comp++)
					  (*j)[comp][comp] = static_cast<MFieldType>(1);
				}
			  (*f)[i.index()] = 0;
			  continue;
			}

		  // insert dirichlet ans processor boundary conditions
		  if (!this->hanging[i.index()])
			{
			  for (int icomp=0; icomp<m; icomp++)
				if (essential[i.index()][icomp]!=BoundaryConditions::neumann)
				  {
					coliterator endj=(*i).end();
					for (coliterator j=(*i).begin(); j!=endj; ++j)
					  if (j.index()==i.index())
						{
						  for (int jcomp=0; jcomp<m; jcomp++)
							if (jcomp==dirichletIndexGlobal[i.index()][icomp])
							//if (icomp==jcomp)
							  (*j)[icomp][jcomp] = static_cast<MFieldType>(1);
							else
							  (*j)[icomp][jcomp] = MFieldType();							  
						}
					  else
						{
						  for (int jcomp=0; jcomp<m; jcomp++)
							(*j)[icomp][jcomp] = MFieldType();
						}
					(*u)[i.index()][icomp] = (*f)[i.index()][icomp];
				  }
			}
		}
	} 

	//! assemble operator, rhs and Dirichlet boundary conditions
	void interpolateHangingNodes (P1FunctionExtended<GV,RT,LC,m>& u)
	{
	  // allocate flag vector to note hanging nodes whose row has been assembled
	  std::vector<unsigned char> treated(this->vertexmapper.size());
	  for (std::size_t i=0; i<treated.size(); i++) treated[i] = false;

	  // run over all leaf elements
	  Iterator eendit = this->gv.template end<0>();
	  for (Iterator it = this->gv.template begin<0>(); it!=eendit; ++it)
		{
		  // get access to shape functions for P1 elements
		  Dune::GeometryType gt = it->type();
		  const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type& 
			sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1);

		  // determine the interpolation factors WITH RESPECT TO NODES OF FATHER
		  for (int k=0; k<sfs.size(); k++) // loop over rows, i.e. test functions
			{
			  int alpha = this->vertexmapper.template map<n>(*it,sfs[k].entity());
			  if (this->hanging[alpha] && !treated[alpha])
				{
				  // determine position of hanging node in father
				  EEntityPointer father=it->father(); // the father element
				  GeometryType gtf = father->type(); // fathers type
				  assert(gtf==gt); // in hanging node refinement the element type is preserved
				  const FieldVector<DT,n>& cpos=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1)[k].position();
				  FieldVector<DT,n> pos = it->geometryInFather().global(cpos); // map corner to father element

				  // evaluate interpolation factors and local to global mapping for father
				  VBlockType value; value=0;
				  for (int i=0; i<Dune::LagrangeShapeFunctions<DT,RT,n>::general(gtf,1).size(); ++i)
					{
					  int beta = this->vertexmapper.template map<n>(*father,i);
					  value.axpy(Dune::LagrangeShapeFunctions<DT,RT,n>::general(gtf,1)[i].evaluateFunction(0,pos),(*u)[beta]);
					}
				  (*u)[alpha] = value;
				  treated[alpha] = true;
				}
			}
		}
	}

	void preMark ()
	{
	  marked.resize(this->vertexmapper.size());
	  for (std::size_t i=0; i<marked.size(); i++) marked[i] = false;
	  return;
	}

	void postMark (G& g)
	{
	  // run over all leaf elements
	  int extra=0;
	  Iterator eendit = this->gv.template end<0>();
	  for (Iterator it = this->gv.template begin<0>(); it!=eendit; ++it)
		{
		  // get access to shape functions for P1 elements
		  Dune::GeometryType gt = it->type();
		  //		  if (gt!=Dune::simplex && gt!=Dune::triangle && gt!=Dune::tetrahedron) continue;

		  const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type& 
			sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1);

		  // count nodes with mark
		  int count=0;
		  for (int k=0; k<sfs.size(); k++)
			{
			  int alpha = this->vertexmapper.template map<n>(*it,sfs[k].entity());
			  if (marked[alpha]) count++;
			}

		  // refine if a marked edge exists
		  if (count>0) {
			extra++;
			g.mark(1,it);
		  }
		}

	  marked.clear();
	  return;
	}

	void mark (G& g, EEntityPointer& it)
	{
	  // refine this element
	  g.mark(1,it);

	  // check geom type, exit if not simplex
	  Dune::GeometryType gt = it->type();
	  //if (gt!=Dune::simplex && gt!=Dune::triangle && gt!=Dune::tetrahedron) return;
	  assert(it->isLeaf());

	  // determine if element has hanging nodes
	  bool hasHangingNodes = false;
	  const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type& 
		sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1);
	  for (int k=0; k<sfs.size(); k++) // loop over rows, i.e. test functions
		if (this->hanging[this->vertexmapper.template map<n>(*it,sfs[k].entity())])
		  {
			hasHangingNodes = true;
			break;
		  }

	  // if no hanging nodes we are done
	  if (!hasHangingNodes) return;

	  // mark all corners of father and element
	  EEntityPointer father=it->father();
	  for (int k=0; k<sfs.size(); k++) // same geometry type ...
		{
		  int alpha=this->vertexmapper.template map<n>(*father,k);
		  marked[alpha] = true;
		}
	  for (int k=0; k<sfs.size(); k++) // same geometry type ...
		{
		  int alpha=this->vertexmapper.template map<n>(*it,k);
		  marked[alpha] = true;
		}

	  return;
	}

	
  private:
	std::vector<bool> marked;

	template<class M, class K>
	void accumulate (M& A, const K& alpha, const M& B)
	{
	  for (std::size_t i=0; i<A.N(); i++)
		for (std::size_t j=0; j<A.M(); j++)
		  A[i][j] += alpha*B[i][j];
	}
  };
  

  /*! @brief Leafwise assembler

  This class serves as a base class for local assemblers. It provides
  space and access to the local stiffness matrix. The actual assembling is done
  in a derived class via the virtual assemble method.

  The template parameters are:

  - G    A grid type
  - RT   The field type used in the elements of the stiffness matrix
  - m    number of degrees of freedom per node (system size)
   */
  template<class G, class RT, int m>
  class LeafP1OperatorAssembler : public P1OperatorAssembler<LeafTag,G,RT,typename G::LeafGridView,LeafCommunicate<G>,m>
  {
  public:
	LeafP1OperatorAssembler (const G& grid, bool extendoverlap=false) 
	  : P1OperatorAssembler<LeafTag,G,RT,typename G::LeafGridView,LeafCommunicate<G>,m>(grid,grid.leafView(),LeafCommunicate<G>(grid),extendoverlap)
	{}
  };


  /*! @brief Levelwise assembler

  This class serves as a base class for local assemblers. It provides
  space and access to the local stiffness matrix. The actual assembling is done
  in a derived class via the virtual assemble method.

  The template parameters are:

  - G    A grid type
  - RT   The field type used in the elements of the stiffness matrix
  - m    number of degrees of freedom per node (system size)
   */
  template<class G, class RT, int m>
  class LevelP1OperatorAssembler : public P1OperatorAssembler<LevelTag,G,RT,typename G::LevelGridView,LevelCommunicate<G>,m>
  {
  public:
	LevelP1OperatorAssembler (const G& grid, int level, bool extendoverlap=false) 
	  : P1OperatorAssembler<LevelTag,G,RT,typename G::LevelGridView,LevelCommunicate<G>,m>(grid,grid.levelView(level),LevelCommunicate<G>(grid,level),extendoverlap)
	{}
  };


  /** @} */

}
#endif
