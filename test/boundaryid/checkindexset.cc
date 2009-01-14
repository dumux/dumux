#ifndef DUNE_CHECKINDEXSET_CC
#define DUNE_CHECKINDEXSET_CC

#include "config.h"
#include <iostream>
#include <algorithm>

#include <dune/common/fvector.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/referenceelements.hh>


#include <map>
#include <set>
#include <vector>
#include <limits>

/** @file
  @author Robert Kloefkorn
  @brief Provides a check of the grids index set.
*/

namespace Dune {

// compare 2 FieldVectors
template <typename ctype, int dim>
bool compareVec(const FieldVector<ctype,dim> & vx1 , const FieldVector<ctype,dim> & vx2 )
{
  const ctype eps = 1e5 * std::numeric_limits<ctype>::epsilon();
  bool comp = true;
  for(int i=0;i<dim; i++)
  {
    if(std::abs( vx1[i] - vx2[i] ) > eps ) comp = false;
  }
  return comp;
}

// check som functionality of grid
template <int codim, class GridType, class EntityType,
  class IndexSetType, class OutputStreamImp,
  class MapType1 , class MapType2 , class MapType3 >
void checkSubEntity ( const GridType & grid,
    const EntityType &en , const IndexSetType & lset,
    OutputStreamImp & sout , MapType1 & subEntities , MapType2 & vertices ,
    MapType3 & vertexCoordsMap )
{
  GeometryType type = en.geometry().type();
  enum { dim = EntityType::dimension };
  typedef typename EntityType :: ctype coordType;

    const ReferenceElement< coordType, dim > & refElem =
      ReferenceElements< coordType, dim >::general(type);


      // check all subEntities of codimension  codim
      if(en.template count<codim>() != refElem.size(0,0,codim))
      {
        std::cerr << "entity index = " << lset.index(en)
                  << ", type = " << type
                  << std::endl
                  << "codim = " << codim
                  << std::endl
                  << "count<codim>() = " << en.template count<codim>()
                  << std::endl
                  << "refElem.size(codim) = " << refElem.size(0,0,codim)
                  << std::endl;
        DUNE_THROW(GridError,
                   "wrong number of subEntities of codim " << codim);
      }

      for(int subEntity = 0; subEntity < refElem.size(0,0,codim); subEntity++)
      {
        typedef std::pair < int , GeometryType > SubEntityKeyType;
        {
          int numSubEntities = refElem.size(subEntity,codim,dim);
          // every entity have at least one vertex
          assert( numSubEntities > 0 );

          // create vectors of number of vertices on sub entity
          std::vector<int> local (numSubEntities,-1);
          std::vector<int> global(numSubEntities,-1);

          for(int j=0 ;j<numSubEntities; j++ )
            local[j] = refElem.subEntity(subEntity , codim , j , dim );


          sout << numSubEntities << " Vertices on subEntity<codim=" << codim << ">\n";
          sout << "check suben [";
          for(int j=0 ;j<numSubEntities-1; j++ )
            sout << local[j] <<  ", ";
          sout << local[numSubEntities-1] << "]\n";

          for(int j=0 ;j<numSubEntities; j++ )
          {
            global[j] = lset.template subIndex<dim> ( en, local[j]);
          }

          SubEntityKeyType globalSubEntity =
            SubEntityKeyType ( lset.template subIndex<codim>(en,subEntity),
                               (en.template entity<codim> (subEntity))->geometry().type());
          assert( globalSubEntity.first >= 0 );
          sout << "local subentity " << subEntity << " consider subentity with global key (" << globalSubEntity.first << "," << globalSubEntity.second << ") on en = " << lset.index(en) << "\n";

          sout << "Found global numbers of entity [ ";
          for(int j=0 ;j<numSubEntities; j++ )
          {
            sout << global[j] << " ";
          }
          sout << "]\n";

          for(int j=0; j<numSubEntities; j++)
          {
            {
              // get entity pointer of sub entity codim=dim (Vertex)
              typedef typename GridType :: template Codim<dim> :: EntityPointer VertexPointerType;
              VertexPointerType vxp = en.template entity<dim> (local[j]);

              FieldVector<coordType,dim> vx ( vxp->geometry().corner(0));
              if(vertexCoordsMap.find(global[j]) != vertexCoordsMap.end())
              {
                FieldVector<coordType,dim> vxcheck ( vertexCoordsMap[global[j]] );
                if( ! compareVec( vxcheck, vx ) )
                {
                  std::cerr << "ERROR map global vertex [" << global[j] << "] vx " << vxcheck << " is not " << vx << "\n";
                  assert( compareVec( vxcheck, vx ) );
                }
              }
            }

            typedef typename GridType :: template Codim<codim> :: EntityPointer SubEnPointerType;
            SubEnPointerType subenp = en.template entity<codim> (subEntity);

            // assert that all sub entities have the same level
            // otherwise one of the theoretical conditions is violated
            assert( subenp.level() == en.level() );

            FieldVector<coordType,dim> vx ( subenp->geometry().corner(j));
            if(vertexCoordsMap.find(global[j]) != vertexCoordsMap.end())
            {
              FieldVector<coordType,dim> vxcheck ( vertexCoordsMap[global[j]] );
              if( ! compareVec( vxcheck, vx ) )
              {
                std::cerr << "Error map global vertex [" << global[j] << "] vx " << vxcheck << " is not " << vx << "\n";
                assert( compareVec( vxcheck, vx ) );
              }
            }
            sout << "vx[" << global[j] << "] = "  << vx << "\n";
          }
            sout << "sort vector of global vertex\n";

          // sort vector of global vertex number for storage in map
          // the smallest entry is the first entry
          std::sort( global.begin(), global.end() );

          // check whether vertex key is already stored in map
          if(vertices.find(global) == vertices.end())
          {
              vertices[global] = globalSubEntity;
          }
          else
          {
              SubEntityKeyType otherSubEntity = vertices[global];
              assert( globalSubEntity == otherSubEntity );
          }

          // check whether subEntity is already stored in map
          if(subEntities.find(globalSubEntity) == subEntities.end() )
          {
            subEntities[globalSubEntity] = global;
          }
          else
          {
            std::vector<int> globalcheck = subEntities[globalSubEntity];
            if(! (global == globalcheck ))
            {
              std::cerr << "For subEntity key (" << globalSubEntity.first << "," << globalSubEntity.second << ") \n";
              std::cerr << "Got ";
              for(int j=0 ;j<numSubEntities; j++ )
              {
                std::cerr << global[j] << " ";
              }
              std::cerr << "\n";
              std::cerr << "Found ";
              for(int j=0 ;j<numSubEntities; j++ )
              {
                std::cerr << globalcheck [j] << " ";
              }
              std::cerr << "\n";
              DUNE_THROW(Dune::GridError, "global != globalcheck");
            }
          }
        }
      } // end check sub entities
      sout << "end check sub entities\n";
}

// check some functionality of grid
template <int codim, class GridType,
          class IndexSetType, class OutputStreamImp >
void checkIndexSetForCodim ( const GridType &grid , const IndexSetType & lset,
    OutputStreamImp & sout , bool levelIndex )
{
  enum { dim = GridType :: dimension };

  sout <<"\n\nStart consistency check of index set \n\n";
  typedef typename IndexSetType :: template Codim<0>::template Partition<All_Partition> :: Iterator Iterator;
  typedef typename GridType :: ctype coordType;
  typedef typename GridType :: template Codim<0>:: Entity EntityCodim0Type;

  // ////////////////////////////////////////////////////////////////
  //   Check whether geomTypes() returns correct result
  // ////////////////////////////////////////////////////////////////
  typedef typename IndexSetType :: template Codim<codim>::
      template Partition<All_Partition> :: Iterator IteratorType;

  IteratorType endit  = lset.template end<codim,All_Partition>   ();
  IteratorType it = lset.template begin<codim,All_Partition> ();

  std::set<GeometryType> geometryTypes;

  if (it == endit) return;

  for (; it!=endit; ++it)
      geometryTypes.insert(it->geometry().type());

  bool geomTypesError = false;
  // Check whether all entries in the official geometry types list are contained in our self-computed one
  for (size_t i=0; i<lset.geomTypes(codim).size(); i++)
      if (geometryTypes.find(lset.geomTypes(codim)[i])==geometryTypes.end())
          geomTypesError = true;


  // And vice versa
  for (std::set<GeometryType>::iterator it = geometryTypes.begin(); it!=geometryTypes.end(); ++it) {
      bool found = false;
      for (size_t i=0; i<lset.geomTypes(codim).size(); i++)
          if (*it == lset.geomTypes(codim)[i]) {
              found = true;
              break;
          }

      if (!found)
          geomTypesError = true;

  }

  if (geomTypesError) {

      std::cerr << "There is a mismatch in the list of geometry types of codim " << codim << "." << std::endl;
      std::cerr << "Geometry types present in the grid are:" << std::endl;
      for (std::set<GeometryType>::iterator it = geometryTypes.begin(); it!=geometryTypes.end(); ++it)
          std::cerr << "  " << *it << std::endl;

      std::cerr << std::endl << "but the method geomTypes() returned:" << std::endl;
      for (size_t j=0; j<lset.geomTypes(codim).size(); j++)
          std::cerr << "  " << lset.geomTypes(codim)[j] << std::endl;

      DUNE_THROW(GridError, "!");
  }

  //*****************************************************************
  // check size of index set
  int gridsize = 0;
  {
    typedef typename IndexSetType :: template Codim<codim>::
      template Partition<All_Partition> :: Iterator IteratorType;

    int count = 0;
    IteratorType endit  = lset.template end<codim,All_Partition>   ();
    for(IteratorType it = lset.template begin<codim,All_Partition> ();
        it != endit ; ++it )
    {
      count ++ ;
    }

    int lsetsize = lset.size(codim);
    if( count != lsetsize)
    {
      derr << "WARNING: walk = "<< count << " entities | set = "
        << lsetsize << " for codim " << codim << std::endl;
    }
    gridsize = count;
    // lsetsize should be at least the size of iterated entities
    assert( count <= gridsize );
  }

  {
    typedef typename GridType :: Traits :: LocalIdSet :: IdType IdType;
    std::map < IdType , bool > entityfound;
    int mycount = 0;
    Iterator endit  = lset.template end  <0,All_Partition> ();
    if (lset.template begin<0,All_Partition> () == endit)
      return;
    for(Iterator it = lset.template begin<0,All_Partition> ();
        it != endit ; ++it )
    {
      assert( lset.contains ( *it ) );
      int suben = it->template count<codim> ();
      for(int i=0 ;i<suben; i++)
      {
        IdType id = grid.localIdSet().id ( *(it->template entity<codim>(i) ) );
        if( entityfound.find(id) == entityfound.end())
        {
          mycount ++ ;
          entityfound[id] = true;
        }
      }
    }

    if ( gridsize != (int)entityfound.size() )
    {
      derr << "WARNING: gridsize = "<< gridsize << " entities | map of entities = "
        << entityfound.size() << " for codim " << codim << std::endl;
    }

    // gridsize should be at least the size of found entities
    //assert( gridsize <= (int) entityfound.size() );
  }

  //******************************************************************

  typedef std::pair < int , GeometryType > SubEntityKeyType;
  typedef std::map < int , std::pair<int,int> > subEntitymapType;
  std::map < SubEntityKeyType , std::vector<int> > subEntities;
  std::map < std::vector<int> , SubEntityKeyType > vertices;

  std::map < int , FieldVector<coordType,dim> > vertexCoordsMap;
  // setup vertex map , store vertex coords for vertex number
  {
    unsigned int count = 0;
    typedef typename IndexSetType :: template Codim<dim>::template Partition<All_Partition> :: Iterator VxIterator;
    VxIterator end = lset.template end <dim,All_Partition>();
    for(VxIterator it = lset.template begin <dim,All_Partition>();
        it != end; ++it )
    {
      count ++ ;
      // get coordinates of vertex
      FieldVector<coordType,dim> vx ( it->geometry().corner(0) );

      // get index of vertex
      sout << "Vertex " << vx << "\n";
      assert( lset.contains ( *it ) );
      int idx = lset.index( *it );

      sout << "Vertex " << idx << " = [" << vx << "]\n";

      // if vertex not in map insert it
      if( vertexCoordsMap.find(idx) == vertexCoordsMap.end())
        vertexCoordsMap[idx] = vx;
    }
    sout << "Found " << vertexCoordsMap.size() << " vertices for that index set!\n\n";

    // check whether size of map equals all found vertices
    assert( vertexCoordsMap.size() == count );

    // check whether size of vertices of set equals all found vertices
    sout << "Checking size of vertices "
     << count
     << " equals all found vertices "
     << (unsigned int)lset.size(Dune::GeometryType(0))
     << "\n";
    // assertion goes wrong for parallel grid since no iteration over ghost
    // subentities
    assert( count == (unsigned int)lset.size(Dune::GeometryType(0)) );
  }

  {
    // choose the right reference element
    Iterator refend = lset.template end  <0,All_Partition>();
    Iterator refit  = lset.template begin<0,All_Partition>();
    assert( refit != refend );

    GeometryType type = refit->geometry().type();

    const ReferenceElement< coordType, dim > & refElem =
      ReferenceElements< coordType, dim >::general(type);

    // print dune reference element
    sout << "Dune reference element provides: \n";
    for(int i=0; i<refElem.size(codim); i++)
    {
      sout << i << " subEntity [";
      int s = refElem.size(i,codim,dim);
      for(int j=0; j<s; j++)
      {
        sout << refElem.subEntity(i , codim , j , dim );
        if(j != s-1) sout << ",";
      }
      sout << "]\n";
    }
  }

  {
    Iterator endit  = lset.template end  <0,All_Partition>();
    for(Iterator it = lset.template begin<0,All_Partition>();
        it != endit; ++it)
    {
      // if (it->partitionType()==4) continue;
      sout << "****************************************\n";
      sout << "Element = " << lset.index(*it) << " on level " << it->level () << "\n";
      sout << "Vertices      = [";
      int svx = it->template count<dim>();

      // print all vertex numbers
      for(int i=0; i<svx; i++)
      {
        int idx = lset.template subIndex<dim> (*it,i);
        if(i == svx-1) sout << idx << "]\n";
        else sout << idx << ", ";
      }

      // print all vertex coordinates
      sout << "Vertex Coords = [";
      for(int i=0; i<svx; i++)
      {
        // get entity pointer of sub entity codim=dim (Vertex)
        typedef typename GridType :: template Codim<dim> :: EntityPointer VertexPointerType;
        VertexPointerType vxp = it->template entity<dim> (i);

        // get coordinates of entity pointer
        FieldVector<coordType,dim> vx (vxp->geometry().corner(0));

        // output vertex coordinates
        if(i<svx-1) sout << vx << " , ";
        else sout << vx << "]\n";

        int vxidx = lset.template subIndex<dim> (*it,i);
        int realidx = lset.index( *vxp );

        // the subIndex and the index for subEntity must be the same
        assert( vxidx == realidx );

        // check whether the coordinates are the same
          assert(vertexCoordsMap.find(vxidx)!=vertexCoordsMap.end());
        FieldVector<coordType,dim> vxcheck ( vertexCoordsMap[vxidx] );
        if( ! compareVec( vxcheck, vx ) )
        {
          sout << "ERROR: map global vertex " << vxidx << " vx " << vxcheck << " is not " << vx << " type:" << it->partitionType() << "\n";
          assert( compareVec( vxcheck, vx ) );
        }
      }

      ////////////////////////////////////////////////////////////
      // check sub entities
      ////////////////////////////////////////////////////////////
      checkSubEntity<codim> (grid, *it, lset, sout,
              subEntities, vertices, vertexCoordsMap);

      // check neighbors
      if(codim == 1)
      {
        std::string name = grid.name();
        if( name != "AlbertaGrid" )
        {
          if( levelIndex )
          {
            typedef typename EntityCodim0Type :: LevelIntersectionIterator IntersectionIterator;
            IntersectionIterator endnit  = it->ilevelend();
            for(IntersectionIterator nit = it->ilevelbegin(); nit != endnit; ++nit)
            {
              if(nit.neighbor())
              {
                typedef typename GridType :: template Codim<0> :: EntityPointer EnPointer;
                EnPointer ep = nit.outside();

                checkSubEntity<codim> (grid, *ep, lset, sout,
                        subEntities, vertices, vertexCoordsMap);
              }
            }
          }
        }
        else
        {
          static bool called = false;
          if( !called )
          {
            //std::cerr << "WARNING: skip indices test using LevelIntersectionIterator for AlbertaGrid!\n";
            called = true;
          }
        }

        if( !levelIndex )
        {
          typedef typename EntityCodim0Type :: LeafIntersectionIterator IntersectionIterator;
          IntersectionIterator endnit  = it->ileafend();
          for(IntersectionIterator nit = it->ileafbegin(); nit != endnit; ++nit)
          {
            if(nit.neighbor())
            {
              typedef typename GridType :: template Codim<0> :: EntityPointer EnPointer;
              EnPointer ep = nit.outside();

              checkSubEntity<codim> (grid, *ep, lset, sout,
                      subEntities, vertices, vertexCoordsMap);
            }
          }
        }
      }
    }
  }
}


template <class GridType, class IndexSetType, class OutputStreamImp,
          int codim, bool hasCodim>
struct CheckIndexSet
{
  static void checkIndexSet( const GridType &grid ,
        const IndexSetType & iset, OutputStreamImp & sout, bool levelIndex )
  {
    checkIndexSetForCodim<codim> (grid,iset,sout,levelIndex);
    CheckIndexSet<GridType,IndexSetType,OutputStreamImp,
      codim-1, Dune::Capabilities::hasEntity<GridType, codim-1>::v > ::
      checkIndexSet( grid, iset, sout,levelIndex );
  }
};

template <class GridType, class IndexSetType, class OutputStreamImp,
          int codim>
struct CheckIndexSet<GridType,IndexSetType,OutputStreamImp,codim,false>
{
  static void checkIndexSet( const GridType &grid ,
        const IndexSetType & iset, OutputStreamImp & sout , bool levelIndex )
  {
    //derr << "WARNING: entities for codim " << codim << " are not being tested!" << std::endl;
    CheckIndexSet<GridType,IndexSetType,OutputStreamImp,
      codim-1, Dune::Capabilities::hasEntity<GridType, codim-1>::v > ::
      checkIndexSet( grid, iset, sout , levelIndex );
  }
};

// end loop over codim by specialisation
template <class GridType, class IndexSetType, class OutputStreamImp>
struct CheckIndexSet<GridType,IndexSetType,OutputStreamImp,0,true>
{
  static void checkIndexSet( const GridType &grid ,
        const IndexSetType & iset, OutputStreamImp & sout , bool levelIndex )
  {
    checkIndexSetForCodim<0> (grid,iset,sout,levelIndex);
  }
};

template <class GridType, class IndexSetType, class OutputStreamImp>
void checkIndexSet( const GridType &grid , const IndexSetType & iset,
    OutputStreamImp & sout ,  bool levelIndex = false )
{
  CheckIndexSet<GridType,IndexSetType,OutputStreamImp,
    GridType::dimension, true> ::
    checkIndexSet (grid,iset,sout,levelIndex);
}

} // end namespace Dune
#endif
