#ifndef DUNE_CHECK_INTERSECTIONITERATOR_CC
#define DUNE_CHECK_INTERSECTIONITERATOR_CC

#include <cmath>

#include <dune/grid/common/quadraturerules.hh>
#include <dune/grid/common/gridpart.hh>

/** \file
    \brief Tests for the IntersectionIterator
*/

/** \brief Helper routine: Test a general geometry
 */

template <class GeometryImp>
void checkGeometry(const GeometryImp& geometry)
{
    using namespace Dune;

    // Get the dimensions and types
    const int dim      = GeometryImp::mydimension;
    const int dimworld = GeometryImp::coorddimension;

    typedef typename GeometryImp::ctype ctype;

    // Get the corresponding reference element
    const ReferenceElement<double,dim>& refElement
        = ReferenceElements<double, dim>::general(geometry.type());

    // Check whether the number of corners is correct
    if (geometry.corners() != refElement.size(dim))
        DUNE_THROW(GridError, "Geometry has wrong number of corners!");

    // check consistency between operator[] and global()
    for (int i=0; i<refElement.size(dim); i++) {
        FieldVector<double,dim> localPos = refElement.position(i,dim);
        if ( (geometry[i] - geometry.global(localPos)).infinity_norm() > 1e-6)
            DUNE_THROW(GridError, "Methods operator[] and global() are inconsistent!");
    }

    // Use a quadrature rule to create a few test points for the following checks
    const QuadratureRule<double, dim>& quad 
        = QuadratureRules<double, dim>::rule(geometry.type(), 2);

    for (size_t i=0; i<quad.size(); i++) {

        const FieldVector<double,dim>& testPoint = quad[i].position();

        // Check whether point is within the intersection
        if (!geometry.checkInside(testPoint))
            DUNE_THROW(GridError, "Test point is not within geometry!");
    
        // Transform to global coordinates
        FieldVector<ctype, dimworld> global = geometry.global(testPoint);

        // The back to local coordinates
        FieldVector<ctype, dim> local = geometry.local(global);
    
        // check for correctness
        if ((testPoint-local).infinity_norm() > 1e-6)
            DUNE_THROW(GridError, "local() and global() are not inverse to each other!");
    
        // The integration element at the element center
        ctype intElement = geometry.integrationElement(testPoint);
        if (intElement <=0)
            DUNE_THROW(GridError, "nonpositive integration element found!");
    
#if 0
    // This method exists in the interface, but it is not expected to work
    // unless dim==dimworld
    const FieldMatrix<ctype, dim, dim> jacobi
        = intersectionGlobal.jacobianInverseTransposed(testPoint);
#endif
    }
}

/** \brief Test the IntersectionIterator
*/
template <class GridPartType>
void checkIntersectionIterator(const GridPartType& gridPart,
                               const typename GridPartType::Traits::template Codim<0>::IteratorType& eIt)
{
  using namespace Dune;

  typedef typename GridPartType::GridType GridType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIterator;

  const GridType& grid = gridPart.grid();
  const bool checkOutside = (grid.name() != "AlbertaGrid");
  const typename GridPartType::IndexSetType& indexSet = gridPart.indexSet();

  const int dimworld = GridType::dimensionworld;

  FieldVector<double,dimworld> sumNormal(0.0); 

  // /////////////////////////////////////////////////////////
  //   Check the types defined by the iterator
  // /////////////////////////////////////////////////////////
//  IsTrue< is_same<
//      typename IntersectionIterator::ctype,
//      typename GridType::ctype>::value == true >::yes();

  IsTrue<static_cast<int>(IntersectionIterator::dimension)
    == static_cast<int>(GridType::dimension)>::yes();

  IsTrue<static_cast<int>(IntersectionIterator::dimensionworld)
    == static_cast<int>(GridType::dimensionworld)>::yes();

  IntersectionIterator iIt    = gridPart.ibegin(*eIt);
  IntersectionIterator iEndIt = gridPart.iend(*eIt);
  
  for (;iIt!=iEndIt; ++iIt) 
  { 
      // //////////////////////////////////////////////////////////////////////
      //   Compute the integral of the outer normal over the whole element.
      //   This has to be zero.
      // //////////////////////////////////////////////////////////////////////
      const int interDim = IntersectionIterator::LocalGeometry::mydimension;
      const QuadratureRule<double, interDim>& quad 
          = QuadratureRules<double, interDim>::rule(iIt.intersectionSelfLocal().type(), interDim);

      for (size_t i=0; i<quad.size(); i++)
          sumNormal.axpy(quad[i].weight(), iIt.integrationOuterNormal(quad[i].position()));
      
      typedef typename IntersectionIterator::Entity EntityType; 
      typedef typename EntityType::EntityPointer EntityPointer;

      assert(eIt == iIt.inside());

      // check that boundary id has positive value 
      if( iIt.boundary() )
      {
        if( iIt.boundaryId() <= 0 )
        {
          DUNE_THROW(GridError, "boundary id has non-positive value (" << iIt.boundaryId() << ") !");
        }
      }

      // //////////////////////////////////////////////////////////////////////
      //   Check whether the 'has-intersection-with'-relation is symmetric
      // //////////////////////////////////////////////////////////////////////

      if (iIt.neighbor() && checkOutside ) 
      {
          EntityPointer outside = iIt.outside();
          bool insideFound = false;

          IntersectionIterator outsideIIt    = gridPart.ibegin(*outside);
          IntersectionIterator outsideIEndIt = gridPart.iend(*outside);

          for (; outsideIIt!=outsideIEndIt; ++outsideIIt) {

              if (outsideIIt.neighbor() && outsideIIt.outside() == iIt.inside()) {

                  if (outsideIIt.numberInSelf() != iIt.numberInNeighbor())
                      DUNE_THROW(GridError, "outside()->outside() == inside(), but with incorrect numbering!");
                  else
                      insideFound = true;

              }

          }

          if (!insideFound)
              DUNE_THROW(GridError, "Could not find inside() through intersection iterator of outside()!");

      }
      else if (!checkOutside) 
      {
        static bool called = false;
        if(!called)
        {
          derr << "WARNING: skip reverse intersection iterator test for " << grid.name() << "!"<< std::endl;
          called = true;
        }
      }

      // /////////////////////////////////////////////////////////////
      //   Check the consistency of numberInSelf, numberInNeighbor
      //   and the indices of the subface between.
      // /////////////////////////////////////////////////////////////
      if ( GridPartType::conforming && iIt.neighbor() ) 
      {
          EntityPointer outside = iIt.outside();
          int numberInSelf     = iIt.numberInSelf();
          int numberInNeighbor = iIt.numberInNeighbor();

          assert(indexSet.template subIndex<1>(*eIt, numberInSelf)
                 == indexSet.template subIndex<1>(*outside, numberInNeighbor));

          assert(grid.localIdSet().template subId<1>(*eIt, numberInSelf)
                 == grid.localIdSet().template subId<1>(*outside, numberInNeighbor));
          
          assert(grid.globalIdSet().template subId<1>(*eIt, numberInSelf)
                 == grid.globalIdSet().template subId<1>(*outside, numberInNeighbor));

      }

      // //////////////////////////////////////////////////////////
      //   Check the geometry returned by intersectionGlobal()
      // //////////////////////////////////////////////////////////
      typedef typename IntersectionIterator::Geometry Geometry;
      const Geometry& intersectionGlobal = iIt.intersectionGlobal();

      checkGeometry(intersectionGlobal);

      // //////////////////////////////////////////////////////////
      //   Check the geometry returned by intersectionSelfLocal()
      // //////////////////////////////////////////////////////////

      const typename IntersectionIterator::LocalGeometry& intersectionSelfLocal = iIt.intersectionSelfLocal();
      checkGeometry(intersectionSelfLocal);

      //  Check the consistency of intersectionSelfLocal() and intersectionGlobal
      
      if (intersectionSelfLocal.corners() != intersectionGlobal.corners())
          DUNE_THROW(GridError, "Geometry of intersection is inconsistent from left hand side and global view!");
      
      // (Ab)use a quadrature rule as a set of test points
      for (size_t i=0; i<quad.size(); i++) 
      {
          // check integrationOuterNormal 
          double det = intersectionGlobal.integrationElement(quad[i].position());
          det -= iIt.integrationOuterNormal(quad[i].position()).two_norm();
          if( std::abs( det ) > 1e-8 )
          {
            DUNE_THROW(GridError, "integrationElement and length of integrationOuterNormal do no match!");
          }
          
          FieldVector<double,dimworld> globalPos = intersectionGlobal.global(quad[i].position());
          FieldVector<double,dimworld> localPos  = eIt->geometry().global(intersectionSelfLocal.global(quad[i].position()));

          if ( (globalPos - localPos).infinity_norm() > 1e-6)
              DUNE_THROW(GridError, "global( intersectionSelfLocal(global() ) is not the same as intersectionGlobal.global() at " << quad[i].position() << "!");

      }

      // ////////////////////////////////////////////////////////////////
      //   Check the geometry returned by intersectionNeighborLocal()
      // ////////////////////////////////////////////////////////////////

      if (iIt.neighbor() ) 
      {

          const typename IntersectionIterator::LocalGeometry& intersectionNeighborLocal = iIt.intersectionNeighborLocal();
          
          checkGeometry(intersectionNeighborLocal);

          if (intersectionSelfLocal.corners() != intersectionNeighborLocal.corners())
              DUNE_THROW(GridError, "Geometry of intersection is inconsistent from left and right hand side!");

          // (Ab)use a quadrature rule as a set of test points
          const int interDim = IntersectionIterator::LocalGeometry::mydimension;
          const QuadratureRule<double, interDim>& quad 
              = QuadratureRules<double, interDim>::rule(intersectionNeighborLocal.type(), 2);

          for (size_t i=0; i<quad.size(); i++) 
          {

              FieldVector<double,dimworld> globalPos = intersectionGlobal.global(quad[i].position());
              FieldVector<double,dimworld> localPos  = iIt.outside()->geometry().global(intersectionNeighborLocal.global(quad[i].position()));

              if ( (globalPos - localPos).infinity_norm() > 1e-6)
                  DUNE_THROW(GridError, "global( intersectionNeighborLocal(global() ) is not the same as intersectionGlobal.global() at " << quad[i].position() << "!");
              
          }

      }

  }

  // ////////////////////////////////////////////////////////////////////////
  //   Check whether the integral over the outer normal really is zero
  // ////////////////////////////////////////////////////////////////////////
  if( sumNormal.two_norm() > 1e-8 ) 
        DUNE_THROW(GridError,"Sum of integrationOuterNormals is " << sumNormal.two_norm() << " but it should be Zero!");

}

/** \brief Test both IntersectionIterators 
 */
template <class GridType>
void checkIntersectionIterator(const GridType& grid, bool skipLevelIntersectionTest = false) {

    using namespace Dune;

    // Loop over all levels
    if(skipLevelIntersectionTest) 
    {
      std::cerr<<"WARNING: skip test of LevelIntersectionIterator! \n";
    }
    else   
    {
      for (int i=0; i<=grid.maxLevel(); i++) 
      {
              
          typedef typename GridType::template Codim<0>::LevelIterator ElementIterator;

          LevelGridPart<const GridType, All_Partition> levelGridPart(grid, i);

          ElementIterator eIt    = grid.template lbegin<0>(i);
          ElementIterator eEndIt = grid.template lend<0>(i);

          for (; eIt!=eEndIt; ++eIt) 
            checkIntersectionIterator(levelGridPart, eIt);

      }
    }

    // test leaf intersection iterator 
    {
        typedef typename GridType::template Codim<0>::LeafIterator ElementIterator;

        LeafGridPart<const GridType, All_Partition> leafGridPart(grid);

        ElementIterator eEndIt = grid.template leafend<0>();
        for (ElementIterator eIt = grid.template leafbegin<0>(); eIt!=eEndIt; ++eIt) 
            checkIntersectionIterator(leafGridPart, eIt);

    }

}

#endif
