#ifndef DUNE_CHECK_GEOMETRYINFATHER_CC
#define DUNE_CHECK_GEOMETRYINFATHER_CC

#include "config.h"
#include <dune/common/typetraits.hh>

/** \file
    \brief A test for the Method Geometry::geometryInFather()
*/

/** \brief Test the Method Geometry::geometryInFather()

    This test works by comparing the output of geometryInFather with the vertex positions
    obtained by directly expressing the son vertex positions in local coordinates of the
    father.  That should work for all grid implementations that realize truly nested
    grids.  One exception is UGGrid with parametrized boundaries.
*/
template <class GridType>
void checkGeometryInFather(const GridType& grid) {

    using namespace Dune;

    typedef typename GridType::ctype ctype;
    const int dim      = GridType::dimension;

    // We need at least two levels to do any checking
    if (grid.maxLevel()==0)
        {
            dwarn << "SKIPPING check geometryInFather(), because grid has only one level! \n";
            return;
        }

    // Loop over all levels except the lowest one
    for (int i=1; i<=grid.maxLevel(); i++) {

        typedef typename GridType::template Codim<0>::LevelIterator ElementIterator;
        ElementIterator eIt    = grid.template lbegin<0>(i);
        ElementIterator eEndIt = grid.template lend<0>(i);

        for (; eIt!=eEndIt; ++eIt)
            {
                // check the father method
                if( eIt->level () > 0 )
                    {
                        typedef typename GridType::template Codim<0>::EntityPointer EntityPointer;
                        EntityPointer father = eIt->father();

                        while ( father->level() > 0 )
                            {
                                EntityPointer grandPa = father->father();
                                typedef typename GridType :: Traits :: HierarchicIterator HierarchicIterator;

                                const int mxl = grandPa->level() + 1;
                                bool foundChild = false;

                                HierarchicIterator end = grandPa->hend(mxl);
                                for(HierarchicIterator sons = grandPa->hbegin(mxl);
                                    sons != end; ++sons)
                                    {
                                        if(father == sons) foundChild = true;
                                    }

                                if(!foundChild)
                                    DUNE_THROW(GridError, "father-child error while iterating over childs of father!");
                                father = grandPa;
                            }
                    }

                // hierarchy check
                {
                    typedef typename GridType :: Traits :: HierarchicIterator HierarchicIterator;
                    const int mxl = grid.maxLevel();
                    int countChildren = 0;

                    HierarchicIterator end = eIt->hend(mxl);
                    for(HierarchicIterator sons = eIt->hbegin(mxl);
                        sons != end; ++sons)
                        {
                            ++countChildren;
                            int count = sons->level();
                            if( sons->level() > 0 )
                                {
                                    typedef typename GridType::template Codim<0>::EntityPointer EntityPointer;
                                    EntityPointer father = sons->father();
                                    --count;
                                    while ( father->level() > 0 )
                                        {
                                            father = father->father();
                                            --count;
                                        }
                                }
                            assert( count == 0 );
                        }

                    if( eIt->isLeaf () && countChildren > 0 )
                        DUNE_THROW(GridError, "leaf entity has children ==> entity is not leaf");
                }

                // Get geometry in father
                typedef typename GridType::template Codim<0>::Entity::Geometry Geometry;
                typedef typename GridType::template Codim<0>::Entity::LocalGeometry LocalGeometry;

                const LocalGeometry& geometryInFather = eIt->geometryInFather();

                // //////////////////////////////////////////////////////
                //   Check for types and constants
                // //////////////////////////////////////////////////////

                //            IsTrue< is_same<
                //              typename Geometry::ctype,
                //              typename GridType::ctype>::value == true >::yes();


                // ///////////////////////////////////////////////////////
                //   Check the different methods
                // ///////////////////////////////////////////////////////
                if (geometryInFather.type() != eIt->geometry().type())
                    DUNE_THROW(GridError, "Type of geometry and geometryInFather differ!");

                if (geometryInFather.corners() != eIt->geometry().corners())
                    DUNE_THROW(GridError, "entity and geometryInFather have different number of corners!");

                // Compute the element center just to have an argument for the following methods
                FieldVector<ctype, dim> center(0);
                for (int j=0; j<geometryInFather.corners(); j++)
                    center += geometryInFather[j];

                if (geometryInFather.integrationElement(center) <=0)
                    DUNE_THROW(GridError, "nonpositive integration element found!");

                /** \todo Missing local() */
                /** \todo Missing global() */
                /** \todo Missing jacobianInverse() */
                /** \todo Missing checkInside() */

                // /////////////////////////////////////////////////////////////////////////////////////
                // Check whether the positions of the vertices of geometryInFather coincide
                // with the ones computed 'by hand'.  This only works if the grids really are nested!
                // /////////////////////////////////////////////////////////////////////////////////////
                for (int j=0; j<geometryInFather.corners(); j++) {

                    FieldVector<ctype, dim> localPos = eIt->father()->geometry().local(eIt->geometry().corner(j));

                    if ( (localPos-geometryInFather[j]).infinity_norm() > 1e-7)
                        {
                            std::cerr << localPos << " lp[" << j << "] | gp " << geometryInFather[j] << "\n";
                            DUNE_THROW(GridError, "geometryInFather yields wrong vertex position!");
                        }

                }

            }

    }

}

#endif
