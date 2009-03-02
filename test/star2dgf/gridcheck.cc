// $Id: gridcheck.cc 3548 2007-02-20 12:53:12Z christi $

/**

   Implements a generic grid check


   \todo check return types

*/

#include "config.h"
#include <dune/grid/common/capabilities.hh>
#include <dune/common/helpertemplates.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/stdstreams.hh>
#include <dune/grid/common/referenceelements.hh>

#include "checkindexset.cc"

#include <limits>

// machine epsilon is multiplied by this factor
static double factorEpsilon = 1.e8;

class CheckError : public Dune::Exception {};

// --- compile-time check of element-interface

template <class Geometry, bool doCheck>
struct JacobianInverse
{
    static void check(const Geometry &e)
    {
        typedef typename Geometry::ctype ctype;
        Dune::FieldVector<ctype, Geometry::mydimension> v;
        e.jacobianInverseTransposed(v);
    }
    JacobianInverse()
    {
        c = check;
    };
    void (*c)(const Geometry&);
};

template <class Geometry>
struct JacobianInverse<Geometry, false>
{
    static void check(const Geometry &e)
    {
    }
    JacobianInverse()
    {
        c = check;
    };
    void (*c)(const Geometry&);
};

template <class Geometry, int codim, int dim>
struct GeometryInterface
{
    static void check(const Geometry &e)
    {
        IsTrue<dim-codim == Geometry::mydimension>::yes();
        IsTrue<dim == Geometry::dimension>::yes();

        typedef typename Geometry::ctype ctype;

        e.type();
        e.corners();
        e[0];

        Dune::FieldVector<ctype, Geometry::mydimension> v;
        e.global(v);
        Dune::FieldVector<ctype, Geometry::coorddimension> g;
        e.local(g);
        e.checkInside(v);
        e.integrationElement(v);
        JacobianInverse<Geometry,
            (int)Geometry::coorddimension == (int)Geometry::mydimension>();
    }
    GeometryInterface()
    {
        c = check;
    };
    void (*c)(const Geometry&);
};

// reduced test on vertices
template <class Geometry, int dim>
struct GeometryInterface <Geometry, dim, dim>
{
    static void check(const Geometry &e)
    {
        IsTrue<0 == Geometry::mydimension>::yes();
        IsTrue<dim == Geometry::dimension>::yes();

        // vertices have only a subset of functionality
        e.type();
        e.corners();
        e[0];
    }
    GeometryInterface()
    {
        c = check;
    };
    void (*c)(const Geometry&);
};

// --- compile-time check of entity-interface

// tests that should work on entities of all codimensions
template <class Entity>
void DoEntityInterfaceCheck (Entity &e)
{
    // exported types
    typedef typename Entity::ctype ctype;

    // methods on each entity
    e.level();
    e.partitionType();
    e.geometry();

    // check interface of attached element-interface
    GeometryInterface<typename Entity::Geometry, Entity::codimension, Entity::dimension>();
}

// recursive check of codim-0-entity methods count(), entity()
template <class Grid, int cd, bool hasEntity>
struct ZeroEntityMethodCheck
{
    typedef typename Grid::template Codim<0>::Entity Entity;
    static void check(Entity &e)
    {
        // check types
        typedef typename Entity::IntersectionIterator IntersectionIterator;
        typedef typename Entity::HierarchicIterator HierarchicIterator;
        typedef typename Entity::EntityPointer EntityPointer;

        e.template count<cd>();
        e.template entity<cd>(0);

        // recursively check on
        ZeroEntityMethodCheck<Grid, cd - 1,
            Dune::Capabilities::hasEntity<Grid, cd - 1>::v >();
    }
    ZeroEntityMethodCheck ()
    {
        c = check;
    }
    void (*c)(Entity &e);
};

// just the recursion if the grid does not know about this codim-entity
template<class Grid, int cd>
struct ZeroEntityMethodCheck<Grid, cd, false>
{
    typedef typename Grid::template Codim<0>::Entity Entity;
    static void check(Entity &e)
    {
        // check types
        typedef typename Entity::IntersectionIterator IntersectionIterator;
        typedef typename Entity::HierarchicIterator HierarchicIterator;
        typedef typename Entity::EntityPointer EntityPointer;

        // recursively check on
        ZeroEntityMethodCheck<Grid, cd - 1,
            Dune::Capabilities::hasEntity<Grid, cd - 1>::v >();
    }
    ZeroEntityMethodCheck ()
    {
        c = check;
    }
    void (*c)(Entity &e);
};

// end recursive checking
template <class Grid>
struct ZeroEntityMethodCheck<Grid, 0, true>
{
    typedef typename Grid::template Codim<0>::Entity Entity;
    static void check(Entity &e)
    {
        // check types
        typedef typename Entity::IntersectionIterator IntersectionIterator;
        typedef typename Entity::HierarchicIterator HierarchicIterator;
        typedef typename Entity::EntityPointer EntityPointer;

        e.template count<0>();
        e.template entity<0>(0);

    }
    ZeroEntityMethodCheck ()
    {
        c = check;
    }
    void (*c)(Entity &e);
};

// end recursive checking - same as true
// ... codim 0 is always needed
template <class Grid>
struct ZeroEntityMethodCheck<Grid, 0, false>
{
    typedef typename Grid::template Codim<0>::Entity Entity;
    static void check(Entity &e)
    {
        // check types
        typedef typename Entity::IntersectionIterator IntersectionIterator;
        typedef typename Entity::HierarchicIterator HierarchicIterator;
        typedef typename Entity::EntityPointer EntityPointer;

        e.template count<0>();
        e.template entity<0>(0);
    }
    ZeroEntityMethodCheck ()
    {
        c = check;
    }
    void (*c)(Entity &e);
};

// IntersectionIterator interface check
template <class Grid>
struct IntersectionIteratorInterface
{
    typedef typename Grid::template Codim<0>::IntersectionIterator IntersectionIterator;
    enum { dim = Grid::dimension };
    typedef typename Grid::ctype ct;

    static void check (IntersectionIterator &i)
    {
        // increment / equality / ...
        IntersectionIterator j = i;
        j++;
        i == j;
        i != j;
        j = i;

        // state
        i.boundary();
        i.neighbor();

        // id of boundary segment
        i.boundaryId();

        // neighbouring elements
        i.inside();
        if(i.neighbor()) i.outside();

        // geometry
        i.intersectionSelfLocal();
        if(i.neighbor()) i.intersectionNeighborLocal();
        i.intersectionGlobal();

        i.numberInSelf();
        if(i.neighbor()) i.numberInNeighbor();

        Dune::FieldVector<ct, dim-1> v(0);
        i.outerNormal(v);
        i.integrationOuterNormal(v);
        i.unitOuterNormal(v);
    }
    IntersectionIteratorInterface ()
    {
        c = check;
    }
    void (*c)(IntersectionIterator&);
};

// check codim-entity and pass on to codim + 1
template <class Grid, int codim, int dim, bool hasEntity>
struct EntityInterface
{
    typedef typename Grid::template Codim<codim>::Entity Entity;

    static void check (Entity &e)
    {
        // consistent?
        IsTrue<codim == Entity::codimension>::yes();
        IsTrue<dim == Entity::dimension>::yes();

        // do the checking
        DoEntityInterfaceCheck(e);

        // recursively check sub-entities
        EntityInterface<Grid, codim + 1, dim,
            Dune::Capabilities::hasEntity<Grid, codim + 1>::v >();
    }
    EntityInterface ()
    {
        c = check;
    }
    void (*c)(Entity&);
};

// just the recursion if the grid does not know about this codim-entity
template <class Grid, int codim, int dim>
struct EntityInterface<Grid, codim, dim, false>
{
    typedef typename Grid::template Codim<codim>::Entity Entity;

    static void check (Entity &e)
    {
        // recursively check sub-entities
        EntityInterface<Grid, codim + 1, dim,
            Dune::Capabilities::hasEntity<Grid, codim + 1>::v >();
    }
    EntityInterface ()
    {
        c = check;
    }
    void (*c)(Entity&);
};

// codim-0 entities have different interface
template <class Grid, int dim>
struct EntityInterface<Grid, 0, dim, true>
{
    typedef typename Grid::template Codim<0>::Entity Entity;

    static void check (Entity &e,bool checkLevelIter=true)
    {
        // consistent?
        IsTrue<0 == Entity::codimension>::yes();
        IsTrue<dim == Entity::dimension>::yes();

        // do the common checking
        DoEntityInterfaceCheck(e);

        // special codim-0-entity methods which are parametrized by a codimension
        ZeroEntityMethodCheck
            <Grid, dim, Dune::Capabilities::hasEntity<Grid, dim>::v >();

        // grid hierarchy
        e.father();
        e.geometryInFather();


        // intersection iterator
        if (checkLevelIter) {
            e.ilevelbegin();
            e.ilevelend();
            IntersectionIteratorInterface<Grid>(e.ilevelbegin());
        }
        e.ileafbegin();
        e.ileafend();

        if(e.isLeaf())
            IntersectionIteratorInterface<Grid>(e.ileafbegin());

        // hierarchic iterator
        e.hbegin(0);
        e.hend(0);

        // adaption
        e.state();

        // recursively check sub-entities
        EntityInterface<Grid, 1, dim,
            Dune::Capabilities::hasEntity<Grid, 1>::v >();
    }
    EntityInterface ()
    {
        c = check;
    }
    void (*c)(Entity&);
};

// non existinng codim-0 entity
template <class Grid, int dim>
struct EntityInterface<Grid, 0, dim, false>
{
    typedef typename Grid::template Codim<0>::Entity Entity;

    static void check (Entity &e)
    {
        // recursively check sub-entities
        EntityInterface<Grid, 1, dim,
            Dune::Capabilities::hasEntity<Grid, 1>::v >();
    }
    EntityInterface ()
    {
        c = check;
    }
    void (*c)(Entity&);
};

// end the recursion over entity-codimensions
template <class Grid, int dim>
struct EntityInterface<Grid, dim, dim, true>
{
    typedef typename Grid::template Codim<dim>::Entity Entity;

    // end recursion
    static void check (Entity &e)
    {
        // consistent?
        IsTrue<dim == Entity::codimension>::yes();
        IsTrue<dim == Entity::dimension>::yes();

        // run common test
        DoEntityInterfaceCheck(e);
    }

    EntityInterface()
    {
        c = check;
    }
    void (*c)(Entity&);
};

// end the recursion over entity-codimensions
// ... codim dim entity does not exist
template <class Grid, int dim>
struct EntityInterface<Grid, dim, dim, false>
{
    typedef typename Grid::template Codim<dim>::Entity Entity;

    // end recursion
    static void check (Entity &e)
    {
    }

    EntityInterface()
    {
        c = check;
    }
    void (*c)(Entity&);
};

template<class Grid>
struct LeafInterface
{
    static void check(Grid &g)
    {
        g.template leafbegin<0>();
        g.template leafend<0>();
    }
    LeafInterface()
    {
        c = check;
    }
    void (*c)(Grid&);
};

template <class Grid>
struct GridInterface
{
    static void check (const Grid &g)
    {
        // check for exported types
        typedef typename Grid::ctype ctype;
        typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
        typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;
        typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;

        // check for member functions
        g.maxLevel();
        // number of grid entities of a given codim on a given level
        g.size(0,0);
        // number of leaf entities per codim in this process
        g.size(0);
        // number of entities per level and geometry type in this process
        g.size(0, Dune::GeometryType(Dune::GeometryType::cube,Grid::dimension));
        // number of leaf entities per geometry type in this process
        g.size(Dune::GeometryType(Dune::GeometryType::cube,Grid::dimension));

        // Check overlap and ghost size on level 0
        g.overlapSize(0,0);
        g.ghostSize(0,0);

        // Check overlap and ghost size on the leaf level
        g.overlapSize(0);
        g.ghostSize(0);

        // check for iterator functions
        g.template lbegin<0>(0);
        g.template lend<0>(0);

        LeafInterface< Grid>();

        // Check for index sets
        typedef typename Grid::template Codim<0>::LevelIndexSet LevelIndexSet;
        typedef typename Grid::template Codim<0>::LeafIndexSet LeafIndexSet;
        typedef typename Grid::template Codim<0>::LocalIdSet LocalIdSet;
        typedef typename Grid::template Codim<0>::GlobalIdSet GlobalIdSet;

        g.levelIndexSet(0);
        if (g.template lbegin<0>(0) !=g.template lend<0>(0) ) {
            // Instantiate all methods of LevelIndexSet
            g.levelIndexSet(0).index(*g.template lbegin<0>(0));
            /** \todo Test for subindex is missing, because I don't know yet
                how to test for the existence of certain codims */
        }
        g.levelIndexSet(0).
            size(Dune::GeometryType(Dune::GeometryType::simplex,Grid::dimension));
        for (int codim = 0; codim < Grid::dimension; codim++)
            g.levelIndexSet(0).geomTypes(codim);

        if (g.template lbegin<0>(0) !=g.template lend<0>(0) ) {
            // Instantiate all methods of LeafIndexSet
            g.leafIndexSet().index(*g.template leafbegin<0>());
        }
        /** \todo Test for subindex is missing, because I don't know yet
            how to test for the existence of certain codims */
        g.leafIndexSet().size(Dune::GeometryType(Dune::GeometryType::simplex,Grid::dimension));
        for (int codim = 0; codim < Grid::dimension; codim++)
            g.leafIndexSet().geomTypes(codim);

        if (g.template lbegin<0>(0) !=g.template lend<0>(0) ) {
            // Instantiate all methods of LocalIdSet
            /** \todo Test for subindex is missing, because I don't know yet
                how to test for the existence of certain codims */
            g.localIdSet().id(*g.template lbegin<0>(0));
            // Instantiate all methods of GlobalIdSet
            /** \todo Test for subindex is missing, because I don't know yet
                how to test for the existence of certain codims */
            g.globalIdSet().id(*g.template lbegin<0>(0));
        }
        // recursively check entity-interface
        // ... we only allow grids with codim 0 zero entites
        IsTrue<Dune::Capabilities::hasEntity<Grid, 0>::v>::yes();
        IsTrue<Dune::Capabilities::hasEntity<const Grid, 0>::v>::yes();

        //EntityInterface<Grid, 0, Grid::dimension,
        //  Dune::Capabilities::hasEntity<Grid, 0>::v >();

        // !!! check for parallel grid?
        /*
          g.template lbegin<0, Dune::Ghost_Partition>(0);
          g.template lend<0, Dune::Ghost_Partition>(0);
        */
    }
    GridInterface()
    {
        c = check;
    }
    // member just to avoid "unused variable"-warning in constructor
    void (*c)(const Grid&);
};

// check
// Entity::geometry().corner(c] == Entity::entity<dim>.geometry()[0)
// for codim=cd
template <int cd, class Grid, class Entity, bool doCheck>
struct subIndexCheck
{
    subIndexCheck (const Grid & g, const Entity & e)
    {
        const int imax = e.template count<cd>();
        for (int i=0; i<imax; ++i)
            {
                if( g.levelIndexSet(e.level()).index( *(e.template entity<cd>(i)) )
                    != g.levelIndexSet(e.level()).template subIndex<cd>(e,i) )
                    {
                        int id_e =
                            g.levelIndexSet(e.level()).index(e);
                        int id_e_i =
                            g.levelIndexSet(e.level()).index( *(e.template entity<cd>(i)) );
                        int subid_e_i =
                            g.levelIndexSet(e.level()).template subIndex<cd>(e,i);
                        DUNE_THROW(CheckError,
                                   "g.levelIndexSet.index( *(e.template entity<cd>(i)) ) "
                                   << "== g.levelIndexSet.template subIndex<cd>(e,i) failed "
                                   << "[with cd=" << cd << ", i=" << i << "]"
                                   << " ... index(e)=" << id_e
                                   << " ... index(e.entity<cd>(i))=" << id_e_i
                                   << " ... subIndex(e,i)=" << subid_e_i
                                   );
                    }
            }
        subIndexCheck<cd-1,Grid,Entity,
            Dune::Capabilities::hasEntity<Grid,cd-1>::v> sick(g,e);
    }
};
// end recursion of subIndexCheck
template <class Grid, class Entity, bool doCheck>
struct subIndexCheck<-1, Grid, Entity, doCheck>
{
    subIndexCheck (const Grid & g, const Entity & e)
    {
        return;
    }
};
// do nothing if doCheck==false
template <int cd, class Grid, class Entity>
struct subIndexCheck<cd, Grid, Entity, false>
{
    subIndexCheck (const Grid & g, const Entity & e)
    {
        subIndexCheck<cd-1,Grid,Entity,
            Dune::Capabilities::hasEntity<Grid,cd-1>::v> sick(g,e);
    }
};
template <class Grid, class Entity>
struct subIndexCheck<-1, Grid, Entity, false>
{
    subIndexCheck (const Grid & g, const Entity & e)
    {
        return;
    }
};

// name says all
template <class Grid>
void zeroEntityConsistency (Grid &g)
{
    typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
    typedef typename Grid::template Codim<0>::Geometry Geometry;
    typedef typename Grid::template Codim<0>::Entity Entity;
    LevelIterator it = g.template lbegin<0>(g.maxLevel());
    const LevelIterator endit = g.template lend<0>(g.maxLevel());

    for (; it!=endit; ++it)
        {
            // Entity::entity<0>(0) == Entity
            assert( g.levelIndexSet(g.maxLevel()).index( *(it->template entity<0>(0)) )
                    == g.levelIndexSet(g.maxLevel()).index( *it ) );
            assert( g.leafIndexSet().index( *(it->template entity<0>(0)) )
                    == g.leafIndexSet().index( *it ) );

            assert( g.globalIdSet().id( *(it->template entity<0>(0)) )
                    == g.globalIdSet().id( *it ) );

            assert( g.localIdSet().id( *(it->template entity<0>(0)) )
                    == g.localIdSet().id( *it ) );
            assert( it->template entity<0>(0)->level() == it->level() );
            // Entity::count<dim>() == Entity::geometry().corners();
            assert( it->template count<Grid::dimension>() == it->geometry().corners() );
            // Entity::geometry().corner(c] == Entity::entity<dim>.geometry()[0);
            const int cmax = it->template count<Grid::dimension>();
            for (int c=0; c<cmax; ++c)
                {
                    Dune::FieldVector<typename Grid::ctype, Grid::dimensionworld> c1(it->geometry().corner(c));
                    Dune::FieldVector<typename Grid::ctype, Grid::dimensionworld> c2(it->template entity<Grid::dimension>(c)->geometry().corner(0));
                    if( (c2-c1).two_norm() > 10 * std::numeric_limits<typename Grid::ctype>::epsilon() )
                        {
                            DUNE_THROW(CheckError, "geometry[i] == entity<dim>(i) failed: || c1-c2 || = || " <<
                                       c1 << " - " << c2 << " || = " << (c2-c1).two_norm() << " [ with i = " << c << " ]");
                        }
                }
            subIndexCheck<Grid::dimension, Grid, Entity, true> sick(g,*it);
        }
}

/*
 * search the LevelIterator for each IntersectionIterator
 */
template <class Grid>
void assertNeighbor (Grid &g)
{
    typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
    typedef typename Grid::template Codim<0>::LevelIntersectionIterator IntersectionIterator;
    enum { dim = Grid::dimension };
    typedef typename Grid::ctype ct;

    typedef typename Grid::template Codim<0>::GlobalIdSet GlobalIdSet;
    const GlobalIdSet & globalid = g.globalIdSet();

    LevelIterator e = g.template lbegin<0>(0);
    const LevelIterator eend = g.template lend<0>(0);

    LevelIterator next = e;
    if (next != eend)
        {
            ++next;
            if (g.name()=="AlbertaGrid") {
                ;//std::cerr << "WARNING: skip indices test using LevelIntersectionIterator for AlbertaGrid!\n";
            } else {
                for (;e != eend; ++e)
                    {
                        // flag vector for elements faces
                        std::vector<bool> visited(e->template count<1>(), false);
                        // loop over intersections
                        IntersectionIterator endit = e->ilevelend();
                        IntersectionIterator it = e->ilevelbegin();
                        // state
                        it.boundary();
                        it.neighbor();
                        // id of boundary segment
                        it.boundaryId();
                        // check id
                        //assert(globalid.id(*e) >= 0);
                        assert(it != endit);

                        if(! e->isLeaf() )
                            {
                                if( e->ileafbegin() != e->ileafend())
                                    {
                                        DUNE_THROW(CheckError, "On non-leaf entities ileafbegin should be equal to ileafend!");
                                    }
                            }

                        // for all intersections
                        for(; it != endit; ++it)
                            {
                                // mark visited face
                                visited[it.numberInSelf()] = true;
                                // check id
                                assert(globalid.id(*(it.inside())) ==
                                       globalid.id(*e));

                                // numbering
                                int num = it.numberInSelf();
                                assert( num >= 0 && num < e->template count<1> () );

                                if(it.neighbor())
                                    {
                                        // geometry
                                        it.intersectionNeighborLocal();
                                        // numbering
                                        num = it.numberInNeighbor();
                                        assert( num >= 0 && num < it.outside()->template count<1> () );
                                    }

                                // geometry
                                it.intersectionSelfLocal();
                                it.intersectionGlobal();

                                // normal vectors
                                Dune::FieldVector<ct, dim-1> v(0);
                                it.outerNormal(v);
                                it.integrationOuterNormal(v);
                                it.unitOuterNormal(v);
                                // search neighbouring cell
                                if (it.neighbor())
                                    {
                                        //assert(globalid.id(*(it.outside())) >= 0);
                                        assert(globalid.id(*(it.outside())) !=
                                               globalid.id(*e));

                                        LevelIterator n    = g.template lbegin<0>(e->level());
                                        LevelIterator nend = g.template lend<0>  (e->level());

                                        while (n != it.outside() && n != nend)
                                            {
                                                assert(globalid.id(*(it.outside())) !=
                                                       globalid.id(*n));
                                                ++n;
                                            }
                                    }
                            }
                        // check that all faces were visited
                        for (size_t i=0; i<visited.size(); i++) assert(visited[i] == true);
                    }
            }
        }
}

template <class GridType, bool c>
struct CheckMark
{
    template <class IteratorType>
    static void check(GridType & grid, IteratorType & it)
    {
        // last marker is 0, so the grid is not changed after this check
        const int refCount[4] = {1,0,-1,0};
        for(int k=0; k<4; ++k)
            {
                // mark entity
                bool marked = grid.mark( refCount[k] , it);
                // if element was marked, check that the marker was set correctly
                if(marked)
                    {
                        // now getMark should return the mark we just set, otherwise error
                        if( grid.getMark(it) != refCount[k] )
                            DUNE_THROW(CheckError,"mark/getMark method not working correctly!");
                    }
            }
    }
};

template <class GridType>
struct CheckMark<GridType,false>
{
    template <class IteratorType>
    static void check(const GridType & grid, IteratorType & it )
    {
    }
};

/*
 * Iterate over the grid und do some runtime checks
 */

template <bool checkMark , class Grid>
void iterate(Grid &g)
{
    typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
    typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;
    typedef typename Grid::template Codim<0>::HierarchicIterator HierarchicIterator;
    typedef typename Grid::template Codim<0>::Geometry Geometry;
    int l = g.maxLevel();
    LevelIterator it = g.template lbegin<0>(l);
    const LevelIterator endit = g.template lend<0>(l);

    Dune::FieldVector<typename Grid::ctype, Grid::dimension> origin(1);
    Dune::FieldVector<typename Grid::ctype, Grid::dimension> result;

    for (;it != endit; ++it)
        {
            LevelIterator l1 = it;
            LevelIterator l2 = l1; ++l1;
            assert(l2 == it);
            assert(l1 != it);
            ++l2;
            assert(l1 == l2);

            result = it->geometry().local(it->geometry().global(origin));
            typename Grid::ctype error = (result-origin).two_norm();
            if(error >= factorEpsilon * std::numeric_limits<typename Grid::ctype>::epsilon())
                {
                    DUNE_THROW(CheckError, "|| geom.local(geom.global(" << origin
                               << ")) - origin || != 0 ( || " << result << " - origin || ) = " << error);
                };
            it->geometry().integrationElement(origin);
            if((int)Geometry::coorddimension == (int)Geometry::mydimension)
                it->geometry().jacobianInverseTransposed(origin);

            it->geometry().type();
            it->geometry().corners();
            it->geometry().corner(0);

        }

    typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;
    LeafIterator lit = g.template leafbegin<0>();
    const LeafIterator lend = g.template leafend<0>();

    // if empty grid, do nothing
    if(lit == lend) return;

    for (;lit != lend; ++lit)
        {
            LeafIterator l1 = lit;
            LeafIterator l2 = l1; ++l1;
            assert(l2 == lit);
            assert(l1 != lit);
            ++l2;
            assert(l1 == l2);

            // leaf check
            if( !lit->isLeaf() )
                DUNE_THROW(CheckError,"LeafIterator gives non-leaf entity!");

            // check adaptation mark for leaf entity mark
            CheckMark<Grid,checkMark>::check(g,lit);

            result = lit->geometry().local(lit->geometry().global(origin));
            typename Grid::ctype error = (result-origin).two_norm();
            if(error >= factorEpsilon * std::numeric_limits<typename Grid::ctype>::epsilon())
                {
                    DUNE_THROW(CheckError, "|| geom.local(geom.global(" << origin
                               << ")) - origin || != 0 ( || " << result << " - origin || ) = " << error);
                };
            lit->geometry().integrationElement(origin);
            if((int)Geometry::coorddimension == (int)Geometry::mydimension)
                lit->geometry().jacobianInverseTransposed(origin);

            lit->geometry().type();
            lit->geometry().corners();
            lit->geometry().corner(0);
        }

}

template <class Grid>
void iteratorEquals (Grid &g)
{
    typedef typename Grid::template Codim<0>::LevelIterator LevelIterator;
    typedef typename Grid::template Codim<0>::LeafIterator LeafIterator;
    typedef typename Grid::template Codim<0>::HierarchicIterator HierarchicIterator;
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;
    typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;

    LevelIterator l1 = g.template lbegin<0>(0);
    LevelIterator l2 = g.template lbegin<0>(0);
    LeafIterator L1 = g.template leafbegin<0>();
    LeafIterator L2 = g.template leafbegin<0>();

    if (l1 == g.template lend<0>(0))
        return;

    HierarchicIterator h1 = l1->hbegin(99);
    HierarchicIterator h2 = l2->hbegin(99);
    IntersectionIterator i1 = l1->ileafbegin();
    IntersectionIterator i2 = l2->ileafbegin();
    EntityPointer e1 = l1;
    EntityPointer e2 = h2;

    // assign
    l1 = l2;
    L1 = L2;
    h1 = h2;
    i1 = i2;
    e1 = e2;

    // equals
#define TestEquals(i) {                                                 \
        i == e2;                                                        \
        i == l2;                                                        \
        i == h2;                                                        \
        i == L2;                                                        \
        if (i2 != l2->ileafend()) i == i2.inside();                     \
        if (i2 != l2->ileafend() && i2.neighbor()) i == i2.outside();   \
    }
    TestEquals(e1);
    TestEquals(l1);
    TestEquals(h1);
    TestEquals(L1);
    if (i1 != l1->ileafend()) TestEquals(i1.inside());
    if (i1 != l1->ileafend() && i1.neighbor()) TestEquals(i1.outside());
}

template <class Grid>
void gridcheck (Grid &g)
{
    /*
     * first do the compile-test: this will not produce any code but
     * fails if an interface-component is missing
     */
    GridInterface<Grid>();

    enum { dim      = Grid :: dimension };
    enum { dimworld = Grid :: dimensionworld };
    typedef typename Grid  :: ctype ctype;
    typedef typename Grid  :: GridFamily GridFamily;

    // type of GridInterface == GridDefaultImplementation
    typedef Dune::GridDefaultImplementation<dim,dimworld,ctype,GridFamily> GridIF;
    const GridIF & gridIF = g;
    // check functionality when grid is interpreted as reference to interface
    GridInterface<GridIF>::check(gridIF);
    /*
     * now the runtime-tests
     */
    const Grid & cg = g;
    iteratorEquals(g);
    iteratorEquals(cg);
    iterate<true>(g);
    iterate<false>(cg);
    zeroEntityConsistency(g);
    zeroEntityConsistency(cg);
    assertNeighbor(g);
    assertNeighbor(cg);
    // note that for some grid this might fail
    // then un comment this test
    Dune::checkIndexSet (g,g.leafIndexSet(),Dune::dvverb);
    for(int lvl = 0; lvl <= g.maxLevel () ; lvl ++ )
        Dune::checkIndexSet (g,g.levelIndexSet(lvl), Dune::dvverb,true);
}
