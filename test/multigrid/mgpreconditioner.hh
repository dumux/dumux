#include <stdio.h>

#include <map>

#include <dune/disc/operators/p1mgtransfer.hh>

#include <dune/common/geometrytype.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/disc/functions/p1function.hh>
#include <dune/disc/operators/p1operator.hh>
#include <dune/disc/operators/p1mgtransfer.hh>

#include <dune/istl/istlexception.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/io.hh>
#include <dune/istl/gsetc.hh>
#include <dune/istl/ilu.hh>

namespace Dune {

/*!
  \brief A geometric multigrid preconditioner for P1 elements on the leaf grid.
*/
template<class Matrix, class G, class LS, class X, class Y, int m=1>
class SeqP1GeomMG : public Preconditioner<X,Y> {
public:

    typedef typename G::ctype DT;
    enum {n=G::dimension};

    template<int dim>
    struct P1Layout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.isVertex();
        }
    };

    typedef typename G::template Codim<0>::LevelIterator ElementLevelIterator;
    typedef typename G::template Codim<0>::LeafIterator ElementLeafIterator;
    typedef typename G::template Codim<0>::LevelIndexSet LevelIndexSet;
    typedef typename G::template Codim<0>::LeafIndexSet LeafIndexSet;
    typedef MultipleCodimMultipleGeomTypeMapper<G,LevelIndexSet,P1Layout> LevelMapper;
    typedef MultipleCodimMultipleGeomTypeMapper<G,LeafIndexSet,P1Layout> LeafMapper;

    typedef    std::map<int,std::pair<int,int> > LeafToLevelMap;

    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;
    //! \brief The matrix type
    typedef typename LevelP1OperatorAssembler<G,field_type,m>::RepresentationType M;

    typedef typename P1MGTransfer<G,field_type,m>::RepresentationType PM;


    // define the category
    enum {
        //! \brief The category the precondtioner is part of.
        category=SolverCategory::sequential};

    /*! \brief Constructor.

      constructor gets all parameters to operate the prec.
      \param grid_   The grid to operate on.
      \param loc_    Local assembler that treats level boundary as Dirichlet
      \param level0_ The coarsest level where exact solver is called
    */
    SeqP1GeomMG (const Matrix& leafmat_, const G& grid_, LS& loc_,
                 int level0_, int n1_, int n2_, bool additive_)
        : leafmat(leafmat_), grid(grid_), loc(loc_), level0(level0_), n1(n1_), n2(n2_), additive(additive_)
    {
        level0 = std::max(0,level0);
        assert(m==1); // needed for the leaftolevelmap

        // create hierarchy
        for (int level = level0; level<=grid.maxLevel(); level++)
            {
                // allocate vectors and matrices
                vHierarchy[level] = new Dune::LevelP1Function<G,field_type>(grid,level);
                dHierarchy[level] = new Dune::LevelP1Function<G,field_type>(grid,level);
                AHierarchy[level] = new Dune::LevelP1OperatorAssembler<G,field_type,m>(grid,level);

                // assemble matrix
                AHierarchy[level]->assemble(loc,*(vHierarchy[level]),*(dHierarchy[level]));

                // assemble grid transfer operators
                if (level>level0)
                    {
                        pHierarchy[level] = new Dune::P1MGTransfer<G,field_type,m>(grid,level);
                        pHierarchy[level]->assemble(loc);
                    }

                // create smoothers
                if (level>level0)
                    {
                        presmoother[level] = new SeqSSOR<M,X,Y>(**AHierarchy[level],n1,1.0);
                        postsmoother[level] = new SeqSSOR<M,X,Y>(**AHierarchy[level],n2,1.0);
                        bpxsmoother[level] = new SeqJac<M,X,Y>(**AHierarchy[level],n1,1.0);
                    }
            }
        alloclevel = grid.maxLevel();

        // create the coarse grid solver objects
        op0 = new MatrixAdapter<M,X,Y>(**AHierarchy[level0]);
        prec0 = new SeqSSOR<M,X,Y>(**AHierarchy[level0],1,1.0);
        //solver0 = new BiCGSTABSolver<X>(*op0,*prec0,1E-4,800,0);
        solver0 = new CGSolver<X>(*op0,*prec0,1E-3,800,0);

        // build leaf to level map to speed up transfer
        // make mapper for the leaf grid
        LeafMapper leafmapper(grid,grid.leafIndexSet());

        // transfer input defect to levels
        for (int level=level0; level<=grid.maxLevel(); level++)
            {
                // make a mapper for this level
                LevelMapper levelmapper(grid,grid.levelIndexSet(level));

                // allocate a flag vector to handle each vertex exactly once
                std::vector<bool> treated(levelmapper.size(),false);

                // traverse elements
                for (ElementLevelIterator it = grid.template lbegin<0>(level);
                     it!=grid.template lend<0>(level); ++it)
                    {
                        // need only consider leaf elements
                        if (!it->isLeaf()) continue;

                        // match level and leaf index for all vertices
                        for (int i=0; i<it->template count<n>(); i++)
                            {
                                // compute index of vertex in level grid
                                int levelindex = levelmapper.template map<n>(*it,i);

                                // compute index of vertex in leaf grid
                                int leafindex = leafmapper.template map<n>(*it,i);

                                // handle each vertex only once
                                if (treated[levelindex]) continue;
                                treated[levelindex] = true;

                                // on level 0 just copy
                                if (level==level0)
                                    {
                                        leaftolevelmap[leafindex] = std::pair<int,int>(level,levelindex);
                                        continue;
                                    }

                                // on higher levels consider the skip flags; THIS WORKS ONLY FOR  m=1 !!
                                if (!pHierarchy[level]->skipFlag(0,levelindex))
                                    //                     {
                                    //                       // find entry
                                    //                       for (typename PM::ColIterator cit=(**pHierarchy[level])[levelindex].begin();
                                    //                            cit!=(**pHierarchy[level])[levelindex].end(); ++cit)
                                    //                         if ( std::abs((*cit)[0][0]-1)<1E-6 && pHierarchy[level]->coarseSkipFlag(0,cit.index())==false )
                                    //                           {
                                    //                             std::cout << "Bingo: "
                                    //                                       << " i_leaf=" << leafindex
                                    //                                       << " i_fine=" << levelindex
                                    //                                       << " i_coarse=" << cit.index()
                                    //                                       << std::endl;
                                    //                             // entries in leaf matrix
                                    //                             for (typename Matrix::ConstColIterator lit=leafmat[leafindex].begin();
                                    //                                  lit!=leafmat[leafindex].end(); ++lit)
                                    //                               if (leaftolevelmap.find(lit.index())!=leaftolevelmap.end())
                                    //                                 std::cout << "BingoBingo: "
                                    //                                           << " i_leaf=" << leafindex
                                    //                                           << " j_leaf=" << lit.index()
                                    //                                           << " j_level=" << leaftolevelmap.find(lit.index())->second.first
                                    //                                           << std::endl;
                                    //                           }
                                    //                     }
                                    //                   else
                                    {
                                        if (leaftolevelmap.find(leafindex)!=leaftolevelmap.end())
                                            {
                                                std::cout << "leafindex alread in map: " << leafindex << " -> "
                                                          << leaftolevelmap.find(leafindex)->second.first << ","
                                                          << leaftolevelmap.find(leafindex)->second.second
                                                          << std::endl;
                                            }
                                        leaftolevelmap[leafindex] = std::pair<int,int>(level,levelindex);
                                    }
                            }
                    }
            }
    }

    ~SeqP1GeomMG ()
    {
        // destroy hierarchy
        for (int level = level0; level<=alloclevel; level++)
            {
                delete AHierarchy[level];
                delete vHierarchy[level];
                delete dHierarchy[level];
                if (level>level0) delete pHierarchy[level];
                if (level>level0) delete presmoother[level];
                if (level>level0) delete postsmoother[level];
                if (level>level0) delete bpxsmoother[level];
            }
        delete op0;
        delete prec0;
        delete solver0;
    }

    /*!
      \brief Prepare the preconditioner.

      \copydoc Preconditioner::pre(X&,Y&)
    */
    virtual void pre (X& x, Y& b)
    {
        for (int level = level0+1; level<=grid.maxLevel(); level++)
            {
                presmoother[level]->pre(**vHierarchy[level],**dHierarchy[level]);
                postsmoother[level]->pre(**vHierarchy[level],**dHierarchy[level]);
                bpxsmoother[level]->pre(**vHierarchy[level],**dHierarchy[level]);
            }
    }

    /*!
      \brief Apply the precondtioner

      \copydoc Preconditioner::apply(X&,Y&)
    */
    virtual void apply (X& v, const Y& d)
    {
        // clear defect and solution on all levels
        for (int level=level0; level<=grid.maxLevel(); level++)
            {
                **vHierarchy[level] = 0;
                **dHierarchy[level] = 0;
            }

        // transfer input defect to levels
        for (LeafToLevelMap::iterator i=leaftolevelmap.begin(); i!=leaftolevelmap.end(); ++i)
            (**dHierarchy[i->second.first])[i->second.second] = d[i->first];

        // call multigrid cycle on hierarchy
        if (additive)
            bpx(grid.maxLevel());
        else
            mgc(grid.maxLevel());

        // transfer correction from levels to leaf output
        for (LeafToLevelMap::iterator i=leaftolevelmap.begin(); i!=leaftolevelmap.end(); ++i)
            v[i->first] = (**vHierarchy[i->second.first])[i->second.second];
    }

    /*!
      \brief Clean up.

      \copydoc Preconditioner::post(X&)
    */
    virtual void post (X& x)
    {
        for (int level = level0+1; level<=grid.maxLevel(); level++)
            {
                presmoother[level]->post(**vHierarchy[level]);
                postsmoother[level]->post(**vHierarchy[level]);
                bpxsmoother[level]->post(**vHierarchy[level]);
            }
    }

private:
    // the multigrid cycle
    void mgc (int level)
    {
        if (level==level0)
            {
                // coarse grid solve
                InverseOperatorResult r;
                solver0->apply(**vHierarchy[level],**dHierarchy[level],r);
            }
        else
            {
                // pre smoothing
                presmoother[level]->apply(**vHierarchy[level],**dHierarchy[level]);

                // recompute defect in a copy
                Y d(**dHierarchy[level]);
                (**AHierarchy[level]).mmv(**vHierarchy[level],d);

                // restrict defect; coarse defect has been cleared, update is OK
                (**pHierarchy[level]).umtv(d,**dHierarchy[level-1]);

                // recursive call; coarse grid solution has already been cleard in apply
                mgc(level-1);

                // prolongate correction; fine correction has been cleared above
                (**pHierarchy[level]).umv(**vHierarchy[level-1],**vHierarchy[level]);

                // fix dirichlet nodes to interpolated value
                for (int i=0; i<(**vHierarchy[level]).size(); ++i)
                    if (pHierarchy[level]->skipFlag(0,i))
                        (**dHierarchy[level])[i] = (**vHierarchy[level])[i];

                // post smoothing
                postsmoother[level]->apply(**vHierarchy[level],**dHierarchy[level]);
            }
    }

    // the multigrid cycle
    void bpx (int level)
    {
        if (level==level0)
            {
                // coarse grid solve
                InverseOperatorResult r;
                solver0->apply(**(vHierarchy)[level],**(dHierarchy)[level],r);
            }
        else
            {
                // restrict defect; coarse defect has been cleared, update is OK
                (**pHierarchy[level]).umtv(**dHierarchy[level],**dHierarchy[level-1]);

                // smoothing
                bpxsmoother[level]->apply(**vHierarchy[level],**dHierarchy[level]);

                // recursive call; coarse grid solution has already cleared been cleard in apply
                bpx(level-1);

                // prolongate correction; fine correction has been cleared above
                (**pHierarchy[level]).umv(**vHierarchy[level-1],**vHierarchy[level]);
            }
    }

    // parameters
    const Matrix& leafmat;
    const G& grid;
    LS& loc;
    int level0;
    int n1,n2;
    bool additive;

    // local variables
    int alloclevel;

    // grid hierarchy
    LevelP1OperatorAssembler<G,field_type,m>* AHierarchy[64];
    LevelP1Function<G,field_type>* vHierarchy[64];
    LevelP1Function<G,field_type>* dHierarchy[64];
    P1MGTransfer<G,field_type>* pHierarchy[64];

    // map leafindex to (level,levelindex)
    LeafToLevelMap leaftolevelmap;

    // smoothers
    SeqSSOR<M,X,Y>* presmoother[64];
    SeqSSOR<M,X,Y>* postsmoother[64];
    SeqJac<M,X,Y>* bpxsmoother[64];

    // coarse grid solver
    MatrixAdapter<M,X,Y>* op0;
    SeqSSOR<M,X,Y>* prec0;
    InverseOperator<X,X>* solver0;
};

}
