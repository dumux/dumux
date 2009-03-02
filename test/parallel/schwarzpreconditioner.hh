// $Id: schwarz.hh 120 2006-05-16 08:48:05Z christi $
#ifndef __DUNE_SCHWARZPREC_HH__
#define __DUNE_SCHWARZPREC_HH__

#include<iostream>
#include<vector>
#include<math.h>
#include<stdio.h>
#include<dune/common/helpertemplates.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/scalarproducts.hh>
#include<dune/istl/io.hh>

/*! \brief Single grid overlapping Schwarz for arbitrary subdomain solver
 */
template<class X, class Y, class S, class C>
class AdditiveOverlappingSchwarzPreconditioner : public Dune::Preconditioner<X,Y> {
public:
    //! export types, they come from the derived class
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;

    //! define the category
    enum {category=Dune::SolverCategory::overlapping};

    /*! \brief Constructor.

      constructor gets all parameters to operate the prec.
    */
    AdditiveOverlappingSchwarzPreconditioner (S& solver_, const C& communication_)
        : solver(solver_), communication(communication_)
    {
        k = 0;
    }

    //! prepare: nothing to do here
    virtual void pre (X& x, Y& b)
    {
        communication.copyOwnerToAll(x,x); // make dirichlet values consistent
    }

    //! just calls the istl functions
    virtual void apply (X& v, const Y& d)
    {
        k++;

        // copy rhs because solver overwrites rhs argument
        Y rhs(d);
        communication.addOwnerOverlapToAll(rhs,rhs);

        // solve subdomain problem
        //     char title[32], rowtitle[32];
        //     sprintf(rowtitle,"[%3d]",communication.communicator().rank());
        //     sprintf(title,"rhs before rank=%d k=%d",communication.communicator().rank(),k);
        //     communication.communicator().barrier();
        //     printvector(std::cout,rhs,title,rowtitle,4,9,1);
        //     sprintf(title,"correction before rank=%d k=%d",communication.communicator().rank(),k);
        //     communication.communicator().barrier();
        //     printvector(std::cout,v,title,rowtitle,4,9,1);

        Dune::InverseOperatorResult r;
        solver.apply(v,rhs,r);

        //     sprintf(title,"correction after solve rank=%d k=%d",communication.communicator().rank(),k);
        //     communication.communicator().barrier();
        //     printvector(std::cout,v,title,rowtitle,4,9,1);

        // add up all corrections
        communication.addOwnerOverlapToAll(v,v);

        //      sprintf(title,"correction after communicate rank=%d k=%d",communication.communicator().rank(),k);
        //     communication.communicator().barrier();
        //     printvector(std::cout,v,title,rowtitle,4,9,1);
    }

    // nothing to do here
    virtual void post (X& x)
    {
    }

private:
    S& solver;        // subdomain solver
    const C& communication; // communication object
    int k;
};

#endif
