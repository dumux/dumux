// $Id$

#ifndef DUNE_NONLINEARMODEL_HH
#define DUNE_NONLINEARMODEL_HH

#include <dune/disc/operators/p1operator.hh>
#include"dumux/operators/mixedoperator.hh"

namespace Dune
{

/** \todo Please doc me! */

template<class G, class RT, class ProblemType, class LocalJacobian, class FunctionType, class OperatorAssembler>
class NonlinearModel {
public:
    typedef typename FunctionType::RepresentationType VectorType;
    typedef typename OperatorAssembler::RepresentationType MatrixType;
    ProblemType& problem;
    FunctionType u;
    FunctionType f;
    OperatorAssembler A;

protected:
    LocalJacobian localJacobian_;
public:

    LocalJacobian &localJacobian()
    { return localJacobian_; };

    const LocalJacobian &localJacobian() const
    { return localJacobian_; };

    //! return const reference to solution vector
    const FunctionType& operator* () const
    {
        return u;
    }

    //! return reference to solution vector
    FunctionType& operator* ()
    {
        return u;
    }

    virtual double getDt()
    {
        return localJacobian().getDt();
    }

    virtual void setDt(double& dt)
    {
        localJacobian().setDt(dt);
    }

    virtual void assemble()
    {
        *f = 0;
        localJacobian().clearVisited();
        A.assemble(localJacobian(), u, f);
    }

    virtual MatrixType& matrix()
    {
        return *A;
    }

    virtual VectorType& rhs()
    {
        return *f;
    }

    virtual VectorType& sol()
    {
        return *u;
    }
    virtual void writerestartfile()
    {
        return;
    }
    virtual void restart()
    {
        return;
    }

    //! always define virtual destructor in abstract base class
    virtual ~NonlinearModel () {}

    NonlinearModel(const G& g, ProblemType& prob)
        : problem(prob), u(g, g.overlapSize(0)==0), f(g, g.overlapSize(0)==0), A(g, g.overlapSize(0)==0),
          localJacobian_(prob, false, g, u, g.overlapSize(0)>0)
          //: problem(prob), u(g), f(g), A(g),
          //localJacobian_(prob, false, g, u)
    { }

    NonlinearModel(const G& g, ProblemType& prob, int level)
        : problem(prob), u(g, level, g.overlapSize(0)==0), f(g, level, g.overlapSize(0)==0), A(g, level, g.overlapSize(0)==0),
          localJacobian_(prob, false, g, u, g.overlapSize(0)>0)
    { }
};

/** \todo Please doc me! */

template<class G, class RT, class ProblemType, class LocalJacobian, int m=1>
class LeafP1NonlinearModel
    : public NonlinearModel<G, RT, ProblemType, LocalJacobian, LeafP1Function<G, RT, m>, LeafP1OperatorAssembler<G, RT, m> >
{
public:
    // define the function type:
    typedef LeafP1Function<G, RT> FunctionType;

    // define the operator assembler type:
    typedef LeafP1OperatorAssembler<G, RT, m> OperatorAssembler;

    typedef NonlinearModel<G, RT, ProblemType, LocalJacobian, FunctionType, OperatorAssembler> ThisNonlinearModel;

    LeafP1NonlinearModel (const G& g, ProblemType& prob)
        : ThisNonlinearModel(g, prob)
    { }
};

/** \todo Please doc me! */

template<class G, class RT, class ProblemType, class LocalJacobian, int m=1>
class MixedNonlinearModel
    : public NonlinearModel<G, RT, ProblemType, LocalJacobian, LeafMixedFunction<G, RT, m>, LeafMixedOperatorAssembler<G, RT, m> >
{
public:
    // define the function type:
    typedef LeafMixedFunction<G, RT> FunctionType;

    // define the operator assembler type:
    typedef LeafMixedOperatorAssembler<G, RT, m> OperatorAssembler;

    typedef NonlinearModel<G, RT, ProblemType, LocalJacobian, FunctionType, OperatorAssembler> ThisNonlinearModel;

    MixedNonlinearModel (const G& g, ProblemType& prob)
        : ThisNonlinearModel(g, prob)
    { }
};

}
#endif
