#ifndef DUNE_NONLINEARMODEL_HH
#define DUNE_NONLINEARMODEL_HH

#include<dune/disc/functions/p1function.hh>
//#include<dune/disc/operators/p1operator.hh>
#include"dumux/operators/p1operatorextended.hh"

namespace Dune
{
template<class G, class RT, class ProblemType, class LocalJacobian, class FunctionType, class OperatorAssembler>
class NonlinearModel {
public:	
	ProblemType& problem;
	FunctionType u;
	FunctionType f;
	OperatorAssembler A;
	LocalJacobian localJacobian;

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

	//! always define virtual destructor in abstract base class
	virtual ~NonlinearModel () {}

	NonlinearModel(const G& g, ProblemType& prob)
: problem(prob), u(g, g.overlapSize(0)==0), f(g, g.overlapSize(0)==0), A(g, g.overlapSize(0)==0), 
localJacobian(prob, false, g, u, g.overlapSize(0)>0)
//: problem(prob), u(g), f(g), A(g), 
//localJacobian(prob, false, g, u)
	{ }

	NonlinearModel(const G& g, ProblemType& prob, int level)
	: problem(prob), u(g, level, g.overlapSize(0)==0), f(g, level, g.overlapSize(0)==0), A(g, level, g.overlapSize(0)==0), 
	localJacobian(prob, false, g, u, g.overlapSize(0)>0)
	{ }
};

template<class G, class RT, class ProblemType, class LocalJacobian, int m=1>
class LeafP1NonlinearModel 
: public NonlinearModel<G, RT, ProblemType, LocalJacobian, LeafP1Function<G, RT, m>, LeafP1OperatorAssembler<G, RT, m> >
{
public:
	// define the function type:
	typedef LeafP1Function<G, RT> FunctionType;

	// define the operator assembler type:
	typedef LeafP1OperatorAssembler<G, RT, m> OperatorAssembler;

	typedef NonlinearModel<G, RT, ProblemType, LocalJacobian, 
	FunctionType, OperatorAssembler> NonlinearModel;

	LeafP1NonlinearModel (const G& g, ProblemType& prob) 
	: NonlinearModel(g, prob)
	{ }
};

}
#endif
