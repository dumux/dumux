// $Id$ 

#ifndef DUNE_COUPLEDMODEL_HH
#define DUNE_COUPLEDMODEL_HH

namespace Dune
{

template<class FirstModel, class SecondModel>
class CoupledModel {
public:
	typedef typename FirstModel::ProblemType FirstProblem;
	typedef typename SecondModel::ProblemType SecondProblem;
	typedef typename FirstModel::GridType FirstGrid;
	typedef typename SecondModel::GridType SecondGrid;
	
	void assemble() 
	{
		
	}
	
	FirstModel& firstModel()
	{
		return firstModel_;
	}
	
	SecondModel& secondModel()
	{
		return secondModel_;
	}

	CoupledModel(const FirstGrid& firstGrid, FirstModel& firstModel, const SecondGrid& secondGrid, SecondModel& secondModel)
	: firstGrid_(firstGrid), secondGrid_(secondGrid), firstModel_(firstModel), secondModel_(secondModel)
	{}
	
protected: 
	const FirstGrid& firstGrid_;
	const SecondGrid& secondGrid_;
	FirstModel& firstModel_;
	SecondModel& secondModel_;
};

}

#endif
