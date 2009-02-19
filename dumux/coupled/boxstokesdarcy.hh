// $Id: coupledfirstsecond.hh 1128 2009-02-06 09:07:10Z anneli $

#ifndef DUNE_BOXSTOKESDARCY_HH
#define DUNE_BOXSTOKESDARCY_HH

#include "coupledbox.hh"

namespace Dune
{
/** \todo Please doc me! */
template<class StokesModel, class DarcyModel>
class BoxStokesDarcy : public CoupledBox<StokesModel, DarcyModel, BoxStokesDarcy<StokesModel, DarcyModel> > 
{
public:
    typedef CoupledBox<StokesModel, DarcyModel, BoxStokesDarcy<StokesModel, DarcyModel> > BaseType;
    typedef typename BaseType::FirstGrid StokesGrid;
    typedef typename BaseType::SecondGrid DarcyGrid;
    enum {dim = StokesGrid::dimension};
    
    template <class FirstFV, class SecondFV>
    void localCoupling12(FirstFV& firstSol, SecondFV& secondSol, int firstIndex, int secondIndex, 
                            FieldVector<double, dim> qGlobal, FirstFV& result)
    {
        return;
    }

    template <class FirstFV, class SecondFV>
    void localCoupling21(FirstFV& firstSol, SecondFV& secondSol, int firstIndex, int secondIndex, 
                            FieldVector<double, dim> qGlobal, SecondFV& result)
    {
        return;
    }

//    template <class A12Type, class A21Type>
//    void assembleCoupling(A12Type& A_12, A21Type& A_21)
//    {
//        BaseType::template assembleCoupling<A12Type, A21Type>(A_12, A_21);
//    }

    BoxStokesDarcy(const StokesGrid& stokesGrid, StokesModel& stokesModel,
            const DarcyGrid& darcyGrid, DarcyModel& darcyModel, bool assembleGlobalSystem)
    : BaseType(stokesGrid, stokesModel, darcyGrid, darcyModel, assembleGlobalSystem)
    {}
};

}

#endif
