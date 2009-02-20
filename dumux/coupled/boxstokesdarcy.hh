// $Id: coupledfirstsecond.hh 1128 2009-02-06 09:07:10Z anneli $

#ifndef DUNE_BOXSTOKESDARCY_HH
#define DUNE_BOXSTOKESDARCY_HH

#include "coupledbox.hh"
#include "dumux/nonlinear/newtonmethodmatrix.hh"

namespace Dune
{
/** \todo Please doc me! */
template<class StokesModel, class DarcyModel, class Scalar>
class BoxStokesDarcy : public CoupledBox<StokesModel, DarcyModel, BoxStokesDarcy<StokesModel, DarcyModel, Scalar> >
{
public:
    typedef CoupledBox<StokesModel, DarcyModel, BoxStokesDarcy<StokesModel, DarcyModel, Scalar> > BaseType;
    typedef BoxStokesDarcy<StokesModel, DarcyModel, Scalar> ThisType;
    typedef typename BaseType::FirstGrid StokesGrid;
    typedef typename BaseType::SecondGrid DarcyGrid;
    typedef typename StokesGrid::HostGridType HostGrid;
    typedef BlockVector<FieldVector<Scalar, 1> > VectorType;
    typedef BCRSMatrix<FieldMatrix<Scalar, 1, 1> > MatrixType;
    enum {dim = StokesGrid::dimension};

    template <class FirstFV, class SecondFV>
    void localCoupling12(FirstFV& stokesSol, SecondFV& darcySol, int stokesIndex, int darcyIndex,
                            FieldVector<double, dim> qGlobal, FieldVector<double, dim> normal, FirstFV& result)
    {
        return;
    }

    template <class FirstFV, class SecondFV>
    void localCoupling21(FirstFV& stokesSol, SecondFV& darcySol, int stokesIndex, int darcyIndex,
            FieldVector<double, dim> qGlobal, FieldVector<double, dim> normal, SecondFV& result)
    {
        FieldVector<double, dim> stokesVel(0);
        for (int i = 0; i < dim; i++)
            stokesVel[i] = stokesSol[i];

        //std::cout << "q = " << qGlobal << ", stokesVel = " << stokesVel << ", normal = " << normal<< std::endl;
        result = (stokesVel*normal);
    }

//    template <class A12Type, class A21Type>
//    void assembleCoupling(A12Type& A_12, A21Type& A_21)
//    {
//        BaseType::template assembleCoupling<A12Type, A21Type>(A_12, A_21);
//    }

    virtual void update(double& dt)
    {
        this->firstModel().localJacobian().setDt(dt);
        this->firstModel().localJacobian().setOldSolution(this->firstModel().uOldTimeStep);
        this->secondModel().localJacobian().setDt(dt);
        this->secondModel().localJacobian().setOldSolution(this->secondModel().uOldTimeStep);
        NewtonMethodMatrix<HostGrid, ThisType> newtonMethod(this->firstGrid().getHostGrid(), *this);
        newtonMethod.execute();
        dt = this->firstModel().localJacobian().getDt();
        this->uOld = this->u;

        return;
    }

    double getDt()
    {
        return this->firstModel().localJacobian().getDt();
    }

    void setDt(double& dt)
    {
        this->firstModel().localJacobian().setDt(dt);
    }

    BoxStokesDarcy(const StokesGrid& stokesGrid, StokesModel& stokesModel,
            const DarcyGrid& darcyGrid, DarcyModel& darcyModel, bool assembleGlobalSystem)
    : BaseType(stokesGrid, stokesModel, darcyGrid, darcyModel, assembleGlobalSystem)
    {}
};

}

#endif
