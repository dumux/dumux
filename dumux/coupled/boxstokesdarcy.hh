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
    void localCoupling12(FirstFV& stokesSol, SecondFV& darcySol,
                         int stokesIndex, int darcyIndex,
                         FieldVector<double, dim> globalPos,
                         FieldVector<double, dim> normal, FirstFV& result)
    {
        result = 0;

        for (int comp = 0; comp < dim; comp++)
            result[comp] = normal[comp];

        result *= -darcySol;

        return;
    }

    template <class FirstFV, class SecondFV>
    void localCoupling21(FirstFV& stokesSol, SecondFV& darcySol,
                         int stokesIndex, int darcyIndex,
                         FieldVector<double, dim> globalPos,
                         FieldVector<double, dim> normal, SecondFV& result)
    {
        result = 0;

        FieldVector<double, dim> stokesVel(0);
        for (int comp = 0; comp < dim; comp++)
            stokesVel[comp] = stokesSol[comp];

        result = -(stokesVel*normal);

        return;
    }

    virtual void update(double& dt)
    {
        this->firstModel().localJacobian().setDt(dt);
        this->firstModel().localJacobian().setOldSolution(this->firstModel().uOldTimeStep);
        this->secondModel().localJacobian().setDt(dt);
        this->secondModel().localJacobian().setOldSolution(this->secondModel().uOldTimeStep);
        double dtol = 1e-1;
        double rtol = 1e7;
        int maxIt = 10;
        double mindt = 1e-5;
        int goodIt = 4;
        int maxInc = 2;
        NewtonMethodMatrix<HostGrid, ThisType> newtonMethod(this->firstGrid().getHostGrid(), *this,
                                                            dtol, rtol, maxIt, mindt, goodIt, maxInc);
        newtonMethod.execute();
        dt = this->firstModel().localJacobian().getDt();
        this->uOldTimeStep = this->u;

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
