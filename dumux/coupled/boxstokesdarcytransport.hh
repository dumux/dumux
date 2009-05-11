// $Id: coupledfirstsecond.hh 1128 2009-02-06 09:07:10Z anneli $

#ifndef DUNE_BOXSTOKESDARCYTRANSPORT_HH
#define DUNE_BOXSTOKESDARCYTRANSPORT_HH

#include "coupledbox.hh"
#include "dumux/nonlinear/newtonmethodmatrix.hh"

namespace Dune
{
/** \todo Please doc me! */
template<class StokesModel, class DarcyModel, class Scalar>
class BoxStokesDarcyTransport : public CoupledBox<StokesModel, DarcyModel, BoxStokesDarcyTransport<StokesModel, DarcyModel, Scalar> >
{
public:
    typedef CoupledBox<StokesModel, DarcyModel, BoxStokesDarcyTransport<StokesModel, DarcyModel, Scalar> > BaseType;
    typedef BoxStokesDarcyTransport<StokesModel, DarcyModel, Scalar> ThisType;
    typedef typename BaseType::FirstGrid StokesGrid;
    typedef typename BaseType::SecondGrid DarcyGrid;
    typedef typename StokesGrid::HostGridType HostGrid;
    typedef BlockVector<FieldVector<Scalar, 1> > VectorType;
    typedef BCRSMatrix<FieldMatrix<Scalar, 1, 1> > MatrixType;
    enum {dim = StokesGrid::dimension};
    // indices for Stokes domain (1)
    enum {velocityXIdx1=0, velocityYIdx1=1, velocityZIdx1=2,
    	  massFracIdx1=dim, pressureIdx1=dim+1};
    // indices for Darcy domain (2)
    enum {pressureIdx2 = 0, switchIdx2 = 1};

    /** \todo Please doc me! */
    template <class FirstFV, class SecondFV>
    void localCoupling12(FirstFV& stokesSol, SecondFV& darcySol,
                         int stokesIndex, int darcyIndex,
                         FieldVector<double, dim> globalPos,
                         FieldVector<double, dim> normal, FirstFV& result)
    {
        result = 0;

        for (int comp = 0; comp < dim; comp++)
            result[comp] = normal[comp];

        // pressure times normal, continuity of normal traction
        // this is used as NEUMANN condition for the Stokes domain (momentum balance)
        result *= -darcySol[pressureIdx2];

        return;
    }

    /** \todo Please doc me! */
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

        Scalar xWG = stokesSol[massFracIdx1];
        Scalar xAG = 1-xWG;

        //TODO: multiply by density, fluid properties?
        // velocity times normal
        result[pressureIdx2] = -xWG*(stokesVel*normal);
        result[switchIdx2] = -xAG*(stokesVel*normal);

        return;
    }

    virtual void update(double& dt)
    {
        this->firstModel().localJacobian().setDt(dt);
        this->firstModel().localJacobian().setOldSolution(this->firstModel().uOldTimeStep);
        this->secondModel().localJacobian().setDt(dt);
        this->secondModel().localJacobian().setOldSolution(this->secondModel().uOldTimeStep);
        double dtol = 1e3;
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

    BoxStokesDarcyTransport(const StokesGrid& stokesGrid, StokesModel& stokesModel,
                            const DarcyGrid& darcyGrid, DarcyModel& darcyModel, bool assembleGlobalSystem)
        : BaseType(stokesGrid, stokesModel, darcyGrid, darcyModel, assembleGlobalSystem)
    {}
};

}

#endif
