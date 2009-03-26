// $Id$

#ifndef NUMERICALFLUX_HH
#define NUMERICALFLUX_HH

namespace Dune
{
template<class Grid, class Scalar> class NumericalFlux
{
public:

    enum
    {   dim=Grid::dimension};
    enum
    {   dimWorld = Grid::dimensionworld};
    typedef Dune::FieldVector<Scalar, dim> VelType;
    typedef Dune::FieldVector<Scalar, dim+1> SystemType;
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;

    virtual SystemType operator()(VelType velFaceL, VelType velFaceR, Scalar wDepthFaceL,
            Scalar wDepthFaceR, VelType nVec)=0;

    virtual ~NumericalFlux()
    {
    }
};

//First Order Upwind Method
template<class Grid, class Scalar> class FirstOrderUpwindModified :
    public NumericalFlux<Grid,Scalar>
{
    enum
    {   dim=Grid::dimension};
    typedef Dune::FieldVector<Scalar, dim> VelType;
    typedef Dune::FieldVector<Scalar, dim+1> SystemType; //System has 3 equations -> 3 positions in a vector
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

    VelType velocityI;
    VelType velocityJ;
    Scalar waterDepthI;
    Scalar waterDepthJ;

private:
    Scalar eps;
    Scalar gravity;

public:

    FirstOrderUpwindModified() :
        NumericalFlux<Grid, Scalar>(), gravity(9.81), eps(1e-8)
    {

    }
    SystemType operator()(VelType velocityI, VelType velocityJ, Scalar waterDepthI,
            Scalar waterDepthJ, VelType nVec)
    {
        Scalar direction =0;
        VelType unityVector=1;
        Scalar contiI =0;
        Scalar contiJ =0;
        Scalar momentumI=0;
        Scalar momentumJ=0;

        Scalar fluxConti=0;
        Scalar fluxMomentumXI=0;
        Scalar fluxMomentum=0;
        SystemType fluxUpwindFace(0);

        direction= nVec;
        direction*= unityVector;
        Scalar help=0;

        contiI = velocityI*waterDepthI;
        contiJ = velocityJ*waterDepthJ;
        momentumI = waterDepthI*pow(velocityI, 2)+0.5*gravity*pow(waterDepthI,
                2);
        momentumJ = waterDepthJ*pow(velocityJ, 2)+0.5*gravity*pow(waterDepthJ,
                2);

        if (direction >0)
        {
            fluxConti = std::max(contiI, 0.0)+std::max(contiJ*(-1), 0.0);

            if (momentumI==momentumJ)
            {
                fluxMomentum = momentumI;
            }
            else
            {
                fluxMomentum = std::max(momentumI, 0.0) +std::max(momentumJ,
                        0.0);
            }

        }
        if (direction <0)
        {
            fluxConti = std::max(contiJ, 0.0)+std::max(contiI*(-1), 0.0);

            if (momentumI==momentumJ)
            {
                fluxMomentum = momentumI;
            }
            else
            {
                fluxMomentum = std::max(momentumI, 0.0) +std::max(momentumJ,
                        0.0);
            }
        }

        fluxUpwindFace[0]=fluxConti*(-1)*nVec;
        fluxUpwindFace[1]=fluxMomentum*(-1)*nVec;

        //   std::cout<<fluxUpwindFace<<std::endl;
        return fluxUpwindFace;

    }

};

//First Order Upwind Method
template<class Grid, class Scalar> class FirstOrderUpwind :
    public NumericalFlux<Grid,Scalar>
{
    enum
    {   dim=Grid::dimension};
    typedef Dune::FieldVector<Scalar, dim> VelType;
    typedef Dune::FieldVector<Scalar, dim+1> SystemType; //System has 3 euqations -> 3 positions in a vector
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

    Scalar eps;
    Scalar gravity;

public:

    FirstOrderUpwind() :
        NumericalFlux<Grid, Scalar>(), gravity(9.81), eps(1e-8)
    {

    }

    SystemType operator()(VelType velocityI, VelType velocityJ, Scalar waterDepthI,
            Scalar waterDepthJ, VelType nVec)
    {
        Scalar direction =0;
        VelType unityVector=1;
        Scalar contiI =0;
        Scalar contiJ =0;
        Scalar momentumI=0;
        Scalar momentumJ=0;

        Scalar fluxConti=0;
        Scalar fluxMomentumXI=0;
        Scalar fluxMomentum=0;
        SystemType fluxUpwindFace(0);

        direction = nVec*unityVector;

        contiI = velocityI*waterDepthI;
        contiJ=velocityJ*waterDepthJ;
        momentumI = waterDepthI*pow(velocityI, 2)+0.5*gravity*pow(waterDepthI,
                2);
        momentumJ = waterDepthJ*pow(velocityJ, 2)+0.5*gravity*pow(waterDepthJ,
                2);

        if (velocityI*nVec > 0 && velocityJ*nVec >= 0)
        {
            fluxUpwindFace[0]=contiI*nVec;
        }
        if (velocityJ*nVec*(-1) > 0 && velocityI*nVec*(-1) >= 0)
        {
            fluxUpwindFace[0]=contiJ*nVec;
        }
        else
        {
            Scalar contiSum = contiI+contiJ;
            if (contiSum *nVec > 0)
            {
                fluxUpwindFace[0]= contiSum*nVec;
            }
            if (contiSum *nVec < 0)
            {
                fluxUpwindFace[0]= contiSum*nVec;
            }
            else
            {
                fluxUpwindFace[0]=0;
            }
        }

        if (momentumI > momentumJ)
        {
            fluxUpwindFace[1]= momentumI*nVec;
        }
        if (momentumI==momentumJ)
        {
            fluxUpwindFace[1]=0;
        }
        else
        {
            fluxUpwindFace[1]= momentumJ*(-1)*nVec;
        }

        //   std::cout<<fluxUpwindFace<<std::endl;
        return fluxUpwindFace;

    }
};

//HLL numerical flux: calculates flux vector through current intersection by means of reconstructed data
template<class Grid, class Scalar> class HllFlux :
    public NumericalFlux<Grid,Scalar>
{
    enum
    {   dim=Grid::dimension};
    typedef Dune::FieldVector<Scalar, dim> VelType;
    typedef Dune::FieldVector<Scalar, dim+1> SystemType; //System has 3 equations -> 3 positions in a vector
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

    Scalar gravity;
    Scalar eps;

public:

    HllFlux() :
        NumericalFlux<Grid, Scalar>(), gravity(9.81), eps(1e-8)
    {

    }

    SystemType operator()(VelType velocityI, VelType velocityJ, Scalar waterDepthI,
            Scalar waterDepthJ, VelType nVec)
    {
        Scalar sigVelL=0;
        Scalar sigVelR=0;
        Scalar waterDepthFaceL=0;
        Scalar waterDepthFaceR=0;
        Scalar velocityFaceR=0;
        Scalar velocityFaceL=0;
        Scalar fluxConti=0;
        Scalar fluxMomentum = 0;
        Scalar partConti = 0;
        Scalar partMomentum = 0;
        Scalar partContiL=0;
        Scalar partMomentumL =0;
        Scalar partContiR=0;
        Scalar partMomentumR =0;
        Scalar velAvg = 0;
        Scalar gravWDepthAvg = 0;
        Scalar help1 = 0;
        Scalar help2 =0;

        VelType unityVector(1);
        Scalar direction;

        direction = nVec*unityVector;

        //determine left and right side value of a face
        if (direction > 0)
        {
            waterDepthFaceL = waterDepthI;
            waterDepthFaceR = waterDepthJ;
            velocityFaceL = velocityI;
            velocityFaceR = velocityJ;
        }

        if (direction < 0)
        {
            waterDepthFaceL = waterDepthJ;
            waterDepthFaceR = waterDepthI;
            velocityFaceL = velocityJ;
            velocityFaceR = velocityI;
        }

        //Calculate HLL Average Values for velocity and waterdepth

        velAvg=(velocityFaceL+velocityFaceR)*nVec;
        velAvg*=0.5;
        velAvg += sqrt(gravity*waterDepthFaceL);
        velAvg -= sqrt(gravity*waterDepthFaceR);

        //std::cout<<"velAvg = "<<velAvg<<std::endl;

        help1 =(velocityFaceL-velocityFaceR)*nVec;
        help1 *=0.25;
        help2 = sqrt(gravity*waterDepthFaceL)+sqrt(gravity *waterDepthFaceR);
        help2 *=0.5;
        gravWDepthAvg=(help1+help2);

        //    std::cout<<"VelAvg "<<velAvg<<std::endl;
        //     std::cout<<"gravWDepthAvg "<<gravWDepthAvg<<std::endl;


        //Calculate signal velocities used in HLL-flux calculation

        if (waterDepthFaceL == 0)
        {
            sigVelR=velocityFaceR*nVec;
            sigVelR+=sqrt(gravity*waterDepthFaceR);

            sigVelL = sqrt(gravity*waterDepthFaceR);
            sigVelL*=(-2);
            sigVelL+=velocityFaceR*nVec;

        }
        if (waterDepthFaceR == 0)
        {
            sigVelR=sqrt(gravity*waterDepthFaceL);
            sigVelR *=(2);
            sigVelR +=velocityFaceL*nVec;

            sigVelL = velocityFaceL*nVec;
            sigVelL-=sqrt(gravity*waterDepthFaceL);
        }

        else
        {
            sigVelL = std::min((velocityFaceL*nVec -sqrt(gravity
                    *waterDepthFaceL)), (velAvg-gravWDepthAvg));
            sigVelR = std::max((velocityFaceR*nVec +sqrt(gravity
                    *waterDepthFaceR)), (velAvg+gravWDepthAvg));
        }

        std::cout<<"sigVelL "<<sigVelL<<std::endl;
        std::cout<<"sigVelR "<<sigVelR<<std::endl;

        //Compute the flux vector components Conti, XMomentum, YMomentum for the left and right state
        //Initialize vectors and fill with flux data


        if (sigVelL >= 0 && sigVelR ==0)
        {
            partConti = velocityFaceL*waterDepthFaceL;
            fluxConti=partContiL*nVec;

            partMomentum = waterDepthFaceL*pow(velocityFaceL, 2) +0.5 *gravity
                    *pow(waterDepthFaceL, 2);
            fluxMomentum = partMomentum*nVec;
        }

        if (sigVelL <= 0 <= sigVelR)
        {
            partContiL = velocityFaceL*waterDepthFaceL;

            partContiR = velocityFaceR*waterDepthFaceR;

            partMomentumL = waterDepthFaceL*pow(velocityFaceL, 2) +0.5 *gravity
                    *pow(waterDepthFaceL, 2);

            partMomentumR =waterDepthFaceR*pow(velocityFaceR, 2)+0.5 *gravity
                    *pow(waterDepthFaceR, 2);

            fluxConti=sigVelR*(partContiL*nVec);
            fluxConti -= sigVelL*(partContiR*nVec);
            fluxConti += (sigVelR*sigVelL*(waterDepthFaceR-waterDepthFaceL)*nVec);
            fluxConti/= (sigVelR-sigVelL);

            fluxMomentum=sigVelR*(partMomentumL*nVec);
            fluxMomentum-= sigVelL*(partMomentumR*nVec);
            fluxMomentum+= (sigVelR*sigVelL*(waterDepthFaceR *velocityFaceR
                    -waterDepthFaceL *velocityFaceL)*nVec);
            fluxMomentum/=(sigVelR-sigVelL);

        }
        if (sigVelR <= 0 && sigVelL == 0)
        {
            partConti = velocityFaceR *waterDepthFaceR;
            fluxConti=partContiR*nVec;

            partMomentum = waterDepthFaceR*pow(velocityFaceR, 2) +0.5 *gravity
                    *pow(waterDepthFaceR, 2);
            fluxMomentum = partMomentumR*nVec;
        }

        //      std::cout<<"partContiL "<<partContiL<<std::endl;
        //     std::cout<<"partContiR "<<partContiR<<std::endl;
        //    std::cout<<"fluxConti "<<fluxConti<<std::endl;
        //   std::cout<<"fluxMomentum "<<fluxMomentum<<std::endl;


        SystemType fluxHllFace(0);
        fluxHllFace[0]=fluxConti*(-1);
        fluxHllFace[1]=fluxMomentum*(-1);

        std::cout<<"Hll flux="<<fluxHllFace<<std::endl;

        return fluxHllFace;
    }
}

;
}
#endif

