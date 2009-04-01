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


//First Order Upwind (not working well)
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

//HLL Flux (it seems to work in 1D
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
        Scalar signalVelocityL=0;
        Scalar signalVelocityR=0;
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
        Scalar starVelocity = 0;
        Scalar starWaveCelerity = 0;
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

        //Calculate HLL Values for the star region for velocity and wavecelerity

        starVelocity=(velocityFaceL+velocityFaceR);
        starVelocity*=0.5;
        starVelocity += sqrt(gravity*waterDepthFaceL);
        starVelocity -= sqrt(gravity*waterDepthFaceR);

        //std::cout<<"starVelocity = "<<starVelocity<<std::endl;

        help1 =(velocityFaceL-velocityFaceR);
        help1 *=0.25;
        help2 = sqrt(gravity*waterDepthFaceL)+sqrt(gravity*waterDepthFaceR);
        help2 *=0.5;
        starWaveCelerity=(help1+help2);

        std::cout<<"starVelocity "<<starVelocity<<std::endl;
        std::cout<<"starWaveCelerity "<<starWaveCelerity<<std::endl;
        //   std::cout<<"waterDepthFaceR "<<waterDepthFaceR<<std::endl;
        //   std::cout<<"waterDepthFaceL "<<waterDepthFaceL<<std::endl;


        //Calculate signal velocities used in HLL-flux calculation

        if (waterDepthFaceL == 0)
        {
            signalVelocityR =velocityFaceR;
            signalVelocityR +=sqrt(gravity*waterDepthFaceR);

            signalVelocityL = sqrt(gravity*waterDepthFaceR);
            signalVelocityL*=(-2);
            signalVelocityL+=velocityFaceR;

        }
        if (waterDepthFaceR == 0)
        {
            signalVelocityR=sqrt(gravity*waterDepthFaceL);
            signalVelocityR *=(2);
            signalVelocityR +=velocityFaceL;

            signalVelocityL = velocityFaceL;
            signalVelocityL-=sqrt(gravity*waterDepthFaceL);
        }

        else
        {
            signalVelocityL = std::min((velocityFaceL* -sqrt(gravity
                    *waterDepthFaceL)), (starVelocity-starWaveCelerity));
            signalVelocityR = std::max((velocityFaceR* +sqrt(gravity
                    *waterDepthFaceR)), (starVelocity+starWaveCelerity));
        }

        std::cout<<"signalVelocityL "<<signalVelocityL<<std::endl;
        std::cout<<"signalVelocityR "<<signalVelocityR<<std::endl;

        //Compute the flux vector components Conti, XMomentum, YMomentum for the left and right state


        if (signalVelocityL >= 0)
        {
            partConti = velocityFaceL*waterDepthFaceL;
            fluxConti=partConti;

            partMomentum = waterDepthFaceL*(velocityFaceL*velocityFaceL) + 0.5 * gravity
                    * (waterDepthFaceL*waterDepthFaceL);
            fluxMomentum = partMomentum;
        }
        else if (signalVelocityL <= 0 && signalVelocityR >= 0)
        {
            partContiL = velocityFaceL*waterDepthFaceL;

            partContiR = velocityFaceR*waterDepthFaceR;

            partMomentumL = waterDepthFaceL*(velocityFaceL*velocityFaceL) +0.5 *gravity
                    *(waterDepthFaceL*waterDepthFaceL);

            partMomentumR =waterDepthFaceR*(velocityFaceR*velocityFaceR)+0.5 *gravity
                    *(waterDepthFaceR*waterDepthFaceR);

            Scalar denominator = (signalVelocityR-signalVelocityL);
            
            fluxConti=signalVelocityR*(partContiL);
            fluxConti -= signalVelocityL*(partContiR);
            fluxConti += (signalVelocityR*signalVelocityL*(waterDepthFaceR
                    -waterDepthFaceL));
            fluxConti/= denominator;

            fluxMomentum=signalVelocityR*(partMomentumL);
            fluxMomentum-= signalVelocityL*(partMomentumR);
            fluxMomentum+= (signalVelocityR*signalVelocityL*(waterDepthFaceR
                    *velocityFaceR - waterDepthFaceL *velocityFaceL));
            fluxMomentum/= denominator;

        }
        else if (signalVelocityR <= 0)
        {
            partConti = velocityFaceR *waterDepthFaceR;
            fluxConti=partConti;

            partMomentum = waterDepthFaceR*(velocityFaceR*velocityFaceR) +0.5 *gravity
                    *(waterDepthFaceR*waterDepthFaceR);
            fluxMomentum = partMomentumR;
        }
        else
        {
            DUNE_THROW(NotImplemented, "Flux Range");
        }

        //      std::cout<<"partContiL "<<partContiL<<std::endl;
        //     std::cout<<"partContiR "<<partContiR<<std::endl;
        //    std::cout<<"fluxConti "<<fluxConti<<std::endl;
        //    std::cout<<"fluxMomentum "<<fluxMomentum<<std::endl;
        //   std::cout<<"partMomentumL "<<partMomentumL<<std::endl;
        //   std::cout<<"partMomentumR "<<partMomentumR<<std::endl;


        SystemType fluxHllFace(0);
        fluxHllFace[0]=fluxConti*nVec;
        fluxHllFace[1]=fluxMomentum*nVec;

        std::cout<<"Hll flux="<<fluxHllFace<<std::endl;

        return fluxHllFace;
    }
}

;
}
#endif

