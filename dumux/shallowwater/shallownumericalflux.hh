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

//Hll Flux 2d implementation
template<class Grid, class Scalar> class HllFlux :
    public NumericalFlux<Grid,Scalar>
{
    enum
    {   dim=Grid::dimension};
    typedef Dune::FieldVector<Scalar, dim> VelType;
    typedef Dune::FieldVector<Scalar, dim+1> SystemType; //System has 3 equations -> 3 positions in a vector
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

    Scalar gravity_;
    Scalar eps_;

public:

    HllFlux() :
        NumericalFlux<Grid, Scalar>(), gravity_(9.81), eps_(1e-8)
    {

    }

    SystemType operator()(VelType velocityI, VelType velocityJ, Scalar waterDepthI,
            Scalar waterDepthJ, VelType nVec)
    {

        //Hll is an approximate Riemann solver 
        //It needs the states of the variables on the left and right side of the interface. 
        //If there is no reconstruction, these values coincide with the I, J values of teh cell centers
        //L = left side, R = right side

        Scalar waterDepthFaceL = 0;
        Scalar waterDepthFaceR = 0;
        Scalar signalVelocityL = 0;
        Scalar signalVelocityR = 0;
        VelType velocityFaceL(0);
        VelType velocityFaceR(0);
        Scalar primVelocityFaceL = 0;
        Scalar primVelocityFaceR = 0;
        Scalar secVelocityFaceL = 0;
        Scalar secVelocityFaceR = 0;
        Scalar partConti = 0;
        Scalar partContiL = 0;
        Scalar partContiR = 0;
        Scalar fluxConti = 0;
        Scalar partMomentumX = 0;
        Scalar partMomentumY = 0;
        Scalar partMomentumXL = 0;
        Scalar partMomentumYL = 0;
        Scalar partMomentumXR = 0;
        Scalar partMomentumYR = 0;
        Scalar fluxMomentumX = 0;
        Scalar fluxMomentumY = 0;
        Scalar starVelocity = 0;
        Scalar starWaveCelerity = 0;
        Scalar help1 = 0;
        Scalar help2 =0;
        VelType unityVector(1);
        Scalar direction;

        direction = nVec*unityVector;

        //determine left and right side value of a face, 
        //depends on direction of the normal vector
        if (direction > 0)
        {
            waterDepthFaceL = waterDepthI;
            waterDepthFaceR = waterDepthJ;
            velocityFaceL = velocityI;
            velocityFaceR = velocityJ;
        }

        else if (direction < 0)
        {
            waterDepthFaceL = waterDepthJ;
            waterDepthFaceR = waterDepthI;
            velocityFaceL = velocityJ;
            velocityFaceR = velocityI;
        }

        switch (dim)
        {
        case 1:
            primVelocityFaceL = velocityFaceL[0];
            primVelocityFaceR = velocityFaceR[0];
            break;

        case 2:
            if (nVec[0] != 0)
            {
                primVelocityFaceL = velocityFaceL[0]; //dominating velocity (e.g. in x-direction u)
                secVelocityFaceL = velocityFaceL[1]; //secondary velocity (e.g. in y-direction v)
                primVelocityFaceR = velocityFaceR[0];
                secVelocityFaceR = velocityFaceR[1];

            }
            else if (nVec[1]!= 0)
            {
                primVelocityFaceL = velocityFaceL[1]; //dominating velocity (e.g. in x-direction u)
                secVelocityFaceL = velocityFaceL[0]; //secondary velocity (e.g. in y-direction v)
                primVelocityFaceR = velocityFaceR[1];
                secVelocityFaceR = velocityFaceR[0];
            }

            else
            {
                DUNE_THROW(NotImplemented, "Primary and Secondary Velocities");
            }
            break;
        }

        //     std::cout<<"primVelocityFaceL "<<primVelocityFaceL<<std::endl;
        //     std::cout<<"primVelocityFaceR "<<primVelocityFaceR<<std::endl;
        //     std::cout<<"waterDepthFaceR "<<waterDepthFaceR<<std::endl;
        //   std::cout<<"waterDepthFaceL "<<waterDepthFaceL<<std::endl;
        //   std::cout<<"velocityFaceR "<<velocityFaceR<<std::endl;
        // std::cout<<"velocityFaceL "<<velocityFaceL<<std::endl;

        //Calculate HLL Values for the star region for velocity and wavecelerity

        starVelocity=(primVelocityFaceL+primVelocityFaceR);
        starVelocity*=0.5;
        starVelocity += sqrt(gravity_*waterDepthFaceL);
        starVelocity -= sqrt(gravity_*waterDepthFaceR);

        help1 =(primVelocityFaceL-primVelocityFaceR);
        help1 *=0.25;
        help2 = sqrt(gravity_*waterDepthFaceL)+sqrt(gravity_*waterDepthFaceR);
        help2 *=0.5;
        starWaveCelerity=(help1+help2);

        // std::cout<<"starVelocity "<<starVelocity<<std::endl;
        // std::cout<<"starWaveCelerity "<<starWaveCelerity<<std::endl;
        // std::cout<<"waterDepthFaceR "<<waterDepthFaceR<<std::endl;
        // std::cout<<"waterDepthFaceL "<<waterDepthFaceL<<std::endl;


        //Calculate signal velocities used in HLL-flux calculation

        if (waterDepthFaceL == 0 && waterDepthFaceR != 0)
        {
            signalVelocityL = sqrt(gravity_*waterDepthFaceR);
            signalVelocityL *= (-2);
            signalVelocityL += primVelocityFaceR;

            signalVelocityR = primVelocityFaceR;
            signalVelocityR += sqrt(gravity_*waterDepthFaceR);
        }
        else if (waterDepthFaceR == 0 && waterDepthFaceL != 0)
        {
            signalVelocityL = primVelocityFaceL;
            signalVelocityL -= sqrt(gravity_*waterDepthFaceL);

            signalVelocityR = sqrt(gravity_*waterDepthFaceL);
            signalVelocityR *= 2;
            signalVelocityR += primVelocityFaceL;
        }

        else if (waterDepthFaceL != 0 && waterDepthFaceR != 0)
        {
            signalVelocityL = std::min((primVelocityFaceL-sqrt(gravity_
                    *waterDepthFaceL)), (starVelocity-starWaveCelerity));
            signalVelocityR = std::max((primVelocityFaceR+sqrt(gravity_
                    *waterDepthFaceR)), (starVelocity+starWaveCelerity));
        }
        else if (waterDepthFaceL == 0 && waterDepthFaceR == 0)
        {
            signalVelocityL = 0;
            signalVelocityR = 0;
        }
        else
        {
            DUNE_THROW(NotImplemented, "SignalVelocity");
        }

        Scalar denominator = (signalVelocityR-signalVelocityL);

        // std::cout<<"signalVelocityL "<<signalVelocityL<<std::endl;
        // std::cout<<"signalVelocityR "<<signalVelocityR<<std::endl;

        //Compute the flux vector components Conti, XMomentum, YMomentum for the left and right state

        if (signalVelocityL >= 0)
        {
            partConti = primVelocityFaceL*waterDepthFaceL;
            fluxConti = partConti;

            switch (dim)
            {
            case 1:

                partMomentumX = waterDepthFaceL*(primVelocityFaceL
                        *primVelocityFaceL) + 0.5 * gravity_ * (waterDepthFaceL
                        *waterDepthFaceL);
                fluxMomentumX = partMomentumX;
                break;

            case 2:

                if (nVec[0] !=0)
                {
                    partMomentumX = waterDepthFaceL*(primVelocityFaceL
                            *primVelocityFaceL) + 0.5 * gravity_
                            * (waterDepthFaceL *waterDepthFaceL);
                    fluxMomentumX = partMomentumX;

                    partMomentumY = waterDepthFaceL*primVelocityFaceL
                            *secVelocityFaceL;
                    fluxMomentumY = partMomentumY;
                }
                else if (nVec[1] !=0)
                {
                    partMomentumX = waterDepthFaceL*primVelocityFaceL
                            *secVelocityFaceL;
                    fluxMomentumX = partMomentumX;

                    partMomentumY = waterDepthFaceL*(primVelocityFaceL
                            *primVelocityFaceL) + 0.5 * gravity_
                            * (waterDepthFaceL *waterDepthFaceL);
                    fluxMomentumY = partMomentumY;
                }
                break;
            }
        }
        else if (signalVelocityL <= 0 && signalVelocityR >= 0)
        {
            partContiL = primVelocityFaceL*waterDepthFaceL;

            partContiR = primVelocityFaceR*waterDepthFaceR;

            switch (dim)
            {
            case 1:
                partMomentumXL = waterDepthFaceL*(primVelocityFaceL
                        *primVelocityFaceL) +0.5 *gravity_ *(waterDepthFaceL
                        *waterDepthFaceL);

                partMomentumXR =waterDepthFaceR*(primVelocityFaceR
                        *primVelocityFaceR) +0.5 *gravity_ *(waterDepthFaceR
                        *waterDepthFaceR);

                fluxConti = signalVelocityR*(partContiL);
                fluxConti -= signalVelocityL*(partContiR);
                fluxConti += signalVelocityR*signalVelocityL*(waterDepthFaceR
                        -waterDepthFaceL);
                fluxConti /= denominator;

                fluxMomentumX = signalVelocityR*(partMomentumXL);
                fluxMomentumX -= signalVelocityL*(partMomentumXR);
                fluxMomentumX += signalVelocityR*signalVelocityL
                        *(waterDepthFaceR *primVelocityFaceR - waterDepthFaceL
                                *primVelocityFaceL);
                fluxMomentumX /= denominator;
                break;

            case 2:
                if (nVec[0] !=0)
                {
                    //std::cout<<"x-direction"<<std::endl;

                    partMomentumXL = waterDepthFaceL*(primVelocityFaceL
                            *primVelocityFaceL) + 0.5 * gravity_
                            * (waterDepthFaceL *waterDepthFaceL);

                    partMomentumYL = waterDepthFaceL*primVelocityFaceL
                            *secVelocityFaceL;

                    partMomentumXR = waterDepthFaceR*(primVelocityFaceR
                            *primVelocityFaceR) + 0.5 * gravity_
                            * (waterDepthFaceR *waterDepthFaceR);

                    partMomentumYR = waterDepthFaceR*primVelocityFaceR
                            *secVelocityFaceR;

                }
                else if (nVec[1] !=0)
                {
                    //std::cout<<"y-direction"<<std::endl;

                    partMomentumXL = waterDepthFaceL*primVelocityFaceL
                            *secVelocityFaceL;

                    partMomentumYL = waterDepthFaceL*(primVelocityFaceL
                            *primVelocityFaceL) + 0.5 * gravity_
                            * (waterDepthFaceL *waterDepthFaceL);

                    partMomentumXR = waterDepthFaceR*primVelocityFaceR
                            *secVelocityFaceR;

                    partMomentumYR = waterDepthFaceR*(primVelocityFaceR
                            *primVelocityFaceR) + 0.5 * gravity_
                            * (waterDepthFaceR *waterDepthFaceR);
                }

                else
                {
                    DUNE_THROW(NotImplemented, "check normal vector");
                }

                /*     std::cout<<"partMomentumXL "<<partMomentumXL<<std::endl;
                 std::cout<<"partMomentumXR "<<partMomentumXR<<std::endl;
                 std::cout<<"partMomentumYL "<<partMomentumYL<<std::endl;
                 std::cout<<"partMomentumYR "<<partMomentumYR<<std::endl;
                 std::cout<<"velocityFaceR "<<velocityFaceR<<std::endl;
                 std::cout<<"velocityFaceL "<<velocityFaceL<<std::endl;
                 */
                fluxConti = signalVelocityR*(partContiL);
                fluxConti -= signalVelocityL*(partContiR);
                fluxConti += signalVelocityR*signalVelocityL*(waterDepthFaceR
                        -waterDepthFaceL);
                fluxConti /= denominator;

                fluxMomentumX = signalVelocityR*(partMomentumXL);
                fluxMomentumX -= signalVelocityL*(partMomentumXR);
                fluxMomentumX += signalVelocityR*signalVelocityL
                        *(waterDepthFaceR *velocityFaceR[0] - waterDepthFaceL
                                *velocityFaceL[0]);
                fluxMomentumX /= denominator;

                fluxMomentumY = signalVelocityR*(partMomentumYL);
                fluxMomentumY -= signalVelocityL*(partMomentumYR);
                fluxMomentumY += signalVelocityR*signalVelocityL
                        *(waterDepthFaceR *velocityFaceR[1] - waterDepthFaceL
                                *velocityFaceL[1]);
                fluxMomentumY /= denominator;

                // std::cout<<"fluxMomentumY "<<fluxMomentumY<<std::endl;

                break;
            }

        }
        else if (signalVelocityR <= 0)
        {
            partConti = primVelocityFaceR*waterDepthFaceR;
            fluxConti = partConti;

            switch (dim)
            {
            case 1:
                partMomentumX = waterDepthFaceR*(primVelocityFaceR
                        *primVelocityFaceR) +0.5 *gravity_ *(waterDepthFaceR
                        *waterDepthFaceR);
                fluxMomentumX = partMomentumX;
                break;

            case 2:
                if (nVec[0] !=0)
                {
                    partMomentumX = waterDepthFaceR*(primVelocityFaceR
                            *primVelocityFaceR) + 0.5 * gravity_
                            * (waterDepthFaceR *waterDepthFaceR);
                    fluxMomentumX = partMomentumX;

                    partMomentumY = waterDepthFaceR*primVelocityFaceR
                            *secVelocityFaceR;
                    fluxMomentumY = partMomentumY;
                }
                else if (nVec[1] !=0)
                {
                    partMomentumX = waterDepthFaceR*primVelocityFaceR
                            *secVelocityFaceR;
                    fluxMomentumX = partMomentumX;

                    partMomentumY = waterDepthFaceR*(primVelocityFaceR
                            *primVelocityFaceR) + 0.5 * gravity_
                            * (waterDepthFaceR *waterDepthFaceR);
                    fluxMomentumY = partMomentumY;
                }
                break;
            }
        }
        else
        {
            DUNE_THROW(NotImplemented,
                    "Signal Velocity Combination not implemented");
        }

        SystemType fluxHllFace(0);

        switch (dim)
        {
        case 1:
            fluxHllFace[0]=fluxConti*nVec[0];
            fluxHllFace[1]=fluxMomentumX*nVec[0];
            break;

        case 2:
            fluxHllFace[0]=fluxConti*(nVec[0]+nVec[1]);
            fluxHllFace[1]=fluxMomentumX*(nVec[0]+nVec[1]);
            fluxHllFace[2]=fluxMomentumY*(nVec[0]+nVec[1]);

            break;
        }
        // std::cout<<"Hll flux="<<fluxHllFace<<std::endl;
        return fluxHllFace;
    }
};
}
#endif

