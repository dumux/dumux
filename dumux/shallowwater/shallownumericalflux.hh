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
/*
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
 */
//****************************************************************************************+

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
        NumericalFlux<Grid, Scalar>(), gravity(9.81), eps(1e-6)
    {

    }

    SystemType operator()(VelType velocityI, VelType velocityJ, Scalar waterDepthI,
            Scalar waterDepthJ, VelType nVec)
    {
        Scalar signalVelocityL = 0;
        Scalar signalVelocityR = 0;
        Scalar waterDepthFaceL = 0;
        Scalar waterDepthFaceR = 0;
        Scalar velocityFaceR = 0;
        Scalar velocityFaceL = 0;
        Scalar fluxConti = 0;
        Scalar fluxMomentum = 0;
        Scalar partConti = 0;
        Scalar partMomentum = 0;
        Scalar partContiL = 0;
        Scalar partMomentumL =0;
        Scalar partContiR = 0;
        Scalar partMomentumR = 0;
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

        else if (direction < 0)
        {
            waterDepthFaceL = waterDepthJ;
            waterDepthFaceR = waterDepthI;
            velocityFaceL = velocityJ;
            velocityFaceR = velocityI;
        }

        //     std::cout<<"waterDepthFaceR "<<waterDepthFaceR<<std::endl;
        //   std::cout<<"waterDepthFaceL "<<waterDepthFaceL<<std::endl;
        //   std::cout<<"velocityFaceR "<<velocityFaceR<<std::endl;
        // std::cout<<"velocityFaceL "<<velocityFaceL<<std::endl;

        //Calculate HLL Values for the star region for velocity and wavecelerity

        starVelocity=(velocityFaceL+velocityFaceR);
        starVelocity*=0.5;
        starVelocity += sqrt(gravity*waterDepthFaceL);
        starVelocity -= sqrt(gravity*waterDepthFaceR);

        help1 =(velocityFaceL-velocityFaceR);
        help1 *=0.25;
        help2 = sqrt(gravity*waterDepthFaceL)+sqrt(gravity*waterDepthFaceR);
        help2 *=0.5;
        starWaveCelerity=(help1+help2);

        // std::cout<<"starVelocity "<<starVelocity<<std::endl;
        // std::cout<<"starWaveCelerity "<<starWaveCelerity<<std::endl;
        //  std::cout<<"waterDepthFaceR "<<waterDepthFaceR<<std::endl;
        // std::cout<<"waterDepthFaceL "<<waterDepthFaceL<<std::endl;


        //Calculate signal velocities used in HLL-flux calculation

        if (waterDepthFaceL == 0 && waterDepthFaceR != 0)
        {
            signalVelocityL = sqrt(gravity*waterDepthFaceR);
            signalVelocityL *= (-2);
            signalVelocityL += velocityFaceR;

            signalVelocityR = velocityFaceR;
            signalVelocityR += sqrt(gravity*waterDepthFaceR);
        }
        else if (waterDepthFaceR == 0 && waterDepthFaceL != 0)
        {
            signalVelocityL = velocityFaceL;
            signalVelocityL -= sqrt(gravity*waterDepthFaceL);

            signalVelocityR = sqrt(gravity*waterDepthFaceL);
            signalVelocityR *= 2;
            signalVelocityR += velocityFaceL;
        }

        else if (waterDepthFaceL != 0 && waterDepthFaceR != 0)
        {
            signalVelocityL = std::min((velocityFaceL-sqrt(gravity
                    *waterDepthFaceL)), (starVelocity-starWaveCelerity));
            signalVelocityR = std::max((velocityFaceR+sqrt(gravity
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
        std::cout<<"signalVelocityL "<<signalVelocityL<<std::endl;
        std::cout<<"signalVelocityR "<<signalVelocityR<<std::endl;

        //Compute the flux vector components Conti, XMomentum, YMomentum for the left and right state


        if (signalVelocityL >= 0)
        {
            partConti = velocityFaceL*waterDepthFaceL;
            fluxConti = partConti;

            partMomentum = waterDepthFaceL*(velocityFaceL*velocityFaceL) + 0.5
                    * gravity * (waterDepthFaceL*waterDepthFaceL);
            fluxMomentum = partMomentum;
        }
        else if (signalVelocityL <= 0 && signalVelocityR >= 0)
        {
            partContiL = velocityFaceL*waterDepthFaceL;

            partContiR = velocityFaceR*waterDepthFaceR;

            partMomentumL = waterDepthFaceL*(velocityFaceL*velocityFaceL) +0.5
                    *gravity *(waterDepthFaceL*waterDepthFaceL);

            partMomentumR =waterDepthFaceR*(velocityFaceR*velocityFaceR)+0.5
                    *gravity *(waterDepthFaceR*waterDepthFaceR);

            Scalar denominator = (signalVelocityR-signalVelocityL);

            fluxConti = signalVelocityR*(partContiL);
            fluxConti -= signalVelocityL*(partContiR);
            fluxConti += signalVelocityR*signalVelocityL*(waterDepthFaceR
                    -waterDepthFaceL);
            fluxConti /= denominator;

            fluxMomentum = signalVelocityR*(partMomentumL);
            fluxMomentum -= signalVelocityL*(partMomentumR);
            fluxMomentum += signalVelocityR*signalVelocityL*(waterDepthFaceR
                    *velocityFaceR - waterDepthFaceL *velocityFaceL);
            fluxMomentum /= denominator;
        }
        else if (signalVelocityR <= 0)
        {
            partConti = velocityFaceR *waterDepthFaceR;
            fluxConti = partConti;

            partMomentum = waterDepthFaceR*(velocityFaceR*velocityFaceR)+0.5
                    *gravity *(waterDepthFaceR*waterDepthFaceR);
            fluxMomentum = partMomentum;

            std::cout<<"waterDepthFaceR "<<waterDepthFaceR<<std::endl;
            std::cout<<"velocityFaceR "<<velocityFaceR<<std::endl;

        }
        else
        {
            DUNE_THROW(NotImplemented, "Flux Range");
        }

        //   std::cout<<"partContiL "<<partContiL<<std::endl;
        //   std::cout<<"partContiR "<<partContiR<<std::endl;
        std::cout<<"fluxConti "<<fluxConti<<std::endl;
        std::cout<<"fluxMomentum "<<fluxMomentum<<std::endl;
        //   std::cout<<"partMomentumL "<<partMomentumL<<std::endl;
        //   std::cout<<"partMomentumR "<<partMomentumR<<std::endl;


        SystemType fluxHllFace(0);
        fluxHllFace[0]=fluxConti*nVec;
        fluxHllFace[1]=fluxMomentum*nVec;

        std::cout<<"Hll flux="<<fluxHllFace<<std::endl;

        return fluxHllFace;
    }
};
//*************************************************************************************************

//Hll Flux 2d implementation
template<class Grid, class Scalar> class HllFlux2d :
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

    HllFlux2d() :
        NumericalFlux<Grid, Scalar>(), gravity(9.81), eps(1e-8)
    {

    }

    SystemType operator()(VelType velocityI, VelType velocityJ, Scalar waterDepthI,
            Scalar waterDepthJ, VelType nVec)
    {
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

        //determine left and right side value of a face
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
        starVelocity += sqrt(gravity*waterDepthFaceL);
        starVelocity -= sqrt(gravity*waterDepthFaceR);

        help1 =(primVelocityFaceL-primVelocityFaceR);
        help1 *=0.25;
        help2 = sqrt(gravity*waterDepthFaceL)+sqrt(gravity*waterDepthFaceR);
        help2 *=0.5;
        starWaveCelerity=(help1+help2);

        // std::cout<<"starVelocity "<<starVelocity<<std::endl;
        // std::cout<<"starWaveCelerity "<<starWaveCelerity<<std::endl;
        // std::cout<<"waterDepthFaceR "<<waterDepthFaceR<<std::endl;
        // std::cout<<"waterDepthFaceL "<<waterDepthFaceL<<std::endl;


        //Calculate signal velocities used in HLL-flux calculation

        if (waterDepthFaceL == 0 && waterDepthFaceR != 0)
        {
            signalVelocityL = sqrt(gravity*waterDepthFaceR);
            signalVelocityL *= (-2);
            signalVelocityL += primVelocityFaceR;

            signalVelocityR = primVelocityFaceR;
            signalVelocityR += sqrt(gravity*waterDepthFaceR);
        }
        else if (waterDepthFaceR == 0 && waterDepthFaceL != 0)
        {
            signalVelocityL = primVelocityFaceL;
            signalVelocityL -= sqrt(gravity*waterDepthFaceL);

            signalVelocityR = sqrt(gravity*waterDepthFaceL);
            signalVelocityR *= 2;
            signalVelocityR += primVelocityFaceL;
        }

        else if (waterDepthFaceL != 0 && waterDepthFaceR != 0)
        {
            signalVelocityL = std::min((primVelocityFaceL-sqrt(gravity
                    *waterDepthFaceL)), (starVelocity-starWaveCelerity));
            signalVelocityR = std::max((primVelocityFaceR+sqrt(gravity
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
                        *primVelocityFaceL) + 0.5 * gravity * (waterDepthFaceL
                        *waterDepthFaceL);
                fluxMomentumX = partMomentumX;
                break;

            case 2:

                if (nVec[0] !=0)
                {
                    partMomentumX = waterDepthFaceL*(primVelocityFaceL
                            *primVelocityFaceL) + 0.5 * gravity
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
                            *primVelocityFaceL) + 0.5 * gravity
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
                        *primVelocityFaceL) +0.5 *gravity *(waterDepthFaceL
                        *waterDepthFaceL);

                partMomentumXR =waterDepthFaceR*(primVelocityFaceR
                        *primVelocityFaceR) +0.5 *gravity *(waterDepthFaceR
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
                    partMomentumXL = waterDepthFaceL*(primVelocityFaceL
                            *primVelocityFaceL) + 0.5 * gravity
                            * (waterDepthFaceL *waterDepthFaceL);

                    partMomentumYL = waterDepthFaceL*primVelocityFaceL
                            *secVelocityFaceL;

                    partMomentumXR = waterDepthFaceR*(primVelocityFaceR
                            *primVelocityFaceR) + 0.5 * gravity
                            * (waterDepthFaceR *waterDepthFaceR);

                    partMomentumYR = waterDepthFaceR*primVelocityFaceR
                            *secVelocityFaceR;

                }
                else if (nVec[1] !=0)
                {
                    partMomentumXL = waterDepthFaceL*primVelocityFaceL
                            *secVelocityFaceL;

                    partMomentumYL = waterDepthFaceL*(primVelocityFaceL
                            *primVelocityFaceL) + 0.5 * gravity
                            * (waterDepthFaceL *waterDepthFaceL);

                    partMomentumXR = waterDepthFaceR*primVelocityFaceR
                            *secVelocityFaceR;

                    partMomentumYR = waterDepthFaceR*(primVelocityFaceR
                            *primVelocityFaceR) + 0.5 * gravity
                            * (waterDepthFaceR *waterDepthFaceR);
                }

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
                        *primVelocityFaceR) +0.5 *gravity *(waterDepthFaceR
                        *waterDepthFaceR);
                fluxMomentumX = partMomentumX;
                break;

            case 2:
                if (nVec[0] !=0)
                {
                    partMomentumX = waterDepthFaceR*(primVelocityFaceR
                            *primVelocityFaceR) + 0.5 * gravity
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
                            *primVelocityFaceR) + 0.5 * gravity
                            * (waterDepthFaceR *waterDepthFaceR);
                    fluxMomentumY = partMomentumY;
                }
                break;
            }
        }
        else
        {
            DUNE_THROW(NotImplemented, "Flux Range");
        }

        //   std::cout<<"partContiL "<<partContiL<<std::endl;
        //   std::cout<<"partContiR "<<partContiR<<std::endl;
        //   std::cout<<"fluxConti "<<fluxConti<<std::endl;
        //   std::cout<<"fluxMomentum "<<fluxMomentum<<std::endl;
        //   std::cout<<"partMomentumL "<<partMomentumL<<std::endl;
        //   std::cout<<"partMomentumR "<<partMomentumR<<std::endl;


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

        //std::cout<<"Hll flux="<<fluxHllFace<<std::endl;
        return fluxHllFace;
    }
};
}
#endif

