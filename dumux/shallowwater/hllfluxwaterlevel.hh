// $Id$

#ifndef DUNE_HLLFLUXWATERLEVEL_HH
#define DUNE_HLLFLUXWATERLEVEL_HH

#include "dumux/shallowwater/shallownumericalflux.hh"

namespace Dune
{

//**********************************************************************************************************
//Hll well-balanced Flux 2d implementation (Liang, Marche)
template<class GridView, class Scalar> class HllFluxWaterLevel :
    public NumericalFlux<GridView,Scalar>
{
    enum
    {   dim=GridView::dimension, oneD = 1, twoD = 2};

    typedef Dune::FieldVector<Scalar, dim> VelType;
    typedef Dune::FieldVector<Scalar, dim+1> SystemType; //System has 3 equations -> 3 positions in a vector
    typedef typename GridView::Traits::template Codim<0>::Entity Element;

    Scalar eps_;

public:

    HllFluxWaterLevel() :
        NumericalFlux<GridView, Scalar>(), eps_(1e-6)
    {

    }

    SystemType operator()(VelType velocityI, VelType velocityJ, Scalar waterDepthI,
            Scalar waterDepthJ, VelType nVec, Scalar gravity,
            Scalar bottomElevationI, Scalar bottomElevationJ, Scalar dist,
            Scalar waterLevelI, Scalar waterLevelJ)
    {

        //Hll is an approximate Riemann solver
        //It needs the states of the variables on the left and right side of the interface.
        //If there is no reconstruction, these values coincide with the I, J values of teh cell centers
        //L = left side, R = right side

        //the well-balanced flux computes the fluxes with the waterlevel not the waterdepth

        Scalar waterDepthFaceL = 0;
        Scalar waterDepthFaceR = 0;
        Scalar waterLevelFaceL = 0;
        Scalar waterLevelFaceR = 0;
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

        Scalar direction = 0;
        direction = nVec*unityVector;

        //calculate the Riemann variables
        Scalar bottomElevationFace = std::max(bottomElevationI,
                bottomElevationJ);

        //determine left and right side value of a face,
        //depends on direction of the normal vector
        if (direction> 0)
        {
            waterDepthFaceL = waterDepthI;
            waterDepthFaceR = waterDepthJ;
            waterLevelFaceL = waterLevelI;
            waterLevelFaceR = waterLevelJ;
            velocityFaceL = velocityI;
            velocityFaceR = velocityJ;
        }

        else if (direction < 0)
        {
            waterDepthFaceL = waterDepthJ;
            waterDepthFaceR = waterDepthI;
            waterLevelFaceL = waterLevelJ;
            waterLevelFaceR = waterLevelI;
            velocityFaceL = velocityJ;
            velocityFaceR = velocityI;
        }
        else
        {
            DUNE_THROW(NotImplemented, "Primary and Secondary Velocities");
        }

        Scalar waterDepthFaceLMod = std::max(0.0, waterLevelFaceL
                - bottomElevationFace);
        Scalar waterDepthFaceRMod = std::max(0.0, waterLevelFaceR
                - bottomElevationFace);

        Scalar waterLevelFaceLMod = waterDepthFaceLMod + bottomElevationFace;
        Scalar waterLevelFaceRMod = waterDepthFaceRMod + bottomElevationFace;

        switch (dim)
        {
        case oneD:
            primVelocityFaceL = velocityFaceL[0];
            primVelocityFaceR = velocityFaceR[0];
            break;

        case twoD:
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
        //     std::cout<<"waterDepthFaceL "<<waterDepthFaceL<<std::endl;
        //     std::cout<<"velocityFaceR "<<velocityFaceR<<std::endl;
        //     std::cout<<"velocityFaceL "<<velocityFaceL<<std::endl;

        //Calculate HLL Values for the star region for velocity and wavecelerity

        starVelocity=(primVelocityFaceL+primVelocityFaceR);
        starVelocity*=0.5;
        starVelocity += sqrt(gravity*waterDepthFaceLMod);
        starVelocity -= sqrt(gravity*waterDepthFaceRMod);

        help1 =(primVelocityFaceL-primVelocityFaceR);
        help1 *=0.25;
        help2 = sqrt(gravity*waterDepthFaceLMod)+sqrt(gravity
                *waterDepthFaceRMod);
        help2 *=0.5;
        starWaveCelerity=(help1+help2);

        // std::cout<<"starVelocity "<<starVelocity<<std::endl;
        // std::cout<<"starWaveCelerity "<<starWaveCelerity<<std::endl;
        // std::cout<<"waterDepthFaceRMod "<<waterDepthFaceRMod<<std::endl;
        // std::cout<<"waterDepthFaceLMod "<<waterDepthFaceLMod<<std::endl;


        //Calculate signal velocities used in HLL-flux calculation

        if (waterDepthFaceLMod <= eps_ && waterDepthFaceRMod> eps_)
        {
            signalVelocityL = sqrt(gravity*waterDepthFaceRMod);
            signalVelocityL *= (-2);
            signalVelocityL += primVelocityFaceR;

            signalVelocityR = primVelocityFaceR;
            signalVelocityR += sqrt(gravity*waterDepthFaceRMod);
        }
        else if (waterDepthFaceRMod <= eps_ && waterDepthFaceLMod> eps_)
        {
            signalVelocityL = primVelocityFaceL;
            signalVelocityL -= sqrt(gravity*waterDepthFaceLMod);

            signalVelocityR = sqrt(gravity*waterDepthFaceLMod);
            signalVelocityR *= 2;
            signalVelocityR += primVelocityFaceL;
        }

        else if (waterDepthFaceLMod> eps_ && waterDepthFaceRMod> eps_)
        {
            signalVelocityL = std::min((primVelocityFaceL-sqrt(gravity
                    *waterDepthFaceLMod)), (starVelocity-starWaveCelerity));
            signalVelocityR = std::max((primVelocityFaceR+sqrt(gravity
                    *waterDepthFaceRMod)), (starVelocity+starWaveCelerity));
        }
        else if (waterDepthFaceLMod <= eps_ && waterDepthFaceRMod <= eps_)
        {
            signalVelocityL = 0;
            signalVelocityR = 0;
        }
        else
        {
            DUNE_THROW(NotImplemented, "SignalVelocity");
        }

        Scalar denominator = (signalVelocityR-signalVelocityL);

        //       std::cout<<"signalVelocityL "<<signalVelocityL<<std::endl;
        //       std::cout<<"signalVelocityR "<<signalVelocityR<<std::endl;

        //Compute the flux vector components Conti, XMomentum, YMomentum for the left and right state

        if (signalVelocityL >= 0)
        {
            partConti = primVelocityFaceL*waterDepthFaceLMod;
            fluxConti = partConti;

            switch (dim)
            {
            case oneD:

                partMomentumX = waterDepthFaceLMod*(primVelocityFaceL
                        *primVelocityFaceL) + 0.5 * gravity
                        * (waterLevelFaceLMod *waterLevelFaceLMod - 2
                                *waterLevelFaceLMod*bottomElevationFace);
                fluxMomentumX = partMomentumX;
                break;

            case twoD:

                if (nVec[0] !=0)
                {
                    partMomentumX = waterDepthFaceLMod*(primVelocityFaceL
                            *primVelocityFaceL) + 0.5 * gravity
                            * (waterLevelFaceLMod *waterLevelFaceLMod - 2
                                    *waterLevelFaceLMod*bottomElevationFace);
                    fluxMomentumX = partMomentumX;

                    partMomentumY = waterDepthFaceLMod*primVelocityFaceL
                            *secVelocityFaceL;
                    fluxMomentumY = partMomentumY;
                }
                else if (nVec[1] !=0)
                {
                    partMomentumX = waterDepthFaceLMod*primVelocityFaceL
                            *secVelocityFaceL;
                    fluxMomentumX = partMomentumX;

                    partMomentumY = waterDepthFaceLMod*(primVelocityFaceL
                            *primVelocityFaceL) + 0.5 * gravity
                            * (waterLevelFaceLMod *waterLevelFaceLMod - 2
                                    *waterLevelFaceLMod*bottomElevationFace);
                    fluxMomentumY = partMomentumY;
                }
                break;
            }
        }
        else if (signalVelocityL <= 0 && signalVelocityR >= 0)
        {
            partContiL = primVelocityFaceL*waterDepthFaceLMod;

            partContiR = primVelocityFaceR*waterDepthFaceRMod;

            switch (dim)
            {
            case oneD:
                partMomentumXL = waterDepthFaceLMod*(primVelocityFaceL
                        *primVelocityFaceL) +0.5 * gravity
                        * (waterLevelFaceLMod *waterLevelFaceLMod - 2
                                *waterLevelFaceLMod*bottomElevationFace);

                partMomentumXR =waterDepthFaceRMod*(primVelocityFaceR
                        *primVelocityFaceR) +0.5 * gravity
                        * (waterLevelFaceRMod *waterLevelFaceRMod - 2
                                *waterLevelFaceRMod*bottomElevationFace);

                fluxConti = signalVelocityR*(partContiL);
                fluxConti -= signalVelocityL*(partContiR);
                fluxConti += signalVelocityR*signalVelocityL
                        *(waterDepthFaceRMod -waterDepthFaceLMod);
                fluxConti /= denominator;

                fluxMomentumX = signalVelocityR*(partMomentumXL);
                fluxMomentumX -= signalVelocityL*(partMomentumXR);
                fluxMomentumX += signalVelocityR*signalVelocityL
                        *(waterDepthFaceRMod *primVelocityFaceR
                                - waterDepthFaceLMod *primVelocityFaceL);
                fluxMomentumX /= denominator;
                break;

            case twoD:
                if (nVec[0] !=0)
                {
                    //std::cout<<"x-direction"<<std::endl;

                    partMomentumXL = waterDepthFaceLMod*(primVelocityFaceL
                            *primVelocityFaceL) + 0.5 * gravity
                            * (waterLevelFaceLMod *waterLevelFaceLMod - 2
                                    *waterLevelFaceLMod*bottomElevationFace);

                    partMomentumYL = waterDepthFaceLMod*primVelocityFaceL
                            *secVelocityFaceL;

                    partMomentumXR = waterDepthFaceRMod*(primVelocityFaceR
                            *primVelocityFaceR) + 0.5 * gravity
                            * (waterLevelFaceRMod *waterLevelFaceRMod - 2
                                    *waterLevelFaceRMod*bottomElevationFace);

                    partMomentumYR = waterDepthFaceRMod*primVelocityFaceR
                            *secVelocityFaceR;

                }
                else if (nVec[1] !=0)
                {
                    //std::cout<<"y-direction"<<std::endl;

                    partMomentumXL = waterDepthFaceLMod*primVelocityFaceL
                            *secVelocityFaceL;

                    partMomentumYL = waterDepthFaceLMod*(primVelocityFaceL
                            *primVelocityFaceL) + 0.5 * gravity
                            * (waterLevelFaceLMod *waterLevelFaceLMod - 2
                                    *waterLevelFaceLMod*bottomElevationFace);

                    partMomentumXR = waterDepthFaceRMod*primVelocityFaceR
                            *secVelocityFaceR;

                    partMomentumYR = waterDepthFaceRMod*(primVelocityFaceR
                            *primVelocityFaceR) + 0.5 * gravity
                            * (waterLevelFaceRMod *waterLevelFaceRMod - 2
                                    *waterLevelFaceRMod*bottomElevationFace);
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
                fluxConti += signalVelocityR*signalVelocityL
                        *(waterDepthFaceRMod -waterDepthFaceLMod);
                fluxConti /= denominator;

                fluxMomentumX = signalVelocityR*(partMomentumXL);
                fluxMomentumX -= signalVelocityL*(partMomentumXR);
                fluxMomentumX += signalVelocityR*signalVelocityL
                        *(waterDepthFaceRMod *velocityFaceR[0]
                                - waterDepthFaceLMod *velocityFaceL[0]);
                fluxMomentumX /= denominator;

                fluxMomentumY = signalVelocityR*(partMomentumYL);
                fluxMomentumY -= signalVelocityL*(partMomentumYR);
                fluxMomentumY += signalVelocityR*signalVelocityL
                        *(waterDepthFaceRMod *velocityFaceR[1]
                                - waterDepthFaceLMod *velocityFaceL[1]);
                fluxMomentumY /= denominator;

                // std::cout<<"fluxMomentumY "<<fluxMomentumY<<std::endl;

                break;
            }

        }
        else if (signalVelocityR <= 0)
        {
            partConti = primVelocityFaceR*waterDepthFaceRMod;
            fluxConti = partConti;

            switch (dim)
            {
            case oneD:
                partMomentumX = waterDepthFaceRMod*(primVelocityFaceR
                        *primVelocityFaceR) +0.5 * gravity
                        * (waterLevelFaceRMod *waterLevelFaceRMod - 2
                                *waterLevelFaceRMod*bottomElevationFace);
                fluxMomentumX = partMomentumX;
                break;

            case twoD:
                if (nVec[0] !=0)
                {
                    partMomentumX = waterDepthFaceRMod*(primVelocityFaceR
                            *primVelocityFaceR) + 0.5 * gravity
                            * (waterLevelFaceRMod *waterLevelFaceRMod - 2
                                    *waterLevelFaceRMod*bottomElevationFace);
                    fluxMomentumX = partMomentumX;

                    partMomentumY = waterDepthFaceRMod*primVelocityFaceR
                            *secVelocityFaceR;
                    fluxMomentumY = partMomentumY;
                }
                else if (nVec[1] !=0)
                {
                    partMomentumX = waterDepthFaceRMod*primVelocityFaceR
                            *secVelocityFaceR;
                    fluxMomentumX = partMomentumX;

                    partMomentumY = waterDepthFaceRMod*(primVelocityFaceR
                            *primVelocityFaceR) + 0.5 * gravity
                            * (waterLevelFaceRMod *waterLevelFaceRMod - 2
                                    *waterLevelFaceRMod*bottomElevationFace);
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
        case oneD:
            fluxHllFace[0]=fluxConti*nVec[0];
            fluxHllFace[1]=fluxMomentumX*nVec[0];
            break;

        case twoD:
            fluxHllFace[0]=fluxConti*(nVec[0]+nVec[1]);
            fluxHllFace[1]=fluxMomentumX*(nVec[0]+nVec[1]);
            fluxHllFace[2]=fluxMomentumY*(nVec[0]+nVec[1]);

            break;
        }
        //    std::cout<<"Hll flux="<<fluxHllFace<<std::endl;

        // Compute Source Term
       
        VelType bedSlopeTerm = 0;

        SystemType totalFlux = fluxHllFace - bedSlopeTerm;

        return totalFlux;
    }
};

}
#endif

