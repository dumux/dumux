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
                                  Scalar wDepthFaceR, VelType nVecScaled, VelType nVec)=0;

    virtual ~NumericalFlux()
    {
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
    VelType avgVel;
    VelType faceVel;
    VelType faceWDepth;
    VelType eps;
    VelType gravity;
    Scalar  fluxConti;
    VelType fluxMomentum;

public:

    FirstOrderUpwind() :
        NumericalFlux<Grid, Scalar>()
    {

    }

    SystemType operator()(VelType velI, VelType velJ, Scalar wDepthI, Scalar wDepthJ,
                          VelType nVecScaled, VelType nVec)
    {

        avgVel = velI*nVecScaled;
        avgVel +=velJ*nVecScaled;
        avgVel*= 0.5;

        if (avgVel >= eps)
        {
            faceVel = velI;
            faceWDepth = wDepthI;
        }
        if (avgVel < eps)
        {
            faceVel = velJ;
            faceWDepth = wDepthJ;
        }

        fluxConti = faceVel * faceWDepth;
        fluxMomentum[0] = faceVel * faceWDepth*faceVel;
        fluxMomentum[0] += gravity*faceWDepth*0.5;
        fluxMomentum[0] *= nVecScaled[0];
        
        fluxMomentum[1] = faceVel * faceWDepth*faceVel;
        fluxMomentum[1] += gravity*faceWDepth*0.5;
        fluxMomentum[1] *= nVecScaled[1];

    

        SystemType fluxUpwindFace(0);
        fluxUpwindFace[0]=fluxConti;
        fluxUpwindFace[1]=fluxMomentum;

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
    typedef Dune::FieldVector<Scalar, dim+1> SystemType; //System has 3 euqations -> 3 positions in a vector
    typedef typename Grid::Traits::template Codim<0>::Element Element;

    Scalar gravity;
    Scalar eps;

public:

    HllFlux() :
        NumericalFlux<Grid, Scalar>(), gravity(9.81), eps(1e-8)
    {

    }

    Scalar sigVelJ;
    Scalar sigVelI;

    SystemType operator()(VelType velFaceL, VelType velFaceR, Scalar wDepthFaceL,
                          Scalar wDepthFaceR, VelType nVecScaled, VelType nVec)
    {

        if (wDepthFaceL < eps && wDepthFaceR < eps)
        {
            SystemType fluxHllFace(0);
            return fluxHllFace;
        }

        else
        {

            //Calculate HLL Average Values for velocity and waterdepth

            Scalar velAvg=(velFaceR+velFaceL)*nVec;
            velAvg*=0.5;
            velAvg += sqrt(gravity*wDepthFaceR);
            velAvg -= sqrt(gravity*wDepthFaceL);

            //std::cout<<"velAvg = "<<velAvg<<std::endl;

            Scalar help1 = (velFaceR+velFaceL)*nVec;
            help1 *=0.25;
            Scalar help2 = sqrt(gravity*wDepthFaceR)+sqrt(gravity*wDepthFaceL);
            help2 *=0.5;
            Scalar gravWDepthAvg=(help1+help2);

            //Calculate signal velocities used in HLL-flux calculation


            if (wDepthFaceL < eps)
            {
                sigVelI=velFaceR*nVec;
                sigVelI-=sqrt(gravity*wDepthFaceR);

                sigVelJ = sqrt(gravity*wDepthFaceR);
                sigVelJ*=2;
                sigVelJ=velFaceR*nVec;

            }
            if (wDepthFaceR < eps)
            {
                sigVelI=sqrt(gravity*wDepthFaceL);
                sigVelI *=(-2);
                sigVelI =velFaceL*nVec;

                sigVelJ = velFaceL*nVec;
                sigVelJ+=sqrt(gravity*wDepthFaceL);
            }

            else
            {

                sigVelI = std::min((velFaceR*nVec-sqrt(gravity *wDepthFaceR)),
                                   (velAvg-gravWDepthAvg));
                sigVelJ = std::max((velFaceL*nVec+sqrt(gravity *wDepthFaceL)),
                                   (velAvg+gravWDepthAvg));

            }

            /*std::cout<<"velFaceL = "<<velFaceL<<" velFaceR = "<<velFaceR<<std::endl;
              std::cout<<"sigVelI = "<<sigVelI<<" sigVelJ = "<<sigVelJ<<std::endl;
              std::cout<<"velAvg = "<<velAvg<<std::endl;
            */

            //Compute the flux vector components Conti, XMomentum, YMomentum for the left and right state
            //Initialize vectors and fill with flux data

            if (sigVelI == sigVelJ)
            {
                SystemType fluxHllFace(0);
                return fluxHllFace;
            }

            else
            {
                VelType fluxContiI = velFaceL; //I means the left side of face, J, teh right side of face
                fluxContiI *=wDepthFaceL;

                VelType fluxContiJ = velFaceR;
                fluxContiJ*=wDepthFaceR;

                VelType fluxXMomentumI(0);
                fluxXMomentumI[0] = wDepthFaceL*pow(velFaceL[0], 2)+0.5*gravity
                    *pow(wDepthFaceL, 2);
                fluxXMomentumI[1] = wDepthFaceL*velFaceL[0]*velFaceL[1];

                VelType fluxXMomentumJ(0);
                fluxXMomentumJ[0] =wDepthFaceR*pow(velFaceR[0], 2)+0.5*gravity
                    *pow(wDepthFaceR, 2);
                fluxXMomentumJ[1] =wDepthFaceR*velFaceR[0]*velFaceR[1];

                VelType fluxYMomentumI(0);
                fluxYMomentumI[0] = wDepthFaceL*velFaceL[0]*velFaceL[1];
                fluxYMomentumI[1] = wDepthFaceL*pow(velFaceL[1], 2)+0.5*gravity
                    *pow(wDepthFaceL, 2);

                VelType fluxYMomentumJ(0);
                fluxYMomentumJ[0] =wDepthFaceR*velFaceR[0]*velFaceR[1];
                fluxYMomentumJ[1] =wDepthFaceR*pow(velFaceR[1], 2)+0.5*gravity
                    *pow(wDepthFaceR, 2);

                //std::cout<<"Y_Mom_I ="<<fluxYMomentumI<<" Y_Mom_J = "<<fluxYMomentumJ<<std::endl;

                //Compute HLL-Flux Components Conti, XMomentum, YMomentum

                Scalar fluxConti=sigVelJ*(fluxContiI*nVecScaled);
                fluxConti -= sigVelI*(fluxContiJ*nVecScaled);
                fluxConti += (sigVelJ*sigVelI*(wDepthFaceR-wDepthFaceL));
                fluxConti/= (sigVelJ-sigVelI);

                Scalar fluxXMomentum= (sigVelJ*(fluxXMomentumI*nVecScaled));
                fluxXMomentum-= sigVelI*(fluxXMomentumJ*nVecScaled);
                fluxXMomentum+= (sigVelJ*sigVelI*(wDepthFaceR*velFaceR[0]
                                                  -wDepthFaceL *velFaceL[0]));
                fluxXMomentum/=(sigVelJ-sigVelI);

                Scalar fluxYMomentum(0);

                if (velFaceL[1] < eps && velFaceR[1]< eps)
                {
                    fluxYMomentum = 0;
                }
                else
                {
                    fluxYMomentum = sigVelJ*(fluxYMomentumI*nVecScaled);
                    fluxYMomentum-= sigVelI*(fluxYMomentumJ*nVecScaled);
                    fluxYMomentum+= (sigVelJ*sigVelI*(wDepthFaceR*velFaceR[1]
                                                      -wDepthFaceL *velFaceL[1]));
                    fluxYMomentum/=(sigVelJ-sigVelI);
                }
                //Define Vector containing the components and return it to shallowfvtransport

                SystemType fluxHllFace(0);
                fluxHllFace[0]=fluxConti;
                fluxHllFace[1]=fluxXMomentum;
                fluxHllFace[2]=fluxYMomentum;

                //std::cout<<"Hll flux="<<fluxHllFace<<std::endl;

                return fluxHllFace;
            }

        }

    }

};
}
#endif

