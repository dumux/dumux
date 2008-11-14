// $Id$

#ifndef DUNE_BOX2P2CJACOBIAN_HH
#define DUNE_BOX2P2CJACOBIAN_HH

#include<map>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<sstream>

#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/common/geometrytype.hh>
#include<dune/grid/common/quadraturerules.hh>
#include<dune/grid/utility/intersectiongetter.hh>
#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include"dumux/operators/boxjacobian.hh"
#include"dumux/2p2c/2p2cproblem.hh"
//#include "varswitch.hh"
#include"dumux/io/vtkmultiwriter.hh"

/*!
 * \file
 * \brief  compute local jacobian matrix for box scheme for two-phase two-component flow equation
 * \author Bernd Flemisch, Klaus Mosthaf, Andreas Lauser
 */



namespace Dune
{
    /*! 
     * \addtogroup DISC_Disc
     *
     * @{
     */
    /*!
     * \brief compute local jacobian matrix for the boxfile for two-phase two-component flow equation
     *
     */
    
    
    //! Derived class for computing local jacobian matrices
    /*! A class for computing local jacobian matrix for the two-phase two-component flow equation
      
      div j = q; j = -K grad u; in Omega
      
      u = g on Gamma1; j*n = J on Gamma2.
      
      Uses box scheme with the Lagrange shape functions.  It should
      work for all dimensions and element types.  All the numbering is
      with respect to the reference element and the Lagrange shape
      functions

      Template parameters are:
      
      - ProblemT    The problem specific stuff
    */
    template <class ProblemT>
    class _TwoPTwoCBoxJacobianImp
    {
        typedef ProblemT                                Problem;
        typedef typename Problem::DomainTraits          DomTraits;
        typedef typename Problem::BoxTraits             BoxTraits;
        
        typedef typename DomTraits::Scalar             Scalar;
        typedef typename DomTraits::CoordScalar        CoordScalar;
        typedef typename DomTraits::Grid               Grid;
        typedef typename DomTraits::Cell               Cell;
        typedef typename DomTraits::LocalCoord         LocalCoord;
        typedef typename DomTraits::WorldCoord         WorldCoord;
        
        enum {
            GridDim     = DomTraits::GridDim,
            WorldDim    = DomTraits::WorldDim,
            
            PrimaryVariables   = BoxTraits::PrimaryVariables,
            NumPhases          = BoxTraits::NumPhases,
            NumComponents      = BoxTraits::NumComponents,

            // Phase indices
            WettingIndex    = BoxTraits::WettingIndex,
            NonwettingIndex = BoxTraits::NonettingIndex,
            
            // Phase state
            WettingComponent    = BoxTraits::WettingComponent,
            NonwettingComponent = BoxTraits::NonwettingComponent,
            BothComponents      = BoxTraits::BothComponents,
        };

        typedef typename BoxTraits::UnknownsVector       UnknownsVector;
        typedef typename BoxTraits::FVElementGeometry    FVElementGeometry;
        typedef typename BoxTraits::CachedCellData       CachedCellData;
        typedef typename BoxTraits::CachedSubContVolData CachedSubContVolData;
        typedef typename BoxTraits::LocalFunction        LocalFunction;

#if 0
        typedef typename G::ctype DT;
        typedef typename G::Traits::template Codim<0>::Entity Entity;
        typedef typename Entity::Geometry Geometry;
        typedef Box2P2CJacobian<G,Scalar,BoxFunction> ThisType;
        typedef typename LocalJacobian<ThisType,G,Scalar,2>::VBlockType VBlockType;
        typedef typename LocalJacobian<ThisType,G,Scalar,2>::MBlockType MBlockType;

 	enum {pWIdx = 0, switchIdx = 1, numberOfComponents = 2};	// Solution vector index
	enum {wPhase = 0, nPhase = 1};					// Phase index
	enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};		// Phase state
	enum {water = 0, air = 1};                                      // Component index


    public:
        // define the number of phases (m) and components (c) of your system, this is used outside
        // to allocate the correct size of (dense) blocks with a FieldMatrix
        enum {dim=G::dimension};
        enum {m=2, c=2};
        enum {SIZE=LagrangeShapeFunctionSetContainer<DT,Scalar,dim>::maxsize};
        struct VariableNodeData;

        typedef FieldMatrix<Scalar,dim,dim> FMatrix;
        typedef FieldVector<Scalar,dim> FVector;

        //! Constructor
        Box2P2CJacobian (TwoPTwoCProblem<G,Scalar>& params, bool levelBoundaryAsDirichlet_, const G& grid,
                         BoxFunction& sol, bool procBoundaryAsDirichlet_=true)
            : BoxJacobian<ThisType,G,Scalar,2,BoxFunction>(levelBoundaryAsDirichlet_, grid, sol, procBoundaryAsDirichlet_),
              problem(params), sNDat(this->vertexMapper.size()), vNDat(SIZE), oldVNDat(SIZE)
            {
                this->analytic = false;
            }
#endif // 0

	/*! 
         * \brief compute time dependent term (storage), loop over nodes / subcontrol volumes
	 * \param e entity
	 * \param sol solution vector
	 * \param node local node id
	 * \return storage term
	 */
//        virtual VBlockType computeM (const Entity& e, const VBlockType* sol,
//                                     int node, std::vector<VariableNodeData>& varData)

        static void evalLocalRate(UnknownsVector &result,
                                  const Problem &problem,
                                  const Cell &cell,
                                  const FVElementGeometry &dualCell,
                                  const LocalFunction &localSol,
                                  const CachedCellData &cachedData,
                                  const int scvId)
            {
                int globalIdx = problem.vertexIndex(cell, scvId);
                
                UnknownsVector result;
                Scalar satN = cachedData[scvId].satN;
                Scalar satW = cachedData[scvId].satW;
                
                Scalar densityW = cachedData[scvId].density[WettingIndex];
                Scalar densityN = cachedData[scvId].density[NonwettingIndex];

                // porosity within the FV cell
                Scalar porosity = problem.porosity(cell, scvId, globalIdx);

                // Mass fractions of a component in a phase
                Scalar massfracWinW = cachedData[scvId].massfrac[WettingIndex][WettingIndex];
                Scalar massfracNinW = cachedData[scvId].massfrac[WettingIndex][NonwettingIndex];
                Scalar massfracWinN = cachedData[scvId].massfrac[NonwettingIndex][WettingIndex];
                Scalar massfracNinN = cachedData[scvId].massfrac[NonwettingIndex][NonwettingIndex];

                // storage of component water
                result[WettingIndex] =
                    porosity*(densityW*satW*massfracWinW + densityN*satN*massfracNinW);
                // storage of component air
                result[NonwettingIndex] =
                    porosity*(densityN*satN*massfracNinN + densityW*satW*massfracNinW);
                
                //std::cout << result << " " << scvId << std::endl;
//                return result;
            };

#if 0
        virtual UnknownsVector computeM (const Entity& e, const UnknownsVector* sol, int node, bool old = false)
            {
                if (old)
                    return computeM(e, sol, node, oldVNDat);
                else
                    return computeM(e, sol, node, vNDat);
            }
#endif

#if 0
        /*! \brief compute diffusive/advective fluxes, loop over subcontrol volume faces
	 *  \param e entity
	 *  \param sol solution vector
	 *  \param face face id
	 *  \return flux term
         */
        virtual UnknownsVector computeA (const Entity& e, const UnknownsVector* sol, int face)
#endif
        /*!
         * \brief Evaluates the mass flux over a face of a subcontrol
         *        volume.
         */
        static void evalFluxRate(UnknownsVector &flux,
                                 const Problem &problem,
                                 const Cell& cell,
                                 const FVElementGeometry &fvGeom,
                                 const LocalFunction &localSol,
                                 const CachedCellData &cachedData,
                                 const int faceId)
            {
                int i = this->fvGeom.subContVolFace[face].i;
                int j = this->fvGeom.subContVolFace[face].j;

                // normal vector, value of the area of the scvf
                const WorldCoord &normal(this->fvGeom.subContVolFace[face].normal);

                // get global coordinates of nodes i,j
                const WorldCoord &global_i = this->fvGeom.subContVol[i].global;
                const WorldCoord &global_j = this->fvGeom.subContVol[j].global;
                const WorldCoord &local_i = this->fvGeom.subContVol[i].local;
                const WorldCoord &local_j = this->fvGeom.subContVol[j].local;

                WorldCoord pGrad[NumPhases], xGrad[NumPhases];
                for (int i = 0; i < NumPhases; ++i) {
                    pGrad[i] = Scalar(0.0);
                    xGrad[i] = Scalar(0.0);
                }
                
                WorldCoord temp(0.);
                UnknownsVector flux(0.);
                FMatrix Ki(0), Kj(0);

                //	FieldVector<Scalar,GridDim> Kij(0);
                // effective permeability in edge direction
                // Scalar Kij = sIPDat[global_j].K_eff[face];
                //const FMatrix K;
                //harmonicMeanK(K, e, face);
                //K.umv(normal, Kij);  // Kij=K*n

                // calculate harmonic mean of permeabilities of nodes i and j
                Ki = this->problem.soil().K(global_i,e,local_i);
                Kj = this->problem.soil().K(global_j,e,local_j);
                const FMatrix  = harmonicMeanK(K, Ki, Kj);

                // calculate FE gradient (grad p for each phase)
                for (int k = 0; k < this->fvGeom.nNodes; k++) // loop over adjacent nodes
                {
                    // FEGradient at node k
                    const WorldCoord &feGrad = this->fvGeom.subContVolFace[face].grad[k];
                    
                    UnknownsVector pressure(0.0), massfrac(0.0);

                    pressure[WettingIndex] = cachedData[k].pW;
                    pressure[NonwettingIndex] = cachedData[k].pN;

                    // compute sum of pressure gradients for each phase
                    for (int phase = 0; phase < m; phase++)
                    {
                        temp = feGrad;
                        temp *= pressure[phase];
                        pGrad[phase] += temp;
                    }
                    // for diffusion of air in wetting phase
                    temp = feGrad;
                    temp *= vNDat[k].massfrac[air][WettingIndex];
                    xGrad[WettingIndex] += temp;

                    // for diffusion of water in nonwetting phase
                    temp = feGrad;
                    temp *= vNDat[k].massfrac[WettingIndex][NonwettingIndex];
                    xGrad[NonwettingIndex] += temp;
                }

                // deduce gravity*density of each phase
                WorldCoord contribComp[NumPhases];
                for (int phase=0; phase < NumPhases; phase++)
                {
                    contribComp[phase]  = problem.gravity();
                    contribComp[phase] *= vNDat[i].density[phase];
                    pGrad[phase] -= contribComp[phase]; // grad p - rho*g
                }

                UnknownsVector outward(0);  // Darcy velocity of each phase

                // calculate the advective flux using upwind: K*n(grad p -rho*g)
                for (int phase=0; phase < NumPhases; phase++)
	 	{
                    WorldCoord v_tilde(0);
                    K.mv(pGrad[phase], v_tilde);  // v_tilde=K*gradP
                    outward[phase] = v_tilde*normal;
	 	}

                // evaluate upwind nodes
                int up_w, dn_w, up_n, dn_n;
                if (outward[WettingIndex] <= 0) {
                    up_w = i; 
                    dn_w = j;
                }
                else {
                    up_w = j;
                    dn_w = i;
                };
                
                if (outward[NonwettingIndex] <= 0) {
                    up_n = i;
                    dn_n = j;
                }
                else {
                    up_n = j;
                    dn_n = i;
                };


                Scalar alpha = 1.0;  // Upwind parameter

                // Water conservation
                flux[WettingIndex] =   (alpha* cachedData[up_w].density[WettingIndex]*cachedData[up_w].mobility[WettingIndex]
                                 * cachedData[up_w].massfrac[WettingIndex][WettingIndex]
                                 + (1-alpha)* cachedData[dn_w].density[WettingIndex]*cachedData[dn_w].mobility[WettingIndex]
                                 * cachedData[dn_w].massfrac[WettingIndex][WettingIndex])
                    * outward[WettingIndex];
                flux[WettingIndex] +=  (alpha* cachedData[up_n].density[NonwettingIndex]*cachedData[up_n].mobility[NonwettingIndex]
                                 * cachedData[up_n].massfrac[WettingIndex][NonwettingIndex]
                                 + (1-alpha)* cachedData[dn_n].density[NonwettingIndex]*cachedData[dn_n].mobility[NonwettingIndex]
                                 * cachedData[dn_n].massfrac[WettingIndex][NonwettingIndex])
                    * outward[NonwettingIndex];
                // Air conservation
                flux[air]   =   (alpha* cachedData[up_n].density[NonwettingIndex]*cachedData[up_n].mobility[NonwettingIndex]
                                 * cachedData[up_n].massfrac[air][NonwettingIndex]
                                 + (1-alpha)* cachedData[dn_n].density[NonwettingIndex]*cachedData[dn_n].mobility[NonwettingIndex]
                                 * cachedData[dn_n].massfrac[air][NonwettingIndex])
                    * outward[NonwettingIndex];
                flux[air]  +=   (alpha* cachedData[up_w].density[WettingIndex]*cachedData[up_w].mobility[WettingIndex]
                                 * cachedData[up_w].massfrac[air][WettingIndex]
                                 + (1-alpha)* cachedData[dn_w].density[WettingIndex]*cachedData[dn_w].mobility[WettingIndex]
                                 * cachedData[dn_w].massfrac[air][WettingIndex])
                    * outward[WettingIndex];

                // DIFFUSION
                UnknownsVector normDiffGrad;

                // get local to global id map
                int state_i = cachedData[i].phasestate;
                int state_j = cachedData[j].phasestate;

                Scalar diffusionWW(0.0), diffusionWN(0.0); // diffusion of water
                Scalar diffusionAW(0.0), diffusionAN(0.0); // diffusion of air
                UnknownsVector avgDensity, avgDpm;

                avgDpm[WettingIndex]=2e-9; // needs to be changed !!!
                avgDpm[NonwettingIndex]=2.25e-5; // water in the gasphase

                normDiffGrad[WettingIndex] = xGrad[WettingIndex]*normal;
                normDiffGrad[NonwettingIndex] = xGrad[NonwettingIndex]*normal;

                // calculate the arithmetic mean of densities
                avgDensity[WettingIndex] = 0.5*(cachedData[i].density[WettingIndex] + cachedData[j].density[WettingIndex]);
                avgDensity[NonwettingIndex] = 0.5*(cachedData[i].density[NonwettingIndex] + cachedData[j].density[NonwettingIndex]);

                if (state_i==2 && state_j==2)
                {
                    diffusionAW = avgDpm[WettingIndex] * avgDensity[WettingIndex] * normDiffGrad[WettingIndex];
                    diffusionWW = - diffusionAW;
                    diffusionWN = avgDpm[NonwettingIndex] * avgDensity[NonwettingIndex] * normDiffGrad[NonwettingIndex];
                    diffusionAN = - diffusionWN;
                }
                else if ((state_i == 1 || state_j == 1) || (state_i == 1 && state_j == 1))
                {
                    diffusionAW = avgDpm[WettingIndex] * avgDensity[WettingIndex] * normDiffGrad[WettingIndex];
                    diffusionWW = - diffusionAW;
                }
                else if ((state_i == 0 || state_j == 0) || (state_i == 0 && state_j == 0))
                {
                    diffusionWN = avgDpm[NonwettingIndex] * avgDensity[NonwettingIndex] * normDiffGrad[NonwettingIndex];
                    diffusionAN = - diffusionWN;
                }

                // add diffusion of water to flux
                flux[WettingIndex] += (diffusionWW + diffusionWN);
                //	std::cout << "Water Flux: " << flux[WettingIndex] << std::endl;

                // add diffusion of air to flux
                flux[air] += (diffusionAN + diffusionAW);
                // std::cout << "Air Flux: " << flux[air] << std::endl;


                return flux;
            };

  	/*! \brief integrate sources / sinks
  	 *  \param e entity
	 *  \param sol solution vector
	 *  \param node local node id
	 *  \return source/sink term
	 */
   	virtual UnknownsVector computeQ (const Entity& e, const UnknownsVector* sol, const int node)
            {
   		// ASSUME problem.q already contains \rho.q
   		return problem.q(this->fvGeom.subContVol[node].global, e, this->fvGeom.subContVol[node].local);
            }

#if 0 // -> move this stuff to the problem or to the model
  	/*! \brief perform variable switch
  	 *  \param global global node id
	 *  \param sol solution vector
	 *  \param local local node id
	 */
   	virtual void primaryVarSwitch (const Entity& e, int globalIdx, UnknownsVector* sol, int localIdx)
            {

                bool switched = false;
   		int state = sNDat[globalIdx].phaseState;
//   		Scalar eps = 1e-15;

                Scalar pW = sol[localIdx][pWIdx];
//                if (pW < 0.) pW = (-1)*pW;

                Scalar satW = 0.0;
                if ( state == BothComponents) 
                    satW = 1.0-sol[localIdx][switchIdx];
  		else if (state == WettingComponent)
                    satW = 1.0;
  		else if (state == NonwettingComponent)
                    satW = 0.0;

  		const FVector global = this->fvGeom.subContVol[localIdx].global;
  		const FVector local = this->fvGeom.subContVol[localIdx].local;

                Scalar pC = problem.materialLaw().pC(satW, global, e, local);
  		Scalar pN = pW + pC;

                switch(state)
                {
                    case BothComponents:
                        Scalar xWNmass, xWNmolar, pwn, pWSat;
                        xWNmass = sol[localIdx][switchIdx];
                        xWNmolar = problem.multicomp().convertMassToMoleFraction(xWNmass, gasPhase);
                        pwn = xWNmolar * pN;
                        pWSat = problem.multicomp().vaporPressure(cachedData[localIdx].temperature);

                        if (pwn > 1.01*pWSat && switched == false)
                        {
                            // appearance of water phase
                            std::cout << "Water appears at node " << globalIdx << "  Coordinates: " << global << std::endl;
                            sNDat[globalIdx].phaseState = bothPhases;
                            sol[localIdx][switchIdx] = 1.0 - 1e-6; // initialize solution vector
                            switched = true;
                        }
                        break;

                    case WettingComponent:
                        Scalar pbub, xAWmass, xAWmolar, henryInv;
                        xAWmass = sol[localIdx][switchIdx];
                        xAWmolar = problem.multicomp().convertMassToMoleFraction(xAWmass, waterPhase);
                        henryInv = problem.multicomp().henry(cachedData[localIdx].temperature);
                        pWSat = problem.multicomp().vaporPressure(cachedData[localIdx].temperature);
                        pbub = pWSat + xAWmolar/henryInv;

                        if (pbub > pN && switched == false && pN > 0.0)
                        {
                            // appearance of gas phase
                            std::cout << "Gas appears at node " << globalIdx << "  Coordinates: " << global << std::endl;
                            sNDat[globalIdx].phaseState = bothPhases;
                            sol[localIdx][switchIdx] = 1e-6; // initialize solution vector
                            switched = true;
                        }
                        break;

                    case NonwettingComponent:
                        Scalar satN = sol[localIdx][switchIdx];

                        if (satN < 0.0  && switched == false)
                        {
                            // disappearance of gas phase
                            std::cout << "Gas disappears at node " << globalIdx << "  Coordinates: " << global << std::endl;
                            sNDat[globalIdx].phaseState = waterPhase;
                            sol[localIdx][switchIdx] = 1e-6; // initialize solution vector
                            switched = true;
                        }
                        else if (satW < 0.0  && switched == false)
                        {
                            // disappearance of water phase
                            std::cout << "Water disappears at node " << globalIdx << "  Coordinates: " << global << std::endl;
                            sNDat[globalIdx].phaseState = gasPhase;
                            sol[localIdx][switchIdx] = 1e-6; // initialize solution vector
                            switched = true;
                        }
                        break;

                }
                if (switched == true) updateVariableData(e, sol, localIdx, cachedData, sNDat[globalIdx].phaseState);

                return;
            }

        // harmonic mean computed directly
        static void harmonicMeanK(FMatrix &dest, FMatrix& Ki, const FMatrix& Kj) const
            {
                double eps = 1e-20;

                for (int kx=0; kx<GridDim; kx++){
                    for (int ky=0; ky<GridDim; ky++){
    			if (Ki[kx][ky] != Kj[kx][ky])
    			{
                            Ki[kx][ky] = 2 / (1/(Ki[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));
    			}
                    }
                }
                return Ki;
            }

        virtual void clearVisited ()
            {
                for (int i = 0; i < this->vertexMapper.size(); i++){
                    sNDat[i].visited = false;
//   	 	sNDat[i].switched = false;
                }
                return;
            }

        // updates old phase state after each time step
        virtual void updatePhaseState ()
            {
                for (int i = 0; i < this->vertexMapper.size(); i++){
                    sNDat[i].oldPhaseState = sNDat[i].phaseState;
                }
                return;
            }

        //*********************************************************
        //*														*
        //*	Calculation of Data at Elements (elData) 			*
        //*						 								*
        //*														*
        //*********************************************************

        virtual void computeElementData (const Entity& e)
            {
//  	 // ASSUME element-wise constant parameters for the material law
// 		 elData.parameters = problem.materialLawParameters
// 		 (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
//
//		 // ASSUMING element-wise constant permeability, evaluate K at the cell center
// 		 elData.K = problem.K(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
//
//		 // ASSUMING element-wise constant porosity
// 		 elData.porosity = problem.porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);
                return;
            }


        //*********************************************************
        //*														*
        //*	Calculation of Data at Nodes that has to be			*
        //*	determined only once	(sNDat)						*
        //*														*
        //*********************************************************

        // analog to EvalStaticData in MUFTE
        virtual void updateStaticData (const Entity& e, UnknownsVector* sol)
            {
                // size of the sNDat vector is determined in the constructor

                // get access to shape functions for P1 elements
                GeometryType gt = e.geometry().type();
                const typename LagrangeShapeFunctionSetContainer<DT,Scalar,GridDim>::value_type&
                    sfs=LagrangeShapeFunctions<DT,Scalar,GridDim>::general(gt,1);

                // get local to global id map
                for (int k = 0; k < sfs.size(); k++) {
                    const int globalIdx = this->vertexMapper.template map<GridDim>(e, sfs[k].entity());

                    // if nodes are not already visited
                    if (!sNDat[globalIdx].visited)
                    {
                        // ASSUME porosity defined at elements!
                        sNDat[globalIdx].porosity = problem.soil().porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

                        // global coordinates
                        FieldVector<DT,GridDim> global_i = this->fvGeom.subContVol[k].global;

                        // evaluate primary variable switch
                        primaryVarSwitch(e, globalIdx, sol, k);

                        // mark elements that were already visited
                        sNDat[globalIdx].visited = true;
                    }
                }

                return;
            }

        // for initialization of the Static Data (sets porosity)
        virtual void initiateStaticData (const Entity& e)
            {
                // get access to shape functions for P1 elements
                GeometryType gt = e.geometry().type();
                const typename LagrangeShapeFunctionSetContainer<DT,Scalar,GridDim>::value_type&
                    sfs=LagrangeShapeFunctions<DT,Scalar,GridDim>::general(gt,1);

                // get local to global id map
                for (int k = 0; k < sfs.size(); k++) {
                    const int globalIdx = this->vertexMapper.template map<GridDim>(e, sfs[k].entity());

                    // if nodes are not already visited
                    if (!sNDat[globalIdx].visited)
                    {
                        // ASSUME porosity defined at nodes
                        sNDat[globalIdx].porosity = problem.soil().porosity(this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal);

                        // mark elements that were already visited
                        sNDat[globalIdx].visited = true;
//  			for (int h = 0; h<1600; h++)
//  			{
//  				(*outPermeability)[h] = Knew[h];
//  			}
                    }
                }

                return;
            }

#endif // move to problem
        
        // analog to EvalPrimaryData in MUFTE, uses members of cachedData
        void updateVariableData(const Entity& e,
                                const UnknownsVector* sol,
                                int i,
                                std::vector<VariableNodeData>& varData, 
                                int state)
            {
                const int globalIdx = this->vertexMapper.template map<GridDim>(e, i);
                FVector& global = this->fvGeom.subContVol[i].global;
                FVector& local = this->fvGeom.subContVol[i].local;
                
                varData[i].pW = sol[i][pWIdx];
                if (state == bothPhases) varData[i].satN = sol[i][switchIdx];
                if (state == waterPhase) varData[i].satN = 0.0;
                if (state == gasPhase) varData[i].satN = 1.0;

                varData[i].satW = 1.0 - varData[i].satN;

                varData[i].pC = problem.materialLaw().pC(varData[i].satW, global, e, local);
                varData[i].pN = varData[i].pW + varData[i].pC;
                varData[i].temperature = 283.15; // in [K], constant

                // Solubilities of components in phases
                if (state == bothPhases){
                    varData[i].massfrac[air][WettingIndex] = problem.multicomp().xAW(varData[i].pN, varData[i].temperature);
                    varData[i].massfrac[WettingIndex][NonwettingIndex] = problem.multicomp().xWN(varData[i].pN, varData[i].temperature);
                }
                if (state == waterPhase){
                    varData[i].massfrac[WettingIndex][NonwettingIndex] = 0.0;
                    varData[i].massfrac[air][WettingIndex] =  sol[i][switchIdx];
                }
                if (state == gasPhase){
                    varData[i].massfrac[WettingIndex][NonwettingIndex] = sol[i][switchIdx];
                    varData[i].massfrac[air][WettingIndex] = 0.0;
                }
                varData[i].massfrac[WettingIndex][WettingIndex] = 1.0 - varData[i].massfrac[air][WettingIndex];
                varData[i].massfrac[air][NonwettingIndex] = 1.0 - varData[i].massfrac[WettingIndex][NonwettingIndex];
                varData[i].phasestate = state;

                // Mobilities & densities
                varData[i].mobility[WettingIndex] = problem.materialLaw().mobW(varData[i].satW, global, e, local, varData[i].temperature, varData[i].pW);
                varData[i].mobility[NonwettingIndex] = problem.materialLaw().mobN(varData[i].satN, global, e, local, varData[i].temperature, varData[i].pN);
                // Density of Water is set constant here!
                varData[i].density[WettingIndex] = 1000;//problem.wettingPhase().density(varData[i].temperature, varData[i].pN);
                varData[i].density[NonwettingIndex] = problem.nonwettingPhase().density(varData[i].temperature, varData[i].pN,
                                                                               varData[i].massfrac[WettingIndex][NonwettingIndex]);

                // CONSTANT solubility (for comparison with twophase)
//         varData[i].massfrac[air][WettingIndex] = 0.0; varData[i].massfrac[WettingIndex][WettingIndex] = 1.0;
//         varData[i].massfrac[WettingIndex][NonwettingIndex] = 0.0; varData[i].massfrac[air][NonwettingIndex] = 1.0;

                //std::cout << "water in gasphase: " << varData[i].massfrac[WettingIndex][NonwettingIndex] << std::endl;
                //std::cout << "air in waterphase: " << varData[i].massfrac[air][WettingIndex] << std::endl;

#if 0 // TODO: move to model!
                // for output
                (*outPressureN)[globalIdx] = varData[i].pN;
                (*outCapillaryP)[globalIdx] = varData[i].pC;
                (*outSaturationW)[globalIdx] = varData[i].satW;
                (*outSaturationN)[globalIdx] = varData[i].satN;
                (*outMassFracAir)[globalIdx] = varData[i].massfrac[air][WettingIndex];
                (*outMassFracWater)[globalIdx] = varData[i].massfrac[WettingIndex][NonwettingIndex];
                (*outDensityW)[globalIdx] = varData[i].density[WettingIndex];
                (*outDensityN)[globalIdx] = varData[i].density[NonwettingIndex];
                (*outMobilityW)[globalIdx] = varData[i].mobility[WettingIndex];
                (*outMobilityN)[globalIdx] = varData[i].mobility[NonwettingIndex];
                (*outPhaseState)[globalIdx] = varData[i].phasestate;
#endif 
            }

#if 0
	virtual void updateVariableData(const Entity& e, const UnknownsVector* sol, int i, bool old = false)
            {
		int state;
		const int global = this->vertexMapper.template map<GridDim>(e, i);
		if (old)
		{
                    state = sNDat[global].oldPhaseState;
                    updateVariableData(e, sol, i, oldVNDat, state);
		}
		else
		{
		    state = sNDat[global].phaseState;
                    updateVariableData(e, sol, i, cachedData, state);
		}
            }

	void updateVariableData(const Entity& e, const UnknownsVector* sol, bool old = false)
            {
		int size = this->fvGeom.nNodes;

		for (int i = 0; i < size; i++)
                    updateVariableData(e, sol, i, old);
            }

#endif

#if 0 // TODO: Move to problem or model
        struct StaticNodeData
        {
            bool visited;
//   	 bool switched;
            int phaseState;
            int oldPhaseState;
            Scalar cellVolume;
            Scalar porosity;
            FMatrix K;
        };

        struct ElementData {
//   	 Scalar cellVolume;
//     	 Scalar porosity;
//   	 Scalar gravity;
//   	 FieldVector<Scalar, 4> parameters;
//   	 FieldMatrix<Scalar,GridDim,GridDim> K;
        } elData;


        // parameters given in constructor
        TwoPTwoCProblem<G,Scalar>& problem;
        CWaterAir multicomp;
        std::vector<StaticNodeData> sNDat;
        std::vector<VariableNodeData> vNDat;
        std::vector<VariableNodeData> oldVNDat;

        // for output files
        BlockVector<FieldVector<Scalar, 1> > *outPressureN;
        BlockVector<FieldVector<Scalar, 1> > *outCapillaryP;
        BlockVector<FieldVector<Scalar, 1> > *outSaturationN;
        BlockVector<FieldVector<Scalar, 1> > *outSaturationW;
        BlockVector<FieldVector<Scalar, 1> > *outMassFracAir;
        BlockVector<FieldVector<Scalar, 1> > *outMassFracWater;
        BlockVector<FieldVector<Scalar, 1> > *outDensityW;
        BlockVector<FieldVector<Scalar, 1> > *outDensityN;
        BlockVector<FieldVector<Scalar, 1> > *outMobilityW;
        BlockVector<FieldVector<Scalar, 1> > *outMobilityN;
        BlockVector<FieldVector<Scalar, 1> > *outPhaseState;
	BlockVector<FieldVector<Scalar, 1> > *outPermeability;
#endif
    };

}
#endif

//    // average permeability from the staticNode Vector
//    void harmonicMeanK(FMatrix &Kdest, const Entity& e, int k) const
//    {
//	 FMatrix Kj;
//	 const Scalar eps = 1e-20;
//
//    GeometryType gt = e.geometry().type();
//    const typename LagrangeShapeFunctionSetContainer<DT,Scalar,GridDim>::value_type&
//    sfs=LagrangeShapeFunctions<DT,Scalar,GridDim>::general(gt,1);
//
//     int i = this->fvGeom.subContVolFace[k].i;
//     int j = this->fvGeom.subContVolFace[k].j;
//
//     int global_i = this->vertexMapper.template map<GridDim>(e, sfs[i].entity());
//     int global_j = this->vertexMapper.template map<GridDim>(e, sfs[j].entity());
//
//
//     	Kdest = sNDat[global_i].K;
//     	Kj = sNDat[global_j].K;
//
//     	for (int kx=0; kx<GridDim; kx++){
//     		for (int ky=0; ky<GridDim; ky++){
//     			if (Kdest[kx][ky] != Kj[kx][ky])
//     			{
//     				Kdest[kx][ky] = 2 / (1/(Kdest[kx][ky]+eps) + (1/(Kj[kx][ky]+eps)));
//     			}
//     		}
//     	}
//    }

