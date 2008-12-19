// $Id$

#ifndef DUNE_CAPILLARYDIFFUSION_HH
#define DUNE_CAPILLARYDIFFUSION_HH

#include "dumux/transport/fv/diffusivepart.hh"
#include "dumux/diffusion/diffusionproblem.hh"

//! \ingroup transport
//! \defgroup diffPart Diffusive transport
/**
 * @file
 * @brief  Base class for defining the diffusive part of an advection-diffusion equation
 * @author Bernd Flemisch, last changed by Markus Wolff
 */
namespace Dune
{
  /*!\ingroup diffPart
   * @brief  Base class for defining the diffusive part of an advection-diffusion equation
   */
  template<class Grid, class Scalar, class VC>
  class CapillaryDiffusion : public DiffusivePart<Grid,Scalar>
  {
    enum{dim = Grid::dimension,dimWorld = Grid::dimensionworld};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;
    typedef typename Grid::LevelGridView GridView;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef Dune::FieldVector<Scalar, dim> FieldVector;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    typedef VC::ScalarVectorType SatType;

  public:
    virtual FieldVector operator() (const Element& element, const int numberInSelf,
				    const Scalar satIntersection, const FieldVector& satGradient, const Scalar time,
				    const Scalar satI, const Scalar satJ) const
    {
      // cell geometry type
      GeometryType gt = element.geometry().type();

      // cell center in reference element
      const LocalPosition& localPos = ReferenceElements<Scalar,dim>::general(gt).position(0,0);

      // get global coordinate of cell center
      const GlobalPosition& globalPos = element.geometry().global(localPos);

      // get absolute permeability of cell
      FieldMatrix K(problem_.K(globalPos,element,localPos));

      IntersectionIterator isItEnd = element.ilevelend();
      IntersectionIterator isIt = element.ilevelbegin();
      for (; isIt != isItEnd; ++isIt)
	{
	  if(is->numberInSelf() == numberInSelf)
	    break;
	}

      // get geometry type of face
      GeometryType faceGT = isIt->intersectionSelfLocal().type();

      // center in face's reference element
      const Dune::FieldVector<Scalar,dim-1>& faceLocal = ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

      FieldVector unitOuterNormal = isIt->unitOuterNormal(faceLocal);
      //std::cout<<"unitOuterNormaldiff"<<unitOuterNormal<<std::endl;

      if (isIt->neighbor()) {
	// access neighbor
	ElementPointer neighborPointer = isIt->outside();

	// compute factor in neighbor
	GeometryType neighborGT = neighborPointer->geometry().type();
	const LocalPosition& local = ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

	// neighbor cell center in global coordinates
	const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

	// take arithmetic average of absolute permeability
	K += problem_.K(globalPosNeighbor, *neighborPointer, localPosNeighbor);
	K *= 0.5;
      }

      // set result to grad(S)
      FieldVector helpResult(satGradient);

      //get capillary pressure gradients
      Scalar dPdSI=constRel_.dPdS(satI);
      Scalar dPdSJ=constRel_.dPdS(satJ);

      // set result to (dp_c/dS)*grad(S)
      helpResult *= (dPdSI+dPdSJ)*0.5;


      // add gravitational effects (rho_w - rho_n)*g
      helpResult += gravity;

      // set result to K*((dp_c/dS)*grad(S) + (rho_w - rho_n)*g)
      FieldVector result(0);
      K.umv(helpResult, result);

      //get lambda_bar = lambda_n*f_w
      Scalar mobBarI=constRel_.mobN(1-satI)*constRel.fractionalW(satI);
      Scalar mobBarJ=constRel_.mobN(1-satJ)*constRel.fractionalW(satJ);

      // set result to f_w*lambda_n*K*((dp_c/dS)*grad(S) + (rho_w - rho_n)*g)
      result *= (mobBarI+mobBarJ)*0.5;

      return result;
    }

    CapillaryDiffusion (DiffusionProblem<Grid, Scalar, VC>& problem)
      : problem_(problem), constRel_(problem.materialLaw), wettingPhase_(constRel_.wettingPhase),
	nonwettingPhase_(constRel_.nonwettingPhase)
    {
      Scalar rhoDiff = wettingPhase_.density() - nonwettingPhase_.density();
      gravity_ = problem_.gravity();
      gravity_ *= rhoDiff;
    }

  private:
    DiffusionProblem<Grid, Scalar, VC>& problem_;
    TwoPhaseRelations& constRel_;
    const Medium& wettingPhase_;
    const Medium& nonwettingPhase_;
    FieldVector gravity_;
  };
}

#endif
