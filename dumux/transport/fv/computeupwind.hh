// $Id$

#ifndef DUNE_COMPUTEUPWIND_HH
#define DUNE_COMPUTEUPWIND_HH

#include "dumux/transport/fv/computenumflux.hh"
#include "dumux/transport/transportproblem.hh"

//! \ingroup transport
//! \defgroup numerical flux transport
/**
 * @file
 * @brief  Base class for defining the numerical flux of an advection-diffusion equation
 * @author Yufei Cao
 */

namespace Dune
{
/*!\ingroup transport
 * @brief  Base class for defining the numerical flux of an advection-diffusion equation
 */
template<class GridView, class Scalar, class VC, class Problem = TransportProblem<GridView, Scalar, VC> >
class ComputeUpwind : public ComputeNumFlux<GridView,Scalar> {
    enum{dim = GridView::dimension,dimWorld = GridView::dimensionworld};
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::IntersectionIterator IntersectionIterator;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * Upwind numerical flux.
     */
    virtual Scalar operator() (const Element& element, const int faceNumber, const Scalar satI, const Scalar satJ) const
    {
        // cell geometry type
        GeometryType gt = element.geometry().type();

        // cell center in reference element
        const LocalPosition& localPos = ReferenceElements<Scalar,dim>::general(gt).position(0,0);

        // get global coordinate of cell center
        const GlobalPosition& globalPos = element.geometry().global(localPos);

        // compute fractional flow function
        Scalar fI = problem_.materialLaw().fractionalW(satI, globalPos, element,localPos);

        IntersectionIterator isItEnd = element.ilevelend();
        IntersectionIterator isIt = element.ilevelbegin();
        for (; isIt != isItEnd; ++isIt)
        {
            if(isIt->indexInInside() == faceNumber)
                break;
        }

        // get geometry type of face
        Dune::GeometryType faceGT = isIt->geometryInInside().type();

        // center in face's reference element
        const Dune::FieldVector<Scalar,dim-1>&
        faceLocal = Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

        // get normal vector scaled with volume
        Dune::FieldVector<Scalar,dimWorld> integrationOuterNormal = isIt->integrationOuterNormal(faceLocal);
        integrationOuterNormal *= Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).volume();

        // get the flux through the faceIJ
        Scalar velocityIJ = std::max(problem_.variables().vTotalElementFace(element, faceNumber)*integrationOuterNormal, 0.0);
        Scalar velocityJI = std::max(-(problem_.variables().vTotalElementFace(element, faceNumber)*integrationOuterNormal), 0.0);

        Scalar fJ;

        if (isIt->neighbor()) {
            // access neighbor
            ElementPointer neighborPointer = isIt->outside();

            // compute factor in neighbor
            GeometryType neighborGT = neighborPointer->geometry().type();
            const LocalPosition& localPosNeighbor = ReferenceElements<Scalar,dim>::general(neighborGT).position(0,0);

            // neighbor cell center in global coordinates
            const GlobalPosition& globalPosNeighbor = neighborPointer->geometry().global(localPosNeighbor);

            // compute fractional flow function
            fJ = problem_.materialLaw().fractionalW(satJ, globalPosNeighbor, *neighborPointer, localPosNeighbor);
        }
        else
        {
            // get geometry type of face
            Dune::GeometryType faceGT = isIt->geometryInInside().type();

            // center in face's reference element
            const Dune::FieldVector<Scalar,dim-1>&
            faceLocal = Dune::ReferenceElements<Scalar,dim-1>::general(faceGT).position(0,0);

            // center of face inside volume reference element
            const LocalPosition&
            localPosFace = Dune::ReferenceElements<Scalar,dim>::general(faceGT).position(faceNumber,1);

            // center of face in global coordinates
            GlobalPosition globalPosFace = isIt->geometry().global(faceLocal);

            // compute fractional flow function
            fJ = problem_.materialLaw().fractionalW(satJ, globalPosFace, element, localPosFace);
        }

        return (fI * velocityIJ - fJ * velocityJI);
    }

    ComputeUpwind (Problem& problem)
        : problem_(problem)
    { }

private:
    Problem& problem_;
};
}

#endif
