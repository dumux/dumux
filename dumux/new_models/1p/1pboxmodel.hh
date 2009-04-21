/*****************************************************************************
 *   Copyright (C) 2007 by Peter Bastian                                     *
 *   Institute of Parallel and Distributed System                            *
 *   Department Simulation of Large Systems                                  *
 *   University of Stuttgart, Germany                                        *
 *                                                                           *
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUMUX_1P_BOX_MODEL_HH
#define DUMUX_1P_BOX_MODEL_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>

#include <dumux/auxiliary/math.hh>

namespace Dune
{
///////////////////////////////////////////////////////////////////////////
// traits for the single phase isothermal model
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief The single phase specific traits.
 */
class OnePTraits
{
public:
    enum {
        numEq     = 1,        //!< Number of primary variables
        numPhases = 1,    //!< Number of fluid phases
    };
    enum {
        pressureIdx = 0  //!< Index for the fluid pressure in a field vector
    };
};

/*!
 * \brief Data which is attached to each vert of the and can
 *        be shared between multiple calculations and should
 *        thus be cached in order to increase efficency.
 */
template <class OnePTraits, 
          class Problem>
class OnePVertexData
{
    typedef OnePTraits Tr;
    typedef typename Problem::DomainTraits::Scalar Scalar;
    typedef typename Problem::DomainTraits::Grid Grid;

    typedef typename Grid::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<Scalar, Tr::numEq>      SolutionVector;
    typedef Dune::FieldVector<Scalar, Tr::numPhases>  PhasesVector;

    typedef Dune::FieldVector<Scalar, Grid::dimensionworld>  GlobalPosition;
    typedef Dune::FieldVector<Scalar, Grid::dimension>       LocalPosition;

public:
    /*!
     * \brief Update all quantities for a given control volume.
     */
    template <class JacobianImp>
    void update(const SolutionVector   &sol,
                const Element          &element,
                int                     vertIdx,
                bool                    isOldSol,
                JacobianImp            &jac) 
    {
        const GlobalPosition &global = element.geometry().corner(vertIdx);
        const LocalPosition   &local =
            Problem::DomainTraits::referenceElement(element.type()).position(vertIdx,
                                                                             Grid::dimension);

        pressure = sol[Tr::pressureIdx];
        density = jac.problem().fluid().density(jac.temperature(sol),
                                                pressure);
        viscosity = jac.problem().fluid().viscosity(jac.temperature(sol),
                                                    pressure);
        // porosity
        porosity = jac.problem().soil().porosity(global,
                                                 element,
                                                 local);
    };

    Scalar pressure;
    Scalar density;
    Scalar viscosity;
    Scalar porosity;
};


///////////////////////////////////////////////////////////////////////////
// OnePBoxJacobian (evaluate the local jacobian for the newton method.)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief Local Jacobian for the single phase isothermal model
 */
template<class ProblemT, class BoxTraitsT, class OnePTraitsT>
class OnePBoxJacobian : public BoxJacobian<ProblemT,
                                           BoxTraitsT,
                                           OnePBoxJacobian<ProblemT,
                                                           BoxTraitsT,
                                                           OnePTraitsT>,
                                           OnePVertexData<OnePTraitsT, ProblemT> >
{
private:
    typedef OnePBoxJacobian<ProblemT, BoxTraitsT, OnePTraitsT>  ThisType;
    typedef OnePVertexData<OnePTraitsT, ProblemT>               VertexData;
    typedef BoxJacobian<ProblemT, 
                        BoxTraitsT,
                        ThisType,
                        VertexData>                 ParentType;

    typedef ProblemT                                Problem;
    typedef typename Problem::DomainTraits          DomTraits;
    typedef BoxTraitsT                              BoxTraits;
    typedef OnePTraitsT                             OnePTraits;

    enum {
        dim            = DomTraits::dim,
        dimWorld       = DomTraits::dimWorld,

        numEq          = BoxTraits::numEq,

        numPhases      = OnePTraits::numPhases,
        pressureIdx    = OnePTraits::pressureIdx,
    };

    typedef typename DomTraits::Scalar              Scalar;
    typedef typename DomTraits::CoordScalar         CoordScalar;
    typedef typename DomTraits::Grid                Grid;
    typedef typename DomTraits::Element             Element;
    typedef typename DomTraits::ElementIterator     ElementIterator;
    typedef typename Element::EntityPointer         ElementPointer;
    typedef typename DomTraits::LocalPosition       LocalPosition;
    typedef typename DomTraits::GlobalPosition      GlobalPosition;

    typedef typename BoxTraits::FVElementGeometry   FVElementGeometry;
    typedef typename BoxTraits::SpatialFunction     SpatialFunction;
    typedef typename BoxTraits::LocalFunction       LocalFunction;
    typedef typename BoxTraits::SolutionVector      SolutionVector;


    typedef std::vector<VertexData>        VertexDataArray;
    typedef FieldMatrix<Scalar, dim, dim>  Tensor;

public:
    OnePBoxJacobian(ProblemT &problem)
        : ParentType(problem)
    {};

    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub control
     *        volume of a finite volume element for the OneP
     *        model.
     *
     * This function should not include the source and sink terms.
     */
    void computeStorage(SolutionVector &result, int scvIdx, bool usePrevSol) const
    {
        // if flag usePrevSol is set, the solution from the previous
        // time step is used, otherwise the current solution is
        // used. The secondary variables are used accordingly.  This
        // is required to compute the derivative of the storage term
        // using the implicit euler method.
        const VertexDataArray &elemDat = usePrevSol ? this->prevElemDat_  : this->curElemDat_;
        const VertexData  &vertDat = elemDat[scvIdx];
        
        // partial time derivative of the wetting phase mass
        result[pressureIdx] =  vertDat.density * vertDat.porosity;
    }


    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     */
    void computeFlux(SolutionVector &flux, int faceId) const
    {
        const typename FVElementGeometry::SubControlVolumeFace
            &face = ParentType::curElementGeom_.subContVolFace[faceId];

        const int i = face.i;
        const int j = face.j;

        // normal vector, value of the area of the scvf
        const GlobalPosition &normal(this->curElementGeom_.subContVolFace[faceId].normal);

        // get global coordinates of verts i,j
        const GlobalPosition &global_i = this->curElementGeom_.subContVol[i].global;
        const GlobalPosition &global_j = this->curElementGeom_.subContVol[j].global;

        // get local coordinates of verts i,j
        const LocalPosition &local_i = this->curElementGeom_.subContVol[i].local;
        const LocalPosition &local_j = this->curElementGeom_.subContVol[j].local;

        // calculate FE gradient
        Scalar densityIJ = 0;
        Scalar viscosityIJ = 0;
        LocalPosition pGrad(0);
        for (int k = 0; k < ParentType::curElementGeom_.numVertices; k++) {
            LocalPosition grad(face.grad[k]);
            grad *= this->curElemDat_[k].pressure;
            pGrad += grad;

            densityIJ += this->curElemDat_[k].density*face.shapeValue[k];
            viscosityIJ += this->curElemDat_[k].viscosity*face.shapeValue[k];
        }

        // adjust pressure gradient by gravity force
        LocalPosition gravity = ParentType::problem_.gravity();
        gravity *= densityIJ;
        pGrad   -= gravity;

        // calculate darcy velocity
        const Tensor &Ki = this->problem_.soil().K(global_i, ParentType::curElement_(), local_i);
        const Tensor &Kj = this->problem_.soil().K(global_j, ParentType::curElement_(), local_j);
        Tensor K;
        harmonicMeanMatrix(K, Ki, Kj);

        // temporary vector for the Darcy velocity
        GlobalPosition vDarcy;

        K.mv(pGrad, vDarcy);  // vDarcy = K * grad p
        vDarcy /= viscosityIJ;

        flux[pressureIdx] = densityIJ * (vDarcy * normal);
    }

    /*!
     * \brief Calculate the source term of the equation
     */
    void computeSource(SolutionVector &q, int localVertexIdx)
    {
        this->problem_.source(q,
                              this->curElement_(),
                              this->curElementGeom_,
                              localVertexIdx);
    }

    /*!
     * \brief Return the temperature given the solution vector of a
     *        finite volume.
     */
    template <class SolutionVector>
    Scalar temperature(const SolutionVector &sol)
    { return this->problem_.temperature(); /* constant temperature */ }

    /*!
     * \brief All relevant primary and secondary of a given
     *        solution to an ouput writer.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer, const SpatialFunction &globalSol)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_.numVertices();
        ScalarField *p = writer.template createField<Scalar, 1>(numVertices);

        LocalFunction   tmpSol;
        VertexDataArray elemDat(BoxTraits::ShapeFunctionSetContainer::maxsize);

        ElementIterator elementIt = this->problem_.elementBegin();
        ElementIterator endit = this->problem_.elementEnd();
        for (; elementIt != endit; ++elementIt) {
            int numLocalVerts = elementIt->template count<dim>();
            tmpSol.resize(numLocalVerts);

            setCurrentElement(*elementIt);
            this->restrictToElement(tmpSol, globalSol);
            updateElementData_(elemDat, tmpSol, false);

            for (int i = 0; i < elementIt->template count<dim>(); ++i) {
                int globalIdx = this->problem_.vertexIdx(*elementIt, i);

                (*p)[globalIdx] = elemDat[i].pressure;
            };
        }

        writer.addVertexData(p, "p");
    }

private:
    Scalar temperature_() const
    { return this->problem_.temperature(); }

    ThisType &asImp_()
    { return *static_cast<ThisType *>(this); }

    const ThisType &asImp_() const
    { return *static_cast<const ThisType *>(this); }
};


///////////////////////////////////////////////////////////////////////////
// OnePBoxModel (The actual numerical model.)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief Adaption of the BOX scheme to the single phase isothermal flow model.
 */
template<class ProblemT>
class OnePBoxModel : public BoxScheme< // The implementation of the model
    OnePBoxModel<ProblemT>,

    // The Traits for the BOX method
    P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                typename ProblemT::DomainTraits::Grid,
                OnePTraits::numEq>,

    // The actual problem we would like to solve
    ProblemT,

    // The local jacobian operator
    OnePBoxJacobian<ProblemT,
                    P1BoxTraits<typename ProblemT::DomainTraits::Scalar,
                                typename ProblemT::DomainTraits::Grid,
                                OnePTraits::numEq>,
                    OnePTraits> >
{
    typedef typename ProblemT::DomainTraits::Grid   Grid;
    typedef typename ProblemT::DomainTraits::Scalar Scalar;
    typedef OnePBoxModel<ProblemT>              ThisType;

public:
    typedef P1BoxTraits<Scalar, Grid, OnePTraits::numEq> BoxTraits;
    typedef Dune::OnePTraits                             OnePTraits;

private:
    typedef OnePBoxJacobian<ProblemT, BoxTraits, OnePTraits>  OnePLocalJacobian;
    typedef BoxScheme<ThisType,
                      BoxTraits,
                      ProblemT,
                      OnePLocalJacobian>  ParentType;

public:
    typedef NewNewtonMethod<ThisType> NewtonMethod;

    OnePBoxModel(ProblemT &prob)
        : ParentType(prob, richardsLocalJacobian_),
          richardsLocalJacobian_(prob)
    {
        Api::require<Api::BasicDomainTraits, typename ProblemT::DomainTraits>();
    }

    /*!
     * \brief All relevant primary and secondary of the current
     *        solution to an ouput writer.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer) 
    {
        richardsLocalJacobian_.addVtkFields(writer, this->currentSolution());
    }

private:
    // calculates the jacobian matrix at a given position
    OnePLocalJacobian  richardsLocalJacobian_;
};
}

#endif
