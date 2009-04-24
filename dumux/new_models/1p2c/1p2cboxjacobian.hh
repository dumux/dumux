// $Id$

#ifndef DUMUX_NEW_1P2C_BOX_JACOBIAN_HH
#define DUMUX_NEW_1P2C_BOX_JACOBIAN_HH

#include <dumux/new_models/boxscheme/boxscheme.hh>
#include <dumux/new_models/boxscheme/p1boxtraits.hh>
#include <dumux/new_models/1p2c/1p2ctraits.hh>

#include <dumux/auxiliary/apis.hh>
#include <dune/common/collectivecommunication.hh>
#include <vector>
#include <iostream>

namespace Dune
{

// forward declaration of the 1p2c box model
//template<class ProblemT, class OnePTwoCTraitsT>
//class OnePTwoCBoxModel;

/*!
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the single-phase, two-component model.
 */
template <class OnePTwoCTraits, 
          class Problem>
class OnePTwoCVertexData
{
    typedef OnePTwoCTraits Tr;
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

        double T = 273;
        double p = 1e-5;

        porosity = jac.problem.soil().porosity(global,element,local);
        viscosity = jac.problem.phase().viscosity(T,p);
        tortuosity = jac.problem.soil().tortuosity(global,element,local);
        diffCoeff = jac.problem.phase().diffCoeff();
        molefraction = vertSol[transport];
        pressure = vertSol[konti];
    }

    
    Scalar porosity;
        Scalar viscosity;
        Scalar tortuosity;
        Scalar diffCoeff;
        Scalar molefraction;
        Scalar pressure;
    };


///////////////////////////////////////////////////////////////////////////
// OnePTwoCBoxJacobian (evaluate the local jacobian for the newton method.)
///////////////////////////////////////////////////////////////////////////
/*!
 * \brief 1P-2C specific details needed to approximately calculate
 *        the local jacobian in the BOX scheme.
 *
 * This class is used to fill the gaps in BoxJacobian for the 1P-2C onephase flow.
 */
template<class ProblemT,
         class BoxTraitsT,
         class OnePTwoCTraitsT>
class OnePTwoCBoxJacobian : public BoxJacobian<ProblemT,
                                                   BoxTraitsT,
                                                   OnePTwoCVertexData<OnePTwoCTraitsT, ProblemT>,
                                                   OnePTwoCBoxJacobian<ProblemT,BoxTraitsT,OnePTwoCTraitsT> >
{
protected:
//    friend class OnePTwoCBoxModel<ProblemT, OnePTwoCTraitsT>;

    typedef OnePTwoCBoxJacobian<ProblemT,
                                    BoxTraitsT,
                                    OnePTwoCTraitsT>             Implementation;
    typedef BoxJacobian<ProblemT, BoxTraitsT, OnePTwoCBoxJacobian>    ParentType;
    typedef ProblemT                                Problem;
    typedef typename Problem::DomainTraits          DomTraits;
    typedef BoxTraitsT                              BoxTraits;
    typedef OnePTwoCTraitsT                         OnePTwoCTraits;

    enum {
        dim              = DomTraits::dim,
        dimWorld         = DomTraits::dimWorld,

        numEq            = BoxTraits::numEq,
        numPhases        = OnePTwoCTraits::numPhases,
        numComponents    = OnePTwoCTraits::numComponents,

        konti      		 = OnePTwoCTraits::konti,
        transport        = OnePTwoCTraits::transport
    };

    typedef typename DomTraits::Scalar                Scalar;
    typedef typename DomTraits::CoordScalar           CoordScalar;
    typedef typename DomTraits::Grid                  Grid;
    typedef typename DomTraits::Element               Element;
    typedef typename DomTraits::ElementIterator       ElementIterator;
    typedef typename Element::EntityPointer           ElementPointer;
    typedef typename DomTraits::LocalPosition         LocalPosition;
    typedef typename DomTraits::GlobalPosition        GlobalPosition;
    typedef typename DomTraits::VertexIterator        VertexIterator;


    typedef typename BoxTraits::SolutionVector      SolutionVector;
    typedef typename BoxTraits::FVElementGeometry   FVElementGeometry;
    typedef typename BoxTraits::SpatialFunction     SpatialFunction;
    typedef typename BoxTraits::LocalFunction       LocalFunction;
    typedef typename Grid::CollectiveCommunication  CollectiveCommunication;

    typedef typename OnePTwoCTraits::PhasesVector        PhasesVector;
    typedef typename OnePTwoCTraits::VariableVertexData  VariableVertexData;
    typedef FieldMatrix<Scalar, dim, dim>  Tensor;
    typedef FieldVector<Scalar,dim> FieldVector;


    /*!
     * \brief Data which is attached to each vertex and is not only
     *        stored locally.
     */
    struct StaticVertexData {

    };


public:
    OnePTwoCBoxJacobian(ProblemT &problem)
        : ParentType(problem),
          staticVertexDat_(problem.numVertices())
    {

    };


    /*!
     * \brief Evaluate the rate of change of all conservation
     *        quantites (e.g. phase mass) within a sub control
     *        volume of a finite volume element in the 2P-2C
     *        model.
     *
     * This function should not include the source and sink terms.
     */
    void computeStorage(SolutionVector &result, int scvId, bool usePrevSol) const
    {
    	result = Scalar(0);

        // if flag usePrevSol is set, the solution from the previous time step is used,
        // otherwise the current solution is used. The secondary variables are used accordingly.
        // This computes the derivative of the storage term.
        //const LocalFunction &sol   = usePrevSol ? this->prevSol_ : this->curSol_;
        const ElementData &elementCache = usePrevSol ? prevElemDat_  : curElemDat_;

        // storage term of continuity equation
        result[konti] = 0;

        // storage term of the transport equation
        result[transport] = elementCache.vertex[scvId].porosity * elementCache.vertex[scvId].molefraction;

    }

    /*!
     * \brief Evaluates the mass flux over a face of a subcontrol
     *        volume.
     */
    void computeFlux(SolutionVector &flux, int faceId) const
    {
        // set flux vector to zero
        int i = this->curElementGeom_.subContVolFace[faceId].i;
        int j = this->curElementGeom_.subContVolFace[faceId].j;

        // normal vector, value of the area of the scvf
        const GlobalPosition &normal(this->curElementGeom_.subContVolFace[faceId].normal);

        // get global coordinates of verts i,j
        const GlobalPosition &global_i = this->curElementGeom_.subContVol[i].global;
        const GlobalPosition &global_j = this->curElementGeom_.subContVol[j].global;

        // get local coordinates of verts i,j
        const LocalPosition &local_i = this->curElementGeom_.subContVol[i].local;
        const LocalPosition &local_j = this->curElementGeom_.subContVol[j].local;

        const ElementData &elemDat = this->curElemDat_;
        const VariableVertexData &vDat_i = elemDat.vertex[i];
        const VariableVertexData &vDat_j = elemDat.vertex[j];

        GlobalPosition pGrad;
        GlobalPosition xGrad;
        GlobalPosition KpGrad;

        pGrad = Scalar(0);
        xGrad = Scalar(0);
        KpGrad = Scalar(0);


        GlobalPosition tmp(0.0);
        PhasesVector pressure(0.0), molefraction(0.0);


        // calculate FE gradient (grad p)
        for (int idx = 0; idx < this->curElementGeom_.numVertices; idx++) // loop over adjacent vertices
            {
                // FEGradient at vertex idx
                const LocalPosition &feGrad = this->curElementGeom_.subContVolFace[faceId].grad[idx];


                // the pressure gradient
                tmp = feGrad;
                tmp *= elemDat.vertex[idx].pressure;
                pGrad += tmp;


                // the concentration gradient
                tmp = feGrad;
                tmp *= elemDat.vertex[idx].molefraction;
                xGrad += tmp;

            }


        // calculate the permeability tensor
        Tensor K         = this->problem_.soil().K(global_i, ParentType::curElement_(), local_i);
        const Tensor &Kj = this->problem_.soil().K(global_j, ParentType::curElement_(), local_j);
        harmonicMeanK_(K, Kj);

        ///////////////////////
        // calculation of the flux of the continuity equation
        //////////////////////

        K.mv(pGrad, KpGrad);  // KpGrad = K * grad p
        flux[konti] = (KpGrad*normal)/elemDat.vertex[i].viscosity;



        /////////////////////////////
        // calculation of the flux of the transport equation
        /////////////////////////////

		//arithmetic mean of the tortuosity
		Scalar avgTortuosity = 0.5 * (vDat_i.tortuosity + vDat_j.tortuosity);

		//arithmetic mean of the porosity
		Scalar avgPorosity = 0.5 * (vDat_i.porosity + vDat_j.porosity);

		Scalar outward = KpGrad * normal;

		//upwinding of the advective flux
		if(outward<0)
			flux[transport] = (outward * vDat_i.molefraction / vDat_i.viscosity) +
			(avgTortuosity * avgPorosity * vDat_i.diffCoeff * (xGrad * normal));
		else
			flux[transport] = (outward * vDat_j.molefraction /vDat_j.viscosity) +
			(avgTortuosity * avgPorosity * vDat_i.diffCoeff * (xGrad * normal));

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
     * \brief Add the mass fraction of air in water to VTK output of
     *        the current timestep.
     */
    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer, const SpatialFunction &globalSol)
    {
        typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > ScalarField;

        // create the required scalar fields
        unsigned numVertices = this->problem_.numVertices();
        //unsigned numElements = this->problem_.numElements();

        ScalarField *pressure = writer.template createField<Scalar, 1>(numVertices);
        ScalarField *molefraction = writer.template createField<Scalar, 1>(numVertices);



        LocalFunction curSol(numVertices);
        ElementData   elemDat;
        //VariableVertexData tmp;

        ElementIterator elementIt = this->problem_.elementBegin();
        ElementIterator endit = this->problem_.elementEnd();

        for (; elementIt != endit; ++elementIt)
            {
                int numLocalVerts = elementIt->template count<dim>();

                setCurrentElement(*elementIt);
                this->restrictToElement(curSol, globalSol);
                updateElementData_(elemDat, curSol, false);

                for (int i = 0; i < numLocalVerts; ++i)
                    {
                        int globalIdx = this->problem_.vertexIdx(*elementIt, i);

                        (*pressure)[globalIdx] = elemDat.vertex[i].pressure;
                        (*molefraction)[globalIdx] = elemDat.vertex[i].molefraction;

                    };

            }

        writer.addVertexData(pressure, "pressure");
        writer.addVertexData(molefraction, "molefraction");
    }


protected:
    Implementation *asImp_()
    { return static_cast<Implementation *>(this); }
    const Implementation *asImp_() const
    { return static_cast<const Implementation *>(this); }
};


} // end namepace

#endif
