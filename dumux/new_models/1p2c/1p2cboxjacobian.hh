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
     * \brief Cached data for the each vertex of the element.
     */
    struct ElementData
    {
        VariableVertexData vertex[BoxTraits::ShapeFunctionSetContainer::maxsize];
    };

    /*!
     * \brief Data which is attached to each vertex and is not only
     *        stored locally.
     */
    struct StaticVertexData {

    };

    /*!
     * \brief Function to update variable data of the vertices of the
     *        the current element (essentially secondary variables)
     */
    void updateVarVertexData_(VariableVertexData &vertDat,
                              const SolutionVector &vertSol,
                              const Element &element,
                              int localIdx,
                              Problem &problem) const
    {
        const GlobalPosition &global = element.geometry().corner(localIdx);
        const LocalPosition &local =
            DomTraits::referenceElement(element.type()).position(localIdx,
                                                                 dim);

        double T = 273;
        double p = 1e-5;

        vertDat.porosity = problem.soil().porosity(global,element,local);
        vertDat.viscosity = problem.phase().viscosity(T,p);
        vertDat.tortuosity = problem.soil().tortuosity(global,element,local);
        vertDat.diffCoeff = problem.phase().diffCoeff();
        vertDat.molefraction = vertSol[transport];
        vertDat.pressure = vertSol[konti];

    }

public:
    OnePTwoCBoxJacobian(ProblemT &problem)
        : ParentType(problem),
          staticVertexDat_(problem.numVertices())
    {

    };

    /*!
     * \brief Set the current grid element.
     */
    void setCurrentElement(const Element &element)
    {
        ParentType::setCurrentElement_(element);
    };

    /*!
     * \brief Set the parameters for the calls to the remaining
     *        members.
     */
    void setParams(const Element &element, LocalFunction &curSol, LocalFunction &prevSol)
    {
        setCurrentElement(element);

        // TODO: scheme which allows not to copy curSol and
        // prevSol all the time
        curSol_ = curSol;
        updateElementData_(curElemDat_, curSol_, false);
        curSolDeflected_ = false;

        prevSol_ = prevSol;
        updateElementData_(prevElemDat_, prevSol_, true);
    };

    /*!
     * \brief Vary a single component of a single vertex of the
     *        local solution for the current element.
     *
     * This method is a optimization, since if varying a single
     * component at a degree of freedom not the whole element cache
     * needs to be recalculated. (Updating the element cache is very
     * expensive since material laws need to be evaluated.)
     */
    void deflectCurSolution(int vert, int component, Scalar value)
    {
        // make sure that the original state can be restored
        if (!curSolDeflected_) {
            curSolDeflected_ = true;

            curSolOrigValue_ = curSol_[vert][component];
            curSolOrigVarData_ = curElemDat_.vertex[vert];
        }


        curSol_[vert][component] = value;
        asImp_()->updateVarVertexData_(curElemDat_.vertex[vert],
                                               curSol_[vert],
                                               this->curElement_(),
                                               vert,
                                               this->problem_);
    }

    /*!
     * \brief Restore the local jacobian to the state before
     *        deflectCurSolution() was called.
     *
     * This only works if deflectSolution was only called with
     * (vertex, component) as arguments.
     */
    void restoreCurSolution(int vert, int component)
    {
        curSolDeflected_ = false;
        curSol_[vert][component] = curSolOrigValue_;
        curElemDat_.vertex[vert] = curSolOrigVarData_;
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
     * \brief Initialize the static data with the initial solution.
     *
     * Called by TwoPTwoCBoxModel::initial()
     */
    void initStaticData()
    {

    }

    /*!
     * \brief Update the static data of a single vert and do a
     *        variable switch if necessary.
     */
    void updateStaticData(SpatialFunction &curGlobalSol, SpatialFunction &oldGlobalSol)
    {

    }

    /*!
     * \brief Set the old phase of all verts state to the current one.
     */
    void updateOldPhaseState()
    {

    }

    /*!
     * \brief Reset the current phase state of all verts to the old one after an update failed
     */
    void resetPhaseState()
    {

    }

    /*!
     * \brief Return true if the primary variables were switched
     *        after the last timestep.
     */
    bool switched() const
    {
        return;
    }

    /*!
     * \brief Set whether there was a primary variable switch after in the last
     *        timestep.
     */
    void setSwitched(bool yesno)
    {

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


    void updateElementData_(ElementData &dest, const LocalFunction &sol, bool isOldSol)
    {

        int numVertices = this->curElement_().template count<dim>();
        for (int i = 0; i < numVertices; i++) {

            asImp_()->updateVarVertexData_(dest.vertex[i],
                                           sol[i],
                                           this->curElement_(),
                                           i,
                                           this->problem_);
        }
    }


    //  perform variable switch at a vertex; Returns true if a
    //  variable switch was performed.
    bool primaryVarSwitch_(SpatialFunction &globalSol,
                           int globalIdx,
                           const GlobalPosition &globalPos)
    {
    	return;
    }

    // harmonic mean of the permeability computed directly.  the
    // first parameter is used to store the result.
    static void harmonicMeanK_(Tensor &Ki, const Tensor &Kj)
    {
        for (int kx=0; kx < Tensor::rows; kx++){
            for (int ky=0; ky< Tensor::cols; ky++){
                if (Ki[kx][ky] != Kj[kx][ky]) {
                    Ki[kx][ky] = harmonicMean_(Ki[kx][ky], Kj[kx][ky]);
                }
            }
        }
    }

    // returns the harmonic mean of two scalars
    static Scalar harmonicMean_(Scalar x, Scalar y)
    {
        if (x == 0 || y == 0)
            return 0;
        return (2*x*y)/(x + y);
    };

    // parameters given in constructor
    std::vector<StaticVertexData> staticVertexDat_;

    // current solution
    LocalFunction      curSol_;
    ElementData        curElemDat_;

    // needed for restoreCurSolution()
    bool               curSolDeflected_;
    Scalar             curSolOrigValue_;
    VariableVertexData curSolOrigVarData_;

    // previous solution
    LocalFunction      prevSol_;
    ElementData        prevElemDat_;
    CollectiveCommunication collectiveCom_;
};


} // end namepace

#endif
