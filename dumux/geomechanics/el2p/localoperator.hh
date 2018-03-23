// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This file contains a local operator for PDELab which
 * wraps the contributions from
 * el2plocalresidual (box discretized mass balances)
 * and alphaMomentum (FE discretized momentum balance).
 */
#ifndef DUMUX_EL2P_LOCAL_OPERATOR_HH
#define DUMUX_EL2P_LOCAL_OPERATOR_HH

#include<dune/common/version.hh>
#include<dune/geometry/quadraturerules.hh>

#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/gridfunctionspace/localvector.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include "properties.hh"

namespace Dumux {

namespace PDELab {

/*!
 * \brief A local operator for PDELab which wraps the contributions from
 * el2plocalresidual (box discretized mass balances)
 * and alphaMomentum (FE discretized momentum balance).
 */
template<class TypeTag>
class El2PLocalOperator
    :
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags
{
    // copying the local operator for PDELab is not a good idea
    El2PLocalOperator(const El2PLocalOperator &);

    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity::Geometry::JacobianInverseTransposed JacobianInverseTransposed;
    typedef typename GridView::Intersection Intersection;
    typedef typename Dune::PDELab::IntersectionGeometry<Intersection>::ctype DT;

    enum{numEq = GET_PROP_VALUE(TypeTag, NumEq)};
    enum{dim = GridView::dimension};
    enum{dimWorld = GridView::dimensionworld};
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    enum {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx
    };

public:
    // pattern assembly flags
    enum { doPatternVolume = true };

    // residual assembly flags
    enum { doAlphaVolume = true };

    /*!
     * \param model The physical model for the box scheme.
     */
    El2PLocalOperator(Model &model)
        : model_(model)
    {}

    /*!
     * \brief Volume integral depending on test and ansatz functions
     *
     * \tparam EG The entity geometry type from PDELab
     * \tparam LFSU The type of the local function space  of the ansatz functions
     * \tparam X The type of the container for the coefficients for the ansatz functions
     * \tparam LFSV The type of the local function space of the test functions
     * \tparam R The range type (usually FieldVector<double>)
     *
     * \param eg The entity geometry object
     * \param lfsu The local function space object of the ansatz functions
     * \param x The object of the container for the coefficients for the ansatz functions
     * \param lfsv The local function space object of the test functions
     * \param r The object storing the volume integral
     */
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x,
                       const LFSV& lfsv, R& r) const
    {
        typedef typename LFSU::Traits::SizeType size_type;

        // evaluate the local residual of the box mass balance equation for the current element
        model_.localResidual().eval(eg.entity());

        // pressure and saturation local function space (mass balance equations)
        typedef typename LFSU::template Child<0>::Type PressSatLFS;
        // local function space for pressure
        typedef typename PressSatLFS::template Child<0>::Type PressLFS;
        const PressSatLFS& pressSatLFS = lfsu.template child<0>();
        const PressLFS& pressLFS = pressSatLFS.template child<0>();
        // local function space for saturation
        typedef typename PressSatLFS::template Child<1>::Type SatLFS;
        const SatLFS& satLFS = pressSatLFS.template child<1>();

        unsigned int numScv = eg.entity().subEntities(dim);

        for (size_type i = 0; i < (numEq-dim) * numScv; i++)
        {
            // retrieve the local residual value for vertex=i%Vertices and equation i/numScv (here 0 or 1)
            Scalar tmp = model_.localResidual().residual(i%numScv)[i/numScv];
            // get residual for brine phase mass balance equation
            if(i < numScv)
                r.rawAccumulate(pressLFS, i, tmp);
            // get residual for CO2 phase mass balance equation
            else
                r.rawAccumulate(satLFS,i-numScv, tmp);
        }
        // get residual for momentum balance equation
        alphaMomentum(eg, lfsu, x, lfsv, r);
    }


    /*!
     * \brief Calculate the local residual of the momentum balance equation
     *             with the finite element method. This requires numerical
     *             integration which is done via a quadrature rule.
     *
     * \tparam EG The entity geometry type from PDELab
     * \tparam LFSU The type of the local function space  of the ansatz functions
     * \tparam X The type of the container for the coefficients for the ansatz functions
     * \tparam LFSV The type of the local function space of the test functions
     * \tparam R The range type (usually FieldVector<double>)
     *
     * \param eg The entity geometry object
     * \param lfsu The local function space object of the ansatz functions
     * \param x The object of the container for the coefficients for the ansatz functions
     * \param lfsv The local function space object of the test functions
     * \param r The object storing the volume integral
     *
     *
     */
    template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
    void alphaMomentum (const EG& eg, const LFSU& lfsu, const X& x,
                    const LFSV& lfsv, R& r) const
    {
        FVElementGeometry fvGeometry;
        fvGeometry.update(model_.problem().gridView(), eg.entity());
        // retrieve lame parameters for calculation of effective stresses
        const Dune::FieldVector<Scalar,2> lameParams = model_.problem().spatialParams().lameParams(eg.entity(), fvGeometry, 0);
        Scalar lambda = lameParams[0];
        Scalar mu = lameParams[1];
        // retrieve materialParams for calculate of capillary pressure
        const MaterialLawParams& materialParams =
            model_.problem().spatialParams().materialLawParams(eg.entity(), fvGeometry, 0);
        // retrieve initial porosity
        Scalar porosity = model_.problem().spatialParams().porosity(eg.entity(), fvGeometry, 0);

        // order of quadrature rule
        const int qorder = 3;

        // extract local function spaces
        // pressure and saturation local function space (mass balance equations)
        typedef typename LFSU::template Child<0>::Type PressSatLFS;
        const PressSatLFS& pressSatLFS = lfsu.template child<0>();
        // local function space for pressure
        typedef typename PressSatLFS::template Child<0>::Type PressLFS;
        const PressLFS& pressLFS = pressSatLFS.template child<0>();
        const unsigned int pressSize = pressLFS.size();
        // local function space for saturation
        typedef typename PressSatLFS::template Child<1>::Type SatLFS;
        const SatLFS& satLFS = pressSatLFS.template child<1>();
        // local function space for solid displacement
        typedef typename LFSU::template Child<1>::Type DisplacementLFS;
        typedef typename DisplacementLFS::template Child<0>::Type DisplacementScalarLFS;
        const DisplacementLFS& displacementLFS = lfsu.template child<1>();
        const DisplacementScalarLFS& uScalarLFS = displacementLFS.template child<0>();
        const unsigned int dispSize = displacementLFS.template child<0>().size();

        // domain and range field type
        typedef typename DisplacementScalarLFS::Traits::FiniteElementType::
                        Traits::LocalBasisType::Traits::RangeFieldType RF;
        typedef typename DisplacementScalarLFS::Traits::FiniteElementType::
                        Traits::LocalBasisType::Traits::RangeType RT_V;
        typedef typename DisplacementScalarLFS::Traits::FiniteElementType::
                        Traits::LocalBasisType::Traits::JacobianType JacobianType_V;
        typedef typename PressLFS::Traits::FiniteElementType::
                        Traits::LocalBasisType::Traits::DomainFieldType DF;
        typedef typename PressLFS::Traits::FiniteElementType::
                        Traits::LocalBasisType::Traits::RangeType RT_P;

        // select quadrature rule for the element geometry type and with the order=qorder
        const auto geometry = eg.geometry();
        Dune::GeometryType geomType = geometry.type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(geomType,qorder);

        // loop over quadrature points
        for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
        {
            // evaluate reference element gradients of shape functions at quadrature point
            // (we assume Galerkin method lfsu=lfsv)
            std::vector<JacobianType_V> vGradRef(dispSize);
            uScalarLFS.finiteElement().localBasis().evaluateJacobian(it->position(),vGradRef);


             // get inverse transposed jacobian for quadrature point
             const JacobianInverseTransposed jacobian = geometry.jacobianInverseTransposed(it->position());

             // calculate shape function gradients at the quadrature point in global coordinates. This is done
             // by multiplying the reference element shape functions with the inverse transposed jacobian
             std::vector<Dune::FieldVector<RF,dim> > vGrad(dispSize);
             for (size_t i = 0; i < dispSize; i++)
             {
                vGrad[i] = 0.0;
                jacobian.umv(vGradRef[i][0],vGrad[i]);
             }

             // calculate the gradient of the solid displacement vector uGrad
             // x(uLFS,i) is the solid displacement entry of the solution vector
             // for element vertex i and coordinate direction coordDir

             Dune::FieldMatrix<RF,dim,dim> uGrad(0.0);
             for(int coordDir = 0; coordDir < dim; ++coordDir) {
                const DisplacementScalarLFS& uLFS = displacementLFS.child(coordDir);
                // compute gradient of u
                for (size_t i = 0; i < dispSize; i++)
                    uGrad[coordDir].axpy(x(uLFS,i),vGrad[i]);
             }
             // calculate the strain tensor epsilon
             Dune::FieldMatrix<RF,dim,dim> epsilon;
             for(int i = 0; i < dim; ++i)
                for(int j = 0; j < dim; ++j)
                    epsilon[i][j] = 0.5*(uGrad[i][j] + uGrad[j][i]);

             RF traceEpsilon = 0;
             for(int i = 0; i < dim; ++i)
                traceEpsilon += epsilon[i][i];

             // calculate the effective stress tensor effStress
             Dune::FieldMatrix<RF,dim,dim> effStress(0.0);
             for(int i = 0; i < dim; ++i)
             {
                effStress[i][i] = lambda*traceEpsilon;
                for(int j = 0; j < dim; ++j)
                    effStress[i][j] += 2.0*mu*epsilon[i][j];
             }

             // retrieve the shape functions for interpolating the primary variables at the
             // current quadrature point
             std::vector<RT_P> q(pressSize);
             pressLFS.finiteElement().localBasis().evaluateFunction(it->position(),q);

             RT_P pw(0.0);
             RT_P sn(0.0);
             RT_P ux(0.0);
             RT_P uy(0.0);
             RT_P uz(0.0);

             // interpolate primary variables at current quadrature point
             for (size_t i = 0; i < pressLFS.size(); i++)
             {
                pw += x(pressLFS,i) * q[i];
                sn += x(satLFS,i) * q[i];
                ux += x(displacementLFS.child(0),i) * q[i];
                if (dim > 1)
                    uy += x(displacementLFS.child(1),i) * q[i];
                if (dim > 2)
                    uz += x(displacementLFS.child(2),i) * q[i];
             }
             RT_P sw = 1.0 - sn;
             RT_P pn = pw + MaterialLaw::pc(materialParams, sw);
             RT_P pEff;

             const GlobalPosition& globalPos = geometry.global(it->position());

             // calculate change in effective pressure with respect to initial conditions pInit (pInit is negativ)
             pEff = pw*sw + pn*sn + model_.problem().pInit(globalPos, it->position(), eg.entity());
             RF uDiv = traceEpsilon;

             // fill primary variable vector for current quadrature point
             PrimaryVariables primVars;

             primVars[wPhaseIdx] = pw;
             primVars[nPhaseIdx] = sn;
             primVars[Indices::uxIdx] = ux;
             if (dim > 1)
                 primVars[Indices::uyIdx] = uy;
             if (dim > 2)
                 primVars[Indices::uzIdx] = uz;

             VolumeVariables volVars;
             // evaluate volume variables for this quadrature point
             // NOTE: this overwrites the entries of the volumevariables of node 0
             //       and can cause errors
             volVars.update(primVars, model_.problem(), eg.entity(), fvGeometry, 0, false);

             ElementVolumeVariables elemVolVars;
             elemVolVars.update(model_.problem(), eg.entity(), fvGeometry, false);

             // calculate the density difference for the gravity term
             RF rhoDiff = volVars.density(nPhaseIdx) - volVars.density(wPhaseIdx);

             // geometric weight need for quadrature rule evaluation (numerical integration)
             RF qWeight = it->weight() * geometry.integrationElement(it->position());

             // evaluate basis functions
             std::vector<RT_V> vBasis(dispSize);
             displacementLFS.child(0).finiteElement().localBasis().evaluateFunction(it->position(), vBasis);

             for(int coordDir = 0; coordDir < dim; ++coordDir) {
                const DisplacementScalarLFS& uLFS = displacementLFS.child(coordDir);
                // assemble momentum balance equation
                for (size_t i = 0; i < dispSize; i++){
                    // multiply effective stress with gradient of weighting function and geometric weight of quadrature rule
                    Scalar tmp = (effStress[coordDir] * vGrad[i]) * qWeight;
                    r.rawAccumulate(uLFS,i,tmp);

                    // subtract effective pressure change contribution multiplied with gradient of weighting function
                    // and geometric weight of quadrature rule (soil mechanics sign conventions, compressive stresses are negative)
                    tmp = -(pEff * vGrad[i][coordDir]) * qWeight;
                    r.rawAccumulate(uLFS,i,tmp);

                    // evaluate gravity term (soil mechanics sign conventions, compressive stresses are negative)
                    // multiplied with weighting function and geometric weight of quadrature rule.
                    // This assumes that the solid phase density remains constant, that the changes in porosity are very small,
                    // and that the density of the brine phase remains constant
                    tmp = sn*elemVolVars[i].effPorosity*rhoDiff*model_.problem().gravity()[coordDir]*vBasis[i]* qWeight;
                    r.rawAccumulate(uLFS,i,tmp);
                }
            }
        }
        // include boundary conditions
        // iterate over element intersections of codim dim-1
        for (const auto& intersection : intersections(model_.problem().gridView(), eg.entity()))
        {
            // handle only faces on the boundary
            if (!intersection.boundary())
                continue;

            // select quadrature rule for intersection faces (dim-1)
            Dune::GeometryType gtface = intersection.geometryInInside().type();
            const Dune::QuadratureRule<DF,dim-1>& faceRule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

            // get face index of this intersection
            int fIdx = intersection.indexInInside();
            // get dimension of face
            const int dimIs = Dune::PDELab::IntersectionGeometry<Intersection>::Entity::Geometry::mydimension;

            // get reference element for intersection geometry (reference element for face if dim = 3)
            const Dune::ReferenceElement<DT,dimIs>& refElement = Dune::ReferenceElements<DT,dimIs>::general(geomType);
            // get reference element for edges of intersection geometry (reference element for edge if dim = 3), needed for Dirichlet BC
            const Dune::ReferenceElement<DT,dimIs-1> &face_refElement =
                Dune::ReferenceElements<DT,dimIs-1>::general(intersection.geometryInInside().type());

            // Treat Neumann boundary conditions
            // loop over quadrature points and integrate normal stress changes (traction changes)
            for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=faceRule.begin(); it!=faceRule.end(); ++it)
            {
                // position of quadrature point in local coordinates of element
                DimVector local = intersection.geometryInInside().global(it->position());

                GlobalPosition globalPos = geometry.global(local);

                // evaluate boundary condition type
                BoundaryTypes boundaryTypes;
                model_.problem().boundaryTypesAtPos(boundaryTypes, globalPos);

                // skip rest if we are on Dirichlet boundary
                if (!boundaryTypes.hasNeumann())
                    continue;

                // evaluate basis functions of all all element vertices for quadrature point location "local"
                std::vector<RT_V> vBasis(dispSize);
                displacementLFS.child(0).finiteElement().localBasis().evaluateFunction(local, vBasis);

                // evaluate stress boundary condition. The stress change is assumed to be in normal direction (i.e. traction)
                PrimaryVariables traction;
                model_.problem().neumannAtPos(traction, globalPos);

                // get quadrature rule weight for intersection
                const RF qWeight = it->weight() * intersection.geometry().integrationElement(it->position());

                for(unsigned int coordDir=0; coordDir<dim; ++coordDir){
                    const DisplacementScalarLFS& uLFS = displacementLFS.child(coordDir);
                    // get the traction values for the current quadrature point,
                    // multiply it with the basis function and the quadrature rule weight
                    // and add it to the residual
                    if (boundaryTypes.isNeumann(Indices::momentum(coordDir)))
                        for (size_t i = 0; i < dispSize; i++){
                            Scalar tmp = -traction[Indices::momentum(coordDir)] * vBasis[i] * qWeight;
                            r.rawAccumulate(uLFS,i,tmp);
                        }

                }
            }

            // Treat Dirichlet boundary conditions, for Dirichlet boundaries we need to check vertices
            // first do loop over degrees of freedom for displacement vector entry, then check codim of this degree of freedom
            // then do loop over the current intersection face for the degrees of freedom with the given codim
            // compare the subentity of the element loop with the subentity of the intersection face loop
            // if the subentities are identical retrieve the coordinates of the intersection face subentity and evaluate the boundary
            // condition type and if it is a Dirichlet boundary condition then retrieve the Dirichlet value.
            // subtract the Dirichlet value from the corresponding solution vector entry (for this the outer element loop is needed)
            // and also subtract the residual value which has already been calculated for this degree of freedom
            // write the result into the residual

            for(unsigned int coordDir=0; coordDir<dim; ++coordDir){
                const DisplacementScalarLFS& uLFS = displacementLFS.child(coordDir);

                // loop over number of element vertices
                for (size_t i = 0; i <  dispSize; i++)
                {
                    // Get the codim to which this degree of freedom is attached to (should be a vertex)
                    unsigned int codim = displacementLFS.child(0).finiteElement().localCoefficients().localKey(i).codim();
                    // if we are within the element do nothing (this could happen if second order approximations are applied)
                    if (codim==0) continue;

                    // iterate over number of degrees of freedom with the given codim which are attached to the current intersection face
                    for (int j = 0; j <  refElement.size(fIdx,1,codim); j++)
                    {   // check if degree of freedom is located on a vertex of the current intersection (boundary face)
                        if (displacementLFS.child(0).finiteElement().localCoefficients().localKey(i).subEntity() ==
                                        refElement.subEntity(fIdx,1,j,codim))
                        {
                            // get local coordinate for this degree of freedom
//                             this doesn't work: DimVector local = intersection.geometryInInside().global(face_refElement.position(j,codim-1));
                            DimVector local = refElement.template geometry<1>(fIdx).global(face_refElement.position(j, codim-1));

                            GlobalPosition globalPos = geometry.global(local);

                            // evaluate boundary condition type
                            BoundaryTypes boundaryTypes;
                            model_.problem().boundaryTypesAtPos(boundaryTypes, globalPos);

                            if (boundaryTypes.isDirichlet(Indices::u(coordDir)))
                            {
                                // set value of dirichlet BC
                                PrimaryVariables dirichletValues;
                                model_.problem().dirichletAtPos(dirichletValues, globalPos);
                                // retrieve residual value which has already been calculated for the given vertex before it
                                // was clear that we are on a Dirichlet boundary
                                Scalar tmpResVal = r.container().base()[(numEq-dim)*dispSize + coordDir*dispSize + i];
                                // subtract the dirichletValue and the stored residual value from the solution vector entry
                                // if the solution vector entry equals the dirichletValue the residual will be zero
                                Scalar tmp = x(uLFS,i) - dirichletValues[Indices::u(coordDir)] - tmpResVal;
                                // write result into the residual vector
                                r.rawAccumulate(uLFS,i,tmp);
                            }
                        }
                    }
                }
            }
        }
    }

    /*!
     * \brief Jacobian of volume term
     *
     * \tparam EG The entity geometry type from PDELab
     * \tparam LFSU The type of the local function space  of the ansatz functions
     * \tparam X The type of the container for the coefficients for the ansatz functions
     * \tparam LFSV The type of the local function space of the test functions
     * \tparam M The matrix type
     *
     * \param eg The entity geometry object
     * \param lfsu The local function space object of the ansatz functions
     * \param x The object of the container for the coefficients for the ansatz functions
     * \param lfsv The local function space object of the test functions
     * \param mat The object containing the local jacobian matrix
     */
    template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
    void jacobian_volume (const EG& eg,
                          const LFSU& lfsu,
                          const X& x,
                          const LFSV& lfsv,
                          M& mat) const
    {
        typedef typename LFSU::Traits::SizeType size_type;

        model_.localJacobian().assemble(eg.entity());
        // pressure and saturation local function space (mass balance equations)
        typedef typename LFSU::template Child<0>::Type PressSatLFS;
        typedef typename PressSatLFS::template Child<0>::Type PressLFS;
        const PressSatLFS& pressSatLFS = lfsu.template child<0>();
        const PressLFS& pressLFS = pressSatLFS.template child<0>();
        typedef typename PressSatLFS::template Child<1>::Type SatLFS;
        const SatLFS& satLFS = pressSatLFS.template child<1>();
        // local function space for solid displacement
        typedef typename LFSU::template Child<1>::Type DisplacementLFS;
        typedef typename DisplacementLFS::template Child<0>::Type DisplacementScalarLFS;
        const DisplacementLFS& displacementLFS = lfsu.template child<1>();

        // type of local residual vector
        typedef typename M::value_type R;
        typedef Dune::PDELab::LocalVector<R> LocalResidualVector;
        typedef Dune::PDELab::WeightedVectorAccumulationView<LocalResidualVector> ResidualView;

        unsigned int numScv = eg.entity().subEntities(dim);

        // loop over all degrees of freedom of the current element
        for (size_type j = 0; j < numScv*numEq; j++)
        {
            // assemble entries for mass balance equations
            for (size_type i = 0; i < (numEq-dim)*numScv; i++)
            {
                // local jacobian value of location idxI=i%numScv, idxJ=j%numScv for equation i/numScv and unknown j/numScv
                Scalar tmp = (model_.localJacobian().mat(i%numScv,j%numScv))[i/numScv][j/numScv];
                // mass balance entries for pressure
                if (j < numScv){
                    if(i < numScv)
                        mat.rawAccumulate(pressLFS,i,pressLFS,j,tmp);
                    else
                        mat.rawAccumulate(satLFS,i-numScv,pressLFS,j,tmp);
                }
                // mass balance entries for saturation
                else if (j < 2*numScv){
                    if(i < numScv)
                        mat.rawAccumulate(pressLFS,i,satLFS,j-numScv,tmp);
                    else
                        mat.rawAccumulate(satLFS,i-numScv,satLFS,j-numScv,tmp);
                }
                // mass balance entries for solid displacement in x-direction
                else if  (j < 3*numScv)
                {
                    const DisplacementScalarLFS& uScalarLFS = displacementLFS.template child<0>();
                    if(i < numScv)
                        mat.rawAccumulate(pressLFS,i,uScalarLFS,j-2*numScv,tmp);
                    else
                        mat.rawAccumulate(satLFS,i-numScv,uScalarLFS,j-2*numScv,tmp);
                }
                // mass balance entries for solid displacement in y-direction
                else if (j < 4*numScv && dim >=2)
                {
                    const DisplacementScalarLFS& uScalarLFS = displacementLFS.template child<1>();
                    if(i < numScv)
                        mat.rawAccumulate(pressLFS,i,uScalarLFS,j-3*numScv,tmp);
                    else
                        mat.rawAccumulate(satLFS,i-numScv,uScalarLFS,j-3*numScv,tmp);
                }
                // mass balance entries for solid displacement in z-direction
                else if(j < 5*numScv && dim >=3)
                {
                    const DisplacementScalarLFS& uScalarLFS = displacementLFS.template child<dim-1>();
                    if(i < numScv)
                        mat.rawAccumulate(pressLFS,i,uScalarLFS,j-(numEq-1)*numScv,tmp);
                    else
                        mat.rawAccumulate(satLFS,i-numScv,uScalarLFS,j-(numEq-1)*numScv,tmp);
                }

            }
        }

        // calculate local jacobian entries and assemble for momentum balance equation
        const int m=lfsv.size();
        const int n=lfsu.size();

        X u(x);
        LocalResidualVector down(mat.nrows(),0);

        // evaluate momentum residual for momentum balance equation
        ResidualView downView = down.weightedAccumulationView(1.0);
        alphaMomentum(eg, lfsu, u, lfsv, downView);

        // loop over all columns (number of element vertices * number of equations)
        using std::abs;
        for (int j = 0; j < n; j++)
        {
          // vary the solution vector entry (lfsu,j) by a small value delta (forward differencing)
          // this comprises presure, saturation, ux, uy and uz
          Scalar delta = 1e-4*(1.0+abs(u(lfsu,j)));
          u(lfsu,j) += delta;

          // evaluate momentum balance residual for the varied solution vector
          LocalResidualVector up(mat.nrows(), 0);
          ResidualView upView = up.weightedAccumulationView(1.0);
          alphaMomentum(eg, lfsu, u, lfsv, upView);

          // calculate partial derivative for momentum balance equations and assemble
          for (int i = (numEq-dim)*numScv; i < m; i++)
          {
              Scalar entry = (up(lfsv, i) - down(lfsv, i))/delta;
              // accumulate resulting partial derivatives into jacobian
              mat.rawAccumulate(lfsv,i, lfsu,j,entry);
          }

          // reset solution
          u(lfsu,j) = x(lfsu,j);
        }
    }

private:
    Model& model_;
};

} // namespace PDELab
} // namespace Dumux

#endif
