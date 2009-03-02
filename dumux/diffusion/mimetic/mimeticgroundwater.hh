// $Id$

#ifndef DUNE_MIMETICGROUNDWATER_HH
#define DUNE_MIMETICGROUNDWATER_HH

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

#include<dune/disc/operators/boundaryconditions.hh>
#include<dune/disc/groundwater/groundwater.hh>
#include"dumux/shapefunctions/CRshapefunctions.hh"
#include"dumux/operators/localstiffnessext.hh"

/**
 * @file
 * @brief  compute local stiffness matrix for conforming finite elements for diffusion equation
 * @author Peter Bastian
 */


namespace Dune
{
/** @addtogroup DISC_Disc
 *
 * @{
 */
/**
 * @brief compute local stiffness matrix for conforming finite elements for diffusion equation
 *
 */


//! A class for computing local stiffness matrices
/*! A class for computing local stiffness matrix for the
  diffusion equation

  div j = q; j = -K grad u; in Omega

  u = g on Gamma1; j*n = J on Gamma2.

  Uses conforming finite elements with the CR shape functions.
  It should work for all dimensions and element types.
  All the numbering is with respect to the reference element and the
  CR shape functions

  Template parameters are:

  - Grid  a DUNE grid type
  - RT    type used for return values
*/
template<class G, class RT, class VC>
class MimeticGroundwaterEquationLocalStiffness
    : public LocalStiffnessExt<MimeticGroundwaterEquationLocalStiffness<G,RT,VC>,G,RT,1>
{
    template<int dim>
    struct ElementLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    // grid types
    enum {n=G::dimension};
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef Dune::LevelP0Function<G,DT,(int)(0.5*n*(n+1))> KType;
    typedef typename G::Traits::LevelIndexSet IS;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;

public:
    // define the number of components of your system, this is used outside
    // to allocate the correct size of (dense) blocks with a FieldMatrix
    enum {m=1};
    enum {SIZE=CRShapeFunctionSetContainer<DT,RT,n>::maxsize};

    //! Constructor
    MimeticGroundwaterEquationLocalStiffness (DeprecatedDiffusionProblem<G,RT,VC>& params,
                                              bool levelBoundaryAsDirichlet_, const G& grid, int level=0,
                                              bool procBoundaryAsDirichlet_=true)
        : problem(params),levelBoundaryAsDirichlet(levelBoundaryAsDirichlet_),
          procBoundaryAsDirichlet(procBoundaryAsDirichlet_),
          elementmapper(grid, grid.levelIndexSet(level))
    {    }


    //! assemble local stiffness matrix for given element and order
    /*! On exit the following things have been done:
      - The stiffness matrix for the given entity and polynomial degree has been assembled and is
      accessible with the mat() method.
      - The boundary conditions have been evaluated and are accessible with the bc() method
      - The right hand side has been assembled. It contains either the value of the essential boundary
      condition or the assembled source term and neumann boundary condition. It is accessible via the rhs() method.
      @param[in]  e    a codim 0 entity reference
      @param[in]  k    order of CR basis
    */
    void assemble (const Entity& e, int k=1)
    {
        // extract some important parameters
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::CRShapeFunctionSetContainer<DT,RT,n>::value_type&
            sfs=Dune::CRShapeFunctions<DT,RT,n>::general(gt,k);
        setcurrentsize(sfs.size());

        // clear assemble data
        for (int i=0; i<sfs.size(); i++)
            {
                this->b[i] = 0;
                this->bctype[i][0] = BoundaryConditions::neumann;
                for (int j=0; j<sfs.size(); j++)
                    this->A[i][j] = 0;
            }

        assembleV(e,k);
        assembleBC(e,k);
    }

    // TODO/FIXME: this is only valid for linear problems where
    // the local stiffness matrix is independend of the current
    // solution. We need to implement this properly, but this
    // should at least make the thing compile...
    typedef Dune::FieldVector<RT, m> VBlockType;
    void assemble(const Entity &cell, const Dune::BlockVector<VBlockType>& localSolution, int orderOfShapeFns = 1)
    {
        assemble(cell, orderOfShapeFns);
    }


    //! assemble only boundary conditions for given element
    /*! On exit the following things have been done:
      - The boundary conditions have been evaluated and are accessible with the bc() method
      - The right hand side contains either the value of the essential boundary
      condition or the assembled neumann boundary condition. It is accessible via the rhs() method.
      @param[in]  e    a codim 0 entity reference
      @param[in]  k    order of CR basis
    */
    void assembleBoundaryCondition (const Entity& e, int k=1)
    {
        // extract some important parameters
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::CRShapeFunctionSetContainer<DT,RT,n>::value_type&
            sfs=Dune::CRShapeFunctions<DT,RT,n>::general(gt,k);
        setcurrentsize(sfs.size());

        // clear assemble data
        for (int i=0; i<sfs.size(); i++)
            {
                this->b[i] = 0;
                this->bctype[i][0] = BoundaryConditions::neumann;
            }

        this->template assembleBC(e,k);
    }

    void assembleElementMatrices(const Entity& e, Dune::FieldVector<DT,2*n>& faceVol,
                                 Dune::FieldMatrix<RT,2*n,2*n>& W, Dune::FieldVector<DT,2*n>& c,
                                 Dune::FieldMatrix<RT,2*n,2*n>& Pi, RT& dinv, Dune::FieldVector<DT,2*n>& F, RT& qmean)
    {
        // extract some important parameters
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::CRShapeFunctionSetContainer<DT,RT,n>::value_type&
            sfs=Dune::CRShapeFunctions<DT,RT,n>::general(gt,1);
        setcurrentsize(sfs.size());

        // cell center in reference element
        const Dune::FieldVector<DT,n>& centerLocal = Dune::ReferenceElements<DT,n>::general(gt).position(0,0);

        // get global coordinate of cell center
        Dune::FieldVector<DT,n> centerGlobal = e.geometry().global(centerLocal);

        // eval diffusion tensor, ASSUMING to be constant over each cell
        Dune::FieldMatrix<DT,n,n> K(0);
        K = problem.K(centerGlobal,e,centerLocal);

        //int elemId = elementmapper.map(e);

        K *= problem.materialLaw.mobTotal(problem.variables.sat(centerGlobal,e,centerLocal));

        // cell volume
        DT volume = e.geometry().volume();

        // build the matrices R and ~N
        Dune::FieldMatrix<DT,2*n,n> R(0), N(0);

        //       std::cout << "element " << elemId << ": center " << centerGlobal << std::endl;;

        typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

        IntersectionIterator endit = IntersectionIteratorGetter<G,LeafTag>::end(e);
        for (IntersectionIterator it = IntersectionIteratorGetter<G,LeafTag>::begin(e); it!=endit; ++it)
            {
                // get geometry type of face
                Dune::GeometryType gtf = it->intersectionSelfLocal().type();

                // local number of facet
                int i = it->numberInSelf();

                const Dune::FieldVector<DT,n>& faceLocal = sfs[i].position();
                Dune::FieldVector<DT,n> faceGlobal = e.geometry().global(faceLocal);
                faceVol[i] = it->intersectionGlobal().volume();

                //           std::cout << "  face " << i << ": local = " << faceLocal << ", global = " << faceGlobal << std::endl;
                //           std::cout << "    boundary = " << it.boundary() << ", neighbor = " << it->neighbor() << std::endl;
                //", outside elemId = " << elementmapper.map(*(it->outside())) << std::endl;

                // center in face's reference element
                const Dune::FieldVector<DT,n-1>&
                    faceLocalNm1 = Dune::ReferenceElements<DT,n-1>::general(gtf).position(0,0);

                // get normal vector
                Dune::FieldVector<DT,n> unitOuterNormal = it->unitOuterNormal(faceLocalNm1);

                N[i] = unitOuterNormal;

                for (int k = 0; k < n; k++)
                    // move origin to the center of gravity
                    R[i][k] = faceVol[i]*(faceGlobal[k] - centerGlobal[k]);
            }

        //      std::cout << "N =\n" << N;
        //      std::cout << "R =\n" << R;

        // proceed along the lines of Algorithm 1 from
        // Brezzi/Lipnikov/Simonicini M3AS 2005
        // (1) orthonormalize columns of the matrix R
        RT norm = R[0][0]*R[0][0];
        for (int i = 1; i < sfs.size(); i++)
            norm += R[i][0]*R[i][0];
        norm = sqrt(norm);
        for (int i = 0; i < sfs.size(); i++)
            R[i][0] /= norm;
        RT weight = R[0][1]*R[0][0];
        for (int i = 1; i < sfs.size(); i++)
            weight += R[i][1]*R[i][0];
        for (int i = 0; i < sfs.size(); i++)
            R[i][1] -= weight*R[i][0];
        norm = R[0][1]*R[0][1];
        for (int i = 1; i < sfs.size(); i++)
            norm += R[i][1]*R[i][1];
        norm = sqrt(norm);
        for (int i = 0; i < sfs.size(); i++)
            R[i][1] /= norm;
        if (n == 3) {
            RT weight1 = R[0][2]*R[0][0];
            RT weight2 = R[0][2]*R[0][1];
            for (int i = 1; i < sfs.size(); i++) {
                weight1 += R[i][2]*R[i][0];
                weight2 += R[i][2]*R[i][1];
            }
            for (int i = 0; i < sfs.size(); i++)
                R[i][1] -= weight1*R[i][0] + weight2*R[i][1];
            norm = R[0][2]*R[0][2];
            for (int i = 1; i < sfs.size(); i++)
                norm += R[i][2]*R[i][2];
            norm = sqrt(norm);
            for (int i = 0; i < sfs.size(); i++)
                R[i][2] /= norm;
        }
        //      std::cout << "~R =\n" << R;

        // (2) Build the matrix ~D
        FieldMatrix<DT,2*n,2*n> D(0);
        for (int s = 0; s < sfs.size(); s++) {
            Dune::FieldVector<DT,2*n> es(0);
            es[s] = 1;
            for (int k = 0; k < sfs.size(); k++) {
                D[k][s] = es[k];
                for (int i = 0; i < n; i++) {
                    D[k][s] -= R[s][i]*R[k][i];
                }
            }
        }

        DT traceK = K[0][0];
        for (int i = 1; i < n; i++)
            traceK += K[i][i];
        D *= traceK/volume;
        //      std::cout << "u~D =\n" << D;

        // (3) Build the matrix W = Minv
        FieldMatrix<DT,2*n,n> NK(N);
        NK.rightmultiply (K);
        for (int i = 0; i < sfs.size(); i++) {
            for (int j = 0; j < sfs.size(); j++) {
                W[i][j] = NK[i][0]*N[j][0];
                for (int k = 1; k < n; k++)
                    W[i][j] += NK[i][k]*N[j][k];
            }
        }

        W /= volume;
        W += D;
        //       std::cout << "W = \n" << W;
        //       std::cout << D[2][2] << ", " << D[2][0] << ", " << D[2][3] << ", " << D[2][1] << std::endl;
        //       std::cout << D[0][2] << ", " << D[0][0] << ", " << D[0][3] << ", " << D[0][1] << std::endl;
        //       std::cout << D[3][2] << ", " << D[3][0] << ", " << D[3][3] << ", " << D[3][1] << std::endl;
        //       std::cout << D[1][2] << ", " << D[1][0] << ", " << D[1][3] << ", " << D[1][1] << std::endl;


        // Now the notation is borrowed from Aarnes/Krogstadt/Lie 2006, Section 3.4.
        // The matrix W developed so far corresponds to one element-associated
        // block of the matrix B^{-1} there.

        // Corresponding to the element under consideration,
        // calculate the part of the matrix C coupling velocities and element pressures.
        // This is just a row vector of size sfs.size().
        // scale with volume
        for (int i = 0; i < sfs.size(); i++)
            c[i] = faceVol[i];

        // Set up the element part of the matrix \Pi coupling velocities
        // and pressure-traces. This is a diagonal matrix with entries given by faceVol.
        for (int i = 0; i < sfs.size(); i++)
            Pi[i][i] = faceVol[i];

        // Calculate the element part of the matrix D^{-1} = (c W c^T)^{-1} which is just a scalar value.
        Dune::FieldVector<DT,2*n> Wc(0);
        W.umv(c, Wc);
        dinv = 1.0/(c*Wc);

        // Calculate the element part of the matrix F = Pi W c^T which is a column vector.
        F = 0;
        Pi.umv(Wc, F);
        //      std::cout << "Pi = \n" << Pi << "c = " << c << ", F = " << F << std::endl;

        // Calculate the source f
        int p = 0;
        qmean = 0;
        for (size_t g=0; g<Dune::QuadratureRules<DT,n>::rule(gt,p).size(); ++g) // run through all quadrature points
            {
                const Dune::FieldVector<DT,n>& local = Dune::QuadratureRules<DT,n>::rule(gt,p)[g].position(); // pos of integration point
                Dune::FieldVector<DT,n> global = e.geometry().global(local);     // ip in global coordinates
                double weight = Dune::QuadratureRules<DT,n>::rule(gt,p)[g].weight();// weight of quadrature point
                DT detjac = e.geometry().integrationElement(local);              // determinant of jacobian
                RT factor = weight*detjac;
                RT q = problem.source(global,e,local);
                qmean += q*factor;
            }

    }

private:
    void assembleV (const Entity& e, int k=1)
    {
        // extract some important parameters
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::CRShapeFunctionSetContainer<DT,RT,n>::value_type&
            sfs=Dune::CRShapeFunctions<DT,RT,n>::general(gt,1);
        setcurrentsize(sfs.size());

        // The notation is borrowed from Aarnes/Krogstadt/Lie 2006, Section 3.4.
        // The matrix W developed here corresponds to one element-associated
        // block of the matrix B^{-1} there.
        Dune::FieldVector<DT,2*n> faceVol(0);
        Dune::FieldMatrix<DT,2*n,2*n> W(0);
        Dune::FieldVector<DT,2*n> c(0);
        Dune::FieldMatrix<DT,2*n,2*n> Pi(0);
        Dune::FieldVector<RT,2*n> F(0);
        RT dinv;
        RT qmean;
        this->assembleElementMatrices(e, faceVol, W, c, Pi, dinv, F, qmean);

        // Calculate the element part of the matrix Pi W Pi^T.
        Dune::FieldMatrix<RT,2*n,2*n> PiWPiT(W);
        PiWPiT.rightmultiply(Pi);
        PiWPiT.leftmultiply(Pi);

        // Calculate the element part of the matrix F D^{-1} F^T.
        Dune::FieldMatrix<RT,2*n,2*n> FDinvFT(0);
        for (int i = 0; i < sfs.size(); i++)
            for (int j = 0; j < sfs.size(); j++)
                FDinvFT[i][j] = dinv*F[i]*F[j];

        // Calculate the element part of the matrix S = Pi W Pi^T - F D^{-1} F^T.
        for (int i = 0; i < sfs.size(); i++)
            for (int j = 0; j < sfs.size(); j++)
                this->A[i][j] = PiWPiT[i][j] - FDinvFT[i][j];

        // Calculate the source term F D^{-1} f
        // NOT WORKING AT THE MOMENT
        RT factor = dinv*qmean;
        for (int i = 0; i < sfs.size(); i++)
            this->b[i] = F[i]*factor;


        //        std::cout << "faceVol = " << faceVol << std::endl << "W = " << std::endl << W << std::endl
        //              << "c = " << c << std::endl << "Pi = " << std::endl << Pi << std::endl
        //              << "dinv = " << dinv << std::endl << "F = " << F << std::endl;
        //        std::cout << "dinvF = " << dinvF << ", q = " << qmean
        //             << ", b = " << this->b[0] << ", " << this->b[1] << ", " << this->b[2] << ", " << this->b[3] << std::endl;
    }

    void assembleBC (const Entity& e, int k=1)
    {
        // extract some important parameters
        Dune::GeometryType gt = e.geometry().type();
        const typename Dune::CRShapeFunctionSetContainer<DT,RT,n>::value_type&
            sfs=Dune::CRShapeFunctions<DT,RT,n>::general(gt,k);
        setcurrentsize(sfs.size());

        // determine quadrature order
        int p=0;

        // evaluate boundary conditions via intersection iterator
        typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator
            IntersectionIterator;

        //std::cout << "new element." << std::endl;
        IntersectionIterator endit = IntersectionIteratorGetter<G,LeafTag>::end(e);
        for (IntersectionIterator it = IntersectionIteratorGetter<G,LeafTag>::begin(e);
             it!=endit; ++it)
            {
                //std::cout << "\tnew intersection iterator." << std::endl;

                // if we have a neighbor then we assume there is no boundary (forget interior boundaries)
                // in level assemble treat non-level neighbors as boundary
                if (it->neighbor())
                    {
                        if (levelBoundaryAsDirichlet && it->outside()->level()==e.level())
                            continue;
                        if (!levelBoundaryAsDirichlet)
                            continue;
                    }

                //std::cout << "\t\tsurvived first if statements" << std::endl;

                // determine boundary condition type for this face, initialize with processor boundary
                typename BoundaryConditions::Flags bctypeface = BoundaryConditions::process;

                // handle face on exterior boundary, this assumes there are no interior boundaries
                //if (it.boundary())
                if (!it->neighbor())
                    {
                        //std::cout << "\t\t\tsurvived second if-statements." << std::endl;
                        Dune::GeometryType gtface = it->intersectionSelfLocal().type();
                        for (size_t g = 0; g < Dune::QuadratureRules<DT,n-1>::rule(gtface,p).size(); ++g)
                            {
                                const Dune::FieldVector<DT,n-1>& faceLocalNm1 = Dune::QuadratureRules<DT,n-1>::rule(gtface,p)[g].position();
                                FieldVector<DT,n> local = it->intersectionSelfLocal().global(faceLocalNm1);
                                FieldVector<DT,n> global = it->intersectionGlobal().global(faceLocalNm1);
                                bctypeface = problem.bctype(global,e,local); // eval bctype
                                //std::cout << "\t\t\tlocal = " << local << ", global = " << global << ", bctypeface = " << bctypeface
                                //        << ", size = " << Dune::QuadratureRules<DT,n-1>::rule(gtface,p).size() << std::endl;


                                if (bctypeface!=BoundaryConditions::neumann) break;

                                RT J = problem.neumannPress(global,e,local);
                                double weightface = Dune::QuadratureRules<DT,n-1>::rule(gtface,p)[g].weight();
                                DT detjacface = it->intersectionGlobal().integrationElement(faceLocalNm1);
                                for (int i=0; i<sfs.size(); i++) // loop over test function number
                                    if (this->bctype[i][0]==BoundaryConditions::neumann)
                                        {
                                            this->b[i] -= J*sfs[i].evaluateFunction(0,local)*weightface*detjacface;
                                        }
                            }
                        if (bctypeface==BoundaryConditions::neumann) continue; // was a neumann face, go to next face
                    }

                // If we are here, then it is
                // (i)   an exterior boundary face with Dirichlet condition, or
                // (ii)  a processor boundary (i.e. neither boundary() nor neighbor() was true), or
                // (iii) a level boundary in case of level-wise assemble
                // How processor boundaries are handled depends on the processor boundary mode
                if (bctypeface==BoundaryConditions::process && procBoundaryAsDirichlet==false
                    && levelBoundaryAsDirichlet==false)
                    continue; // then it acts like homogeneous Neumann

                // now handle exterior or interior Dirichlet boundary
                for (int i=0; i<sfs.size(); i++) // loop over test function number
                    {
                        if (sfs[i].codim()==0) continue; // skip interior dof
                        if (sfs[i].codim()==1) // handle face dofs
                            {
                                if (sfs[i].entity()==it->numberInSelf())
                                    {
                                        if (this->bctype[i][0]<bctypeface)
                                            {
                                                this->bctype[i][0] = bctypeface;
                                                if (bctypeface==BoundaryConditions::process)
                                                    this->b[i] = 0;
                                                if (bctypeface==BoundaryConditions::dirichlet)
                                                    {
                                                        Dune::FieldVector<DT,n> global = e.geometry().global(sfs[i].position());
                                                        //std::cout << "i = " << i << ", loop 1, global = " << global << std::endl;
                                                        this->b[i] = problem.dirichletPress(global,e,sfs[i].position());
                                                    }
                                            }
                                    }
                                continue;
                            }
                        // handle subentities of this face
                        for (int j=0; j<ReferenceElements<DT,n>::general(gt).size(it->numberInSelf(),1,sfs[i].codim()); j++)
                            if (sfs[i].entity()==ReferenceElements<DT,n>::general(gt).subEntity(it->numberInSelf(),1,j,sfs[i].codim()))
                                {
                                    if (this->bctype[i][0]<bctypeface)
                                        {
                                            this->bctype[i][0] = bctypeface;
                                            if (bctypeface==BoundaryConditions::process)
                                                this->b[i] = 0;
                                            if (bctypeface==BoundaryConditions::dirichlet)
                                                {
                                                    Dune::FieldVector<DT,n> global = e.geometry().global(sfs[i].position());
                                                    //std::cout << "loop 2, global = " << global << std::endl;
                                                    this->b[i] = problem.dirichletPress(global,e,sfs[i].position());
                                                }
                                        }
                                }
                    }
            }
    }

    // parameters given in constructor
    DeprecatedDiffusionProblem<G,RT,VC>& problem;
    bool levelBoundaryAsDirichlet;
    bool procBoundaryAsDirichlet;
    EM elementmapper;
};

/** @} */
}
#endif
