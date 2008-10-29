// $Id$

template<class G, int v_order, int p_order>
void DGFiniteElementMethod<G,v_order,p_order>::assembleVolumeTerm(Entity& ent, LocalMatrixBlock& Aee,LocalVectorBlock& Be) const
{
  Gradient grad_phi_ei[dim],grad_phi_ej[dim],temp;
  ctype  psi_ei,psi_ej;
  ctype entry;
  
  //get the shape function set
  //for  velocity
  ShapeFunctionSet vsfs(v_order);
  // for pressure
  ShapeFunctionSet psfs(p_order); 
 
  //shape function size and total dof
  int vdof=vsfs.size()*dim; // dim velocity components 

  //get the geometry type
  Dune::GeometryType gt = ent.type();
  //specify the quadrature order ?
  //  #warning fixed quadrature order 
  int qord=6;
  for (unsigned int nqp=0;nqp<Dune::QuadratureRules<ctype,dim>::rule(gt,qord).size();++nqp) 
	{
	  //local position of quad points
	  const Dune::FieldVector<ctype,dim> & quad_point_loc = Dune::QuadratureRules<ctype,dim>::rule(gt,qord)[nqp].position();
	  //global position
	  Dune::FieldVector<ctype,dim> quad_point_glob = ent.geometry().global(quad_point_loc); 
	  // calculate inv jacobian
	  InverseJacobianMatrix inv_jac=ent.geometry().jacobianInverseTransposed(quad_point_loc);
	  // quadrature weight 
	  double quad_wt=Dune::QuadratureRules<ctype,dim>::rule(gt,qord)[nqp].weight();
	  // get the determinant jacobian
	  ctype detjac=ent.geometry().integrationElement(quad_point_loc);

	  ctype rhsval[dim+1];

	  
	  //================================================//
	  // source term: TERM:14 : f* v
	  // TERM 1 : \int (mu*grad_u*grad_v)
	  //================================================//	  
	  // dim velocity comps. 
	  for(int dm=1;dm<=dim;++dm) 
		{
		  for (int i=0;i<vsfs.size();++i) 
			{
			  //space dimension-sd  
			  for (int sd=0; sd<dim; sd++)
				temp[sd] = vsfs[i].evaluateDerivative(0,sd,quad_point_loc);
			  grad_phi_ei[dm-1] = 0;
			  //matrix vect multiplication
			  // transform gradient to global coordinates by multiplying with inverse jacobian
			  inv_jac.umv(temp,grad_phi_ei[dm-1]);
			  int ii=(dm-1)*vsfs.size()+i;
			  // get the rhs value
			  rhsval[dm-1] = (problem_.q(quad_point_glob, ent, quad_point_loc))[dm-1];
			  Be[ii]+=rhsval[dm-1]*vsfs[i].evaluateFunction(0,quad_point_loc)*detjac*quad_wt;
			  for (int j=0;j<vsfs.size();++j)
				{
				  for (int sd=0; sd<dim; sd++)//space dimension -sd
					temp[sd] = vsfs[j].evaluateDerivative(0,sd,quad_point_loc);
				  grad_phi_ej[dm-1] = 0;
				  inv_jac.umv(temp,grad_phi_ej[dm-1]);
				  int jj=(dm-1)*vsfs.size()+j;
				  entry =parameter.mu*(grad_phi_ei[dm-1]*grad_phi_ej[dm-1])*detjac*quad_wt;
				  Aee[ii][jj]+=entry;
				  
				}
			}
		}
	  //================================================//	
	  // -  p * div v 
	  //================================================//
	  for(int dm=1;dm<=dim;++dm)
		{
		  for (int i=0;i<vsfs.size();++i) 
			{
			  int ii=(dm-1)*vsfs.size()+i;
			  for (int sd=0; sd<dim; sd++)//space dimension -sd
				temp[sd] = vsfs[i].evaluateDerivative(0,sd,quad_point_loc);
			  grad_phi_ei[dm-1] = 0;
			  inv_jac.umv(temp,grad_phi_ei[dm-1]);
			  for (int j=0;j<psfs.size();++j) // pressure shapefns
				{
				  int jj=vdof+j;
				  psi_ej=psfs[j].evaluateFunction(0,quad_point_loc);
				  entry =-(grad_phi_ei[dm-1][dm-1]*psi_ej)*detjac*quad_wt;
				  Aee[ii][jj]+=entry;
				}
			}
		}


	  
	  //================================================//	
	  // 	 -  q* div u 	
	  //================================================//			  	  
	 
		  for (int i=0;i<psfs.size();++i) // pressure shapefns
			{
			  int ii=vdof+i;
			  psi_ei=psfs[i].evaluateFunction(0,quad_point_loc);
			   for(int dm=1;dm<=dim;++dm)
				 {
			  for (int j=0;j<vsfs.size();++j) 
				{
				  int jj=(dm-1)*vsfs.size()+j;
				  for (int sd=0; sd<dim; sd++)//space dimension -sd
					temp[sd] = vsfs[j].evaluateDerivative(0,sd,quad_point_loc);
				  grad_phi_ej[dm-1] = 0;
				  inv_jac.umv(temp,grad_phi_ej[dm-1]);
				  entry =-(grad_phi_ej[dm-1][dm-1]*psi_ei)*detjac*quad_wt;
				  Aee[ii][jj]+=entry;
				}
			}
		}


		  
	  //================================================//	

	} // end of volume term quadrature
  //printmatrix(std::cout,Aee,"Matrix A: ","row");
  //printvector(std::cout,Be,"Volume Be: ","row");
}// end of assemble volume term

 
 
template<class G, int v_order, int p_order>
void DGFiniteElementMethod<G,v_order,p_order>::assembleFaceTerm(Entity& ent, IntersectionIterator& isit,
													 LocalMatrixBlock& Aee, LocalMatrixBlock& Aef,
													 LocalMatrixBlock& Afe, LocalVectorBlock& Be) const
{
  Gradient grad_phi_ei[dim],grad_phi_ej[dim],temp;
  ctype   phi_ei[dim],phi_ej[dim],phi_fi[dim],phi_fj[dim],psi_ei,psi_ej;
  ctype entry;
  //get the shape function set
  //self shape functions
  ShapeFunctionSet vsfs(v_order); //for  velocity
  ShapeFunctionSet psfs(p_order); // for pressure
  //neighbor shape functions
  ShapeFunctionSet nbvsfs(v_order); //for  velocity
 
  //shape function size and total dof
  int vdof=vsfs.size()*dim; // two velocity components and total velocity sfs size
  
  //get parameter
  // DGStokesParameters parameter;
  //get the geometry type of the face
  Dune::GeometryType gtface = isit->intersectionSelfLocal().type();
  //std::cout<<"----gtface: "<<gtface<<std::endl;
  Dune::GeometryType nbgtface = isit->intersectionNeighborLocal().type();
  //specify the quadrature order ?
  // #warning now fixed quadrature order
  int qord=6;
 //  Grid grid;

  
  for (unsigned int qedg=0; qedg<Dune::QuadratureRules<ctype,dim-1>::rule(gtface,qord).size(); ++qedg)
	{	
	  //quadrature position on the edge/face in local=facelocal
	  //note that this is dim entity
	  const Dune::FieldVector<ctype,dim-1>& local = Dune::QuadratureRules<ctype,dim-1>::rule(gtface,qord)[qedg].position();
	  Dune:: FieldVector<ctype,dim> face_self_local = isit->intersectionSelfLocal().global(local);
	  Dune:: FieldVector<ctype,dim> face_neighbor_local = isit->intersectionNeighborLocal().global(local);

	  Dune::FieldVector<ctype,dim> global = isit->intersectionGlobal().global(local);
	  //std::cout<<"local: "<<local<<" face_self_local: "<<face_self_local
 	  //<<"  face_neighbor_local: "<<face_neighbor_local<<"  glob: "<<global<<std::endl;
	  // calculating the inverse jacobian 
	  InverseJacobianMatrix inv_jac= ent.geometry().jacobianInverseTransposed(face_self_local);
	  // get quadrature weight
	  ctype quad_wt_face = Dune::QuadratureRules<ctype,dim-1>::rule(gtface,qord)[qedg].weight();
	  ctype detjacface = isit->intersectionGlobal().integrationElement(local);
	  // get the face normal: unit normal.
	  Dune::FieldVector<ctype,dim> normal = isit->unitOuterNormal(local);
	  ctype norm_e= isit->intersectionGlobal().integrationElement(local);
	  

	  
	  //================================================//	
	  // term to be evaluated : TERM:2
	  //- \mu \int average(\nabla u). normal . jump(v) 
	  //================================================//	
	  // diagonal block
	  // -mu* 0.5 * grad_phi_ei * normal* phi_ej 
	  for(int dm=1;dm<=dim;++dm)
		{
		  for (int i=0;i<vsfs.size();++i) 
			{
			  int ii=(dm-1)*vsfs.size()+i;
			  phi_ei[dm-1] = vsfs[i].evaluateFunction(0,face_self_local);
			  for (int j=0;j<vsfs.size();++j) 
				{
				  int jj=(dm-1)*vsfs.size()+j;
				  for (int sd=0; sd<dim; sd++)
					temp[sd] = vsfs[j].evaluateDerivative(0,sd,face_self_local);
				  grad_phi_ej[dm-1] = 0;
				  // transform gradient to global ooordinates by multiplying with inverse jacobian
				  inv_jac.umv(temp,grad_phi_ej[dm-1]);
				  entry =-0.5 * parameter.mu * ((grad_phi_ej[dm-1]*normal)*phi_ei[dm-1])*detjacface*quad_wt_face;
				  //Aee.add(ii,jj,entry);
				  Aee[ii][jj]+=entry;
				}
			}
		}
	  // offdiagonal entry
	  // mu* 0.5 * grad_phi_ei * normal* phi_fj 
	  for(int dm=1;dm<=dim;++dm)
		{
		  for (int i=0;i<nbvsfs.size();++i) 
			{
			  int ii=(dm-1)*nbvsfs.size()+i; 
			  phi_fi[dm-1] = nbvsfs[i].evaluateFunction(0,face_neighbor_local);
			  for (int j=0;j<vsfs.size();++j) 
				{
				  int jj=(dm-1)*vsfs.size()+j; 
				  for (int sd=0; sd<dim; sd++)
					temp[sd] = vsfs[j].evaluateDerivative(0,sd,face_self_local);
				  grad_phi_ej[dm-1] = 0;
				  inv_jac.umv(temp,grad_phi_ej[dm-1]);
				  entry =  0.5*parameter.mu*((grad_phi_ej[dm-1]*normal)*phi_fi[dm-1])*detjacface*quad_wt_face;
				  //Afe.add(ii,jj,entry);
				  Afe[ii][jj]+=entry;
				}
			}
		}
	  //================================================//	
	  // term to be evaluated TERM:4
	  // \mu \parameter.epsilon .\int average(\nabla v). normal . jump(u)
	  //================================================//	
	  // diagonal term 
	  // mu* 0.5 * parameter.epsilon* phi_ei * grad_phi_ej* normal 
	  for(int dm=1;dm<=dim;++dm)
		{
		  for (int i=0;i<vsfs.size();++i) 
			{
			  int ii=(dm-1)*vsfs.size()+i; 
			  for (int sd=0; sd<dim; sd++)
				temp[sd] = vsfs[i].evaluateDerivative(0,sd,face_self_local);
			  grad_phi_ei[dm-1] = 0;
			  inv_jac.umv(temp,grad_phi_ei[dm-1]);
			  for (int j=0;j<vsfs.size();++j) 
				{
				  int jj=(dm-1)*vsfs.size()+j;
				  phi_ej[dm-1] = vsfs[j].evaluateFunction(0,face_self_local); 
				  entry=  0.5*parameter.mu*parameter.epsilon*(phi_ej[dm-1]*(grad_phi_ei[dm-1]*normal))*detjacface*quad_wt_face;
				  Aee[ii][jj]+=entry;
				}
			}
		}
	  // offdiagonal block 
	  // -mu* 0.5 * parameter.epsilon * grad_phi_ej * normal* phi_fi 
	  for(int dm=1;dm<=dim;++dm)
		{
		  for (int i=0;i<vsfs.size();++i) 
			{
			  int ii=(dm-1)*vsfs.size()+i; 
			  for (int sd=0; sd<dim; sd++)
				temp[sd] = vsfs[i].evaluateDerivative(0,sd,face_self_local);
			  grad_phi_ei[dm-1] = 0;
			  inv_jac.umv(temp,grad_phi_ei[dm-1]);
			  // note test fns are now from neighbor element
			  for (int j=0;j<nbvsfs.size();++j) 
				{
				  int jj=(dm-1)*nbvsfs.size()+j;
				  phi_fj[dm-1] = nbvsfs[j].evaluateFunction(0,face_neighbor_local); 
				  entry =-0.5*parameter.mu*parameter.epsilon*( phi_fj[dm-1]*(grad_phi_ei[dm-1]*normal))*detjacface*quad_wt_face;
				  Aef[ii][jj]+=entry;
				}
			}
		}
	  
	  //================================================//	
	  // term to be evaluated TERM:6
	  //  term J0 =  \mu. (\parameter.sigma/norm(e)). jump(u). jump (v)
	  //================================================//	
	  // Diagonalblock :
	  // \mu. (\parameter.sigma/abs(e)).\int phi_ie. phi_je  where abs(e) =  norm (e) ???
	  for(int dm=1;dm<=dim;++dm)
		{
		  for (int i=0;i<vsfs.size();++i) 
			{ 
			  int ii=(dm-1)*vsfs.size()+i; 
			  phi_ei[dm-1] = vsfs[i].evaluateFunction(0,face_self_local); 
			  for (int j=0;j<vsfs.size();++j) 
				{
				  int jj=(dm-1)*vsfs.size()+j;
				  phi_ej[dm-1] = vsfs[j].evaluateFunction(0,face_self_local); 	 
				  entry=parameter.mu*(parameter.sigma/norm_e)*phi_ei[dm-1]*phi_ej[dm-1]*detjacface*quad_wt_face;
				  Aee[ii][jj]+=entry;
				}
			}
		}
	  // offdiagonal block
	  //- mu*(parameter.sigma/norm_e)*phi_ei*phi_fj*detjacface*quad_wt_face;
	  for(int dm=1;dm<=dim;++dm)
		{
		  for (int i=0;i<vsfs.size();++i) 
			{
			  phi_ei[dm-1] = vsfs[i].evaluateFunction(0,face_self_local);
			  int ii=(dm-1)*vsfs.size()+i; 
			  for (int j=0;j<nbvsfs.size();++j) // neighbor basis
				{
				  int jj=(dm-1)*nbvsfs.size()+j;
				  phi_fj[dm-1] = nbvsfs[j].evaluateFunction(0,face_neighbor_local);
				  entry=-parameter.mu*(parameter.sigma/norm_e)*phi_ei[dm-1]*phi_fj[dm-1]*detjacface*quad_wt_face;
				  Aef[ii][jj]+=entry;
				}
			}
		}
	  //================================================//	
	  // term to be evaluated TERM:9
	  //edge  term from B(v,p)
	  // term \int average(p). jump(v).normal
	  //================================================//	
	  //diagonal block
	  // term==  0.5 * psi_ei. phi_ej* normal
	  for(int dm=1;dm<=dim;++dm)
		{
		  for (int i=0;i<vsfs.size();++i) 
			{
			  int ii=(dm-1)*vsfs.size()+i;
			  phi_ei[dm-1] = vsfs[i].evaluateFunction(0,face_self_local);
			  for (int j=0;j<psfs.size();++j) 
				{
				  int jj=vdof+j;
				  psi_ej = psfs[j].evaluateFunction(0,face_self_local); 	
				  entry =0.5*(phi_ei[dm-1]*psi_ej*normal[dm-1])*detjacface*quad_wt_face;
				  Aee[ii][jj]+=entry;
				}
			}
		}
	  
	  //offdiagonal block
	  // term==  -0.5 * psi_ei. phi_fj* normal
	  for(int dm=1;dm<=dim;++dm)
		{
		  for (int i=0;i<nbvsfs.size();++i) 
			{
			  int ii=(dm-1)*nbvsfs.size()+i; 
			  phi_fi[dm-1] = nbvsfs[i].evaluateFunction(0,face_neighbor_local); 
			  for (int j=0;j<psfs.size();++j) 
				{
				  int jj=vdof+j;
				  psi_ej = psfs[j].evaluateFunction(0,face_self_local); 
				  entry = -0.5*(phi_fi[dm-1]*psi_ej*normal[dm-1])*detjacface*quad_wt_face;
				  Afe[ii][jj]+=entry;
				}
			}
		}
	  
	  
	  //================================================//	
	  // term to be evaluated TERM:12
	  //edge  term from B(q,u)
	  // term \int average(q). jump(u).normal
	  // TERM:12 
	  //================================================//	
	  //diagonal block
	  // term==  0.5 * psi_ej. phi_ei* normal
	  
		  for (int i=0;i<psfs.size();++i) 
			{
			  int ii=vdof+i;
			  psi_ei = psfs[i].evaluateFunction(0,face_self_local);
			  for(int dm=1;dm<=dim;++dm)
				{
			  for (int j=0;j<vsfs.size();++j) 
				{
				  int jj=(dm-1)*vsfs.size()+j;
				  phi_ej[dm-1] = vsfs[j].evaluateFunction(0,face_self_local);
				  entry =0.5*(phi_ej[dm-1]*psi_ei*normal[dm-1])*detjacface*quad_wt_face;
				  Aee[ii][jj]+=entry;
				}
			}
		}
	  
	  //offdiagonal block
	  // term==  -0.5 * psi_ej. phi_fi* normal
	  
	 
		  for (int i=0;i<psfs.size();++i) 
			{
			  int ii=vdof+i;
			  psi_ei = psfs[i].evaluateFunction(0,face_self_local);
			   for(int dm=1;dm<=dim;++dm)
				 {
			  for (int j=0;j<nbvsfs.size();++j) // neighbor
				{
				  phi_fj[dm-1] = nbvsfs[j].evaluateFunction(0,face_neighbor_local); 
				  int jj=(dm-1)*nbvsfs.size()+j; 
				  entry = -0.5*(phi_fj[dm-1]*psi_ei*normal[dm-1])*detjacface*quad_wt_face;
				  Aef[ii][jj]+=entry;
				}
			}
		}
	  
	  //================================================//
	  
	}// end of assemble face quadrature loop

 //  printmatrix(std::cout,Aee,"Matrix Aee: ","row");
//   printmatrix(std::cout,Aef,"Matrix Aef: ","row");
//   printmatrix(std::cout,Afe,"Matrix Afe: ","row");
	
  
}// end of assemble face term



template<class G, int v_order, int p_order>
void DGFiniteElementMethod<G,v_order,p_order>::assembleDirichletBoundaryTerm(Entity& ent, IntersectionIterator& isit, 
									     LocalMatrixBlock& Aee, LocalVectorBlock& Be) const
{
  Gradient grad_phi_ei[dim],grad_phi_ej[dim],temp;
  ctype   phi_ei[dim],phi_ej[dim],psi_ei,psi_ej;
  ctype entry;
  ctype dirichlet[dim+1]; // dim velocity and 1 pressure
  //get the shape function set
  //self shape functions
  ShapeFunctionSet vsfs(v_order);; //for  velocity
  ShapeFunctionSet psfs(p_order); // for pressure
  //neighbor shape functions
   
  //shape function size and total dof
  int vdof=vsfs.size()*dim; // two velocity components and total velocity sfs size
  
  //get parameter
  //DGStokesParameters parameter;

  //get the geometry type of the face
  Dune::GeometryType gtboundary = isit->intersectionSelfLocal().type();
  
  //specify the quadrature order ?
  int qord=6;
  for(unsigned int bq=0;bq<Dune::QuadratureRules<ctype,dim-1>::rule(gtboundary,qord).size();++bq)
    {
      const Dune::FieldVector<ctype,dim-1>& boundlocal = Dune::QuadratureRules<ctype,dim-1>::rule(gtboundary,qord)[bq].position();
      Dune:: FieldVector<ctype,dim> blocal = isit->intersectionSelfLocal().global(boundlocal);
      const Dune::FieldVector<ctype,dim> bglobal = isit->intersectionGlobal().global(boundlocal);
      ctype norm_eb=isit->intersectionGlobal().integrationElement(boundlocal);  
      // calculating the inverse jacobian 
      InverseJacobianMatrix inv_jac= ent.geometry().jacobianInverseTransposed(blocal);
      // get quadrature weight
      ctype quad_wt_bound = Dune::QuadratureRules<ctype,dim-1>::rule(gtboundary,qord)[bq].weight();
      ctype detjacbound = isit->intersectionGlobal().integrationElement(boundlocal);
      // get the boundary normal 
      Dune::FieldVector<ctype,dim> boundnormal = isit->unitOuterNormal(boundlocal);
      // velocity boundary condition
      // dirichlet boundary
       for(int i=0;i<dim;++i)
	 dirichlet[i] = (problem_.g(bglobal, ent, isit, blocal))[i];
      
	  
      //================================================//
      // 
      // TERM:3
      //- (\mu \int \nabla u. normal . v)  
      //================================================//

      for(int dm=1;dm<=dim;++dm)
	{
	  
	  for (int i=0;i<vsfs.size();++i) 
	    {
	      int ii=(dm-1)*vsfs.size()+i; 
	      phi_ei[dm-1] = vsfs[i].evaluateFunction(0,blocal);
	      for (int j=0;j<vsfs.size();++j) 
		{
		  int jj=(dm-1)*vsfs.size()+j;
		  for (int sd=0; sd<dim; sd++)
		    temp[sd] = vsfs[j].evaluateDerivative(0,sd,blocal);
		  grad_phi_ej[dm-1] = 0;
		  inv_jac.umv(temp,grad_phi_ej[dm-1]);
		  entry = ( - parameter.mu * ((grad_phi_ej[dm-1]*boundnormal)*phi_ei[dm-1])) * detjacbound*quad_wt_bound;
		  Aee[ii][jj]+=entry;
		}
	    }
	}
      //================================================//				  
      //TERM:5=  \mu parameter.epsilon \nabla v . normal. u
      //  TERM:15
      // rhs entry:  parameter.mu * parameter.epsilon* g * \nabla v * n
      //================================================//				  
      for(int dm=1;dm<=dim;++dm)
	{
	  
	  for (int i=0;i<vsfs.size();++i) 
	    {
	      int ii=(dm-1)*vsfs.size()+i; 
	      for (int sd=0; sd<dim; sd++)
		temp[sd] = vsfs[i].evaluateDerivative(0,sd,blocal);
	      grad_phi_ei[dm-1] = 0;
	      inv_jac.umv(temp,grad_phi_ei[dm-1]);
	      for (int j=0;j<vsfs.size();++j) 
		{
		  int jj=(dm-1)*vsfs.size()+j;
		  phi_ej[dm-1] = vsfs[j].evaluateFunction(0,blocal);
		  //TERM:5 \mu parameter.epsilon \nabla v . normal. u		 
		  entry = parameter.mu *(parameter.epsilon*(grad_phi_ei[dm-1]*boundnormal)*phi_ej[dm-1] ) * detjacbound*quad_wt_bound;
		  Aee[ii][jj]+=entry;
		}
	      //------------------------------------
	      //  TERM:15
	      // rhs entry:  parameter.mu * parameter.epsilon* g * \nabla v * n
	      //------------------------------------
	      
	      
	      Be[ii]+= (parameter.epsilon*parameter.mu*(dirichlet[dm-1])*(grad_phi_ei[dm-1]*boundnormal)) * detjacbound * quad_wt_bound;
	    }
	}
      
      //================================================//
      //  TERM:7
      // + \mu parameter.sigma/norm_e . v . u
      // TERM:16 
      // rhs entry: mu*parameter.sigma/norm_e * g * v 
      //================================================//			  
      for(int dm=1;dm<=dim;++dm)
	{
	  
	  for (int i=0;i<vsfs.size();++i) 
	    {
	      
	      phi_ei[dm-1] =  vsfs[i].evaluateFunction(0,blocal);
	      int ii=(dm-1)*vsfs.size()+i; 
	      for (int j=0;j<vsfs.size();++j) 
		{
		  
		  int jj=(dm-1)*vsfs.size()+j;
		  phi_ej[dm-1] = vsfs[j].evaluateFunction(0,blocal);
		  entry = ((parameter.mu*(parameter.sigma/norm_eb)*phi_ej[dm-1]*phi_ei[dm-1]))* detjacbound*quad_wt_bound;
		  Aee[ii][jj]+=entry;
		}
	      //------------------------------------
	      // TERM:16 
	      // rhs entry: mu*parameter.sigma/norm_e * g * v 
	      //------------------------------------
	      Be[ii]+= (parameter.mu*(parameter.sigma/norm_eb)*(dirichlet[dm-1])*phi_ei[dm-1])* detjacbound * quad_wt_bound;
	      
	    }
	  
	}
      
      //================================================//
      // TERM:10
      // 	  \int p v n
      //================================================//			  
      for(int dm=1;dm<=dim;++dm)
	{
	  for (int i=0;i<vsfs.size();++i) 
	    {
	      int ii=(dm-1)*vsfs.size()+i;
	      phi_ei[dm-1] = vsfs[i].evaluateFunction(0,blocal);
	      for (int j=0;j<psfs.size();++j) 
		{
		  psi_ej = psfs[j].evaluateFunction(0,blocal);
		  int jj=vdof+j;
		  entry= (psi_ej*(phi_ei[dm-1]*boundnormal[dm-1]))* detjacbound * quad_wt_bound;
		  Aee[ii][jj]+=entry;
		}
	    }
	}
      
      //================================================//
      // \int q . u . n  --> TERM:13
      // psi_ej * phi_ei * normal
      //================================================//
      
      
      
      for (int i=0;i<psfs.size();++i) 
	{
	  int ii=vdof+i;
	  psi_ei = psfs[i].evaluateFunction(0,blocal);
	  for(int dm=1;dm<=dim;++dm)
	    {		
	      for (int j=0;j<vsfs.size();++j) 
		{
		  phi_ej[dm-1] = vsfs[j].evaluateFunction(0,blocal);
		  int jj=(dm-1)*vsfs.size()+j;
		  entry= (psi_ei*(phi_ej[dm-1]*boundnormal[dm-1]) )* detjacbound * quad_wt_bound;
		  Aee[ii][jj]+=entry;
		  
		}
	    }
	  
	}
      
      //================================================// 
      //TERM:17 (rhs)
      // \int q . g . n
      //================================================//
      for (int i=0;i<psfs.size();++i) 
	{
	  int ii=vdof+i;
	  psi_ei = psfs[i].evaluateFunction(0,blocal);
	  ctype val=0;
	  for (int dm=0;dm<dim;++dm)
	    {
	      val+=dirichlet[dm]*boundnormal[dm];
	      //std::cout<<"D: "<<dirichlet[dm]<<" BN: "<<boundnormal[dm]<<std::endl;
	    }
	  Be[ii]+=val*psi_ei*detjacbound*quad_wt_bound;
	  //Be[ii]+=(dirichlet[0]*boundnormal[0]+dirichlet[1]*boundnormal[1])*psi_ei*detjacbound*quad_wt_bound;
	  //std::cout<<"ii: "<<ii<<" val: "<<Be[ii]<<std::endl;
	}
    }// end of quadrature loop
  
  //printvector(std::cout,Be,"Vector Be: ","row");
}

template<class G, int v_order, int p_order>
void DGFiniteElementMethod<G,v_order,p_order>::assembleNeumannBoundaryTerm(Entity& ent, IntersectionIterator& isit, 
									   LocalMatrixBlock& Aee, LocalVectorBlock& Be) const
{
  Gradient temp;
  ctype   phi_ei[dim];
  //get the shape function set
  //self shape functions
  ShapeFunctionSet vsfs(v_order);; //for  velocity
  ShapeFunctionSet psfs(p_order); // for pressure
  //neighbor shape functions
   
  //get the geometry type of the face
  Dune::GeometryType gtboundary = isit->intersectionSelfLocal().type();
  
  //specify the quadrature order ?
  int qord=2;
  for(unsigned int bq=0;bq<Dune::QuadratureRules<ctype,dim-1>::rule(gtboundary,qord).size();++bq)
    {
      const Dune::FieldVector<ctype,dim-1>& boundlocal = Dune::QuadratureRules<ctype,dim-1>::rule(gtboundary,qord)[bq].position();
      Dune:: FieldVector<ctype,dim> blocal = isit->intersectionSelfLocal().global(boundlocal);
      const Dune::FieldVector<ctype,dim> bglobal = isit->intersectionGlobal().global(boundlocal);
      // calculating the inverse jacobian 
      InverseJacobianMatrix inv_jac= ent.geometry().jacobianInverseTransposed(blocal);
      // get quadrature weight
      ctype quad_wt_bound = Dune::QuadratureRules<ctype,dim-1>::rule(gtboundary,qord)[bq].weight();
      ctype detjacbound = isit->intersectionGlobal().integrationElement(boundlocal);
      // get the boundary normal 
      Dune::FieldVector<ctype,dim> boundnormal = isit->unitOuterNormal(boundlocal);
      // normal traction BC 
      ctype normalTraction = problem_.Jn(bglobal, ent, isit, blocal);
      // tangential traction BC 
      Gradient tangentialTraction =  problem_.Jt(bglobal, ent, isit, blocal);
//       std::cout << "x = " << bglobal << ", normalF = " << normalTraction << ", tangentialF = " << tangentialTraction << std::endl;

      //================================================//
      // RHS: - \int p_D v n
      //================================================//			  
      for(int dm=1;dm<=dim;++dm)
	{
	  for (int i=0;i<vsfs.size();++i) 
	    {
	      int ii=(dm-1)*vsfs.size()+i;
	      phi_ei[dm-1] = vsfs[i].evaluateFunction(0,blocal);
	      Be[ii] -= normalTraction*phi_ei[dm-1]*boundnormal[dm-1]* detjacbound * quad_wt_bound;
	    }
	}

      //================================================//
      // RHS: \int g_t.v 
      //================================================//			  
      for(int dm=1;dm<=dim;++dm)
	{
	  for (int i=0;i<vsfs.size();++i) 
	    {
	      phi_ei[dm-1] =  vsfs[i].evaluateFunction(0,blocal);
	      int ii=(dm-1)*vsfs.size()+i; 
	      Be[ii] += (tangentialTraction[dm-1]*phi_ei[dm-1])* detjacbound * quad_wt_bound;
	    }
	}
    }
}

template<class G, int v_order, int p_order>
void DGFiniteElementMethod<G,v_order,p_order>::assembleInterfaceTerm(Entity& ent, IntersectionIterator& isit, 
									   LocalMatrixBlock& Aee, LocalVectorBlock& Be) const
{
  Gradient temp;
  ctype   phi_ei[dim], phi_ej[dim], entry;
  //get the shape function set
  //self shape functions
  ShapeFunctionSet vsfs(v_order);; //for  velocity
  ShapeFunctionSet psfs(p_order); // for pressure
  //neighbor shape functions
   
  //get the geometry type of the face
  Dune::GeometryType gtboundary = isit->intersectionSelfLocal().type();
  
  //specify the quadrature order ?
  int qord=2;
  for(unsigned int bq=0;bq<Dune::QuadratureRules<ctype,dim-1>::rule(gtboundary,qord).size();++bq)
    {
      const Dune::FieldVector<ctype,dim-1>& boundlocal = Dune::QuadratureRules<ctype,dim-1>::rule(gtboundary,qord)[bq].position();
      Dune:: FieldVector<ctype,dim> blocal = isit->intersectionSelfLocal().global(boundlocal);
      const Dune::FieldVector<ctype,dim> bglobal = isit->intersectionGlobal().global(boundlocal);
      // calculating the inverse jacobian 
      InverseJacobianMatrix inv_jac= ent.geometry().jacobianInverseTransposed(blocal);
      // get quadrature weight
      ctype quad_wt_bound = Dune::QuadratureRules<ctype,dim-1>::rule(gtboundary,qord)[bq].weight();
      ctype detjacbound = isit->intersectionGlobal().integrationElement(boundlocal);
      // get the boundary normal 
      Dune::FieldVector<ctype,dim> boundnormal = isit->unitOuterNormal(boundlocal);
      // normal traction BC 
      ctype normalTraction = problem_.Jn(bglobal, ent, isit, blocal);
      // Beavers-Joseph proportionality constant c = sqrt(k)/alpha such that u_t = - c (grad u . n)_t 
      ctype beaversJosephC = problem_.beaversJosephC(bglobal, ent, isit, blocal);

      //================================================//
      // RHS: - \int p_D v n
      //================================================//			  
      for(int dm=1;dm<=dim;++dm)
	{
	  for (int i=0;i<vsfs.size();++i) 
	    {
	      int ii=(dm-1)*vsfs.size()+i;
	      phi_ei[dm-1] = vsfs[i].evaluateFunction(0,blocal);
	      Be[ii] -= normalTraction*phi_ei[dm-1]*boundnormal[dm-1]* detjacbound * quad_wt_bound;
	    }
	}

      //================================================//
      // Beavers-Joseph interface condition 
      // \int 1/c u_t . v
      //================================================//			  
      for(int dm=1;dm<=dim;++dm)
	{
	  for (int i=0;i<vsfs.size();++i) 
	    {
	      phi_ei[dm-1] =  vsfs[i].evaluateFunction(0,blocal);
	      int ii=(dm-1)*vsfs.size()+i; 
	      for (int j=0;j<vsfs.size();++j) 
		{
		  int jj=(dm-1)*vsfs.size()+j;
		  phi_ej[dm-1] = vsfs[j].evaluateFunction(0,blocal);
		  ctype uN = phi_ej[dm-1]*boundnormal[dm-1];
		  ctype uT = phi_ej[dm-1] - uN*boundnormal[dm-1];
		  entry = 1.0/beaversJosephC * (uT*phi_ei[dm-1])* detjacbound*quad_wt_bound;
		  Aee[ii][jj]+=entry;
		}
	    }
	}
    }
}




template<class G, int v_order, int p_order>
void DGStokes<G,v_order,p_order>::assembleStokesSystem()
{
std::cout << "Assembling the matrix and rhs: \n";
  
  ShapeFunctionSet vsfs(v_order);; //for  velocity
  ShapeFunctionSet psfs(p_order); // for pressure
  
  int vdof=vsfs.size()*dim; // dim velocity components and total velocity sfs size
  //int pdof=psfs.size();
  //int ndof=vdof+pdof; // total dofs per element
  //vsfs.print(std::cout);

  // assembling on the finest level / leaf level
  //level=grid.maxLevel();
    //istl matrix
  // N is = (the no of blocks of size=BlockSize) = no of elements
   //BlockSize is ndof;
  // N is now no of elements and each elements contain a block of ndof = BlockSize
  //int N = grid.size(level, 0); // per level
  int N = grid.size(0); // per leaf
  std::cout<<"leaf size: "<<N<<" fine level size: "<<grid.size(level,0)<<std::endl;
  LMatrix tmp(N,N,LMatrix::row_wise);
  typename LMatrix::CreateIterator mit=tmp.createbegin();

  // build up the matrix structure

    
   ElementLeafIterator eit = grid.template leafbegin<0>();
   ElementLeafIterator eitend = grid.template leafend<0>();
  
    // ElementLevelIterator eit = grid.template lbegin<0>(level);
//     ElementLevelIterator eitend = grid.template lend<0>(level);
  
  for (; eit != eitend; ++eit)
    {
      // insert a non zero entry for myself
	  //mit.insert(grid.levelIndexSet(level).index(*eit));
	   mit.insert(grid.leafIndexSet().index(*eit));      
      assert(mit != tmp.createend());

	  //IntersectionLevelIterator endit = eit->ilevelend();
      //IntersectionLevelIterator iit = eit->ilevelbegin();
	   IntersectionIterator endit = eit->ileafend();
       IntersectionIterator iit = eit->ileafbegin();

	  
	  // insert a non zero entry for each neighbour
      for(; iit != endit; ++iit)
		{
		  
		  if (iit->neighbor())
			{
			  //mit.insert(grid.levelIndexSet(level).index(*iit->outside()));
			  mit.insert(grid.leafIndexSet().index(*iit->outside()));
			}
		}
      ++mit;
    }

  
  
  tmp = 0.0;
  A = tmp;
  LVector tmpv(N);
  b = tmpv;
  // b.resize(N, false);
  b = 0.0;
  solution =tmpv;
  solution =0;

  
  
  // loop over all elements or leaf

 
 //  ElementLevelIterator it = grid.template lbegin<0>(level);
//   ElementLevelIterator itend = grid.template lend<0>(level);
  
  ElementLeafIterator it = grid.template leafbegin<0>();
  ElementLeafIterator itend = grid.template leafend<0>();

  
  for (; it != itend; ++it)
    {
      EntityPointer epointer = it;
      //int eid = grid.levelIndexSet(level).index(*epointer);
      int eid = grid.leafIndexSet().index(*epointer);
      
      dgfem.assembleVolumeTerm(*it,A[eid][eid],b[eid]);
      //IntersectionLevelIterator endis = it->ilevelend();
      //   IntersectionLevelIterator is = it->ilevelbegin();
      IntersectionIterator endis = it->ileafend();
      IntersectionIterator is = it->ileafbegin();
      
      for(; is != endis; ++is)
	{ 
	  if(is->neighbor())
	    {
	      int eid = grid.leafIndexSet().index(*is->inside());
	      int fid = grid.leafIndexSet().index(*is->outside());
	      dgfem.assembleFaceTerm(*it,is,A[eid][eid],A[eid][fid],A[fid][eid],b[eid]);
	      
	    }
	  if (is->boundary())
	    {
	      GeometryType gtf = is->intersectionSelfLocal().type();
	      const FieldVector<ctype,dim-1>& faceLocal = ReferenceElements<ctype,dim-1>::general(gtf).position(0,0);
	      FieldVector<ctype,dim> faceGlobal = is->intersectionGlobal().global(faceLocal);
	      const FieldVector<ctype,dim>& faceLocalDim = ReferenceElements<ctype,dim>::general(gtf).position(is->numberInSelf(),1);
	      BoundaryConditions::Flags bctype = dgfem.problem().bctype(faceGlobal, *it, is, faceLocalDim);

	      if (bctype == BoundaryConditions::dirichlet) 
		dgfem.assembleDirichletBoundaryTerm(*it,is,A[eid][eid],b[eid]);	
	      else if (bctype == BoundaryConditions::neumann) 
		dgfem.assembleNeumannBoundaryTerm(*it,is,A[eid][eid],b[eid]);
	      else // ASSUME that we are on a interface to porous media
		dgfem.assembleInterfaceTerm(*it,is,A[eid][eid],b[eid]);
	    }
	}
    }
  
  




 //istl--------------	 

 //printmatrix(std::cout,A,"Matrix A: ","row");

   //printvector(std::cout,b,"Vector b: ","row");
   //modify matrix for introducing pressure boundary condition
 for (typename LMatrix::RowIterator i=A.begin(); i!=A.end(); ++i)
   for (typename LMatrix::ColIterator j=(*i).begin(); j!=(*i).end(); ++j)
	 {
	   if(i.index()==0) // 0'th block
		 {
		   for(int n=0;n<BlockSize;++n)
			 {
			   A[i.index()][j.index()][vdof][n]=0.0;
			   if((j.index()==0) &(n==vdof))
				 {				
				   A[i.index()][j.index()][vdof][n]=1.0;
				 }
			 }
		 }
	 }
//istl--------------	



 
	//istl--------------	
  // applying pressure boundary condition
 // changing rhs entry
 //printvector(std::cout,b,"Vector b: ","row");
 for(typename LVector::iterator i=b.begin();i!=b.end();++i)
   {
	 if(i.index()==0)
	   {
		 b[i.index()][vdof]=0.0;
	   }
   }
 
 //printvector(std::cout,b,"Vector b: ","row");
 //printmatrix(std::cout,A,"Matrix A: ","row");
 //istl--------------	
  

 std::cout<<"Size of the Matrix(in blocks): "<<A.N()<<" X "<<A.M()<<" with blocksize = "<<BlockSize<<std::endl;


 
}// end of assemble



 template<class G, int v_order, int p_order>
void DGStokes<G,v_order,p_order>::solveStokesSystem()
{
 //------------------ISTL solver--------------------------// 
std::cout << "Solving Stokes System using ISTL solver\n";
 
//printmatrix(std::cout,A,"Matrix A: ","row");
//printvector(std::cout,b,"Vector b: ","row");
 std::cout<<"============================================="<<std::endl;

  Dune::MatrixAdapter<LMatrix,LVector,LVector> op(A);
  int maxIterations=1000;
  double reduction=1E-18;
  //Dune::SeqILUn<LMatrix,LVector,LVector> ilu0(A,0,0.92);
  //Dune::SeqILUn<LMatrix,LVector,LVector> ilu0(A,1,1.0);
  Dune::SeqILUn<LMatrix,LVector,LVector> ilu0(A,1,0.92);
  //Dune::SeqGS<Matrix,Vector,Vector> seqgs(A,1,0.92);
  //Dune::BiCGSTABSolver<Vector> bcgsolver(op,ilu0,1E-14,8000,1);

  Dune::BiCGSTABSolver<LVector> bcgsolver(op,ilu0,reduction,maxIterations,0);
  //Dune::CGSolver<LVector> cgsolver(op,ilu0,1E-10,10000,2);
  
  Dune::InverseOperatorResult r;
  solution = 1.0;
  bcgsolver.apply(solution,b,r);
  //cgsolver.apply(solution,b,r);
  std::cout<<"Iterations: "<<r.iterations<<std::endl;
 while (! r.converged && maxIterations < 8000)
   {
	 maxIterations *= 2;
	 std::cout << "warning ... BiCGStab did not converge..."
			   << " increasing maxIterations to "
			   << maxIterations << std::endl;
	 Dune::BiCGSTABSolver<LVector> bcgs(op,ilu0,reduction,maxIterations,0);
	 bcgs.apply(solution,b,r);
	 //    gs.apply(x,b);
   }
  
  //cgsolver.apply(solution,b,r);
  //printvector(std::cout,solution,"Solution","");
   b=solution;//for l2 error calculation
//   for(typename LVector::iterator i=solution.begin();i!=solution.end();++i)
//	 {
//	   for(int n=0;n<VBlockSize;++n)
//		 (*xv)[i.index()][n]=solution[i.index()][n];
//	   for(int n=0;n<PBlockSize;++n)
//		 (*xp)[i.index()][n]=solution[i.index()][VBlockSize+n];
//	   
//	 }
    //printvector(std::cout,solution,"Solution","");
   //printvector(std::cout,(*xv),"Velocity Coeff","");
   //printvector(std::cout,(*xp),"Pressure Coeff","");
    //------------------ISTL solver--------------------------//
}


template <class G, int v_order, int p_order>
inline const typename DGStokes<G,v_order,p_order>::ShapeFunctionSet &
DGStokes<G,v_order,p_order>::getVelocityShapeFunctionSet(const EntityPointer & ep) const
{
  return dgfem.getVelocityShapeFunctionSet(ep->type());
}

template <class G, int v_order, int p_order>
inline const typename DGStokes<G,v_order,p_order>::ShapeFunctionSet &
DGStokes<G,v_order,p_order>::getPressureShapeFunctionSet(const EntityPointer & ep) const
{
  return dgfem.getPressureShapeFunctionSet(ep->type());
}

// template <class G, int v_order, int p_order>
// inline double DGStokes<G,v_order,p_order>::evaluateSolution(const EntityPointer & e,
// 												 const Dune::FieldVector<ctype, dim> & local
// 												 ) const
// {
//   int eid = grid.levelIndexSet(level).index(*e);
//   return dgfem.evaluateSolution(0,*e, local, b[eid]);  
// }



template <class G, int v_order, int p_order>
inline const typename DGFiniteElementMethod<G,v_order,p_order>::ShapeFunctionSet &
DGFiniteElementMethod<G,v_order,p_order>::getVelocityShapeFunctionSet(Dune::GeometryType gt) const
{
  return vspace(gt, v_order);
  }

template <class G, int v_order, int p_order>
inline const typename DGFiniteElementMethod<G,v_order,p_order>::ShapeFunctionSet &
DGFiniteElementMethod<G,v_order,p_order>::getPressureShapeFunctionSet(Dune::GeometryType gt) const
{
  return pspace(gt, p_order);
}


// evaluate value at local coord in e
template <class G, int v_order, int p_order>
double
DGFiniteElementMethod<G,v_order,p_order>::evaluateSolution(int variable,
							   const Entity& element,
							   const Dune::FieldVector< ctype, dim > & coord,
							   const LocalVectorBlock & xe) const
{
  // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
  const ShapeFunctionSet&  vsfs = getVelocityShapeFunctionSet(element.type());
  const ShapeFunctionSet&  psfs = getPressureShapeFunctionSet(element.type());
  int nvsfs = vsfs.size();
  int npsfs = psfs.size();
  
  ctype value[dim+1];
  value[variable]= 0;
  if (variable<dim)
	for (int i=0; i<nvsfs; ++i)
	  {
		int ii=variable*nvsfs+i;
		value[variable] += xe[ii] * vsfs[i].evaluateFunction(0, coord);
	  }
  else if(variable==dim)
	for (int i=0; i<npsfs; ++i)
	  {
		int ii=(dim*nvsfs)+i;
		//std::cout<<"val: "<<value[variable];
		value[variable] += xe[ii] * psfs[i].evaluateFunction(0, coord);
	  //std::cout<<"  xe: "<<xe[ii]<<"  coord:  "<<coord<<"  psfsvalue: "<<psfs[i].evaluateFunction(0, coord)<<std::endl;
		//std::cout<<"value:  "<<value[variable]<<std::endl;
	  }
  return value[variable];
}


// evaluate gradient at local coord in e
template <class G, int v_order, int p_order>
typename DGFiniteElementMethod<G,v_order,p_order>::Gradient
DGFiniteElementMethod<G,v_order,p_order>::evaluateGradient(int variable,
							   const Entity& element,
							   const Dune::FieldVector< ctype, dim > & coord,
							   const LocalVectorBlock & xe) const
{
  // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
  const ShapeFunctionSet&  vsfs = getVelocityShapeFunctionSet(element.type());
  int nvsfs = vsfs.size();
  Gradient grad[dim+1];
  Gradient grad_phi_ei[dim+1];
  InverseJacobianMatrix invj;
  invj=element.geometry().jacobianInverseTransposed(coord);
  grad[variable]=0.0;
   
  if (variable<dim)
	for (int i=0; i<nvsfs; ++i)
	  {
		Gradient temp;
		int ii=variable*nvsfs+i;
		for (int sd=0; sd<dim; sd++)
		  temp[sd] = vsfs[i].evaluateDerivative(0, sd, coord);
    grad_phi_ei[variable] = 0;
    invj.umv(temp,grad_phi_ei[variable]);
    for (int sd=0; sd<dim; sd++)
      grad[variable][sd] += grad_phi_ei[variable][sd] * xe[ii];
	
	  }

  return grad[variable];
}

template<class G,int v_order, int p_order>
void Dune::DGStokes<G,v_order,p_order>::convertToCellData(int variable, BlockVector<FieldVector<double, 1> >& cellData) 
{
	ElementLevelIterator it = grid.template lbegin<0>(level);
	ElementLevelIterator itend = grid.template lend<0>(level);
	for (; it != itend; ++it)
	{
		GeometryType gt = it->geometry().type();
		const FieldVector<ctype,dim>& local = ReferenceElements<ctype,dim>::general(gt).position(0, 0);

		int eid = grid.levelIndexSet(level).index(*it);
		cellData[eid] = dgfem.evaluateSolution(variable, *it, local, b[eid]);
	}

	return;
}

template<class G,int v_order, int p_order>
void Dune::DGStokes<G,v_order,p_order>::vtkout (const char* name, int k) 
{
	VTKWriter<typename G::LevelGridView> 
		vtkwriter(grid.levelView(level));
	typedef Dune::BlockVector<Dune::FieldVector<double, 1> > BlockVector;
	BlockVector pressureCellData(grid.size(0));
	convertToCellData(dim, pressureCellData);
	vtkwriter.addCellData(pressureCellData, "pressure");
	BlockVector velXCellData(grid.size(0));
	convertToCellData(0, velXCellData);
	vtkwriter.addCellData(velXCellData, "velocity x-comp");
	BlockVector velYCellData(grid.size(0));
	convertToCellData(1, velYCellData);
	vtkwriter.addCellData(velYCellData, "velocity y-comp");
	if (dim > 2) {
		BlockVector velZCellData(grid.size(0));
		convertToCellData(2, velZCellData);
		vtkwriter.addCellData(velZCellData, "velocity z-comp");		
	}
	char fname[128];
	sprintf(fname, "%s-%05d", name, k);
	vtkwriter.write(fname, Dune::VTKOptions::ascii);		
}



