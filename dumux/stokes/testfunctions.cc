// $Id$

// here are some sample problem with known solution

template<int dim, class ct>
class Example : public ExactSolution<ct, dim>
{
  typedef Dune::FieldVector< ct, dim > Point; //global coord
  typedef Dune::FieldVector<ct,dim> Gradient;

public:
   Example(){}

  ct velocity(int comp,const Point & glob) const
    {
      if (comp==0)
    return 2*glob[0]*glob[1]; // x*x
      if (comp==1)
    return -glob[1]*glob[1];// -2xy*

      return 0;// check this 3D prob
    }
  ct pressure(const Point & glob) const
  {
    return glob[0] - 1.0; // x

  }
    ct rhsvalue(int variable, const Point& glob) const
  {
    // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
    if (variable==0)
      return 1.0;
    if (variable==1)
      return 2.0;
    if (variable==2)
      return 0.0;
    return 0.0;

  }
Gradient velocityGradient(int comp,const Point &glob)const
  {
    Gradient result(0);
    if (comp == 0)
      {
    result[0] = 2.0*glob[1];
    result[1] = 2.0*glob[0];
      }
    else
      result[1] = -2.0*glob[1];

    return result;
  }


  virtual ~Example(){}
};


template<int dim, class ct>
class PoiseuilleFlow : public ExactSolution<ct, dim>
{
  typedef Dune::FieldVector< ct, dim > Point;
  typedef Dune::FieldVector<ct,dim> Gradient;

public:
   PoiseuilleFlow(){}

  ct velocity(int comp,const Point & glob) const
    {
      if (comp==0) return glob[1]*(1.0-glob[1]); // y*(1-y)
      if (comp==1) return 0; // 0
          }
  ct pressure(const Point & glob) const
  {
    return -2*glob[0];

  }
  ct rhsvalue(int variable, const Point& glob) const
  {
    // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
    if (variable==0) return 0;
    if (variable==1) return 0;
    if (variable==2) return 0;
  }
Gradient velocityGradient(int comp,const Point &glob)const
  {
    // //DUNE_THROW(NotImplemented, "velocityGradient not implemented yet");
    //return 0;
  }

  virtual ~ PoiseuilleFlow(){}
};


template<int dim, class ct>
class Example1 : public ExactSolution<ct, dim>
{
  typedef Dune::FieldVector< ct, dim > Point;
  typedef Dune::FieldVector<ct,dim> Gradient;

public:
   Example1(){}

  ct velocity(int comp,const Point & glob) const
    {
      if (comp==0) return sin(glob[0]); // sin(x)
      if (comp==1) return -glob[1]*cos(glob[0]);// -y*cos(x)
    }
  ct pressure(const Point & glob) const
  {
    return glob[0]*glob[1]; // x*y

  }
  ct rhsvalue(int variable, const Point& glob) const
  {
    // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
    if (variable==0) return sin(glob[0])+glob[1]; // sin(x)+y
    if (variable==1) return -glob[1]*cos(glob[0])+glob[0];// -y*cos(x)+x
    if (variable==2) return 0.0;
  }
  Gradient velocityGradient(int comp,const Point &glob)const
    {
      // //DUNE_THROW(NotImplemented, "velocityGradient not implemented yet");
      //return 0;
    }

  virtual ~Example1(){}
};


template<int dim, class ct>
class Example2 : public ExactSolution<ct, dim>
{
  typedef Dune::FieldVector< ct, dim > Point;
  typedef Dune::FieldVector<ct,dim> Gradient;

public:
   Example2(){}

  ct velocity(int comp,const Point & glob) const
    {
      if (comp==0) return -std::exp(glob[0])*(glob[1]*cos(glob[1])+sin(glob[1])); // -e^x(y*cos(y) +sin(y))
      if (comp==1) return std::exp(glob[0])*glob[1]*sin(glob[1]);// e^x y*sin(y)
    }
  ct pressure(const Point & glob) const
  {
    return 2*std::exp(glob[0])*sin(glob[1]); //2*e^x*sin(y)

  }
  ct rhsvalue(int variable, const Point& glob) const
  {
    // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
    if (variable==0) return 0.0;
    if (variable==1) return 0.0;
    if (variable==2) return 0.0;
  }
  Gradient velocityGradient(int comp,const Point &glob)const
    {
      // //DUNE_THROW(NotImplemented, "velocityGradient not implemented yet");
      //return 0;
    }

  virtual ~Example2(){}
};



template<int dim,class ct>
class Example2D : public ExactSolution<ct,dim>
{
  typedef Dune::FieldVector<ct,dim> Point;
    typedef Dune::FieldVector<ct,dim> Gradient;
public:
  Example2D(){}
  ct velocity(int comp,const Point &glob)const
  {
    if (comp==0)
      return (Dune::SQR(glob[0])
        *Dune::SQR(1-glob[0])
        *(2*glob[1]-6*Dune::SQR(glob[1])+4*glob[1]*Dune::SQR(glob[1]))); //x^2(1-x)^2(2y-6y^2+4y^3)
    if (comp==1)
      return
        -Dune::SQR(glob[1])
        *Dune::SQR(1-glob[1])
        *(2*glob[0]-6*Dune::SQR(glob[0])+4*glob[0]*Dune::SQR(glob[0]))
        ; // -y^2(1-y)^2(2x-6x^2+4x^3)
    if (comp==2) return 0;
  }
  ct pressure (const Point& glob) const
  {
    return glob[0]*(1-glob[0]) ; //x(1-x)
  }
  ct rhsvalue(int variable, const Point& glob) const
  {
    // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
    if (variable==0)
      return -(
               ((24*glob[1]-12)
               *
                (Dune::SQR(glob[0])*Dune::SQR(glob[0])-2*glob[0]*Dune::SQR(glob[0])+Dune::SQR(glob[0])))
               +
               ((2*glob[1]-6*Dune::SQR(glob[1])+4*glob[1]*Dune::SQR(glob[1]))
               *
                (12*Dune::SQR(glob[0])-12*glob[0]+2))
               )
        +
        (1-2*glob[0]);
    if (variable==1)
      return
        ((24*glob[0]-12)
        *
        (
         Dune::SQR(glob[1])*Dune::SQR(glob[1])-2*glob[1]*Dune::SQR(glob[1])+Dune::SQR(glob[1])
         ))
        +
        ((2*glob[0]-6*Dune::SQR(glob[0])+4*glob[0]*Dune::SQR(glob[0]))
        *
         (12*Dune::SQR(glob[1])-12*glob[1]+2))
        ;
    if (variable==2) return 0;
}
Gradient velocityGradient(int comp,const Point &glob)const
  {
    // //DUNE_THROW(NotImplemented, "velocityGradient not implemented yet");
//         return 0;
  }


  virtual ~ Example2D(){}
};


template<int dim,class ct>
class Example2DNew : public ExactSolution<ct,dim>
{
  typedef Dune::FieldVector<ct,dim> Point;
  typedef Dune::FieldVector<ct,dim> Gradient;

public:
  Example2DNew(){}
  ct velocity(int comp,const Point &glob)const
  {
    if (comp==0)
      return
        (Dune::SQR(glob[0])*Dune::SQR(glob[0])-2*glob[0]*Dune::SQR(glob[0])+Dune::SQR(glob[0]))
        *(4*glob[1]*Dune::SQR(glob[1])-6*Dune::SQR(glob[1])+2*glob[1]);//(x^4-2x^3+x^2)(4y^3-6y2+2y)
    if (comp==1)
      return

        -(4*glob[0]*Dune::SQR(glob[0])-6*Dune::SQR(glob[0])+2*glob[0])
        *(Dune::SQR(glob[1])*Dune::SQR(glob[1])-2*glob[1]*Dune::SQR(glob[1])+Dune::SQR(glob[1]));
         // -(4x^3-6x^2+2x)(y^4-2y^3+y^2)

  }

  Gradient velocityGradient(int comp,const Point &glob)const
  {
    Gradient gV;

    if (comp==0)
      {
        gV[0]=(4*glob[1]*Dune::SQR(glob[1])-6*Dune::SQR(glob[1])+2*glob[1])
          *(4*glob[0]*Dune::SQR(glob[0])-6*Dune::SQR(glob[0])+2*glob[0]); // dudx

         gV[1]=(Dune::SQR(glob[0])*Dune::SQR(glob[0])-2*glob[0]*Dune::SQR(glob[0])+Dune::SQR(glob[0]))
               *(2-12*glob[1]+12*Dune::SQR(glob[1])); // dudy
        return gV;
      }

    if (comp==1)
      {
        gV[0]=-(2-12*glob[0]+12*Dune::SQR(glob[0]))
          *(Dune::SQR(glob[1])*Dune::SQR(glob[1])-2*glob[1]*Dune::SQR(glob[1])+Dune::SQR(glob[1]));// dvdx

         gV[1]=-(4*glob[0]*Dune::SQR(glob[0])-6*Dune::SQR(glob[0])+2*glob[0])
          *(4*glob[1]*Dune::SQR(glob[1])-6*Dune::SQR(glob[1])+2*glob[1]); // dvdy


        return gV;
      }

  }

  ct pressure (const Point& glob) const
  {
    return (glob[0]+glob[1]) ; //x+y
  }
  ct rhsvalue(int variable, const Point& glob) const
  {
    // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
    if (variable==0)
      return -(
               ((24*glob[1]-12)
               *
                (Dune::SQR(glob[0])*Dune::SQR(glob[0])-2*glob[0]*Dune::SQR(glob[0])+Dune::SQR(glob[0])))
               +
               ((2*glob[1]-6*Dune::SQR(glob[1])+4*glob[1]*Dune::SQR(glob[1]))
               *
                (12*Dune::SQR(glob[0])-12*glob[0]+2))
               )
        +1;
    if (variable==1)
      return
        ((24*glob[0]-12)
         *
         (
          Dune::SQR(glob[1])*Dune::SQR(glob[1])-2*glob[1]*Dune::SQR(glob[1])+Dune::SQR(glob[1])
          ))
        +
        ((2*glob[0]-6*Dune::SQR(glob[0])+4*glob[0]*Dune::SQR(glob[0]))
         *
         (12*Dune::SQR(glob[1])-12*glob[1]+2))
        +1;
    if (variable==2) return 0;
  }
  virtual ~ Example2DNew(){}
};



template<int dim,class ct>
class Example3D : public ExactSolution<ct,dim>
{
  typedef Dune::FieldVector<ct,dim> Point;
  typedef Dune::FieldVector<ct,dim> Gradient;
public:
  Example3D(){}
  ct velocity(int comp,const Point &glob)const
  {
    if (comp==0)
      //x^2*(1-x)^2(2y-6y^2+4y^3)(2z-6z^2+4z^3)
      return
        Dune::SQR(glob[0])
        *Dune::SQR(1-glob[0])
        *(2*glob[1]-6*Dune::SQR(glob[1])+4*glob[1]*Dune::SQR(glob[1]))
        *(2*glob[2]-6*Dune::SQR(glob[2])+4*glob[2]*Dune::SQR(glob[2]));
    if (comp==1)
      //y^2*(1-y)^2(2x-6x^2+4x^3)(2z-6z^2+4z^3)
      return
          Dune::SQR(glob[1])
        *Dune::SQR(1-glob[1])
        *(2*glob[0]-6*Dune::SQR(glob[0])+4*glob[0]*Dune::SQR(glob[0]))
        *(2*glob[2]-6*Dune::SQR(glob[2])+4*glob[2]*Dune::SQR(glob[2]));
    if (comp==2)
      //-2^z*(1-z)^2(2x-6x^2+4x^3)(2y-6y^2+4y^3)
      return
        -2*Dune::SQR(glob[2])
        *Dune::SQR(1-glob[2])
        *(2*glob[0]-6*Dune::SQR(glob[0])+4*glob[0]*Dune::SQR(glob[0]))
        *(2*glob[1]-6*Dune::SQR(glob[1])+4*glob[1]*Dune::SQR(glob[1]));

  }
  ct pressure (const Point& glob) const
  {
    return Dune::SQR(glob[0])+Dune::SQR(glob[1])+Dune::SQR(glob[2]); //x^2+y^2+z^2
  }
  ct rhsvalue(int variable, const Point& glob) const
  {
    // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
    if (variable==0)
      return
        -(Dune::SQR(glob[0])*Dune::SQR(glob[0])-2*glob[0]*Dune::SQR(glob[0])+Dune::SQR(glob[0]))
        *
        (
         ((2*glob[1]-6*Dune::SQR(glob[1])+4*glob[1]*Dune::SQR(glob[1]))
         *
         (-12+24*glob[2]))
         +
         ((2*glob[2]-6*Dune::SQR(glob[2])+4*glob[2]*Dune::SQR(glob[2]))*
          (-12+24*glob[1]))
         )
        +(
        (2*glob[2]-6*Dune::SQR(glob[2])+4*glob[2]*Dune::SQR(glob[2]))
        *
        (2*glob[1]-6*Dune::SQR(glob[1])+4*glob[1]*Dune::SQR(glob[1]))
        *
        (2-12*glob[0]+12*Dune::SQR(glob[0])))
        +
        2*glob[0];
    if (variable==1)
      return
        -(Dune::SQR(glob[1])*Dune::SQR(glob[1])-2*glob[1]*Dune::SQR(glob[1])+Dune::SQR(glob[1]))
        *
        (
         ((2*glob[0]-6*Dune::SQR(glob[0])+4*glob[0]*Dune::SQR(glob[0]))
         *
          (-12+24*glob[2]))
         +
         ((2*glob[2]-6*Dune::SQR(glob[2])+4*glob[2]*Dune::SQR(glob[2]))
          *(-12+24*glob[0]))
         )
        +
        ((2*glob[2]-6*Dune::SQR(glob[2])+4*glob[2]*Dune::SQR(glob[2]))
        *
        (2*glob[0]-6*Dune::SQR(glob[0])+4*glob[0]*Dune::SQR(glob[0]))
        *
         (2-12*glob[1]+12*Dune::SQR(glob[1])))
        +
        2*glob[1];
    if (variable==2)
      return
        2*(Dune::SQR(glob[2])*Dune::SQR(glob[2])-2*glob[2]*Dune::SQR(glob[2])+Dune::SQR(glob[2]))
        *(
          ((2*glob[1]-6*Dune::SQR(glob[1])+4*glob[1]*Dune::SQR(glob[1]))
           * (-12+24*glob[0]))
          +
          ((2*glob[0]-6*Dune::SQR(glob[0])+4*glob[0]*Dune::SQR(glob[0]))
           *(-12+24*glob[1]))
         )
        +
        2*(2*glob[0]-6*Dune::SQR(glob[0])+4*glob[0]*Dune::SQR(glob[0]))
        *(2*glob[1]-6*Dune::SQR(glob[1])+4*glob[1]*Dune::SQR(glob[1]))
        *
        (2-12*glob[2]+12*Dune::SQR(glob[2]))

        +
        2*glob[2];
    if (variable==3)
          return 0;

  }

      Gradient velocityGradient(int comp,const Point &glob)const
  {
    // //DUNE_THROW(NotImplemented, "velocityGradient not implemented yet");
//         return 0;
  }
  virtual ~ Example3D(){}
};


template<int dim,class ct>
class DrivenCavity2D  : public ExactSolution<ct,dim>
{
  typedef Dune::FieldVector<ct,dim> Point;
  typedef Dune::FieldVector<ct,dim> Gradient;
public:
  DrivenCavity2D(){}
  ct velocity(int comp,const Point &glob)const
  {
    if (glob[1]>1-1E-8)
      {
        if (comp==0)
          return 1;
        if (comp==1)
          return 0;
      }
    else
      {
        if (comp==0)
          return 0;
        if (comp==1)
          return 0;

      }

  }

Gradient velocityGradient(int comp,const Point &glob)const
  {
    // //DUNE_THROW(NotImplemented, "velocityGradient not implemented yet");
//         return 0;
  }


  ct pressure (const Point& glob) const
  {
    return 0;
  }

 ct rhsvalue(int variable, const Point& glob) const
  {
    // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
    if (variable==0) return 0.0;
    if (variable==1) return 0.0;
    if (variable==2) return 0.0;
  }
  virtual ~ DrivenCavity2D(){}
};

template<int dim,class ct>
class DrivenCavity3D  : public ExactSolution<ct,dim>
{
  typedef Dune::FieldVector<ct,dim> Point;
  typedef Dune::FieldVector<ct,dim> Gradient;
public:
  DrivenCavity3D(){}
  ct velocity(int comp,const Point &glob)const
  {
    if (glob[dim-1]>1.0-1E-8)
      {
        if (comp==0)
          return 1;
        else
          return 0;
      }
    else
      {
        // if (comp==0)
//           return 0;
//         if (comp==1)
          return 0;

      }

  }

Gradient velocityGradient(int comp,const Point &glob)const
  {
    // //DUNE_THROW(NotImplemented, "velocityGradient not implemented yet");
//         return 0;
  }


  ct pressure (const Point& glob) const
  {
    return 0;
  }

 ct rhsvalue(int variable, const Point& glob) const
  {
    // stokes system has dim+1 variables (dim velocity comps and 1 pressure)
    if (variable==0) return 0.0;
    if (variable==1) return 0.0;
    if (variable==2) return 0.0;
    if (variable==3) return 0.0;
  }
  virtual ~ DrivenCavity3D(){}
};
