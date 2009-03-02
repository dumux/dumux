
/** \todo Please doc me! */

class PressureBoundary{
public:
    // boundaryType 1 -> Dirichlet, 2 -> Neumann, 3 -> FreeFlow
    double operator() (int Id, double t, int& boundaryType) const
    {
        // Left Boundary
        if (Id == 1)
        {
            boundaryType = 2;
            return 0; // m3/(m2*s) for Neumann boundary
        } // Right Boundary
        else if (Id == 2)
        {
            boundaryType = 1;
            return 1.0e+5;
        }
        else
        {
            boundaryType = 99;
            return 0;
        }
    }

    PressureBoundary()
    {
    }
};

/** \todo Please doc me! */

class VelocityBoundary{
public:
    // boundaryType 1 -> Dirichlet, 2 -> Neumann, 3 -> FreeFlow
    double operator() (int Id, double t, int& boundaryType) const
    {
        // Left Boundary
        if (Id == 1)
        {
            boundaryType = 1;
            return 0;
        } // Right Boundary
        else if (Id == 1)
        {
            boundaryType = 1;
            return 0;
        }
        else
        {
            boundaryType = 99;
            return 0;
        }
    }

    VelocityBoundary()
    {
    }
};

/** \todo Please doc me! */

template<class FieldVector>
class ICPressurePipe {
public:

    double operator() (FieldVector x) const
    {
        if ( (x[0] > 2.0 - 1e-5 ) && (x[0] < 2.0 + 1e-5)  )
        {
            return 1.0e5;
        }
        else
            return 1.0e5;
    }

    ICPressurePipe ()
    {

    }
};

/** \todo Please doc me! */

template<class FieldVector>
class ICVelocityPipe {
public:

    double operator() (FieldVector x) const
    {
        if ( (x[0] > 2.0 - 1e-5 ) && (x[0] < 2.0 + 1e-5)  )
        {
            return 9.0;
        }
        else
            return 9.0;
    }

    ICVelocityPipe ()
    {
    }
};

/** \todo Please doc me! */

template<class FieldVector>
class sourceSinkPipe {
public:
    double sourceSinkValue;

    double operator() (FieldVector x, double t) const
    {

        if ( (x[0] > 2.0 - 1e-5 ) && (x[0] < 2.0 + 1e-5)  )
        {
            if (t==0.0)
                return sourceSinkValue;
            else
                if(t<0.5)
                    return sourceSinkValue;
                else
                    return sourceSinkValue;
        }
        else
        {
            return 0;
        }

    }

    sourceSinkPipe()
    {
        sourceSinkValue = 0.0; //unit kg/m3
    }
};
