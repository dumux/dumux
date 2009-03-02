
/** \todo Please doc me! */

class PressureBoundary{
public:
    // boundaryType 1 -> Dirichlet, 2 -> Neumann [kg/(m2*s)], 3 -> FreeFlow
    double operator() (int Id, double t, int& boundaryType) const
    {
        // Left Boundary
        if (Id == 1)
            {
                boundaryType = 1;
                return 1.0e+5 ; // kg/(m2*s) for Neumann boundary
            } // Right Boundary
        else if (Id == 2)
            {
                boundaryType = 1;
                return 1.0e+5 ;
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

template<class GlobalPosition>
class ICPressurePipe {
public:

    double operator() (GlobalPosition globalPos) const
    {
        if ( (globalPos[0] > 2.0 - 1e-5 ) && (globalPos[0] < 2.0 + 1e-5)  )
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

template<class GlobalPosition>
class ICVelocityPipe {
public:

    double operator() (GlobalPosition globalPos) const
    {
        if ( (globalPos[0] > 2.0 - 1e-5 ) && (globalPos[0] < 2.0 + 1e-5)  )
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

template<class GlobalPosition>
class sourceSinkPipe {
public:
    double sourceSinkValue;

    double operator() (GlobalPosition globalPos, double t) const
    {

        if ( (globalPos[0] > 2.0 - 1e-5 ) && (globalPos[0] < 2.0 + 1e-5)  )
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
                return sourceSinkValue;
            }

    }

    sourceSinkPipe()
    {
        sourceSinkValue = 0.0; //unit [m3/(m3*s)]
    }
};
