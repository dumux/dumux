
void TimeloopOptsPipe( double& tstart, double& tend, double& max_dt, double& first_dt, double& CFL_factor, int& flag,
                       int& n_iter, double& max_def, int& modulo, int& stages );

/** \todo Please doc me! */

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class LmbdLocal, class Grid, class VertexMapper, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
class PipeFlow {
public:
    typedef double     Scalar;
    const Grid& grid;             //!< pointer to grid
    VertexMapper& vMap;     //!< pointer to mapper
    MapperNodeGlobalIDtoOnOutIndexType& mapperGlobalNodeIdtoOnPipeNodeIndex;
    MapperNodeGlobalIDtoOnOutIndexType& mapperGlobalNodeIdtoOutPipeNodeIndex;
    VertexVectorOnLineType& vertexVectorOnLine;
    VertexVectorOutLineType& vertexVectorOutLine;
    BCP boundaryConditionP; //!< functor for the boundary condition
    BCV boundaryConditionV; //!< functor for the boundary condition
    ICP initialConditionP;  //!< functor for the initial condition pressure
    ICV initialConditionV;  //!< functor for the initial condition velocity
    SST source;            //!< functor for source/sink term
    Lmbd frictionCoef;        //!< functor for pipe friction coefficient
    LmbdLocal lambdaLocal;  //!< functor for pipe local loss coefficient
    Press pressurePorous;    //!< Pressure of porous media
    Press mobility;
    Scalar alphaExchange;            //!< mass exchange factor for coupling
    Scalar temp;            //!< gravity vector
    Scalar density;            //!< fluid property
    Scalar viscosity;        //!< fluid property
    Scalar roughness;        //!< pipe sand grain roughness
    Scalar diameter;        //!< pipe diameter
    Press pressure;           //!< vector for the values of the pressure
    Press velocity;          //!< block vector for the total velocity
    typedef Dune::FieldVector<Scalar,Grid::dimension>  FieldVector;
    FieldVector gravity;

    void PrintVertexVector ();
    //! sets the vector of unknowns to initial values.
    void SetInitialSolution ();
    // sets Boundary Condition for velocity, pressure boundary condition does not ly on the pressure index
    void SetDrichletBoundary (Scalar t);

    template<class Matrix>
    void MassEquation (unsigned& k, Matrix& A, Press& f, Scalar t, Scalar dt);

    template<class FaceVector>
    void KFace (unsigned k, FaceVector& kFace, Scalar t, Scalar dt);

    template<class FaceVector>
    void VelocityFace( unsigned k, FaceVector kFace, FaceVector& vFace, Scalar t, Scalar dt);

    //! \brief Calculate the pressure.
    void IterationStep_Mass (Scalar t, Scalar dt, Press& pressureIt);

    void Iteration (Scalar t, Scalar dt, Scalar maxdef);

    //! calculate the update vector and estimate the time step size.
    void update(void);

    template <class MultiWriter>
    void addVtkFields(MultiWriter &writer)
    {
        // create pressureVtkOutput vector with 0 value , size=sizeOfBigGridNodes and copy pressure values from pipe pressure
        const int dim = Grid::dimension;
        int sizeOfNodes = vMap.size();

        typedef Dune::BlockVector<Dune::FieldVector<Scalar,dim>  > VelocityVectorDIM;
        //        VelocityVectorDIM &velocityVectorVtkOutput = *writer.template createField<Scalar, dim>(sizeOfNodes);
        VelocityVectorDIM velocityVectorVtkOutput(sizeOfNodes);
        velocityVectorVtkOutput = 0;

        typedef Dune::BlockVector<Dune::FieldVector<Scalar,1>  > VelocityScalar;
        VelocityScalar &velocityVectorVtkOutputX = *writer.template createField<Scalar, 1>(sizeOfNodes);;
        VelocityScalar &velocityVectorVtkOutputY = *writer.template createField<Scalar, 1>(sizeOfNodes);;
        VelocityScalar &velocityVectorVtkOutputZ = *writer.template createField<Scalar, 1>(sizeOfNodes);;
        velocityVectorVtkOutputX = 0;
        velocityVectorVtkOutputY = 0;
        velocityVectorVtkOutputZ = 0;

        Press &pressureVtkOutput = *writer.template createField<Scalar, 1>(sizeOfNodes);;
        Press &isPipeVtkOutput = *writer.template createField<Scalar, 1>(sizeOfNodes);;
        pressureVtkOutput = 0;
        isPipeVtkOutput = 0;
        for (unsigned n = 0; n < vertexVectorOnLine.size(); n++)
            {
                int indexi = n;

                pressureVtkOutput[vertexVectorOnLine[indexi].globalId] = pressure[indexi];
                isPipeVtkOutput[vertexVectorOnLine[indexi].globalId] = 1;

                unsigned numNeighbor = vertexVectorOnLine[indexi].indexVertexVectorOnLine.size();
                typedef Dune::BlockVector<Dune::FieldVector<Scalar,1>  > FaceVector;
                FaceVector kFace(numNeighbor), vFace(numNeighbor);
                kFace = 0;
                vFace = 0;

                Scalar t= 0;
                Scalar dt= 0;
                KFace<FaceVector>(indexi, kFace, t, dt);
                VelocityFace<FaceVector>(indexi, kFace, vFace, t, dt);

                //                      VelocityVectorDIM vFaceVector(numNeighbor);
                //                      vFaceVector = 0;
                typedef Dune::FieldVector<Scalar,dim> FieldVectorDIM;
                FieldVectorDIM averageVelocityVector;
                averageVelocityVector = 0;
                for (unsigned m = 0; m < numNeighbor; m++)
                    {
                        int indexj = vertexVectorOnLine[indexi].indexVertexVectorOnLine[m];

                        FieldVectorDIM unitPD = vertexVectorOnLine[indexi].unitPD(vertexVectorOnLine, indexj);
                        unitPD *= vFace[m];
                        averageVelocityVector +=  unitPD;
                        /*                std::cout << "indexi: " << indexi
                                          << " velocity" << m
                                          << ": "<< vFace[m]
                                          << " kFace" << m
                                          << ": "<< kFace[m]
                                          << std::endl;
                        */
                    }

                averageVelocityVector /= numNeighbor;

                //            printvector(std::cout,averageVelocityVector,"avergaeVelocityVector","row",200,1);
                int globalID = vertexVectorOnLine[indexi].globalId;
                velocityVectorVtkOutput[globalID]= averageVelocityVector;

                velocityVectorVtkOutputX[globalID] = velocityVectorVtkOutput[globalID][0];
                velocityVectorVtkOutputY[globalID] = velocityVectorVtkOutput[globalID][1];
                velocityVectorVtkOutputZ[globalID] = velocityVectorVtkOutput[globalID][2];

            }

        writer.addVertexData(&isPipeVtkOutput,"isPipe");
        writer.addVertexData(&pressureVtkOutput,"pPipe");
        writer.addVertexData(&velocityVectorVtkOutputX,"vPipeX");
        writer.addVertexData(&velocityVectorVtkOutputY,"vPipeY");
        writer.addVertexData(&velocityVectorVtkOutputZ,"vPipeZ");

        return;
    }

    /*
     *  Initializes the vectors \a pressure and \a velocity, the boundary condition \a boundaryCondition,
     *  and the initial condition \a initialCondition.
     *  Sets the pointers \a grid and \a mapper to \a g and \a m.
     */
    PipeFlow(int nNodes,
             const Grid& g,
             VertexMapper& mCDim,
             MapperNodeGlobalIDtoOnOutIndexType& mapGlobalNodeIDtoPipeNodeOnlineIndex,
             MapperNodeGlobalIDtoOnOutIndexType& mapGlobalNodeIDtoPipeNodeOutlineIndex,
             VertexVectorOnLineType& vertexVectorOnL,
             VertexVectorOutLineType& vertexVectorOutL,
             Scalar alpEx,
             Scalar tmp,
             Scalar dens,
             Scalar kinematicViscosity,
             Scalar roughn,
             Scalar diam,
             FieldVector grav)
        : grid(g), vMap(mCDim), mapperGlobalNodeIdtoOnPipeNodeIndex(mapGlobalNodeIDtoPipeNodeOnlineIndex),
          mapperGlobalNodeIdtoOutPipeNodeIndex (mapGlobalNodeIDtoPipeNodeOutlineIndex),
          vertexVectorOnLine(vertexVectorOnL), vertexVectorOutLine (vertexVectorOutL),
          boundaryConditionP(), boundaryConditionV(),
          initialConditionP(), initialConditionV(), source(), lambdaLocal(),
          pressure(nNodes), velocity(nNodes)
    {
        mobility.resize(mCDim.size());
        pressurePorous.resize(mCDim.size());

        alphaExchange = alpEx;
        temp = tmp;
        density = dens;
        viscosity = kinematicViscosity;
        roughness = roughn;
        diameter = diam;
        gravity = grav;
    }
};

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class LmbdLocal, class Grid, class VertexMapper, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, LmbdLocal, Grid, VertexMapper, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::PrintVertexVector()
{

    std::cout << "number of vertices on line = "<< vertexVectorOnLine.size() << std::endl;
    std::cout << "number of vertices out line = "<< vertexVectorOutLine.size() << std::endl;

    for (unsigned k = 0; k < vertexVectorOnLine.size(); k++)
        {
            std::cout << "vertice on line coord: " <<vertexVectorOnLine[k].nodePoint->geometry()[0] << std::endl;
            std::cout << "       vertice on line: " << k << std::endl;
            std::cout << "       globalId: " << vertexVectorOnLine[k].globalId << std::endl;
            std::cout << "       vertexOnLineIndex: " << (mapperGlobalNodeIdtoOnPipeNodeIndex.find(vertexVectorOnLine[k].globalId))->second << std::endl;
            for (unsigned n = 0; n < vertexVectorOnLine[k].parameter.size(); n++)
                {
                    std::cout << "       parameters: " << vertexVectorOnLine[k].parameter[n] << std::endl;
                }
            std::cout << "           boundary: " << vertexVectorOnLine[k].boundary() << std::endl;
            if (vertexVectorOnLine[k].boundary())
                {
                    std::cout << "           normalBF: " << vertexVectorOnLine[k].normalBF(vertexVectorOnLine) << std::endl;
                    std::cout << "           unitPDBF: " << vertexVectorOnLine[k].unitPDBF(vertexVectorOnLine) << std::endl;
                }
            std::cout << "             length: " << vertexVectorOnLine[k].length(vertexVectorOnLine) << std::endl;
            std::cout << "             on line size: " << vertexVectorOnLine[k].lineVectorOnLine.size() << std::endl;
            for (unsigned m = 0; m < vertexVectorOnLine[k].indexVertexVectorOnLine.size(); m++)
                {
                    std::cout << "                neighbour vertices on line: " << vertexVectorOnLine[k].indexVertexVectorOnLine[m] << std::endl;
                    std::cout << "                normal vector: " << vertexVectorOnLine[k].normal(vertexVectorOnLine, vertexVectorOnLine[k].indexVertexVectorOnLine[m]) << std::endl;
                    std::cout << "                unitPD vector: " << vertexVectorOnLine[k].unitPD(vertexVectorOnLine, vertexVectorOnLine[k].indexVertexVectorOnLine[m]) << std::endl;
                }
            std::cout << "             out line size: " << vertexVectorOnLine[k].lineVectorOutLine.size() << std::endl;
            for (unsigned m = 0; m < vertexVectorOnLine[k].indexVertexVectorOutLine.size(); m++)
                {
                    std::cout << "                neighbour vertices out line: " << vertexVectorOnLine[k].indexVertexVectorOutLine[m] << std::endl;
                    std::cout << "                normal vector: " << vertexVectorOnLine[k].normal(vertexVectorOutLine, vertexVectorOnLine[k].indexVertexVectorOutLine[m]) << std::endl;
                }
        }

    for (unsigned k = 0; k < vertexVectorOutLine.size(); k++)
        {
            std::cout << "vertice out line coord: " <<vertexVectorOutLine[k].nodePoint->geometry()[0] << std::endl;
            std::cout << "       vertice out line: " << k << std::endl;
            std::cout << "       globalId: " << vertexVectorOutLine[k].globalId << std::endl;
            std::cout << "       vertexOutLineIndex: " << (mapperGlobalNodeIdtoOutPipeNodeIndex.find(vertexVectorOutLine[k].globalId))->second << std::endl;
            std::cout << "          out line size: " << vertexVectorOutLine[k].lineVectorOutLine.size() << std::endl;
            for (unsigned m = 0; m < vertexVectorOutLine[k].indexVertexVectorOnLine.size(); m++)
                {
                    std::cout << "             neighbour vertices on line: " << vertexVectorOutLine[k].indexVertexVectorOnLine[m] << std::endl;
                }
        }
    return;
}

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class LmbdLocal, class Grid, class VertexMapper, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, LmbdLocal, Grid, VertexMapper, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::SetInitialSolution ()
{
    for (unsigned k = 0; k < vertexVectorOnLine.size(); k++)
        {
            typedef Dune::FieldVector<Scalar,Grid::dimension>  FieldVector;
            FieldVector globalPos = vertexVectorOnLine[k].nodePoint->geometry()[0];
            pressure[k] = initialConditionP(globalPos);
        }
    return;
}

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class LmbdLocal, class Grid, class VertexMapper, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, LmbdLocal, Grid, VertexMapper, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::SetDrichletBoundary (Scalar t)
{

    for (unsigned k = 0; k < vertexVectorOnLine.size(); k++)
        {
            if (vertexVectorOnLine[k].boundary() )
                {
                    int boundaryId = (int) vertexVectorOnLine[k].parameter[1];
                    int boundaryType;
                    Scalar boundaryValue = boundaryConditionP(boundaryId, t, boundaryType);
                    if (boundaryType == 1) // if dirichlet boundary
                        {
                            pressure[k]=boundaryValue;
                        }
                }
        }
    return;
}

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class LmbdLocal, class Grid, class VertexMapper, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
template<class FaceVector>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, LmbdLocal, Grid, VertexMapper, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::VelocityFace(unsigned k, FaceVector kFace, FaceVector& vFace, Scalar t, Scalar dt)
{
    // first we extract the dimensions of the grid
    const int dim = Grid::dimension;

    int indexi = k;

    //      Scalar PI = 3.14;
    //      Scalar crossArea = PI* diameter * diameter / 4;
    typedef Dune::FieldVector<Scalar,dim> FieldVector;
    FieldVector globalPos_i = vertexVectorOnLine[indexi].nodePoint->geometry().corner(0);

    for (unsigned m = 0; m < vertexVectorOnLine[indexi].indexVertexVectorOnLine.size(); m++)
        {
            int indexj = vertexVectorOnLine[indexi].indexVertexVectorOnLine[m];
            const FieldVector &globalPos_j = vertexVectorOnLine[indexj].nodePoint->geometry().corner(0);
            FieldVector distVect = globalPos_j;
            distVect -= globalPos_i;
            Scalar distVal = distVect.two_norm();

            //          FieldVector unitOuterNormal = vertexVectorOnLine[k].normal(vertexVectorOnLine, indexj);
            FieldVector unitPD = vertexVectorOnLine[k].unitPD(vertexVectorOnLine, indexj);

            Scalar sign = distVect * unitPD;
            if (sign>0) sign=1.0;
            else sign=-1.0;
            Scalar densityFace = density;

            vFace[m]= kFace[m]* (-1.0) * (pressure[indexj] - pressure[indexi]) * sign/distVal + kFace[m] * densityFace * (gravity * unitPD);
        }

    // treat Neumann boundary
    if (vertexVectorOnLine[k].boundary() )
        {

        }

    return;
}

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class LmbdLocal, class Grid, class VertexMapper, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
template<class FaceVector>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, LmbdLocal, Grid, VertexMapper, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::KFace(unsigned k, FaceVector& kFace, Scalar t, Scalar dt)
{
    // first we extract the dimensions of the grid
    const int dim = Grid::dimension;

    int indexi = k;

    Scalar length_i =  vertexVectorOnLine[indexi].length(vertexVectorOnLine);
    if (vertexVectorOnLine[indexi].boundary() )
        {
            length_i *= 2;
        }
    Scalar crossArea = M_PI* diameter * diameter / 4;
    typedef Dune::FieldVector<Scalar,dim> FieldVector;
    FieldVector globalPos_i = vertexVectorOnLine[indexi].nodePoint->geometry().corner(0);

    Scalar TaoNoVel_i= 8 * viscosity * density /diameter * (M_PI * diameter) / crossArea; // Tao(no velocity inside) * perimeter / crossarea uniform
    Scalar localLoss_i= lambdaLocal(globalPos_i, t) * density / 2 * crossArea;
    Scalar q_i = source(globalPos_i, t);
    Scalar alphaEXCHANGE_i = density * mobility[vertexVectorOnLine[indexi].globalId] * alphaExchange * (M_PI * diameter * length_i) / diameter;
    Scalar alphaEXCHANGE_i_Prime = alphaEXCHANGE_i / (crossArea * length_i);
    Scalar qex_i = alphaEXCHANGE_i_Prime * (pressure[indexi]- pressurePorous[vertexVectorOnLine[indexi].globalId]);

    q_i = 0;
    qex_i = 0;

    Scalar K_i = 1 / (2*q_i*density - 2*qex_i + TaoNoVel_i + localLoss_i);

    for (unsigned m = 0; m < vertexVectorOnLine[indexi].indexVertexVectorOnLine.size(); m++)
        {
            int indexj = vertexVectorOnLine[indexi].indexVertexVectorOnLine[m];

            Scalar length_j =  vertexVectorOnLine[indexj].length(vertexVectorOnLine);
            if (vertexVectorOnLine[indexj].boundary() )
                {
                    length_j *= 2;
                }
            FieldVector globalPos_j = vertexVectorOnLine[indexj].nodePoint->geometry().corner(0);

            Scalar TaoNoVel_j= 8 * viscosity * density /diameter * (M_PI * diameter) / crossArea; // Tao(no velocity inside) * perimeter / crossarea uniform
            Scalar localLoss_j= lambdaLocal(globalPos_j, t) * density / 2 * crossArea;
            Scalar q_j = source(globalPos_j, t);
            Scalar alphaEXCHANGE_j = density * mobility[vertexVectorOnLine[indexj].globalId] * alphaExchange * (M_PI * diameter * length_j) / diameter;
            Scalar alphaEXCHANGE_j_Prime = alphaEXCHANGE_j / (crossArea * length_j);
            Scalar qex_j = alphaEXCHANGE_j_Prime * (pressure[indexj]- pressurePorous[vertexVectorOnLine[indexj].globalId]);

            q_j = 0;
            qex_j = 0;

            Scalar K_j = 1 / (2*q_j*density - 2*qex_j + TaoNoVel_j + localLoss_j);

            kFace[m]= 2*(K_i*K_j)/(K_i+K_j);

        }

    // treat Neumann boundary
    if (vertexVectorOnLine[k].boundary() )
        {

        }

    return;
}


template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class LmbdLocal, class Grid, class VertexMapper, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
template<class Matrix>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, LmbdLocal, Grid, VertexMapper, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::MassEquation(unsigned& k, Matrix& A, Press& f, Scalar t, Scalar dt)
{
    // first we extract the dimensions of the grid
    const int dim = Grid::dimension;

    int indexi = k;

    //      Scalar PI = 3.14;

    Scalar length =  vertexVectorOnLine[k].length(vertexVectorOnLine);

    Scalar crossArea = M_PI* diameter * diameter / 4;

    typedef Dune::FieldVector<Scalar,dim> FieldVector;
    FieldVector globalPos_i = vertexVectorOnLine[k].nodePoint->geometry().corner(0);

    f[indexi] += source(globalPos_i, t) * density *length*crossArea; // source [m3/(m3*s)]

    Scalar alphaEXCHANGE = density * mobility[vertexVectorOnLine[k].globalId] * alphaExchange * (M_PI * diameter * length) / diameter;
    f[indexi] += alphaEXCHANGE * pressurePorous[vertexVectorOnLine[k].globalId];

    //      Scalar prePor = pressurePorous[vertexVectorOnLine[k].globalId];
    //      std::cout<<"globalID "<< vertexVectorOnLine[k].globalId << " pressurePorous " << prePor <<std::endl;

    A[indexi][indexi] += alphaEXCHANGE * 1.0;

    unsigned numNeighbor = vertexVectorOnLine[k].indexVertexVectorOnLine.size();
    typedef Dune::BlockVector<Dune::FieldVector<Scalar,1>  > FaceVector;
    FaceVector kFace(numNeighbor), vFace(numNeighbor);
    KFace<FaceVector>(k, kFace, t, dt);
    VelocityFace<FaceVector>(k, kFace, vFace, t, dt);

    for (unsigned m = 0; m < vertexVectorOnLine[k].indexVertexVectorOnLine.size(); m++)
        {
            int indexj = vertexVectorOnLine[k].indexVertexVectorOnLine[m];

            FieldVector globalPos_j = vertexVectorOnLine[indexj].nodePoint->geometry().corner(0);
            FieldVector distVect = globalPos_j;
            distVect -= globalPos_i;
            Scalar distVal = distVect.two_norm();

            typedef Dune::FieldVector<Scalar,dim> FieldVector;
            FieldVector unitOuterNormal = vertexVectorOnLine[k].normal(vertexVectorOnLine, indexj);
            FieldVector unitPD = vertexVectorOnLine[k].unitPD(vertexVectorOnLine, indexj);

            Scalar sign = distVect * unitOuterNormal;
            if (sign>0) sign=1.0;
            else sign=-1.0;

            A[indexi][indexi] += 1.0 * density * kFace[m] * crossArea * 1.0 * sign / distVal;
            A[indexi][indexj] += -1.0 * density * kFace[m] * crossArea * 1.0 * sign / distVal;

            Scalar gravitySign = unitPD * unitOuterNormal;
            if (gravitySign>0) gravitySign=1.0;
            else gravitySign=-1.0;
            Scalar densityFace = density;
            f[indexi] += -1.0 * density * (kFace[m] * densityFace * (gravity*unitPD)) * crossArea * gravitySign; // effect of gravity Force on velocity
        }

    // treat Neumann boundary
    if (vertexVectorOnLine[k].boundary() )
        {
            int boundaryId = (int) vertexVectorOnLine[k].parameter[1];
            int boundaryType;
            Scalar boundaryValue = boundaryConditionP(boundaryId, t, boundaryType);
            if (boundaryType == 2) // if neumann flow boundary
                {

                    f[indexi] += -1 * crossArea * boundaryValue; // positiv boundaryValue [kg/(m2*s)]--> out flow,  negativ boundaryValue [kg/(m2*s)] --> in flow   []
                }
        }

    return;
}


template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class LmbdLocal, class Grid, class VertexMapper, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, LmbdLocal, Grid, VertexMapper, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::update(void)
{
    Scalar tstart;
    Scalar tend;
    Scalar max_dt;
    Scalar first_dt;
    Scalar CFL_factor;
    int flag;
    int n_iter;
    Scalar max_def;
    int modulo, stages;

    // get timeloop options
    TimeloopOptsPipe( tstart, tend, max_dt, first_dt, CFL_factor, flag, n_iter, max_def, modulo, stages );

    // initialize concentration with initial values
    //    SetInitialSolution(); // if not coupled activate the line

    // printvector(std::cout,pressure,"initial pressure","row",200,4);

    //  vtkout(*grid, pressure, permeability, velocity, "multiscalec", 0);
    //printvector(std::cout,velocity,"velocity","row",4,4);

    // generate one meta vtk-file holding the individual timesteps
    //    char fileName[128] = "pipe";
    //    char multiFileName[128];
    //    char fileNameVTK[128];
    //    sprintf(multiFileName,"multi-%s.pvd",fileName );
    //    std::ofstream multiFile(multiFileName);
    //    multiFile << "<?xml version=\"1.0\"?>" << std::endl
    //    << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" " << std::endl
    //    << "compressor=\"vtkZLibDataCompressor\">"  << std::endl
    //    << " <Collection>" << std::endl;

    // now do the time steps
    Scalar dt = max_dt;
    Scalar t = tstart;
    //    int k = 0;
    //    std::cout.setf (std::ios::scientific, std::ios::floatfield);

    std::cout << "start pipeflow  sp sp sp sp sp sp sp sp sp sp sp sp sp sp sp sp sp sp sp sp sp sp sp sp sp sp sp " <<  std::endl;

    while (t<tend) {

        SetDrichletBoundary (t);

        // outputs
        //        std::cout << "timestep: " << k << "\t t=" << t << "\t dt=" << dt << std::endl;
        //        std::cout.setf (std::ios::scientific, std::ios::floatfield);
        // printvector(std::cout,pressure,"pressure","row",200,1);

        //        printvector(std::cout,pressure,"pressure","row",200,1);
        //        printvector(std::cout,velocity,"velocity","row",200,1);

        //        if (k%modulo==0)
        //        {
        //            // create pressureVtkOutput vector with 0 value , size=sizeOfBigGridNodes and copy pressure values from pipe pressure
        //            int sizeOfNodes = vMap.size();
        //            std::cout <<sizeOfNodes << std::endl;
        //            Press pressureVtkOutput(sizeOfNodes);
        //            pressureVtkOutput = 0;
        //            {
        //                for (unsigned n = 0; n < vertexVectorOnLine.size(); n++)
        //                {
        //                    pressureVtkOutput[vertexVectorOnLine[n].globalId] = pressure[n];
        //                }
        //            }
        //            vtkout_pipeflow(grid, pressureVtkOutput, fileName,k/modulo);
        //            //
        //            sprintf(fileNameVTK,"%s-%05d.vtu", fileName, k/modulo);
        //            multiFile << "   <DataSet timestep=\"" << t << "\" file=\""
        //            << fileNameVTK << "\"/>" << std::endl;
        //        }

        Iteration (t, dt, max_def);

        //        k++;
        t += dt;
    }

    std::cout << "end pipeflow ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep ep" <<  std::endl;
    // finalize output
    //      multiFile << " </Collection>" << std::endl << "</VTKFile>" << std::endl;
    //      multiFile.close();

    return;
}

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class LmbdLocal, class Grid, class VertexMapper, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, LmbdLocal, Grid, VertexMapper, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::Iteration (Scalar t, Scalar dt, Scalar max_def)
{
    //    int nElem =  vertexVectorOnLine.size();
    //    int systemSize = nElem;

    Press pressureIt = pressure;

    for (int n = 0; n < 101; n++)
        {
            Press pressureOldIt = pressureIt;

            IterationStep_Mass(t, dt, pressureIt);
            //        printvector(std::cout, pressureIt, "PressureIt","row",100,1,4);

            Press deltaPressure = pressureIt;
            deltaPressure -= pressureOldIt;

            Scalar deltaPressureTwoNorm = deltaPressure.two_norm();
            Scalar pressureItTwoNorm = pressureIt.two_norm();
            Scalar defectPressure;

            if (pressureItTwoNorm == 0.0)
                {
                    defectPressure = 0.0;
                }
            else
                {
                    defectPressure = deltaPressureTwoNorm/pressureItTwoNorm;
                }

            //        std::cout << "defectPressure=" << defectPressure << std::endl;

            if (defectPressure < max_def)
                {
                    pressure= pressureIt;
                    //            std::cout << "Converged in  n= " << n << "iterations" << std::endl;
                    //            printvector(std::cout,pressure,"Pressure","row",100,1,4);
                    break;
                }
            else if (n >= 100)
                {
                    Dune::Exception exception;
                    exception.message("Diverged !!!!");
                    throw exception;
                }
        }

    return;
}

template<class BCP, class BCV, class ICP, class ICV, class SST, class Press, class Lmbd, class LmbdLocal, class Grid, class VertexMapper, class MapperNodeGlobalIDtoOnOutIndexType, class VertexVectorOnLineType, class VertexVectorOutLineType>
void PipeFlow<BCP, BCV, ICP, ICV, SST, Press, Lmbd, LmbdLocal, Grid, VertexMapper, MapperNodeGlobalIDtoOnOutIndexType, VertexVectorOnLineType, VertexVectorOutLineType>::IterationStep_Mass (Scalar t, Scalar dt, Press& pressureIt)
{
    int nElem = vertexVectorOnLine.size();
    int systemSize = nElem;

    // solution vector SolutionVector and right side vector f
    Press f(systemSize);
    Press SolutionVector(systemSize);
    f=0;
    SolutionVector=0;
    typedef Dune::FieldMatrix<Scalar,1,1> MB;
    Dune::BCRSMatrix<MB> A(systemSize, systemSize, Dune::BCRSMatrix<MB>::random);

    // determine matrix row sizes
    for (unsigned k = 0; k < vertexVectorOnLine.size(); k++)
        {
            // cell index
            int indexi = k;

            // initialize row size = diagonal
            int rowSize = systemSize;

            A.setrowsize(indexi, rowSize);
        }
    A.endrowsizes();

    // determine position of matrix entries
    for (int i = 0;i<systemSize; ++i)
        {
            for (int j = 0;j<systemSize; ++j)
                {
                    A.addindex(i, j);
                }
        }
    A.endindices();

    //initiliaze A with 0
    for (int i = 0;i<systemSize; ++i)
        {
            f[i]= 0.0;
            for (int j = 0;j<systemSize; ++j)
                {
                    A[i][j]= 0.0;
                }
        }

    // fill matrix
    for (unsigned k = 0; k < vertexVectorOnLine.size(); k++)
        {
            MassEquation<Dune::BCRSMatrix<MB> > (k, A, f, t, dt);
        }   // end grid traversal

    // set dirichlet boundary in Matrix A and right hand side f
    for (unsigned k = 0; k < vertexVectorOnLine.size(); k++)
        {
            if (vertexVectorOnLine[k].boundary() )
                {
                    int boundaryId = (int) vertexVectorOnLine[k].parameter[1];
                    int boundaryType;
                    Scalar boundaryConditionValue = boundaryConditionP(boundaryId, t, boundaryType);
                    if (boundaryType == 1) // if dirichlet boundary
                        {
                            int indexi = k;
                            A[indexi]=0;
                            A[indexi][indexi]=1.0;
                            f[indexi]= boundaryConditionValue;
                        }
                }
        }   // end grid traversal

    //    printvector(std::cout,f,"right hand side","row",200,1);
    //    printmatrix(std::cout,A,"matrix","",8,1);
    // set up the high-level solver objects
    typedef Dune::FieldVector<Scalar, 1> VB;
    typedef Dune::BlockVector<VB> Vector;
    typedef Dune::BCRSMatrix<MB> Matrix;
    Dune::MatrixAdapter<Matrix,Vector,Vector> op(A);        // make linear operator from A
    Dune::InverseOperatorResult r;

    Dune::SeqPardiso<Matrix,Vector,Vector> pardiso(A);        // preconditioner object
    Dune::LoopSolver<Vector> loop(op, pardiso, 1E-14, 2, 1);  // an inverse operator
    loop.apply(SolutionVector, f, r);                            // call the solver

    //    Dune::SeqGS<Matrix,Vector,Vector> ilu0(bigMatrix, 1, 1.0);         // preconditioner object
    //    Dune::BiCGSTABSolver<Vector> cg(op, ilu0, 1E-14, 1000, 2); // an inverse operator
    //    cg.apply(bigSolutionVector, fBig, r);                               // call the solver


    //fill deltaVelocity
    for (int i = 0;i<nElem; ++i)
        {
            pressureIt[i]=SolutionVector[i];
        }

    //    printvector(std::cout,deltaVelocity,"deltaVelocity","row",200,1);
    //    printvector(std::cout,deltaPressure,"deltaPressure","row",200,1);
    return;
}
