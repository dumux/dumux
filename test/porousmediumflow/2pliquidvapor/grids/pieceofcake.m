% Cake cake grid
clear;

pathdgf = 'pieceofcake.dgf';
pathvtk = 'pieceofcake.vtk';
R = 5; % Domain radius
A = pi/6; % Domain angle
Z = 7; % Domain length in z direction
nSlice = 4; %no of slices
nZ = 28; %approximate no of elements in Z direction
cellSizeZ = Z/nZ;
wellRadius = 0.135; %Well radius at the tip of the cake slice
numCellsR = 21; %number of Cells in R direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%dr Vector

%the spacing in R direction is generated with the following formula:
%dr(i) = a*exp(lambda*i) - wellRadius 
%--> dr(1) = 0
%--> dr(numCellsR)=R-wellRadius


lambda = log(R/wellRadius)/(numCellsR-1)
a = wellRadius/exp(lambda)

for i=1:numCellsR
    dr(i) = a*exp(lambda*i) - wellRadius;
end



%da Vector
da = 0:A/nSlice:A;

%Z-Vector
dz = 0:Z/nZ:Z;

%Create the dgf file
dgfWriterCake(pathdgf, dr, da, dz, wellRadius);

%Create the vtk file
vtkWriterCake(pathvtk, dr, da, dz, wellRadius);


noNodes = length(dr)*length(da)*length(dz);
display(noNodes);
