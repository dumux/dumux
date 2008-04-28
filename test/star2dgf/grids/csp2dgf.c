#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TRI_3 8
#define TETRA_4 4
#define QUAD_4 14
#define HEXA_8 6
#define BAR_2 2

typedef struct points{
  double x;
  double y;
  double z;
  double pbVal;
  int pbFlag;
} *POINTS;

typedef struct elements{
  int objectType;
  int finiteElementType;
  int anz;
  int *nodes;
  int materialType;
} *ELEMENTS;




int main (int argc, char **argv)
{
  int i;
  int dim = 3;

  if (argc != 2)
    {
      printf("usage: csp2dgf basefilename\n");
      return (1);
    }

  char baseFileName[120];
  strcpy(baseFileName, argv[1]);

  char ascFileName[120];
  strcpy(ascFileName, baseFileName);
  strcat(ascFileName, ".asc");

  // test dimension
  FILE *inputFile; 
  if((inputFile=fopen(ascFileName, "r"))==NULL) {
    printf("ERROR: %s doesn't exist\n",ascFileName);
    return 1;
  }
  
  char testString[120]; 
  int triDim=0;
 
  char dgfFileName[120];
  strcpy(dgfFileName, baseFileName);
  strcat(dgfFileName, ".dgf");

  // input coordinates and boundary values
  char datFileName[120];
  strcpy(datFileName, baseFileName);
  strcat(datFileName, ".dat");
 
  if((inputFile = fopen(datFileName, "r"))==NULL) {
    printf("ERROR: %s doesn't exist\n",datFileName);
    return 1;
  }

  int numberOfNodes;  
  fscanf(inputFile,"%d",&numberOfNodes);
  fgets(testString,120,inputFile);

  POINTS nodes;
  nodes=malloc(numberOfNodes*sizeof(*nodes));

  for(i=0;i<numberOfNodes;++i)
    fscanf(inputFile,"%lf", &nodes[i].x);
  for(i=0;i<numberOfNodes;++i)
    fscanf(inputFile,"%lf",&nodes[i].y);
  for(i=0;i<numberOfNodes;++i){
    fscanf(inputFile,"%lf",&nodes[i].z);
    if(nodes[i].z!=0)
      triDim=1;
  }
 
  fscanf(inputFile,"%s %s",&testString,&testString);
  for(i=0;i<numberOfNodes;++i)
    fscanf(inputFile,"%d",&nodes[i].pbFlag);

  fscanf(inputFile,"%s %s",&testString,&testString);
  for(i=0;i<numberOfNodes;++i)
    fscanf(inputFile,"%lf",&nodes[i].pbVal);

  // output vertex
  // coordinates, pbFlag, pbValue
  FILE *dgfFile = fopen(dgfFileName, "w");
  if(triDim)
    fprintf(dgfFile, "DGF\nVERTEX\nparameters 2\nfirstindex 1\n");
  else
    fprintf(dgfFile,"DGF\nVERTEX\nparameters 2\n");
  printf("nodes: %d\n", numberOfNodes);

  if(triDim)
    for(i=0;i<numberOfNodes;++i)
      fprintf(dgfFile, "%14.7e %14.7e %14.7e %d %14.7e \% % %d\n", nodes[i].x, nodes[i].y, nodes[i].z, nodes[i].pbFlag, nodes[i].pbVal,i+1);
  else 
    for(i=0;i<numberOfNodes;++i)
      fprintf(dgfFile, "%14.7e %14.7e %d %14.7e \% % %d\n", nodes[i].x, nodes[i].y, nodes[i].pbFlag, nodes[i].pbVal, i);
   


  // input element types and node ids
  int numberOfElements;
  fscanf(inputFile,"%d",&numberOfElements);
  fgets(testString,120,inputFile);
  printf("elements: %d\n",numberOfElements);

  ELEMENTS elements;
  elements = malloc(numberOfElements*sizeof(*elements));
  for(i=0;i<numberOfElements;++i)
    fscanf(inputFile,"%d", &elements[i].finiteElementType);
  
  int ifSimplex=0;
  for(i=0;i<numberOfElements;++i){
    switch(elements[i].finiteElementType){
    case HEXA_8:
      elements[i].anz=8;
      break;
    case QUAD_4:
    case TETRA_4:
      elements[i].anz=4;
      break;
    case TRI_3:
      elements[i].anz=3;
      ifSimplex=1;
      break;
    case BAR_2:
      elements[i].anz=2;
      break;
    default:
      printf("unknown finiteElementType\n");
      free(nodes);
      free(elements);
      fclose(dgfFile);
      fclose(inputFile);
      return 1;
    } 
  }

  if(triDim)
    if (ifSimplex)
      fprintf(dgfFile,"#\nSIMPLEX\nparameters 1\n");
    else 
      fprintf(dgfFile,"#\nCUBE\nparameters 1\nmap 0 1 3 2 4 5 7 6\n");
  else 
    if(ifSimplex)
      fprintf(dgfFile,"#\nSIMPLEX\nparameters 1\n");
    else
      fprintf(dgfFile,"#\nCUBE\nparameters 1\n");
    
  for(i=0;i<numberOfElements;++i)
    elements[i].nodes=calloc(elements[i].anz,sizeof(int));

  fgets(testString,120,inputFile);
  fgets(testString,120,inputFile);

  int j;
  for(i=0;i<numberOfElements;++i)
    for(j=0;j<elements[i].anz;++j)
      fscanf(inputFile,"%d",&elements[i].nodes[j]);
      
  do {
    fscanf(inputFile,"%s",testString);
  } while(strcmp(testString,"PMATERIAL"));

  for(i=0;i<numberOfElements;++i)
      fscanf(inputFile,"%d",&elements[i].materialType);

  fclose(inputFile);


  // output elements
  int numberOfFamilies;
  inputFile=fopen(ascFileName,"r");  
  fgets(testString,120,inputFile);
  fgets(testString,120,inputFile);
  fscanf(inputFile,"%d",&numberOfFamilies);
  do {
    fscanf(inputFile,"%s",testString);
  } while(strcmp(testString,"sequence"));
 
  fgets(testString,120,inputFile);

  int index;
  int numberOfUsedElements;
  int k;
  
  //im 3d wird bei 1 angefangen zu zaehlen
  if(triDim) 
    for(i=0;i<numberOfElements;++i)
      for(j=0;j<elements[i].anz;++j)
	elements[i].nodes[j]+=1; 

  fscanf(inputFile,"%s %s %d",&testString,&testString,&numberOfUsedElements);
  for(i=0;i<numberOfUsedElements;++i) {
    fscanf(inputFile,"%d",&index);
    for(j=0;j<elements[index].anz;++j)
      fprintf(dgfFile,"%7d ",elements[index].nodes[j]);
    fprintf(dgfFile,"%7d\n",elements[index].materialType);
  }

  //boundary id, node ids
  //boundary id ?! deshalb testhalber family-nr 
  fprintf(dgfFile,"#\nBOUNDARYSEGMENTS\n");  
  for(i=1;i<numberOfFamilies;++i) {
    fscanf(inputFile,"%s %s %d",&testString,&testString,&numberOfUsedElements);
    for(j=0;j<numberOfUsedElements;++j) {
      fscanf(inputFile,"%d",&index);
       if(elements[index].anz>2){
	 fprintf(dgfFile,"%7d ",i); 
	for(k=0;k<elements[index].anz;++k)
	  fprintf(dgfFile,"%7d ",elements[index].nodes[k]);
	fprintf(dgfFile,"\n");
	}
    }
  }
  
  fprintf(dgfFile,"#\n# %s\n",dgfFileName);

  fclose(inputFile);
  fclose(dgfFile);

  for(i=0;i<numberOfElements;++i)
    free(elements[i].nodes);
  free(elements);
  free(nodes);

  return 0;
}
