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
  int isVolume;
  int isFace;
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
  FILE *ascFile;
  if((ascFile=fopen(ascFileName, "r"))==NULL) {
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

  FILE *datFile;
  if((datFile = fopen(datFileName, "r"))==NULL) {
    printf("ERROR: %s doesn't exist\n",datFileName);
    return 1;
  }

  int numberOfNodes;
  fscanf(datFile,"%d",&numberOfNodes);
  fgets(testString,120,datFile);

  POINTS nodes;
  nodes=malloc(numberOfNodes*sizeof(*nodes));

  for(i=0;i<numberOfNodes;++i)
    fscanf(datFile,"%lf", &nodes[i].x);
  for(i=0;i<numberOfNodes;++i)
    fscanf(datFile,"%lf",&nodes[i].y);
  for(i=0;i<numberOfNodes;++i){
    fscanf(datFile,"%lf",&nodes[i].z);
    if(nodes[i].z!=0)
      triDim=1;
  }

  fscanf(datFile,"%s %s",&testString,&testString);
  for(i=0;i<numberOfNodes;++i)
    fscanf(datFile,"%d",&nodes[i].pbFlag);

  fscanf(datFile,"%s %s",&testString,&testString);
  for(i=0;i<numberOfNodes;++i)
    fscanf(datFile,"%lf",&nodes[i].pbVal);

  // output vertex
  // coordinates, pbFlag, pbValue
  FILE *dgfFile = fopen(dgfFileName, "w");
  fprintf(dgfFile, "DGF\nVERTEX\nparameters 2\n");
  printf("nodes: %d\n", numberOfNodes);

  if(triDim)
    for(i=0;i<numberOfNodes;++i)
      fprintf(dgfFile, "%14.7e %14.7e %14.7e %d %14.7e \% % %d\n", nodes[i].x, nodes[i].y, nodes[i].z, nodes[i].pbFlag, nodes[i].pbVal,i+1);
  else
    for(i=0;i<numberOfNodes;++i)
      fprintf(dgfFile, "%14.7e %14.7e %d %14.7e \% % %d\n", nodes[i].x, nodes[i].y, nodes[i].pbFlag, nodes[i].pbVal, i);



  // input element types and node ids
  int numberOfElements;
  fscanf(datFile,"%d",&numberOfElements);
  fgets(testString,120,datFile);
  printf("elements:\n  specified in CSP-files: %d\n",numberOfElements);

  ELEMENTS elements;
  elements = malloc(numberOfElements*sizeof(*elements));
  for(i=0;i<numberOfElements;++i)
    fscanf(datFile,"%d", &elements[i].finiteElementType);

  int ifSimplex=0;
  for(i=0;i<numberOfElements;++i){
    switch(elements[i].finiteElementType){
    case HEXA_8:
      elements[i].anz=8;
      elements[i].isVolume=1;
      elements[i].isFace=0;
      break;
    case QUAD_4:
      if (triDim) {
	elements[i].isVolume=0;
	elements[i].isFace=1;
      }
      else {
	elements[i].isVolume=1;
	elements[i].isFace=0;
      }
      elements[i].anz=4;
      break;
    case TETRA_4:
      elements[i].isVolume=1;
      elements[i].isFace=0;
      elements[i].anz=4;
      break;
    case TRI_3:
      if (triDim) {
	elements[i].isVolume=0;
	elements[i].isFace=1;
      }
      else {
	elements[i].isVolume=1;
	elements[i].isFace=0;
      }
      elements[i].anz=3;
      ifSimplex=1;
      break;
    case BAR_2:
      if (triDim) {
	elements[i].isVolume=0;
	elements[i].isFace=0;
      }
      else {
	elements[i].isVolume=0;
	elements[i].isFace=1;
      }
      elements[i].anz=2;
      break;
    default:
      printf("unknown finiteElementType\n");
      free(nodes);
      free(elements);
      fclose(dgfFile);
      fclose(ascFile);
      fclose(datFile);
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
      fprintf(dgfFile,"#\nCUBE\nparameters 1\nmap 3 2 0 1\n");

  for(i=0;i<numberOfElements;++i)
    elements[i].nodes=calloc(elements[i].anz,sizeof(int));

  fgets(testString,120,datFile);
  fgets(testString,120,datFile);

  int j;
  for(i=0;i<numberOfElements;++i)
    for(j=0;j<elements[i].anz;++j)
      fscanf(datFile,"%d",&elements[i].nodes[j]);

  do {
    fscanf(datFile,"%s",testString);
  } while(strcmp(testString,"PMATERIAL"));

  for(i=0;i<numberOfElements;++i)
      fscanf(datFile,"%d",&elements[i].materialType);

  fclose(datFile);


  // output elements
  int numberOfFamilies;
  fgets(testString,120,ascFile);
  fgets(testString,120,ascFile);
  fscanf(ascFile,"%d",&numberOfFamilies);
  do {
    fscanf(ascFile,"%s",testString);
  } while(strcmp(testString,"sequence"));

  fgets(testString,120,ascFile);

  int index;
  int numberOfUsedElements;
  int k;
  int boundaryStarted = 0;
  int count;
  for(i=0;i<numberOfFamilies;++i) {
    fscanf(ascFile,"%s %s %d",&testString,&testString,&numberOfUsedElements);
    for(j=0;j<numberOfUsedElements;++j) {
      fscanf(ascFile,"%d",&index);

      if (!boundaryStarted && !elements[index].isVolume) {
	boundaryStarted = 1;
	fprintf(dgfFile,"#\nBOUNDARYSEGMENTS\n");
	printf("  volume elements: %d\n", count);
	count = 0;
      }

      if (boundaryStarted && elements[index].isVolume) {
	printf("CSP File lists some volume elements after some boundary elements! Aborting...\n");
	return 1;
      }

      if (boundaryStarted && elements[index].isFace) {
	fprintf(dgfFile,"%7d ",elements[index].materialType);
	for(k=0;k<elements[index].anz;++k)
	  fprintf(dgfFile,"%7d ",elements[index].nodes[k]);
	fprintf(dgfFile,"\n");
	count++;
      }
      else if (elements[index].isVolume) {
	for(k=0;k<elements[index].anz;++k)
	  fprintf(dgfFile,"%7d ",elements[index].nodes[k]);
	fprintf(dgfFile,"%7d\n",elements[index].materialType);
	count++;
      }

    }
  }

  printf("  face elements: %d\n", count);
  fprintf(dgfFile,"#\n# %s\n",dgfFileName);

  fclose(ascFile);
  fclose(dgfFile);

  for(i=0;i<numberOfElements;++i)
    free(elements[i].nodes);
  free(elements);
  free(nodes);

  printf("DGF file %s generated.\n", dgfFileName);

  return 0;
}
