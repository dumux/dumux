#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main (int argc, char **argv)
{
  int i;
  int dim = 3;

  if (argc != 2)
    {
      printf("usage: star2dgf basefilename\n");
      return (1);
    }

  char baseFileName[120];
  strcpy(baseFileName, argv[1]);

  char dgfFileName[120];
  strcpy(dgfFileName, baseFileName);
  strcat(dgfFileName, ".dgf");
  FILE *dgfFile = fopen(dgfFileName, "w");
  fprintf(dgfFile, "DGF\nVERTEX\nfirstindex 1\n");

  char inputFileName[120];
  strcpy(inputFileName, baseFileName);
  strcat(inputFileName, ".inp");
  FILE *inputFile = fopen(inputFileName, "r");

  char testString[120]; 
  do {
    fscanf(inputFile, "%s", testString);
  }
  while(strcmp(testString, "nodes"));
  int numberOfNodes;
  fscanf(inputFile, "%s %d", testString, &numberOfNodes);
  printf("nodes: %d\n", numberOfNodes);

  do {
    fscanf(inputFile, "%s", testString);
  }
  while(strcmp(testString, "elements"));
  int numberOfElements;
  fscanf(inputFile, "%s %d", testString, &numberOfElements);
  printf("elements: %d\n", numberOfElements);

  fclose(inputFile);

  char vertexFileName[120];
  strcpy(vertexFileName, baseFileName);
  strcat(vertexFileName, ".vrt");
  FILE *vertexFile = fopen(vertexFileName, "r");

  for (i = 1; i <= numberOfNodes; i++) {
    int dummy;
    char xCoord[30], yCoord[30], zCoord[30];
    fscanf(vertexFile, " %d %16s %16s %16s", &dummy, xCoord, yCoord, zCoord);
    fprintf(dgfFile, "%14.7e %14.7e %14.7e \% %", strtod(xCoord, NULL), strtod(yCoord, NULL), strtod(zCoord, NULL));
    fprintf(dgfFile, " %d\n", i);
  }
  fprintf(dgfFile, "#\n");

  char cellFileName[120];
  strcpy(cellFileName, baseFileName);
  strcat(cellFileName, ".cel");
  FILE *cellFile = fopen(cellFileName, "r");

  int isVolumeElement = 1;
  int isSimplex = 0;
  int surfaceStarted = 0;
  for (i = 1; i <= numberOfElements; i++) {
    int dummy;
    int vertices[8];
    int boundaryId;
    int surfaceOrVolume; 
    fscanf(cellFile, " %d %d %d %d %d %d %d %d %d %d %d %d", &dummy, 
	   &vertices[0], &vertices[1], &vertices[2], &vertices[3], 
	   &vertices[4], &vertices[5], &vertices[6], &vertices[7], 
	   &boundaryId, &surfaceOrVolume, &dummy);
    if (i == 1) {
      if (vertices[2] == vertices[3]) {
	isSimplex = 1;
	fprintf(dgfFile, "SIMPLEX\n");
      }
      else {
	fprintf(dgfFile, "CUBE\n");
      }
    }
    if (surfaceOrVolume == isVolumeElement) {
      if (isSimplex)
	fprintf(dgfFile, "%7d %7d %7d %7d\n", vertices[0], vertices[1], vertices[2], vertices[4]);
      else 
	fprintf(dgfFile, "%7d %7d %7d %7d %7d %7d %7d %7d\n", vertices[0], vertices[1], vertices[3], vertices[2], 
		vertices[4], vertices[5], vertices[7], vertices[6]);
    }
    else {
      if (!surfaceStarted) {
	surfaceStarted = 1;
	int numberOfVolumeElements = i - 1;
	printf(" volume: %d\n", numberOfVolumeElements);
	printf(" surface: %d\n", numberOfElements - numberOfVolumeElements);
	fprintf(dgfFile, "#\nBOUNDARYSEGMENTS\n");
      }
      if (isSimplex)
	fprintf(dgfFile, "%d %7d %7d %7d\n", boundaryId, vertices[0], vertices[1], vertices[2]);
      else 
	fprintf(dgfFile, "%d %7d %7d %7d %7d\n", boundaryId, vertices[0], vertices[1], vertices[3], vertices[2]);
    }
  }

  fprintf(dgfFile, "#\n# %s\n", dgfFileName);
  fclose(dgfFile);

  return (0);
}
