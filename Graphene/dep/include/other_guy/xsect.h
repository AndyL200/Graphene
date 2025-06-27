/**********************************************
  Cross section tables header
**********************************************/

#ifndef __xsect_h
#define __xsect_h

struct crossSectStruct
{
  int n;/* number of data pairs */
  float thresholdE; /* threshold energy */
  float* E; /* energy array */
  float* sigma; /* cross section array */
};

typedef struct crossSectStruct CrossSect;

int readXSectionTables(CrossSect** xSect, CrossSect* totalMom, char* fileName);
float interpolateTableLinear(CrossSect *xSect, float E);

#endif
