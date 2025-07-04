/***************************************************
  Cross section table functions
***************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#define GLOBALORIGIN
#include "xsect.h"
#undef GLOBALORIGIN

/****************************************************
Read all cross sections from fileName into a array of cross sections,
and return the address of the first element of the array in xSect, the
number of items in the array is the return value.
****************************************************/

int readXSectionTables(CrossSect** xSectPtr, CrossSect* totalMomPtr, char* fileName)
{
#define READCOMMENTLINE fscanf(xSectFile, "%[^\n]%*c", buffer)
  FILE *xSectFile;
  int i, j;
  int number, rxVib, rxElect, rxOther, rxMomXFer; /* number of reactions */
  int nMomXFer; /* number of data pairs */
  float molWeight;
  char buffer[200];
  float scale;
  CrossSect	*xSect;

  if (!(xSectFile = fopen(fileName, "r")))
  {
    printf("readXSection(): file open failed on %s.\n", fileName);
    exit(-1);
  }
  while (fscanf(xSectFile, "%d %d %d %d %d", &rxVib, &rxElect, &rxOther,
		&rxMomXFer, &nMomXFer) < 5) READCOMMENTLINE;
  number = rxVib + rxElect + rxOther + rxMomXFer;
  while (fscanf(xSectFile, "%g", &molWeight) < 1) READCOMMENTLINE;

  /* allocate memory */
  xSect = (CrossSect*) malloc(sizeof(CrossSect)*number);
	
  /* read total Momentum transfer table */
  totalMomPtr->n = nMomXFer;
  totalMomPtr->thresholdE = 0;
  totalMomPtr->E = (float *)malloc(totalMomPtr->n*sizeof(float));
  totalMomPtr->sigma = (float *) malloc(totalMomPtr->n*sizeof(float));
  while (fscanf(xSectFile, "%g", &scale) < 1) READCOMMENTLINE;
  for (j = 0; j < totalMomPtr->n; j++)
  {
    while (fscanf(xSectFile, "%g %g", &totalMomPtr->E[j], &totalMomPtr->sigma[j]) < 2)
      READCOMMENTLINE;
    totalMomPtr->sigma[j] *= 1E-20;  /* convert from angstrom^2 to m^2 */
  }
  
  /* read the rest of the cross section tables */
  for (i = 1; i < number; i++)
  {
    float statWt, scale2;
    while (fscanf(xSectFile, "%g %g %d %g", &xSect[i].thresholdE, &statWt,
		  &xSect[i].n, &scale2) < 4) READCOMMENTLINE;
    xSect[i].E = (float *)malloc(xSect[i].n*sizeof(float));
    xSect[i].sigma = (float *)malloc(xSect[i].n*sizeof(float));
			
    for (j = 0; j < xSect[i].n; j++)
    {
      while (fscanf(xSectFile, "%g %g", &xSect[i].E[j], &xSect[i].sigma[j])
	     < 2) READCOMMENTLINE;
      xSect[i].sigma[j] *= scale2*1E-20; /* convert from angstrom^2 to m^2 */
    }
  }

  fclose(xSectFile);

  /* Now compute the elastic scattering cross section from the total mom. xfer */
  xSect[0].n = totalMomPtr->n;
  xSect[0].E = (float*)malloc(xSect[0].n*sizeof(float));
  xSect[0].sigma = (float*)malloc(xSect[0].n*sizeof(float));
  xSect[0].thresholdE = totalMomPtr->thresholdE;
  for (j = 0; j < xSect[0].n; j++)
  {
    float E;
    E = xSect[0].E[j] = totalMomPtr->E[j];
    xSect[0].sigma[j] = totalMomPtr->sigma[j];
    for (i = 1; i < number; i++)
    {
      xSect[0].sigma[j] -= interpolateTableLinear(&xSect[i], E);
    }
  }

  /* Point the xSectPtr to table */
  *xSectPtr = xSect;
  return number;
}

/*******************************************
Return the linearly interpolated cross section value for the cross section
table xSect given an energy E.  Both E[] and sigma[] are from 0...n-1.
The interpolation assumes sigma is 0 below the range given, and remains
constant at the last value above the energy range.
*******************************************/

float interpolateTableLinear(CrossSect *xSect, float E)
{
  int ju, jl, jm; /* upper, lower and middle indices into x */
  float *x = xSect->E;
  float *y = xSect->sigma;
  if (x[xSect->n - 1] > x[0])
  {
    jl = 0;
    ju = xSect->n - 1;
  }
  else
  {
    jl = xSect->n - 1;
    ju = 0;
  }
  if (E < x[jl]) return 0;
  else if (E > x[ju]) return y[ju];
  
  while (ju - jl > 1)
  {
    jm = (jl + ju)/2;
    if (x[jm] > E) ju = jm;
    else jl = jm;
  }

/*
printf("E=%G E[%d]=%G E[%d]=%G ", E, jl, x[jl], ju, x[ju]);
printf("sigma=%G\n", ((E - x[jl])*y[ju] + (x[ju] - E)*y[jl])/(x[ju] - x[jl]));
*/
  return ((E - x[jl])*y[ju] + (x[ju] - E)*y[jl])/(x[ju] - x[jl]);
}
