#include "pdp1.h"

/***************************************************************/
void One_2_One(float *ary, int nmax);

void gather(int isp)
{
  register int i, j;
  register float s;
  
  //printf("%E \n",sp_n_0[1][1]);
  for (j=0; j< ng; j++) {
    sp_n_0[isp][j]= sp_n_k[isp][j] +sp_n_mcc[isp][j];
    sp_n_mcc[isp][j]= sp_n_k[isp][j]= 0.0;
  }
  //printf("First Step \n");
  for (i=np[isp]-1; i>=0; i--) {
    //printf("%d \n",i);
    j = x[isp][i];
    //printf("x = %E",x[isp][i]);
    s = x[isp][i] - j;
    //printf("s done\n");
    //printf("j = %d",j);
    sp_n_k[isp][j]  += 1. - s;
    //printf("first sp\n");
    sp_n_k[isp][j+1]+= s;
    //printf("second sp\n");
  }
	sp_n_k[isp][0]  *= 2.;
	sp_n_k[isp][nc] *= 2.;

	//printf("Second Step \n");
  /************************************************/
  /* Smoothing the charge density of each species */
  
  for(i=0; i< nsmoothing; i++) One_2_One(sp_n_k[isp], nc);

  //printf("Third Step \n");
}

/***************************************************************/

void setrho(void)
{
  int j, isp;

  for(isp=0; isp<nsp; isp++) {
    k_count[isp] = 0;
    gather(isp);
    for (j=0; j<ng; j++) {
      sp_n[isp][j]= sp_n_k[isp][j];
      sp_n_mcc[isp][j]= 0.0;
    }
  }
}

/***************************************************************/
/* Smoothing the array using the 1-2-1 method with the proper   */
/* boundary conditions to conserve charge.                     */

void One_2_One(float *ary, int nmax)
{
  register int i;
  static int nlocal=0;
  static float *temp=0;
  
  if(nlocal < nmax) {
    if(temp) free(temp);
    temp = (float *)malloc((nmax+1)*sizeof(float));
    nlocal = nmax;
  }
  
	temp[0]= (ary[0] +ary[1])/2;
	for(i=1; i< nmax; i++)  temp[i]= (ary[i-1] +2*ary[i] +ary[i+1])/4;
	temp[nmax]= (ary[nmax] +ary[nmax-1])/2;

  for(i=0; i<= nmax; i++) ary[i]= temp[i];
}

/***************************************************************/
