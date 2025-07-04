#include<math.h>
#include <stdlib.h> 
#include <stdio.h>

#define  NVTS    3.3

#ifndef M_PI
#define M_PI    3.14159265358979323846
#endif

/**************************************************************/
float F(float);
float frand(void);

void maxwellv(float *vx_local, float *vy_local, float *vz_local, float vth)
{
  static int nvel, init_flag= 1;
  static float *vsave_x, *vsave_y, *vsave_z;
  
  int i, n;
  float vmag, aphi, dv, rr, sintheta, costheta;
  
  if (init_flag) {
    nvel= 1./(1-F((float)NVTS));
    dv = sqrt(M_PI)/(4.0*nvel);
    if(nvel > 100000)
      puts("Warning: Your choice of NVTS has made nvel > 1e5"); 
    
    vsave_x= (float *) malloc(nvel*sizeof(float));
    vsave_y= (float *) malloc(nvel*sizeof(float));
    vsave_z= (float *) malloc(nvel*sizeof(float));
    
    init_flag= 0;
    i=n=0;
    for (n=0; n<nvel; n++) {
      rr=(1.0*n)/nvel;
      while (F(i*dv)< rr) i++;
      vmag=i*dv;
      aphi=2*M_PI*frand();
      costheta = 1-2*frand();
      sintheta = sqrt(1-costheta*costheta);
      vsave_x[n] = vmag*sintheta*cos(aphi);
      vsave_y[n] = vmag*sintheta*sin(aphi);
      vsave_z[n] = vmag*costheta;
    }
  }
  n = (nvel-1)*frand();
  *vx_local = vth*vsave_x[n];
  *vy_local = vth*vsave_y[n];
  *vz_local = vth*vsave_z[n];
}

/**************************************************************/

float F(float v)
{
  return(-2*v*exp(-v*v)/sqrt(M_PI) +erf(v));
}

/**************************************************************/
