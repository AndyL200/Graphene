#include "pdp1.h"
#include "xgrafix.h"

#define NEMAX 500
#define NIMAX 500

// 5 collisions are tracked for nitrogen (including both electron and ion collisions)
float nitsigma1(float),nitsigma2(float),nitsigma3(float),nitsigma4(float),nitsigma5(float);
float gden,vgth,extengy0,ionengy0;
int ecolsp,icolsp,ionsp;

void nitnewvel(float,float,float *, float *, float *, float,int);

/**************************************************************/
/*Monte Carlo Collisions for electron-neutral and ion-neutral */

void nitrogenmcc(int isp)
{
  static int init_flag = 1;
  static float ecol_extra,icol_extra;
  static float max_sigmav_e,max_sigmav_i;
  static float col_prob_e, col_prob_i;

  register int j, k, index;
  int N,nnp,s;
  float random,dum,vel,sigma_total,metal;
  float engy,rengy,del,phi1;
  float cosphi,sinphi,coschi,sinchi;
  float sum_sigma,temp,vneutx,vneuty,vneutz;

  //printf("%s \n","Initialization");
  // Initialization step called only the first time
  if (init_flag) {
    if (gtemp)
      gden = NperTORR*pressure/(gtemp);
    else
      gden = NperTORR*pressure/FLOAT_MIN;
    ecolsp = ecollisional-1;
    icolsp = icollisional-1;
    ionsp = ionspecies-1; /* Fixing the indices into the array of species */

    /***************************************/
    /* Calculating the null collision prob */

    //printf("%s \n","null collision prob ecollisional");
    if (ecollisional) {
      if (ionsp > 0) vgth = sqrt(gtemp/Escale[ionsp]);
      extengy0 = 11.03;
      ionengy0 = 15.58;
     
      max_sigmav_e = 0.0;
      for (engy=0;engy<NEMAX;engy+=0.1){
	//printf("%f \n",engy);
	//printf("%E \n",m[ecolsp]);
	//printf("%E \n",nitsigma1(engy));
	//printf("%E %E \n",engy,nitsigma2(engy));
	//printf("%E \n",nitsigma3(engy));
	max_sigmav_e = max(max_sigmav_e,sqrt(2*1.602E-19*engy/m[ecolsp])*(nitsigma1(engy)+nitsigma2(engy)+nitsigma3(engy)));
      }
      //printf("%s \n","col_prob_e");
      //printf("max_sigmav_e = %E",max_sigmav_e);
      //printf("max_sigmav_e = %E\n",max_sigmav_e);
      col_prob_e = 1-exp(-max_sigmav_e*gden*dt*sp_k[ecolsp]);
      //printf("col_prob_e = %E\n",col_prob_e);
      ecol_extra = 0.5;
    }

    //printf("%s \n","null collision prob icollisional");
    if (icollisional) {
      max_sigmav_i = 0.0;
      for (engy=0;engy<NIMAX;engy+=0.1)
	max_sigmav_i=max(max_sigmav_i,sqrt(2*1.602E-19*engy/m[icolsp])*(nitsigma4(engy)+nitsigma5(engy)));
      
      //printf("max_sigmav_i = %E",max_sigmav_i);
      col_prob_i = 1-exp(-max_sigmav_i*gden*dt*sp_k[icolsp]);
      icol_extra = 0.5;
    }
    init_flag = 0;
  }

  //printf("%s \n","About to collide");
  /*******************************/
  /* Electron collisions with N2 */
  if (ecollisional && isp==ecolsp) {
    /**** First Clear the diagnostics ****/
    if (theRunWithXFlag || theRunWithTec) {
      for (j=0; j<ng; j++)
	mccrate[0][j] = mccrate[1][j] = mccrate[2][j] = 0.0;
    }

    /**** Now do the collisions ****/
    ecol_extra += np[ecolsp]*col_prob_e;
    N = ecol_extra;
    ecol_extra -= N;
    
    nnp = np[ecolsp];
    for (j=0; j<N; j++) {
      index=nnp*frand();
      nnp--;
      
      temp = x[ecolsp][nnp];
      x[ecolsp][nnp] = x[ecolsp][index];
      x[ecolsp][index] = temp;

      temp = vx[ecolsp][nnp];
      vx[ecolsp][nnp] = vx[ecolsp][index];
      vx[ecolsp][index] = temp;

      temp = vy[ecolsp][nnp];
      vy[ecolsp][nnp] = vy[ecolsp][index];
      vy[ecolsp][index] = temp;

      temp = vz[ecolsp][nnp];
      vz[ecolsp][nnp] = vz[ecolsp][index];
      vz[ecolsp][index] = temp;
    }

    for (j=nnp; j<nnp+N; j++) {
      dum = (vx[ecolsp][j]*vx[ecolsp][j] + vy[ecolsp][j]*vy[ecolsp][j]+vz[ecolsp][j]*vz[ecolsp][j]);
      engy=Escale[ecolsp]*dum;
      vel=sqrt(dum);
      sigma_total = max_sigmav_e/(vel*vscale);
      //printf("sigma_total = %E\n",sigma_total);
      random=frand();

      /*********************************************************************/
      /* determine the type of collision and calculate new velocities, etc */
      
      /*******************************/
      /* if the collision is elastic */
      if (random <= (sum_sigma = nitsigma1(engy))/sigma_total) {
	N_elastic_el++;
	/* Normalize velocity */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;

	/* scatter the electron */
	nitnewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], x[ecolsp][j],1);
	/* collect diagnostics */
	if (theRunWithXFlag || theRunWithTec) {
	  s=x[ecolsp][j];
	  del = x[ecolsp][j] - s;
	  mccrate[0][s] += (!s) ? 2*(1-del) : 1-del;
	  mccrate[0][s+1]+= (s==nc-1) ? 2*del : del;
	}
      }


      /**********************************/
      /* if the collision is excitation */
      else if (engy >= extengy0 && random <= (sum_sigma += nitsigma2(engy))/sigma_total) {
	N_excite_el++;
	/* Normalize velocity */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;

	engy -= extengy0;
	vel = sqrt(fabs(engy)/Escale[ecolsp]);

	/* scatter the electron */
	nitnewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], x[ecolsp][j], 0);

	/* collect diagnostics */
	if (theRunWithXFlag || theRunWithTec) {
	  s=x[ecolsp][j];
	  del = x[ecolsp][j]-s;
	  mccrate[1][s] += (!s) ? 2*(1-del) : 1-del;
	  mccrate[1][s+1] += (s==nc-1) ? 2*del : del;
	}
      }
      

      /***************************************************************************/
      /* If the collision is ionization, add the released electron and ion first */
      else if(engy >= ionengy0 && random <= (sum_sigma +=nitsigma3(engy))/sigma_total) {
	N_ionize++;
	/* Normalize velocity */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	

	/*********************************************************************/
	/* subtract the ionization energy and partition the remaining energy */
	engy -= ionengy0;
	/* MIGHT HAVE TO CHANGE THE EXACT NUMBERS BELOW */
	rengy = 10.0*tan(frand()*atan(engy/20.0));
	engy -= rengy;

	/********************************/
	/* scatter the created electron */
	vel = sqrt(fabs(rengy)/Escale[ecolsp]);
	k = np[ecolsp];
	vx[ecolsp][k] = vx[ecolsp][j];
	vy[ecolsp][k] = vy[ecolsp][j];
	vz[ecolsp][k] = vz[ecolsp][j];
	nitnewvel(rengy, vel, &vx[ecolsp][k], &vy[ecolsp][k], &vz[ecolsp][k], x[ecolsp][j], 0);
	x[ecolsp][k] = x[ecolsp][j];


	/******************************************************/
	/* assign velocities to the ion created by ionization */
	k = np[ionsp];
	maxwellv(&vx[ionsp][k], &vy[ionsp][k], &vz[ionsp][k], vgth);
	x[ionsp][k] = x[ecolsp][j];
	s = x[ionsp][k];
	del = x[ionsp][k] - s;
	sp_n_mcc[ionsp][s] += 1. - del;
	sp_n_mcc[ionsp][s+1] += del;


	/*****************************************/
	/* finally scatter the incident electron */
	vel = sqrt(fabs(engy)/Escale[ecolsp]);
	nitnewvel(engy,vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], x[ecolsp][j],0);
	

	/*************************************************************/
	/* Check for number of particles exceeding the maximum value */
	if (++np[ionsp] > maxnp[ionsp] || ++np[ecolsp] > maxnp[ecolsp]) {
	  puts("ADJUST(Ionization) : too many particles. MUST EXIT!");
	  exit(-1);
	}

	/***********************/
	/* collect diagnostics */
	if (theRunWithXFlag || theRunWithTec) {
	  mccrate[2][s] += (!s) ? 2*(1-del) : 1-del;
	  mccrate[2][s+1] += (s== nc-1) ? 2*del : del;
	}
      }
    }
  }
	
  //printf("N2 collisions \n");
    /***************************/
    /* N2+ collisions          */
    /* Ion-Neutral collisions  */
    if (icollisional && isp==icolsp) {
      //printf("%d \n",icollisional);
      /* Clear the diagnostics */
      if (theRunWithXFlag || theRunWithTec) {
	for (j=0 ; j<ng; j++)
	  mccrate[3][j] = mccrate[4][j] = 0;
      }

      /* Now do the collisions */
      icol_extra += np[icolsp]*col_prob_i;
      N = icol_extra;
      icol_extra -= N;
      
      nnp = np[icolsp];
      for(j=0; j<N; j++) {
	index = nnp*frand();
	nnp--;
	
	temp = x[icolsp][nnp];
	x[icolsp][nnp] = x[icolsp][index];
	x[icolsp][index] = temp;

	temp = vx[icolsp][nnp];
	vx[icolsp][nnp] = vx[icolsp][index];
	vx[icolsp][index] = temp;

	temp = vy[icolsp][nnp];
	vy[icolsp][nnp] = vy[icolsp][index];
	vy[icolsp][index] = temp;

	temp = vz[icolsp][nnp];
	vz[icolsp][nnp] = vz[icolsp][index];
	vz[icolsp][index] = temp;

      }

      
      for (j=nnp; j<nnp+N; j++) {
	maxwellv(&vneutx, &vneuty, &vneutz, vgth);
	vx[icolsp][j] -= vneutx;
	vy[icolsp][j] -= vneuty;
	vz[icolsp][j] -= vneutz;
	dum = (vx[icolsp][j]*vx[icolsp][j] + vy[icolsp][j]*vy[icolsp][j] + vz[icolsp][j]*vz[icolsp][j]);
	
	if (dum) {
	  engy = Escale[icolsp]*dum;
	  vel = sqrt(dum);
	  sigma_total = max_sigmav_i/(vel*vscale);
	  random = frand();
	}
	else {
	  /* did not collide */
	  sigma_total = 1e7;
	  random = 0;
	}

	//printf("Doing actual collisions \n");
	/**********************************/
	/* if the collision is scattering */
	if (random <= (sum_sigma = nitsigma5(engy))/sigma_total) {
	  float up1, up2, up3, mag;
	  float r11, r12, r13, r21, r22, r23, r31, r32, r33;
	  N_elastic_ion++;
	  coschi = sqrt(frand());
	  sinchi = sqrt(fabs(1.0-coschi*coschi));
	  
	  phi1 = 2*M_PI*frand();
	  cosphi = cos(phi1);
	  sinphi = sin(phi1);
	  
	  r13 = vx[icolsp][j]/vel;
	  r23 = vy[icolsp][j]/vel;
	  r33 = vz[icolsp][j]/vel;

	  if (r33 == 1.0) { up1=0; up2=1; up3=0; }
	  else {up1=0; up2=0; up3=1;}
	  
	  r12 = r23*up3 - r33*up2;
	  r22 = r33*up1 - r13*up3;
	  r32 = r13*up2 - r23*up1;
	  mag = sqrt(r12*r12 + r22*r22 + r32*r32);
	  
	  r12 /= mag;
	  r22 /= mag;
	  r32 /= mag;

	  r11 = r22*r33 - r32*r23;
	  r22 = r32*r13 - r12*r33;
	  r32 = r12*r23 - r22*r13;

	  vel *= coschi;
	  vx[icolsp][j] = vel*coschi*(r11*sinchi*cosphi + r12*sinchi*sinphi + r13*coschi);
	  vy[icolsp][j] = vel*coschi*(r21*sinchi*cosphi + r22*sinchi*sinphi + r23*coschi);
	  vz[icolsp][j] = vel*coschi*(r31*sinchi*cosphi + r32*sinchi*sinphi + r33*coschi);

	  
	  /* collect diagnostics */
	  if (theRunWithXFlag || theRunWithTec) {
	    s = x[icolsp][j];
	    del = x[icolsp][j] - s;
	    mccrate[3][s] += (!s) ? 2*(1-del) : 1-del;
	    mccrate[3][s+1] += (s==nc-1) ? 2*del : del;
	  }
	}

	/*****************************/
	/* Charge Exchange Collision */
	else if (random <= (sum_sigma +=nitsigma4(engy))/sigma_total) {
	  N_exchange_ion++;
	  vx[icolsp][j] = vy[icolsp][j] = vz[icolsp][j] = 0.0;
	  
	  /* collect diagnostics */
	  if (theRunWithXFlag || theRunWithTec) {
	    s = x[icolsp][j];
	    del=x[icolsp][j] - s;
	    mccrate[4][s] += (!s) ? 2*(1-del) : 1-del;
	    mccrate[4][s+1] += (s==nc-1) ? 2*del : del ;
	  }
	}

	vx[icolsp][j] += vneutx;
	vy[icolsp][j] += vneuty;
	vz[icolsp][j] += vneutz;
      }
    }
}


/*****************************************************************************/

/* Function for elastic scattering of electrons */
void nitnewvel(float energy, float vel, float *vx, float *vy, float *vz, float xi, int e_flag)
{
  float phi1, cosphi, sinphi, coschi, sinchi, up1, up2, up3;
  float mag, r11, r12, r13, r21, r22, r23, r31, r32, r33;

  if (energy < 1E-30) coschi = 1;
  /* Eq (9) of Vahedi et.al 1995 */
  else coschi = (energy+2.0-2.0*pow(energy+1.0,frand()))/energy;
  sinchi = sqrt(fabs(1.0-coschi*coschi));

  /* Eq (10) of Vahedi et.al 1995 */
  phi1 = 2*M_PI*base2();
  cosphi = cos(phi1);
  sinphi = sin(phi1);

  if (e_flag) vel *= sqrt(1.0-2.0*m[ecolsp]*(1-coschi)/m[ionsp]);

  r13 = *vx;
  r23 = *vy;
  r33 = *vz;

  if(r33 == 1.0) { 
    up1 = 0; 
    up2 = 1; 
    up3 = 0;
  }
  else {
    up1 = 0;
    up2 = 0; 
    up3 = 1;
  }
  
  r12 = r23*up3 - r33*up2;
  r22 = r33*up1 - r13*up3;
  r32 = r13*up2 - r23*up1;

  mag = sqrt(r12*r12 + r22*r22 + r32*r32);

  r12 /= mag;
  r22 /= mag;
  r32 /= mag;
  
  r11 = r22*r33 - r32*r23;
  r21 = r32*r13 - r12*r33;
  r31 = r12*r23 - r22*r13;

  *vx = vel*(r11*sinchi*cosphi + r12*sinchi*sinphi + r13*coschi);
  *vy = vel*(r21*sinchi*cosphi + r22*sinchi*sinphi + r23*coschi);
  *vz = vel*(r31*sinchi*cosphi + r32*sinchi*sinphi + r33*coschi);
}

/* Collision Cross Sections */
/* e + Ar -> e + Ar (elastic scattering) */

float nitsigma1 (float energy)
{
  int i;
  int count;
  float els_sigma, alpha;
  static int init_flag=1;
  static float *elastic1, *elastic2;
  static float *data_energy, *data_sigma;
  static int maxcount;
  FILE *XSdatafile;

  /* Variables for bisection search */
  int pos_low, pos_high, pos, dummy;
  float slope;
  
  //printf("%s \n","initialization");
  /* Initialization */
  if (init_flag) {
    float engy;
    elastic1 = (float *) malloc(101*sizeof(float));
    elastic2 = (float *) malloc(NEMAX*sizeof(float));
    //printf("%s \n","elastic1,elastic2");
    /* Read in data from a file */
    XSdatafile = fopen("elastic_N2.dat","r");
    printf("%s \n","file opened");
    fscanf(XSdatafile,"%d",&maxcount);
    //printf("%d \n",maxcount);
    data_energy = (float *) malloc(maxcount*sizeof(float));
    data_sigma = (float *) malloc(maxcount*sizeof(float));
    for (count=0; count < maxcount; count++) {
      //printf("%d \n",count);
      fscanf(XSdatafile,"%g",&data_energy[count]);
      fscanf(XSdatafile,"%g",&data_sigma[count]);
    }
    fclose(XSdatafile);
    init_flag = 0;
  }

  //printf("%s \n","read XS data");
  /* Bisection Search for the given energy */

  /* Check if input energy is greater than MAX value of data */
  if (energy > data_energy[maxcount-1]) {
    slope = (data_sigma[maxcount-1]-data_sigma[maxcount-2])/(data_energy[maxcount-1]-data_energy[maxcount-2]);
    els_sigma = slope*(energy-data_energy[maxcount-1]) + data_sigma[maxcount-1];
    els_sigma = els_sigma*1E-20; /* unit conversion */
  }
  /* Check if input energy is lesser than MIN value of data */
  else if (energy < data_energy[0]) {
    slope = (data_sigma[1]-data_sigma[0])/(data_energy[1]-data_energy[0]);
    els_sigma = data_sigma[0]-slope*(data_energy[1]-energy);

    if (els_sigma < 0.) {
      els_sigma = 0.0;
    }
    else {
      els_sigma = els_sigma*1E-20 ; /* unit conversion */
    }
  }

  /* Do the bisection search after making sure that the input energy is in the data range */
  else {
    pos_low = 0;
    pos_high = maxcount-1;
    
    for (dummy=0; dummy<100; dummy++){
      pos = (pos_low+pos_high)/2;
      if (energy < data_energy[pos]){
	pos_high = pos;
      }
      else {
	pos_low = pos;
      }

      if (pos_low+1 >= pos_high) {
	if (energy < data_energy[pos_high]){
	  pos = pos_low;
	}
	else {
	  pos = pos_high;
	}
	break; /* break from the 'for' loop once you have narrowed down on the bin */
      }
    }
  }
  slope = (data_sigma[pos+1]-data_sigma[pos])/(data_energy[pos+1]-data_energy[pos]);
  els_sigma = slope*(energy-data_energy[pos]) + data_sigma[pos];
  els_sigma = els_sigma*1E-20 ; /* unit conversion */
  //printf("els_sigma = %E \n",els_sigma);
  return(els_sigma);
}
/* Excitation reaction (21 of Zhang et. al paper) */
float nitsigma2 (float energy)
{
  int i;
  int count;
  float exc_sigma, alpha;
  static float exc_threshold;
  static int init_flag=1;
  static float *excit;
  static float *data_energy, *data_sigma;
  static int maxcount;
  FILE *XSdatafile;

  /* Variables for bisection search */
  int pos_low, pos_high, pos, dummy;
  float slope;
  
  /* Initialization */
  if (init_flag) {
    float engy;
    excit = (float *) malloc(101*sizeof(float));
    exc_threshold = 11.03; //initialization of the threshold energy for excitation
    /* Read in data from a file */
    XSdatafile = fopen("excite_N2.dat","r");
    fscanf(XSdatafile,"%d",&maxcount);
    data_energy = (float *) malloc(maxcount*sizeof(float));
    data_sigma = (float *) malloc(maxcount*sizeof(float));
    for (count=0; count < maxcount; count++) {
      fscanf(XSdatafile,"%g",&data_energy[count]);
      fscanf(XSdatafile,"%g",&data_sigma[count]);
    }
    fclose(XSdatafile);
    init_flag = 0;
  }

  /* Bisection Search for the given energy */

  /* Check if input energy is greater than the threshold */
  if (energy < exc_threshold) {
    exc_sigma = 0.0;
  }
  /* Check if input energy is greater than MAX value of data */
  else {
  if (energy > data_energy[maxcount-1]) {
    slope = (data_sigma[maxcount-1]-data_sigma[maxcount-2])/(data_energy[maxcount-1]-data_energy[maxcount-2]);
    exc_sigma = slope*(energy-data_energy[maxcount-1]) + data_sigma[maxcount-1];
    exc_sigma = exc_sigma*1E-20; /* unit conversion */
  }
  /* Check if input energy is lesser than MIN value of data */
  else if (energy < data_energy[0]) {
    slope = (data_sigma[1]-data_sigma[0])/(data_energy[1]-data_energy[0]);
    exc_sigma = data_sigma[0]-slope*(data_energy[1]-energy);

    if (exc_sigma < 0.) {
      exc_sigma = 0.0;
    }
    else {
      exc_sigma = exc_sigma*1E-20 ; /* unit conversion */
    }
  }

  /* Do the bisection search after making sure that the input energy is in the data range */
  else {
    pos_low = 0;
    pos_high = maxcount-1;
    
    for (dummy=0; dummy<100; dummy++){
      pos = (pos_low+pos_high)/2;
      if (energy < data_energy[pos]){
	pos_high = pos;
      }
      else {
	pos_low = pos;
      }

      if (pos_low+1 >= pos_high) {
	if (energy < data_energy[pos_high]){
	  pos = pos_low;
	}
	else {
	  pos = pos_high;
	}
	break; /* break from the 'for' loop once you have narrowed down on the bin */
      }
    }
  }
  slope = (data_sigma[pos+1]-data_sigma[pos])/(data_energy[pos+1]-data_energy[pos]);
  exc_sigma = slope*(energy-data_energy[pos]) + data_sigma[pos];
  }
  exc_sigma = exc_sigma*1E-20 ; /* unit conversion */
  //printf("exc_sigma = %E \n",exc_sigma);
  return(exc_sigma);
}

/* Ionization : e + N2 -> e + e + N2+ */ 
float nitsigma3 (float energy)
{
  int i;
  int count;
  float ion_sigma, alpha;
  static float ion_threshold;
  static int init_flag=1;
  static float *ioniz;
  static float *data_energy, *data_sigma;
  static int maxcount;
  FILE *XSdatafile;

  /* Variables for bisection search */
  int pos_low, pos_high, pos, dummy;
  float slope;
  
  /* Initialization */
  if (init_flag) {
    float engy;
    ioniz = (float *) malloc(101*sizeof(float));
    ion_threshold = 15.58; //initialization of the threshold energy for ionization
    /* Read in data from a file */
    XSdatafile = fopen("ionization_N2.dat","r");
    fscanf(XSdatafile,"%d",&maxcount);
    data_energy = (float *) malloc(maxcount*sizeof(float));
    data_sigma = (float *) malloc(maxcount*sizeof(float));
    for (count=0; count < maxcount; count++) {
      fscanf(XSdatafile,"%g",&data_energy[count]);
      fscanf(XSdatafile,"%g",&data_sigma[count]);
    }
    fclose(XSdatafile);
    init_flag = 0;
  }

  /* Bisection Search for the given energy */

  /* Check if input energy is greater than the threshold */
  if (energy < ion_threshold) {
    ion_sigma = 0.0;
  }
  /* Check if input energy is greater than MAX value of data */
  else {
  if (energy > data_energy[maxcount-1]) {
    slope = (data_sigma[maxcount-1]-data_sigma[maxcount-2])/(data_energy[maxcount-1]-data_energy[maxcount-2]);
    ion_sigma = slope*(energy-data_energy[maxcount-1]) + data_sigma[maxcount-1];
    ion_sigma = ion_sigma*1E-20; /* unit conversion */
  }
  /* Check if input energy is lesser than MIN value of data */
  else if (energy < data_energy[0]) {
    slope = (data_sigma[1]-data_sigma[0])/(data_energy[1]-data_energy[0]);
    ion_sigma = data_sigma[0]-slope*(data_energy[1]-energy);

    if (ion_sigma < 0.) {
      ion_sigma = 0.0;
    }
    else {
      ion_sigma = ion_sigma*1E-20 ; /* unit conversion */
    }
  }

  /* Do the bisection search after making sure that the input energy is in the data range */
  else {
    pos_low = 0;
    pos_high = maxcount-1;
    
    for (dummy=0; dummy<100; dummy++){
      pos = (pos_low+pos_high)/2;
      if (energy < data_energy[pos]){
	pos_high = pos;
      }
      else {
	pos_low = pos;
      }

      if (pos_low+1 >= pos_high) {
	if (energy < data_energy[pos_high]){
	  pos = pos_low;
	}
	else {
	  pos = pos_high;
	}
	break; /* break from the 'for' loop once you have narrowed down on the bin */
      }
    }
  }
  slope = (data_sigma[pos+1]-data_sigma[pos])/(data_energy[pos+1]-data_energy[pos]);
  ion_sigma = slope*(energy-data_energy[pos]) + data_sigma[pos];
  }
  ion_sigma = ion_sigma*1E-20 ; /* unit conversion */
  //printf("ion_sigma = %E \n",ion_sigma);
  return(ion_sigma);
}

/* Charge Exchange : N2 + N2+ -> N2+ + N2 */
float nitsigma4 (float energy)
{
  int i;
  int count;
  float cex_sigma, alpha;
  static float cex_threshold;
  static int init_flag=1;
  static float *data_energy, *data_sigma;
  static int maxcount;
  FILE *XSdatafile;

  /* Variables for bisection search */
  int pos_low, pos_high, pos, dummy;
  float slope;
  
  /* Initialization */
  /*  if (init_flag) {
    cex_threshold = 4.0; //initialization of the threshold energy for excitation
    /* Read in data from a file */
  /*    XSdatafile = fopen("elastic_N2.dat","r");
    fscanf(XSdatafile,"%d",maxcount);
    data_energy = (float *) malloc(maxcount*sizeof(float));
    data_sigma = (float *) malloc(maxcount*sizeof(float));
    for (count=0; count < maxcount; count++) {
      fscanf(XSdatafile,"%g",data_energy[count]);
      fscanf(XSdatafile,"%g",data_sigma[count]);
    }
    fclose(XSdatafile);
    }*/

  /* Bisection Search for the given energy */

  /* Check if input energy is greater than the threshold */
  /*if (energy < cex_threshold) {
    cex_sigma = 0.0;
    }*/
  /* Check if input energy is greater than MAX value of data */
  /*else if (energy > data_energy[maxcount-1]) {
    slope = (data_sigma[maxcount-1]-data_sigma[maxcount-2])/(data_energy[maxcount-1]-data_energy[maxcount-2]);
    cex_sigma = slope*(energy-data_energy[maxcount-1]) + data_sigma[maxcount-1];
    cex_sigma = cex_sigma*1E-20; /* unit conversion */
  /*}
  /* Check if input energy is lesser than MIN value of data */
  /*else if (energy < data_energy[0]) {
    slope = (data_sigma[1]-data_sigma[0])/(data_energy[1]-data_energy[0]);
    /*cex_sigma = data_sigma[0]-slope*(data_energy[1]-energy);

    if (cex_sigma < 0.) {
      cex_sigma = 0.0;
    }
    else {
      cex_sigma = cex_sigma*1E-20 ; /* unit conversion */
      /*}
  }

  /* Do the bisection search after making sure that the input energy is in the data range */
  /*else {
    pos_low = 0;
    /*pos_high = maxcount-1;
    
    for (dummy=0; dummy<100; dummy++){
      pos = (pos_low+pos_high)/2;
      if (energy < data_energy[pos]){
	pos_high = pos;
      }
      else {
	pos_low = pos;
      }

      if (pos_low+1 >= pos_high) {
	if (energy < data_energy[pos_high]){
	  pos = pos_low;
	}
	else {
	  pos = pos_high;
	}
	break; /* break from the 'for' loop once you have narrowed down on the bin */
	/*}
    }
    /*}
  slope = (data_sigma[pos+1]-data_sigma[pos])/(data_energy[pos+1]-data_energy[pos]);
  /*cex_sigma = slope*(energy-data_energy[pos]) + data_sigma[pos];
  cex_sigma = cex_sigma*1E-20 ; /* unit conversion */

  cex_sigma = 3E-19;
  return(cex_sigma);
}

/* Scattering N2 + N2+ -> N2 + N2+ */
float nitsigma5 (float energy)
{
  int i;
  int count;
  float ionels_sigma, alpha;
  static float ionels_threshold;
  static int init_flag=1;
  static float *data_energy, *data_sigma;
  static int maxcount;
  FILE *XSdatafile;

  /* Variables for bisection search */
  int pos_low, pos_high, pos, dummy;
  float slope;
  
  /* Initialization */
  /*if (init_flag) {
    ionels_threshold = 4.0; //initialization of the threshold energy for excitation
    /* Read in data from a file */
  /*  XSdatafile = fopen("ionscattering_N2.dat","r");
    fscanf(XSdatafile,"%d",maxcount);
    data_energy = (float *) malloc(maxcount*sizeof(float));
    data_sigma = (float *) malloc(maxcount*sizeof(float));
    for (count=0; count < maxcount; count++) {
      fscanf(XSdatafile,"%g",data_energy[count]);
      fscanf(XSdatafile,"%g",data_sigma[count]);
    }
    fclose(XSdatafile);
    }*/

  /* Bisection Search for the given energy */
  /* Check if input energy is greater than the threshold */
  /*if (energy < ionels_threshold) {
    ionels_sigma = 0.0;
    }*/
  /* Check if input energy is greater than MAX value of data */
  /*if (energy > data_energy[maxcount-1]) {
    slope = (data_sigma[maxcount-1]-data_sigma[maxcount-2])/(data_energy[maxcount-1]-data_energy[maxcount-2]);
    ionels_sigma = slope*(energy-data_energy[maxcount-1]) + data_sigma[maxcount-1];
    ionels_sigma = ionels_sigma*1E-20; /* unit conversion */
  /*}
  /* Check if input energy is lesser than MIN value of data */
/*else if (energy < data_energy[0]) {
    slope = (data_sigma[1]-data_sigma[0])/(data_energy[1]-data_energy[0]);
    ionels_sigma = data_sigma[0]-slope*(data_energy[1]-energy);

    if (ionels_sigma < 0.) {
      ionels_sigma = 0.0;
    }
    else {
      ionels_sigma = ionels_sigma*1E-20 ; /* unit conversion */
/*}
  }

  /* Do the bisection search after making sure that the input energy is in the data range */
/*else {
    pos_low = 0;
    pos_high = maxcount-1;
    
    for (dummy=0; dummy<100; dummy++){
      pos = (pos_low+pos_high)/2;
      if (energy < data_energy[pos]){
	pos_high = pos;
      }
      else {
	pos_low = pos;
      }

      if (pos_low+1 >= pos_high) {
	if (energy < data_energy[pos_high]){
	  pos = pos_low;
	}
	else {
	  pos = pos_high;
	}
	break; /* break from the 'for' loop once you have narrowed down on the bin */
/*   }
    }
  }
  slope = (data_sigma[pos+1]-data_sigma[pos])/(data_energy[pos+1]-data_energy[pos]);
  ionels_sigma = slope*(energy-data_energy[pos]) + data_sigma[pos];
  ionels_sigma = ionels_sigma*1E-20 ; *//* unit conversion */ 

  ionels_sigma = 3E-19;
  return(ionels_sigma);
}
	    



	  


      
  
