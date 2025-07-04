#include "pdp1.h"
#include "xgrafix.h"

#define   NEMAX    500
#define   NIMAX    500

float asigma1(float), asigma2(float), asigma3(float), asigma4(float), asigma5(float);
float gden, vgth, extengy0, ionengy0;
int ecolsp, icolsp, ionsp;
void anewvel(float, float, float *, float *, float *, float, int);

/*********************************************************/
/* Monte Carlo Collisions for electron- and ion-neutrals */

void argonmcc(int isp)
{
  static int init_flag=1;
  static float ecol_extra, icol_extra;
  static float max_sigmav_e, max_sigmav_i;
  static float col_prob_e, col_prob_i;

  register int j, k, index;
  int N, nnp, s;
  float random, dum, vel, sigma_total,metal;
  float engy, rengy, del, phi1;
  float cosphi, sinphi, coschi, sinchi;
  float sum_sigma, temp, vneutx, vneuty, vneutz;

  FILE * rand_file ;
  
  if (init_flag) {
    if (gtemp)
      gden = NperTORR*pressure/(gtemp);  /* Calculating the gas density */
    else
      gden = NperTORR*pressure/FLOAT_MIN;
    ecolsp= ecollisional-1;
    icolsp= icollisional-1;
    ionsp = ionspecies-1;     /* Fixing the indices into the array of species */
    
    /*****************************************/
    /* Calculating the null collision prob.  */
    
    if(ecollisional) {
      if(ionsp>0) vgth= sqrt(gtemp/Escale[ionsp]);
      extengy0 = 11.55;
      ionengy0 = 15.76;
      
      max_sigmav_e = 0.0;
      for (engy=0; engy<NEMAX; engy += 0.1)
	max_sigmav_e= max(max_sigmav_e,
			  sqrt(2*1.602e-19*engy/m[ecolsp])
			  *(asigma1(engy)+asigma2(engy)+asigma3(engy)));

      //printf("max_sigmav_e = %E\n",max_sigmav_e);
      col_prob_e = 1 -exp(-max_sigmav_e*gden*dt*sp_k[ecolsp]);
      //printf("col_prob_e = %E\n",col_prob_e);
      ecol_extra = 0.5;
    }
    if(icollisional) {
      max_sigmav_i = 0.0;
      for (engy=0; engy<NIMAX; engy += 0.1)
	max_sigmav_i= max(max_sigmav_i,
			  sqrt(2*1.602e-19*engy/m[icolsp])
			  *(asigma4(engy)+asigma5(engy)));
      //printf("max_sigmav_i = %E\n",max_sigmav_i);
      col_prob_i = 1 -exp(-max_sigmav_i*gden*dt*sp_k[icolsp]);
      icol_extra = 0.5;
    }
    init_flag = 0;
  }
  
  /********************************/
  /* Electron collisions with Ar  */
  
  if(ecollisional && isp==ecolsp) {
    /**** First Clear the diagnostics  ******/
    if(theRunWithXFlag || theRunWithTec) {
      for (j=0; j<ng; j++)
	mccrate[0][j] = mccrate[1][j] = mccrate[2][j] = 0.0;
    }
    /**** Now do the collisions *****/
    ecol_extra += np[ecolsp]*col_prob_e;
    N = ecol_extra;
    ecol_extra -= N;
    
    nnp = np[ecolsp];
    for(j=0; j< N; j++) {
      index= nnp*frand();
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
    
    for(j=nnp; j<nnp+N; j++) {
      dum = (vx[ecolsp][j]*vx[ecolsp][j] +vy[ecolsp][j]*vy[ecolsp][j]
						 +vz[ecolsp][j]*vz[ecolsp][j]);
      engy= Escale[ecolsp]*dum;
      vel = sqrt(dum);
      sigma_total = max_sigmav_e/(vel*vscale);
      //printf("sigma_total = %E\n",sigma_total);
      random= frand();
      
      /*********************************************************************/
      /* determine the type of collision and calculate new velocities, etc */
      
      /*******************************/
      /* if the collision is elastic */
      metal = frand(); 
      if (metal > metal_conc){
      	ionengy0 = 15.76;
        //printf("gas loop\n");
      }
      else {
      	ionengy0 = 15.76/5 ; // no particular reason - just wanted it to be smaller than 15.76
        //printf("metal loop\n");
      }
      //ionengy0 = 15.76/5;
      //rand_file = fopen("rand_file.txt","a");
      //fprintf(rand_file,"%e \n",metal);
      //fclose(rand_file);
      //fprintf(rand_file,"%e %e %e\n",metal,metal_conc,ionengy0);
      //fclose(rand_file);
      if (random <= (sum_sigma =asigma1(engy))/sigma_total) {
	                        N_elastic_el++; 
				/* first normalize vel */
				vx[ecolsp][j] /= vel;
				vy[ecolsp][j] /= vel;
				vz[ecolsp][j] /= vel;
	
				/* scatter the electron */
				anewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], x[ecolsp][j], 1);
	
				/* collect diagnostics*/
				if(theRunWithXFlag || theRunWithTec) {
					s= x[ecolsp][j];
					del= x[ecolsp][j] - s;
					mccrate[0][s]  += (!s) ? 2*(1-del) : 1-del;
					mccrate[0][s+1]+= (s== nc-1) ? 2*del : del;
				}
      } 
      /**********************************/
      /* if the collision is excitation */
      
      else if (engy>= extengy0 && random <= (sum_sigma +=asigma2(engy))/sigma_total) {
                                N_excite_el++;
				/* first normalize vel */
				vx[ecolsp][j] /= vel;
				vy[ecolsp][j] /= vel;
				vz[ecolsp][j] /= vel;
	
				engy -= extengy0;
				vel   = sqrt(fabs(engy)/Escale[ecolsp]);
	
				/* scatter the electron */
				anewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], x[ecolsp][j], 0);
	
				/* collect diagnostics */
				if(theRunWithXFlag || theRunWithTec) {
					s= x[ecolsp][j];
					del= x[ecolsp][j] - s;
					mccrate[1][s]  += (!s) ? 2*(1-del) : 1-del;
					mccrate[1][s+1]+= (s== nc-1) ? 2*del : del;
				}
      } 
      /***************************************************************************/
      /* if the collision is ionization, add the released electron and ion first */
    
      else if(engy >= ionengy0 && random <= (sum_sigma +=asigma3(engy))/sigma_total) {
	/* increment number of ionization reactions by one */
	N_ionize++ ;
	/* first normalize vel */
	vx[ecolsp][j] /= vel;
	vy[ecolsp][j] /= vel;
	vz[ecolsp][j] /= vel;
	
	/********************************/
	/* subtract the ion. energy and */
	/* partition the remaining energy */
	
	engy -= ionengy0;
	rengy = 10.0*tan(frand()*atan(engy/20.0));
	engy -= rengy;
	
	/********************************/
	/* scatter the created electron */
	
	vel = sqrt(fabs(rengy)/Escale[ecolsp]);
	k = np[ecolsp];
	vx[ecolsp][k] = vx[ecolsp][j];
	vy[ecolsp][k] = vy[ecolsp][j];
	vz[ecolsp][k] = vz[ecolsp][j];
	anewvel(rengy, vel, &vx[ecolsp][k], &vy[ecolsp][k], &vz[ecolsp][k], x[ecolsp][j], 0);
	x[ecolsp][k] = x[ecolsp][j];
	
	/****************************************/
	/* assign velocities to the created ion */
	
	k = np[ionsp];
	maxwellv(&vx[ionsp][k], &vy[ionsp][k], &vz[ionsp][k], vgth);
	x[ionsp][k] = x[ecolsp][j];
	
	s = x[ionsp][k];
	del = x[ionsp][k] - s;
	sp_n_mcc[ionsp][s]  += 1. - del;
	sp_n_mcc[ionsp][s+1]+= del;
	
	/*****************************************/
	/* finally scatter the incident electron */
	
	vel = sqrt(fabs(engy)/Escale[ecolsp]);
	anewvel(engy, vel, &vx[ecolsp][j], &vy[ecolsp][j], &vz[ecolsp][j], x[ecolsp][j], 0);
	
	if(++np[ionsp] >maxnp[ionsp] || ++np[ecolsp] >maxnp[ecolsp]) {
	  puts("ADJUST(Ionization): too many particles. MUST EXIT!");
	  exit(-1);
	}
	
	/* collect diagnostics */
	if(theRunWithXFlag || theRunWithTec) {
	  mccrate[2][s]  += (!s) ? 2*(1-del) : 1-del;
	  mccrate[2][s+1]+= (s== nc-1) ? 2*del : del;
	}
      }
    } 
  }
  
  /**************************************/
  /* Ar+ + Ar -> .....                  */
  
  if(icollisional && isp==icolsp) {
    /**** First Clear the diagnostics  ******/
    if(theRunWithXFlag || theRunWithTec) {
      for (j=0; j<ng; j++)
				mccrate[3][j] = mccrate[4][j] = 0.0;
    }
    /**** Now do the collisions *****/
    icol_extra += np[icolsp]*col_prob_i;
    N = icol_extra;
    icol_extra -= N;
    
    nnp = np[icolsp];
    for(j=0; j< N; j++) {
      index= nnp*frand();
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
      dum = (vx[icolsp][j]*vx[icolsp][j] +vy[icolsp][j]*vy[icolsp][j]
						 +vz[icolsp][j]*vz[icolsp][j]);
			if (dum){
				engy= Escale[icolsp]*dum;
				vel = sqrt(dum);
				sigma_total = max_sigmav_i/(vel*vscale);
				random= frand();
			}
			else{  /* did collide */
				sigma_total = 1e7;
				random = 0;
			}
      
      /**********************************/
      /* if the collision is scattering */
		
      if (random <= (sum_sigma =asigma5(engy))/sigma_total) {
				float up1, up2, up3, mag;
				float r11, r12, r13, r21, r22, r23, r31, r32, r33;
	
                                N_elastic_ion++;
				coschi= sqrt(frand());
				sinchi= sqrt(fabs(1.0 -coschi*coschi));
	
				phi1  = 2*M_PI*frand();
				cosphi= cos(phi1);
				sinphi= sin(phi1);
	
				r13 = vx[icolsp][j]/vel;
				r23 = vy[icolsp][j]/vel;
				r33 = vz[icolsp][j]/vel;
	
				if(r33 == 1.0) { up1= 0;  up2= 1;  up3= 0; }
				else           { up1= 0;  up2= 0;  up3= 1; }
	
				r12 = r23*up3 -r33*up2;
				r22 = r33*up1 -r13*up3;
				r32 = r13*up2 -r23*up1;
				mag = sqrt(r12*r12 + r22*r22 + r32*r32);
	
				r12/= mag;
				r22/= mag;
				r32/= mag;
	
				r11 = r22*r33 -r32*r23;
				r21 = r32*r13 -r12*r33;
				r31 = r12*r23 -r22*r13;
	
				vel*= coschi;
				vx[icolsp][j]= vel*coschi*(r11*sinchi*cosphi +r12*sinchi*sinphi +r13*coschi);
				vy[icolsp][j]= vel*coschi*(r21*sinchi*cosphi +r22*sinchi*sinphi +r23*coschi); 
				vz[icolsp][j]= vel*coschi*(r31*sinchi*cosphi +r32*sinchi*sinphi +r33*coschi);
	
				/* collect diagnostics */
				if(theRunWithXFlag || theRunWithTec) {
					s= x[icolsp][j];
					del= x[icolsp][j] - s;
					mccrate[3][s]  += (!s) ? 2*(1-del) : 1-del;
					mccrate[3][s+1]+= (s== nc-1) ? 2*del : del;
				}
      }
      /***************************************/
      /* if the collision is charge exchange */
      
      else if (random <= (sum_sigma +=asigma4(engy))/sigma_total) {
                                N_exchange_ion++;
				vx[icolsp][j] = vy[icolsp][j] = vz[icolsp][j] = 0.0;
	
				/* collect diagnostics */
				if(theRunWithXFlag || theRunWithTec) {
					s= x[icolsp][j];
					del= x[icolsp][j] - s;
					mccrate[4][s]  += (!s) ? 2*(1-del) : 1-del;
					mccrate[4][s+1]+= (s== nc-1) ? 2*del : del;
				}
      }
      vx[icolsp][j] += vneutx;
      vy[icolsp][j] += vneuty;
      vz[icolsp][j] += vneutz;
    }
  }
}
/**************************************************************/

void anewvel(float energy, float vel, float *vx, float *vy, float *vz, float xi, int e_flag)
{
  float phi1, cosphi, sinphi, coschi, sinchi, up1, up2, up3;
  float mag, r11, r12, r13, r21, r22, r23, r31, r32, r33;
  
  if(energy < 1e-30)  coschi = 1;
  else  coschi = (energy +2.0 -2.0*pow(energy+1.0, frand()))/energy;
  sinchi= sqrt(fabs(1.0 -coschi*coschi));
  
  phi1  = 2*M_PI*base2();
  cosphi= cos(phi1);
  sinphi= sin(phi1);
  
  if(e_flag) vel *= sqrt(1.0 -2.0*m[ecolsp]*(1-coschi)/m[ionsp]);
  
  r13 = *vx;
  r23 = *vy;
  r33 = *vz;
  
  if(r33 == 1.0) { up1= 0;  up2= 1;  up3= 0; }
  else           { up1= 0;  up2= 0;  up3= 1; }
  
  r12 = r23*up3 -r33*up2;
  r22 = r33*up1 -r13*up3;
  r32 = r13*up2 -r23*up1;
  mag = sqrt(r12*r12 + r22*r22 + r32*r32);
  
  r12/= mag;
  r22/= mag;
  r32/= mag;
  
  r11 = r22*r33 -r32*r23;
  r21 = r32*r13 -r12*r33;
  r31 = r12*r23 -r22*r13;
  
  *vx= vel*(r11*sinchi*cosphi +r12*sinchi*sinphi +r13*coschi);
  *vy= vel*(r21*sinchi*cosphi +r22*sinchi*sinphi +r23*coschi); 
  *vz= vel*(r31*sinchi*cosphi +r32*sinchi*sinphi +r33*coschi);
}

/***********************************/
/*  e + Ar -> e + Ar  Elastic      */

float asigma1(float energy)
{
  int i;
  float els_sigma, alpha;
  static int init_flag=1;
  static float *elastic1, *elastic2;

  /******  Initialization   *********/
  if(init_flag) {
    float engy;
    elastic1 = (float *) malloc(101*sizeof(float));
    elastic2 = (float *) malloc(NEMAX*sizeof(float));
    
    /*****  With no Ramseur min.  ******/
    /*
			for(i=0; i<100; i++) elastic1[i]= 1.4e-19;
			for(i=0; i<NEMAX; i++)
			{
			engy = i;
			if(engy <= 20.) elastic2[i]= 1.4e-19;
			else
			elastic2[i]= 9.07e-19*pow(engy, 1.55)*pow(engy+70.0, 1.10)/pow(14.+engy, 3.25);
			}
			*/
    
    /*****  With the Ramseur min.  *****/
    for(i=0; i<101; i++) {
      engy = .01*i;
      if(engy < 0.2) elastic1[i]= 1./pow(10.0, 19.0 +engy/.11);
      else elastic1[i]= 9.07e-19*pow(engy, 1.55)*pow(engy+70.0, 1.10)/pow(14.+engy, 3.25);
    }
    for(i=0; i<NEMAX; i++) {
      engy = i;
      elastic2[i]= 9.07e-19*pow(engy, 1.55)*pow(engy+70.0, 1.10)/pow(14.+engy, 3.25);
    }
    init_flag =0;
  }
  /****************************/
  
  if(energy < 1.0) {
    i= 100*energy;  
    alpha= 100*energy -i;
    els_sigma = elastic1[i] + alpha*(elastic1[i+1] - elastic1[i]);
  }
  else {
    i= energy; 
    if(i < NEMAX-1) {
      alpha= energy -i;
      els_sigma = elastic2[i] + alpha*(elastic2[i+1] - elastic2[i]);
    }
    else
      els_sigma = elastic2[NEMAX-1];
  }
  return(els_sigma);
}

/************************************/
/*  e + Ar -> e + Ar  Excitation    */

float asigma2(float energy)
{
  int i;
  float exc_sigma, alpha;
  static int init_flag=1;
  static float *excit;

  /******  Initialization   *********/
  if(init_flag)
		{
			float engy;
			excit = (float *) malloc(NEMAX*sizeof(float));

			for(i=0; i<NEMAX; i++) {
				engy = i;
				if(engy < 11.55) excit[i] = 0;
				else
					excit[i] = (3.85116e-19*log(engy/3.4015) -4.85227e-19)/engy;
			}
			init_flag =0;
		}
  /****************************/
  
  i= energy; 
  if(i < NEMAX-1) {
    alpha= energy -i;
    exc_sigma = excit[i] + alpha*(excit[i+1] - excit[i]);
  }
  else
    exc_sigma = excit[NEMAX-1];
  
  return(exc_sigma);
}

/************************************/
/*  e + Ar -> e + e + Ar+  Ion.     */

float asigma3(float energy)
{
  int i;
  float ion_sigma, alpha;
  static int init_flag=1;
  static float *ioniz;
  
  /******  Initialization   *********/
  if(init_flag) {
    float engy;
    ioniz = (float *) malloc(NEMAX*sizeof(float));

    for(i=0; i<NEMAX; i++) {
      engy = i;
      if(engy < 15.76) ioniz[i] = 0;
      else
				ioniz[i] = (1.3596e-18/engy)*log((engy +120.0/engy)/15.76)*
					(atan((engy*engy -9.76*engy +2.4)/(20.6*engy +206)) +
					 atan((2*engy -80.0)/(10.3*engy +103.0)));
    }
    init_flag =0;
  }
  /****************************/
  
  i= energy; 
  if(i < NEMAX-1) {
    alpha= energy -i;
    ion_sigma = ioniz[i] + alpha*(ioniz[i+1] - ioniz[i]);
  }
  else
    ion_sigma = ioniz[NEMAX-1];
  
  return(ion_sigma);
}

/************************************/
/*  Ar + Ar+ -> Ar+ + Ar  Charge X  */

float asigma4(float energy)
{
  if(energy > 4.0) return(2.0e-19 +5.5e-19/sqrt(energy));
  return(-2.95e-19*sqrt(energy) +10.65e-19);
}

/***********************************/
/*  Ar + Ar+ -> Ar + Ar+   Scat.   */

float asigma5(float energy)
{
  if(energy > 4.0) return(1.8e-19 +4.0e-19/sqrt(energy));
  return(-2.0e-19*sqrt(energy) +7.8e-19);
}

/**************************************************************/

void makethefile(void)
{
  float e;
  FILE *DMPFile;
  
  DMPFile = fopen("arxsects.dat", "w"); 
  
  for(e=0; e<10; e+=0.01)
    fprintf(DMPFile, "%e %e %e %e %e %e\n", e+1e-30, asigma1(e)+1e-30,
						asigma2(e)+1e-30, asigma3(e)+1e-30, asigma4(e)+1e-30, asigma5(e)+1e-30);
  
  for(e=10; e<100; e+=0.1)
    fprintf(DMPFile, "%e %e %e %e %e %e\n", e+1e-30, asigma1(e)+1e-30,
						asigma2(e)+1e-30, asigma3(e)+1e-30, asigma4(e)+1e-30, asigma5(e)+1e-30);
  
  fclose(DMPFile);
}

/**************************************************************/
