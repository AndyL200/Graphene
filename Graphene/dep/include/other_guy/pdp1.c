#include "pdp1.h"
#include "xgrafix.h"

void display_title(void);
void InitWindows(void);
void tecplot_output(void);

void main(int argc, char *argv[])
{
  display_title();        /* Display XPDP1 title           */
  XGInit(argc, argv, &t); /* Initialize XGrafix stuff      */
  FILE * TimeHistory;
  FILE * second ;
  if (WasDumpFileGiven){
    TimeHistory = fopen("history.plt","a");
    fclose (TimeHistory);
    second = fopen("secondary_el.plt","a");
    fclose (second);
    second = fopen("secondary_ion.plt","a");
    fclose(second);
  }
  else{
    TimeHistory = fopen("history.plt","w");
    fclose (TimeHistory);
    second = fopen("secondary_el.plt","w");
    fclose (second);
    second = fopen("secondary_ion.plt","w");
    fclose(second);
  }
	
  //printf("theRunWithXFlag : %d",theRunWithXFlag);
  start();                /* Allocate arrays and initialize */
  //printf("theRunWithXFlag : %d",theRunWithXFlag);
  if(theRunWithXFlag) {
    mccdiag_init();       /* Initialize MCC rate diagnostic */
    InitWindows();        /* Initialize diagnostic windows  */
  }
  //printf("%s \n","Reached just before mccdiag_init()");
  if (theRunWithTec) mccdiag_init(); // changed by Venkat on Jan 6th 2010
  setrho();               /* Set initial charge density    */
  fields();               /* Initialize field arrays       */
  if(theRunWithXFlag || theRunWithTec) // changed by Venkat on Jan 6th 2010
    history();            /* Initialize history arrays     */
  TimeHistory = fopen("history.plt","a");
  fprintf(TimeHistory,"%e\t%d\t%d\t%d\t%e\t%e\t%e\t%e\t%d\t%d\t%d\t%d\t%d\n",t,np[0],np[1],N_ionize,exti,phi[0],e[nc],j_fe,N_fe,N_elastic_el,N_excite_el,N_elastic_ion,N_exchange_ion);
  fclose (TimeHistory);
  XGStart();              /* Start XGrafix main loop       */
}

/***********************************************************/
/*  The main physics loop                                  */

void XGMainLoop(void)
{
  register int j, isp;
  int DiagFlag;
  float frac;
  FILE * TimeHistory;
  t += dt; 
  //printf("%E \n",t);
  //printf("timestep : %d %d\n",it[0],it[1]);
  //fprintf(TimeHistory,"VARIABLES=\n\"Time\"\n\"Ne\"\n\"Ni\"\nZONE I= %d, J=1, F= POINT\n",it[0]/100);
  if ((theRunWithXFlag || theRunWithTec)&&n_ave)
    DiagFlag=1;
  else
    DiagFlag=0;
  //printf("Print %g \n",metal_conc);
  //printf("\n  - Time step: %e  Electrons:  %d  Ions: %d \n",t, np[0],np[1]);
  //printf("\n n_ave = %d \n",n_ave);
  /*no subcycling loop*/
  
  //printf("%s \n","About to enter for loop");
  for(isp=0; isp< nsp; isp++) {
    it[isp]++;                      /* Advance time for the species isp         */
    //printf("About to move %d \n",isp);
    (*moveptr)(isp, DiagFlag);      /* Advance position and velocity            */
    //printf("Adjust %d \n",isp);
    adjust(isp);                    /* Remove particles that cross boundaries   */
    //printf("Adjust done %d \n",isp);
    //printf("About to mccptr %d \n",isp);
    (*mccptr) (isp);                /* Monte Carlo collisions for species isp   */
    //printf("mccptr done %d \n",isp);
  }
  for(isp=0; isp< nsp; isp++) {
    //printf("%d \n",isp);
    gather(isp);                    /* Assign charge densities to the grid      */
    //printf("Gather done \n");
    
    for (j=0; j<ng; j++)
      //printf("Something \n");
      sp_n[isp][j]=sp_n_k[isp][j];
    //printf("Totally done %d\n",isp);
  }
  //printf("Everything done\n");
  //printf("Before OLD_CODE \n");
  /* gather needs to be in a seperate loop because new particles might
     be created in mcc and adjust that might not have been weighted. */
  
  
  /*subcycling loop*/
  /*subcycling has bugs when used with mcc and secondaries*/
  /*the mcc is mostly right except when smoothing is used */
  /*the subcycling method is first order when using the explicit push*/
  /*ie the subcycling is for the implicit push */
  
#ifdef OLD_CODE  
  for(isp=0; isp< nsp; isp++) {
    printf("\n  - Old Code \n");
    if(!(k_count[isp]%sp_k[isp])) {
      it[isp]++;          /* Advance time for the species isp         */
      (*moveptr)(isp, DiagFlag);    /* Advance position and velocity            */
      adjust(isp);        /* Remove particles that cross boundaries   */
      (*mccptr) (isp);    /* Monte Carlo collisions for species isp   */
      gather(isp);        /* Assign charge densities to the grid      */
      k_count[isp]=0;     /* Reset species counter for subcycling     */
    }
    k_count[isp]++;
    frac = ((float)k_count[isp])/sp_k[isp];
    for (j=0; j<ng; j++)
      sp_n[isp][j]= (1- frac)*sp_n_0[isp][j]
				+frac*sp_n_k[isp][j] +sp_n_mcc[isp][j];
  }
#endif
  //printf("Before Fields \n");
  fields(); 
  //printf("Fields \n");

  //if(theRunWithXFlag) history(); original code
  if(theRunWithXFlag || theRunWithTec) history(); // changed by Venkat on Jan 6th 2010 to compute history even when code is run without X server (for tecplot output )
  if (remainder(it[0],1000) == 0){
    tecplot_output();
    //printf("output_written\n");
    TimeHistory = fopen("history.plt","a");
    //printf("history to be written\n");
    fprintf(TimeHistory,"%e\t%d\t%d\t%d\t%e\t%e\t%e\t%e\t%e\t%d\t%d\t%d\t%d\n",t,np[0],np[1],N_ionize,exti,phi[0],e[nc],j_fe,N_fe,N_elastic_el,N_excite_el,N_elastic_ion,N_exchange_ion);
    fclose (TimeHistory);
  }
  //printf("TimeHistory \n\n");
}

/***********************************************************/

void display_title(void)
{
  puts("\n\nXPDP1 Version 4.11");
  puts("Copyright (C) 1988-2002");
  puts("Regents of the University of California");
  puts("Plasma Theory and Simulation Group");
  puts("University of California - Berkeley\n");
}

/***************************************************************/
void tecplot_output(void)
{
  FILE * TecPlot;
  int i;
  TecPlot = fopen("output.plt","w");
  fprintf(TecPlot,"VARIABLES=\n\"x\"\n\"Potential\"\n\"Charge\"\n\"j_x_e\"\n\"j_x_i\"\n\"j_y_e\"\n\"j_y_i\"\n\"j_z_e\"\n\"j_z_i\"\n\"ne\"\n\"ni\"\n\"Uxe\"\n\"Uxi\"\n\"Uye\"\n\"Uyi\"\n\"Uze\"\n\"Uzi\"\n\"KEx_e\"\n\"KEx_i\"\n\"KEy_e\"\n\"KEy_i\"\n\"KEz_e\"\n\"KEz_i\"\n\"E\"\n\"(j.E)_e\"\n\"(j.E)_i\"\nZONE I= %d, J=1, F= POINT\n",ng);
  for (i=0;i<ng;++i){
    fprintf(TecPlot,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",i*dx,phi[i],rho[i],sp_j_x_show[0][i],sp_j_x_show[1][i],sp_j_y_show[0][i],sp_j_y_show[1][i],sp_j_z_show[0][i],sp_j_z_show[1][i],sp_n[0][i]*nc2p/area/dx,sp_n[1][i]*nc2p/area/dx,0.5*dx/dt*sp_u_x_show[0][i],0.5*dx/dt*sp_u_x_show[1][i],0.5*dx/dt*sp_u_y_show[0][i],0.5*dx/dt*sp_u_y_show[1][i],0.5*dx/dt*sp_u_z_show[0][i],0.5*dx/dt*sp_u_z_show[1][i],0.25*sp_ke_x_show[0][i],0.25*sp_ke_x_show[1][i],0.25*sp_ke_y_show[0][i],0.25*sp_ke_y_show[1][i],0.25*sp_ke_z_show[0][i],0.25*sp_ke_z_show[1][i],e[i],jdote_show[0][i],jdote_show[1][i]);
  }
  fclose (TecPlot);

  TecPlot = fopen("output_ave.plt","w");
  fprintf(TecPlot,"VARIABLES=\n\"X\"\n\"Time Ave n_e\"\n\"Time Ave n_i\"\n\"Time Ave Ux_e\"\n\"Time Ave Ux_i\"\n\"Time Ave Uy_e\"\n\"Time Ave Uy_i\"\n\"Time Ave Uz_e\"\n\"Time Ave Uz_i\"\n\"Time Ave Jx_e\"\n\"Time Ave Jx_i\"\n\"Time Ave Jy_e\"\n\"Time Ave Jy_i\"\n\"Time Ave Jz_e\"\n\"Time Ave Jz_i\"\n\"Time Ave KEx_e\"\n\"Time Ave KEx_i\"\n\"Time Ave KEy_e\"\n\"Time Ave KEy_i\"\n\"Time Ave KEz_e\"\n\"Time Ave KEz_i\"\n\"Time Ave T_e\"\n\"Time Ave T_i\"\n\"Time Ave Tx_e\"\n\"Time Ave Tx_i\"\n\"Time Ave Ty_e\"\n\"Time Ave Ty_i\"\n\"Time Ave Tz_e\"\n\"Time Ave Tz_i\"\nZONE I= %d, J=1, F= POINT\n",ng);
  for (i=0;i<ng;++i){
    fprintf(TecPlot,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",x_grid[i]*dx,sp_n_ave_show[0][i]*nc2p/area/dx,nc2p/area/dx*sp_n_ave_show[1][i],0.5*dx/dt*sp_u_x_ave_show[0][i],0.5*dx/dt*sp_u_x_ave_show[1][i],0.5*dx/dt*sp_u_y_ave_show[0][i],0.5*dx/dt*sp_u_y_ave_show[1][i],0.5*dx/dt*sp_u_z_ave_show[0][i],0.5*dx/dt*sp_u_z_ave_show[1][i],sp_j_x_ave_show[0][i],sp_j_x_ave_show[1][i],sp_j_y_ave_show[0][i],sp_j_y_ave_show[1][i],sp_j_z_ave_show[0][i],sp_j_z_ave_show[1][i],0.25*sp_ke_x_ave_show[0][i],0.25*sp_ke_x_ave_show[1][i],0.25*sp_ke_y_ave_show[0][i],0.25*sp_ke_y_ave_show[1][i],0.25*sp_ke_z_ave_show[0][i],0.25*sp_ke_z_ave_show[1][i],0.5*T_ave_show[0][i],0.5*T_ave_show[1][i],0.5*Tx_ave_show[0][i],0.5*Tx_ave_show[1][i],0.5*Ty_ave_show[0][i],0.5*Ty_ave_show[1][i],0.5*Tz_ave_show[0][i],0.5*Tz_ave_show[1][i]);
  }
  fclose (TecPlot);
}
/***************************************************************/
