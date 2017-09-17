/* readinputfile.h -- This file is part of GNU archimedes

   Archimedes is a simulator for Submicron and Nanoscaled
   2D III-V Semiconductor Devices.

   Copyright (C) 2004-2011 Jean Michel D. Sellier
   <jeanmichel.sellier@gmail.com>
   <jsellier@purdue.edu>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULfAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


// ##################################################################
// Created on 10 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 06 september 2011, Carry le Rouet, France, J.M.Sellier
// ##################################################################

// Parse the script describing the device.
//
// For the syntax see the manual

#include "vec.h"
#include <stdio.h>

void read_input_file(FILE *fp) {
    char s[180];
    double num,num0;
    int ini,fin;
    int LXflag=0, LYflag=0;
    int transportflag=0;


    // setting defaults using configuration object
    g_config->optical_phonon_scattering = ON;     // Optical phonons scattering ON by default
    g_config->acoustic_phonon_scattering = ON;    // Acoustic phonons scattering ON by default
    g_config->impurity_scattering = ON;           // Charged impurity scattering ON by default
    g_config->neutral_impurity_scattering = OFF;  // Neutral impurity scattering OFF by default
    g_config->piezoelectric_scattering = ON;      // Piezoelectric scattering OFF by default
    g_config->electron_hole_scattering = OFF;     // Electron-hole scattering OFF by default
    g_config->conduction_band = KANE;
    g_config->qep_alpha = 0.5;
    g_config->qep_gamma = 1.0;
    g_config->qep_model = QEP_BOHM;
    g_config->save_mesh = OFF;
    g_config->simulation_model = MCE; // model_number
    g_config->time = 0.;
    g_config->tf = 5.0e-12;
    g_config->dt = 0.001e-12;
    g_config->tauw = 0.4e-12;
    g_config->faraday_flag = OFF;
    g_config->poisson_flag = ON;
    g_config->photon_energy = 0.;
    g_config->photoexcitation_flag = OFF;
    g_config->impurity_conc = 1e17; // cimp
    g_config->neutral_impurity_conc = 1e15; // cimp
    g_config->thomas_fermi_screening = OFF;
    g_config->lattice_temp = 300.0;
    g_config->particles_per_cell = 2500; // np1 aka statistical weight
    g_config->avg_steps = 500;
    g_config->avg_alpha = 1. / (real)g_config->avg_steps;
    g_config->quantum_flag = OFF;
    g_config->max_min_output = OFF; // maximini
    g_config->save_step_output = OFF; // savealways
    g_config->scattering_output = OFF;
    g_config->tracking_output = OFF;
    g_config->tracking_mod = 1000;
    g_config->output_format = GNUPLOTFORMAT;
    g_config->load_initial_data = OFF; // leid_flag
    g_config->tcad_data = OFF;
    g_config->num_particles = 0;
    g_config->surface_bb_flag = OFF;
    g_config->surface_bb_direction = direction_t.LEFT;
    g_config->surface_bb_delV = 0.;
    g_config->constant_efield_flag = OFF;


    g_mesh->nx = NXM - 1;
    g_mesh->ny = NYM - 1;
    g_mesh->dx = 0.;
    g_mesh->dy = 0.;
    g_mesh->width = 0.;
    g_mesh->height = 0.;

    for(int i = 1; i <= g_mesh->nx + 1; ++i) {
        for(int j = 1; j <= g_mesh->ny + 1; ++j) {
            g_mesh->nodes[i][j].qep = 0.;
            g_mesh->nodes[i][j].potential = 0.;
            g_mesh->nodes[i][j].efield = (Vec2){.x=0., .y=0.};
            g_mesh->nodes[i][j].magnetic_field = 0.;

            g_mesh->nodes[i][j].donor_conc = NI;
            g_mesh->nodes[i][j].acceptor_conc = NI;

            g_mesh->nodes[i][j].e.density = NI;
            g_mesh->nodes[i][j].e.velocity = (Vec2){.x=0., .y=0.};
            g_mesh->nodes[i][j].e.energy = 0.;

            g_mesh->nodes[i][j].h.density = NI;
            g_mesh->nodes[i][j].h.velocity = (Vec2){.x=0., .y=0.};
            g_mesh->nodes[i][j].h.energy = 0.;
        }
    }

    for(int i = 1; i <= g_mesh->nx + 1; ++i) {
        g_mesh->edges[0][i].boundary = 0;
        g_mesh->edges[0][i].potential = 0.;
        g_mesh->edges[0][i].n = 0.;
        g_mesh->edges[0][i].p = 0.;

        g_mesh->edges[2][i].boundary = 0;
        g_mesh->edges[2][i].potential = 0.;
        g_mesh->edges[2][i].n = 0.;
        g_mesh->edges[2][i].p = 0.;
    }
    for(int j = 1; j <= g_mesh->nx + 1; ++j) {
        g_mesh->edges[1][j].boundary = 0;
        g_mesh->edges[1][j].potential = 0.;
        g_mesh->edges[1][j].n = 0.;
        g_mesh->edges[1][j].p = 0.;

        g_mesh->edges[3][j].boundary = 0;
        g_mesh->edges[3][j].potential = 0.;
        g_mesh->edges[3][j].n = 0.;
        g_mesh->edges[3][j].p = 0.;
    }


// Thess are the default values
// i.e. if nothing is specified
// in the input file we have
// the following values.
// ============================
 XVAL[ALXINXSB]=XVAL[ALXIN1XSB]=XVAL[INXGA1XAS]=XVAL[INXAL1XAS]=XVAL[INXGAXXAS]=0.5;
// standard doping concentration
 for(int i=1;i<=g_mesh->nx+1;i++)
   for(int j=1;j<=g_mesh->ny+1;j++){
      u2d[i][j][1]=NI;
      h2d[i][j][1]=NI;
      u2d[i][j][2]=u2d[i][j][3]=0.;
      h2d[i][j][2]=h2d[i][j][3]=0.;
   }

// =====================

// Reading the input file
// ======================
 printf("Processing the input file\n\
        =========================\n");
 do{
// read the current row
  fscanf(fp,"%s",s);
// if row is a comment then read it and ignore it
  if(strcmp(s,"#")==0){
    fgets(s,80,fp);
    printf("COMMENT ---> %s",s);
  }
// read and check if the material specified exists
  else if(strcmp(s,"MATERIAL")==0){
    int type;
    real xi,xf;
    real yi,yf;
    fscanf(fp,"%s",s);
    if(strcmp(s,"X")==0) fscanf(fp,"%lf %lf",&xi,&xf);
    else{
      printf("%s : X material not specified!\n",progname);
      exit(0);
    }
    fscanf(fp,"%s",s);
    if(strcmp(s,"Y")==0) fscanf(fp,"%lf %lf",&yi,&yf);
    else{
      printf("%s : Y material not specified!\n",progname);
      exit(0);
    }
    fscanf(fp,"%s",s);
    if(strcmp(s,"SILICON")==0){
        type=SILICON;
    }
    else if(strcmp(s,"GAAS")==0){
        type=GAAS;
    }
    else if(strcmp(s,"GERMANIUM")==0){
        type=GERMANIUM;
    }
    else if(strcmp(s,"INSB")==0) type=INSB;
    else if(strcmp(s,"ALSB")==0) type=ALSB;
    else if(strcmp(s,"ALAS")==0) type=ALAS;
    else if(strcmp(s,"ALP")==0) type=ALP;
    else if(strcmp(s,"GAP")==0) type=GAP;
    else if(strcmp(s,"GASB")==0) type=GASB;
    else if(strcmp(s,"INAS")==0) type=INAS;
    else if(strcmp(s,"INP")==0) type=INP;
    else if(strcmp(s,"GAN")==0) type=GAN;
    else{
      printf("%s : unknown specified material!\n",progname);
      exit(0);
    }
    for(int i=0;i<=g_mesh->nx+4;i++)
      for(int j=0;j<=g_mesh->ny+4;j++){
        if((i-0.5)*g_mesh->dx>=xi && (i-1.5)*g_mesh->dx<=xf
         &&(j-0.5)*g_mesh->dy>=yi && (j-1.5)*g_mesh->dy<=yf){
           g_mesh->nodes[i][j].material = &(g_materials[type]);
        }
      }
    printf("MATERIAL %s X=[%g,%g] Y=[%g,%g] ---> Ok\n",s,xi,xf,yi,yf);
  }

// here we choose the charge particle transport simulated
  else if(strcmp(s,"TRANSPORT")==0){
    transportflag=1;
    fscanf(fp,"%s",s);
    if(strcmp(s,"MC")==0){
      fscanf(fp,"%s",s);
      if(strcmp(s,"ELECTRONS")==0){
        g_config->simulation_model = MCE;
        printf("TRANSPORT MC %s ---> Ok\n",s);
      }
      else if(strcmp(s,"HOLES")==0){
        g_config->simulation_model = MCH;
        printf("TRANSPORT %s ---> Ok\n",s);
        printf("%s: hole transport not yet implemented\n",progname);
        exit(EXIT_FAILURE);
      }
      else if(strcmp(s,"BIPOLAR")==0){
        g_config->simulation_model = MCEH;
        printf("TRANSPORT MC %s ---> Ok\n",s);
      }
      else{
        printf("%s: unknown particles in MC transport specification\n",progname);
        exit(EXIT_FAILURE);
      }
    }
    else if(strcmp(s,"MEP")==0){
      fscanf(fp,"%s",s);
      if(strcmp(s,"ELECTRONS")==0){
        g_config->simulation_model = MEPE;
        printf("TRANSPORT MEP %s ---> Ok\n",s);
      }
      else if(strcmp(s,"HOLES")==0){
        g_config->simulation_model = MEPH;
        printf("TRANSPORT %s ---> Ok\n",s);
      }
      else if(strcmp(s,"BIPOLAR")==0){
        g_config->simulation_model = MEPEH;
        printf("TRANSPORT %s ---> Ok\n",s);
      }
      else{
        printf("%s: unknown particles in MEP transport specification\n",progname);
        exit(EXIT_FAILURE);
      }
    }
    else{
        printf("%s: unknown transport model specified\n",progname);
        exit(EXIT_FAILURE);
    }
  }
// Specify the number of cells in x direction
  else if(strcmp(s,"XSPATIALSTEP")==0){
    fscanf(fp,"%lf",&num);
    g_mesh->nx = (int)num;
    if(g_mesh->nx>NXM){
      printf("%s: too large x-spatial step\n",progname);
      exit(EXIT_FAILURE);
    }
    if(LXflag==0){
      printf("%s: you have to define the x-length first\n",progname);
      exit(EXIT_FAILURE);
    }
    g_mesh->dx = g_mesh->width / g_mesh->nx;
    printf("XSPATIALSTEP = %d ---> Ok\n",g_mesh->nx);
  }
// Specify the number of cells in y direction
  else if(strcmp(s,"YSPATIALSTEP")==0){
    fscanf(fp,"%lf",&num);
    g_mesh->ny = (int)num;
    if(g_mesh->ny>NYM){
      printf("%s: too large y-spatial step\n",progname);
      exit(EXIT_FAILURE);
    }
    if(LYflag==0){
      printf("%s: you have to define the y-length first\n",progname);
      exit(EXIT_FAILURE);
    }
    g_mesh->dy = g_mesh->height / g_mesh->ny;
    printf("YSPATIALSTEP = %d ---> Ok\n",g_mesh->ny);
  }
// specify the final time of simulation
  else if(strcmp(s,"FINALTIME")==0){
    fscanf(fp,"%lg",&num);
    g_config->tf = num;
    if(g_config->tf<=0.){
      printf("%s: not valid final time\n",progname);
    }
    printf("FINAL TIME = %lg ---> Ok\n",g_config->tf);
  }
// setting the energy relaxation time
  else if(strcmp(s,"TAUW")==0){
    fscanf(fp,"%lf",&num);
    g_config->tauw = num;
    if(g_config->tauw==0.){
      printf("%s: not valid energy relaxation time\n",progname);
      exit(EXIT_FAILURE);
    }
    printf("TAUW = %g ---> Ok\n",g_config->tauw);
  }
// configuration of the time step
  else if(strcmp(s,"TIMESTEP")==0){
    fscanf(fp,"%lf",&num);
    g_config->dt = num;
    if(g_config->dt<=0.){
      printf("%s: not valid time step\n",progname);
      exit(EXIT_FAILURE);
    }
    printf("TIME STEP = %g ---> Ok\n",g_config->dt);
  }
// load electron initial data ___ LEID = Load Electron Initial Data
  else if(strcmp(s,"LEID")==0){
    FILE *dp;
    FILE *pp;
    FILE *ep;
    double dum,dum0,dum1;
// Open the initial electron density file
    dp=fopen("density_start.xyz","r");
// File Control, just in case the file does not exist...
     if(dp==NULL){
      printf("%s: fatal error in opening the density000.xyz input file\n",progname);
      exit(EXIT_FAILURE);
     }
// Open the initial electron energy file
    ep=fopen("energy_start.xyz","r");
// File Control, just in case the file does not exist...
     if(ep==NULL){
      printf("%s: fatal error in opening the energy000.xyz input file\n",progname);
      exit(EXIT_FAILURE);
     }
// Open the initial potential file
    pp=fopen("potential_start.xyz","r");
// File Control, just in case the file does not exist...
     if(pp==NULL){
      printf("%s: fatal error in opening the potential000.xyz input file\n",progname);
      exit(EXIT_FAILURE);
     }
// Load the initial data for electrons in case of Monte Carlo method
    if(g_config->simulation_model==MCE || g_config->simulation_model==MCEH){
      for(int j=1;j<=g_mesh->ny+1;j++)
       for(int i=1;i<=g_mesh->nx+1;i++){
        fscanf(dp,"%lf %lf %lf",&dum0,&dum1,&dum);
        u2d[i][j][1]=(real)(dum);
        g_mesh->nodes[i][j].e.density = (real)dum;
        fscanf(ep,"%lf %lf %lf",&dum0,&dum1,&dum);
        u2d[i][j][4]=(real)(dum*Q*u2d[i][j][1]);
        g_mesh->nodes[i][j].e.energy = (real)(dum * Q * g_mesh->nodes[i][j].e.density);
        fscanf(pp,"%lf %lf %lf",&dum0,&dum1,&dum);
        g_mesh->nodes[i][j].potential = (real)dum;
       }
    }
// Load the initial data for electrons in case of Simplified MEP model
    if(g_config->simulation_model==MEPE || g_config->simulation_model==MEPEH){
      for(int j=1;j<=g_mesh->ny+1;j++)
       for(int i=1;i<=g_mesh->nx+1;i++){
        fscanf(dp,"%lf %lf %lf",&dum0,&dum1,&dum);
        u2d[i+2][j+2][1]=(real)(dum);
        g_mesh->nodes[i+2][j+2].e.density = (real)dum;
        fscanf(ep,"%lf %lf %lf",&dum0,&dum1,&dum);
        u2d[i+2][j+2][4]=(real)(dum*Q*u2d[i+2][j+2][1]);
        g_mesh->nodes[i+2][j+2].e.energy = (real)(dum * Q * g_mesh->nodes[i+2][j+2].e.density);
        fscanf(pp,"%lf %lf %lf",&dum0,&dum1,&dum);
        g_mesh->nodes[i][j].potential = (real)dum;
       }
    }
    fclose(dp);
    fclose(ep);
    fclose(pp);
    g_config->load_initial_data = ON;
    printf("Electron initial data loaded ---> Ok\n");
  }
// choice of the GaAs impurity concentration
  else if(strcmp(s,"CIMP")==0){
    fscanf(fp,"%lf",&num);
    g_config->impurity_conc = num;
    if(g_config->impurity_conc<0.){
      printf("%s: not valid impurity concentration\n",progname);
      exit(EXIT_FAILURE);
    }
    printf("CIMP = %g ---> Ok\n",g_config->impurity_conc);
  }
  else if(strcmp(s,"NIMP")==0){
    fscanf(fp,"%lf",&num);
    g_config->neutral_impurity_conc = num;
    if(g_config->neutral_impurity_conc<0.){
      printf("%s: not valid impurity concentration\n",progname);
      exit(EXIT_FAILURE);
    }
    printf("NIMP = %g ---> Ok\n",g_config->neutral_impurity_conc);
  }

// read the x length of the device
  else if(strcmp(s,"XLENGTH")==0){
    fscanf(fp,"%lf",&num);
    g_mesh->width = num;
    if(g_mesh->width<=0.){
      printf("%s: not valid x-length\n",progname);
      exit(EXIT_FAILURE);
    }
    LXflag=1;
    printf("LENGTH X = %g ---> Ok\n",g_mesh->width);
  }
// read the y length of the device
  else if(strcmp(s,"YLENGTH")==0){
    fscanf(fp,"%lf",&num);
    g_mesh->height = num;
    if(g_mesh->height<=0.){
      printf("%s: not valid y-length\n",progname);
      exit(EXIT_FAILURE);
    }
    LYflag=1;
    printf("LENGTH Y = %g ---> Ok\n",g_mesh->height);
  }
// read the donor doping density
  else if(strcmp(s,"DONORDENSITY")==0){
    real xmin,ymin,xmax,ymax,conc;
    if(transportflag==0){
      printf("%s: you have to specify a transport model first\n",progname);
      exit(EXIT_FAILURE);
    }
// read and check the xmin value
    fscanf(fp,"%lf",&num);
    xmin=num;
    if(xmin<0.){
      printf("%s: not valid xmin value\n",progname);
      exit(EXIT_FAILURE);
    }
// read and check the ymin value
    fscanf(fp,"%lf",&num);
    ymin=num;
    if(ymin<0.){
      printf("%s: not valid ymin value\n",progname);
      exit(EXIT_FAILURE);
    }
// read and check the xmax value
    fscanf(fp,"%lf",&num);
    xmax=num;
    if(xmax>g_mesh->width){
      printf("%s: not valid xmax value\n",progname);
      exit(EXIT_FAILURE);
    }
// read and check the ymax value
    fscanf(fp,"%lf",&num);
    ymax=num;
    if(ymax>g_mesh->height){
      printf("%s: not valid xmin value\n",progname);
      exit(EXIT_FAILURE);
    }
// read and check the donor concentration value
    fscanf(fp,"%lf",&num);
    conc=num;
    if(conc<0.){
      printf("%s: not valid donor density\n",progname);
      exit(EXIT_FAILURE);
    }
    if(g_config->simulation_model==MCE || g_config->simulation_model==MCEH)
    for(int i=1;i<=g_mesh->nx+1;i++)
      for(int j=1;j<=g_mesh->ny+1;j++)
        if((i-0.5)*g_mesh->dx>=xmin && (i-1.5)*g_mesh->dx<=xmax
         &&(j-0.5)*g_mesh->dy>=ymin && (j-1.5)*g_mesh->dy<=ymax){
           u2d[i][j][1]=conc;
           g_mesh->nodes[i][j].donor_conc = conc;
           g_mesh->nodes[i][j].e.density = conc;
        }
    if(g_config->simulation_model==MEPE || g_config->simulation_model==MEPEH || g_config->simulation_model==MEPH)
    for(int i=1;i<=g_mesh->nx+1;i++)
      for(int j=1;j<=g_mesh->ny+1;j++){
        if((i-0.5)*g_mesh->dx>=xmin && (i-1.5)*g_mesh->dx<=xmax
         && (j-0.5)*g_mesh->dy>=ymin && (j-1.5)*g_mesh->dy<=ymax){
           u2d[i+2][j+2][1]=conc;
           g_mesh->nodes[i][j].donor_conc = conc;
           g_mesh->nodes[i+2][j+2].e.density = conc;
           u2d[i+2][j+2][4]=conc*1.5*KB*g_config->lattice_temp;
           g_mesh->nodes[i+2][j+2].e.energy = conc * 1.5 * KB * g_config->lattice_temp;
        }
      }
    printf("DONOR DENSITY %g %g %g %g %g ---> Ok\n",
           xmin,ymin,xmax,ymax,conc);
  }
// read the acceptor density in the n zone
// read the donor doping density
  else if(strcmp(s,"ACCEPTORDENSITY")==0){
    real xmin,ymin,xmax,ymax,conc;
// read and check the xmin value
    fscanf(fp,"%lf",&num);
    xmin=num;
    if(xmin<0.){
      printf("%s: not valid xmin value\n",progname);
      exit(EXIT_FAILURE);
    }
// read and check the ymin value
    fscanf(fp,"%lf",&num);
    ymin=num;
    if(ymin<0.){
      printf("%s: not valid ymin value\n",progname);
      exit(EXIT_FAILURE);
    }
// read and check the xmax value
    fscanf(fp,"%lf",&num);
    xmax=num;
    if(xmax>g_mesh->width){
      printf("%s: not valid xmax value\n",progname);
      exit(EXIT_FAILURE);
    }
// read and check the ymax value
    fscanf(fp,"%lf",&num);
    ymax=num;
    if(ymax>g_mesh->height){
      printf("%s: not valid xmin value\n",progname);
      exit(EXIT_FAILURE);
    }
// read and check the acceptor concentration value
    fscanf(fp,"%lf",&num);
    conc=num;
    if(conc<0.){
      printf("%s: not valid acceptor density\n",progname);
      exit(EXIT_FAILURE);
    }
    if(g_config->simulation_model==MEPEH || g_config->simulation_model==MEPH || g_config->simulation_model==MCEH)
    for(int i=1;i<=g_mesh->nx+1;i++)
      for(int j=1;j<=g_mesh->ny+1;j++){
        if((i-0.5)*g_mesh->dx>=xmin && (i-1.5)*g_mesh->dx<=xmax
         && (j-0.5)*g_mesh->dy>=ymin && (j-1.5)*g_mesh->dy<=ymax){
           h2d[i+2][j+2][1]=conc;
           g_mesh->nodes[i][j].acceptor_conc = conc;
           g_mesh->nodes[i+2][j+2].h.density = conc;
           h2d[i+2][j+2][4]=conc*1.5*KB*g_config->lattice_temp;
           g_mesh->nodes[i+2][j+2].h.energy = conc * 1.5 * KB * g_config->lattice_temp;
        }
      }
  }
// read the lattice temperature
  else if(strcmp(s,"LATTICETEMPERATURE")==0){
    fscanf(fp,"%lf",&num);
    g_config->lattice_temp = num;
    if(g_config->lattice_temp<=0.){
      printf("%s: not valid lattice temperature\n",progname);
      exit(EXIT_FAILURE);
    }
    printf("LATTICE TEMPERATURE = %g ---> Ok\n",g_config->lattice_temp);
  }
// read the statistical weight for electrons and holes
  else if(strcmp(s,"STATISTICALWEIGHT")==0){
    fscanf(fp,"%lf",&num);
    g_config->particles_per_cell = (int)num;
    printf("STATISTICAL WEIGHT = %d ---> Ok\n",g_config->particles_per_cell);
  }
// configuration of the output format
  else if(strcmp(s,"OUTPUTFORMAT")==0){
    fscanf(fp,"%s",s);
    if(strcmp(s,"GNUPLOT")==0){
     g_config->output_format = GNUPLOTFORMAT;
     printf("OUTPUT FORMAT = GNUPLOT/XD3D\n");
    }
    else if(strcmp(s,"MESH")==0){
     g_config->output_format = MESHFORMAT;
     printf("OUTPUT FORMAT = MESH/BB\n");
    }
    else{
     printf("%s: unknown output format\n",progname);
     exit(EXIT_FAILURE);
    }
  }
// Definition of an eventual contact
  else if(strcmp(s,"CONTACT")==0){
    char pos[80],kind[80];
    real ipos,fpos,delt = 0.,dens = 0.,denshole = 0.;
    real potential = 0.;
    int i = 0,j,k;
// read and check the "qualitative" position of the contact
    fscanf(fp,"%s",pos);
    if(strcmp(pos,"UP")!=0 && strcmp(pos,"DOWN")!=0
       && strcmp(pos,"LEFT")!=0 && strcmp(pos,"RIGHT")!=0){
      printf("%s: unknown position of contact\n",progname);
      exit(EXIT_FAILURE);
    }
    if(strcmp(pos,"DOWN")==0){
     i=0;
     delt=g_mesh->width/g_mesh->nx;
    }
    if(strcmp(pos,"RIGHT")==0){
     i=1;
     delt=g_mesh->height/g_mesh->ny;
    }
    if(strcmp(pos,"UP")==0){
     i=2;
     delt=g_mesh->width/g_mesh->nx;
    }
    if(strcmp(pos,"LEFT")==0){
     i=3;
     delt=g_mesh->height/g_mesh->ny;
    }
// read the "quantitative" position of the contact
    fscanf(fp,"%lf %lf",&num0,&num);
    ipos=num0;
    fpos=num;
    if(ipos>=fpos){
      printf("%s: not valid position of contact\n",progname);
      exit(EXIT_FAILURE);
    }
// read and check the kind of contact
// It can be : Insulator, Schottky, Ohmic, Vacuum.
    fscanf(fp,"%s",kind);
    if(strcmp(kind,"INSULATOR") != 0 &&
       strcmp(kind,"SCHOTTKY") != 0 &&
       strcmp(kind,"OHMIC") != 0 &&
       strcmp(kind, "VACUUM") != 0) {
        printf("%s: specified physical contact unknown\n", progname);
        exit(EXIT_FAILURE);
    }
    if(strcmp(kind, "INSULATOR") ==0) { k = 0; }
    if(strcmp(kind, "SCHOTTKY") ==0)  { k = 1; }
    if(strcmp(kind, "OHMIC") ==0)     { k = 2; }
    if(strcmp(kind, "VACUUM") ==0)    { k = 3; }
// read the voltage applied to the contact
    fscanf(fp,"%lf",&num);
    potential=num;
// read the electron density for the reservoir at contact
// if and only if the contact is ohmic
    if(k==2){
      fscanf(fp,"%lf",&num);
      dens=num;
    }
// read the hole density for the reservoir at contact
// if and only if the contact is ohmic
    if(k==2 && (g_config->simulation_model==MCH || g_config->simulation_model==MEPH
                   || g_config->simulation_model==MCEH || g_config->simulation_model==MEPEH)){
      fscanf(fp,"%lf",&num);
      denshole=num;
    }
// internal definition of the boundary conditions
// EDGE[i][j][k] = ref
// i = 0 BOTTOM EDGE
// i = 1 RIGHT EDGE
// i = 2 UPPER EDGE
// i = 3 LEFT EDGE
// j is the j-th cell
// k = 0
// ref = 0 INSULATOR
// ref = 1 SCHOTTKY
// ref = 2 OHMIC
// ref = 3 VACUUM
// k = 1
// ref is the applied potential reference
// k = 2
// ref is the density of electron reservoirs at the contact
    ini=(int)(ipos/delt)+1;
    fin=(int)(fpos/delt)+2;
    for(j=ini;j<=fin;j++){
      EDGE[i][j][0]=k;
      g_mesh->edges[i][j].boundary = k;
      if(k==0 || k==1){
        EDGE[i][j][1]=potential;
        EDGE[i][j][2]=0;
        g_mesh->edges[i][j].potential = potential;
        g_mesh->edges[i][j].n = 0.;
        g_mesh->edges[i][j].p = 0.;
        if(k==1) {
          EDGE[i][j][2]=NGATE;
          g_mesh->edges[i][j].n = NGATE;
        }
        if(k==1 && (g_config->simulation_model==MCH || g_config->simulation_model==MEPH
                   || g_config->simulation_model==MCEH || g_config->simulation_model==MEPEH)) {
          EDGE[i][j][3]=NI*NI/NGATE;
          g_mesh->edges[i][j].p = NI * NI / NGATE;
        }
      }
      else if(k==2){
        EDGE[i][j][1]=potential;
        EDGE[i][j][2]=dens;
        g_mesh->edges[i][j].potential = potential;
        g_mesh->edges[i][j].n = dens;
        g_mesh->edges[i][j].p = 0.;
        if(g_config->simulation_model==MCH || g_config->simulation_model==MEPH
          || g_config->simulation_model==MCEH || g_config->simulation_model==MEPEH)
         EDGE[i][j][3]=denshole;
          g_mesh->edges[i][j].p = denshole;
      }
      else if(k == 3) {
        EDGE[i][j][1] = potential;
        EDGE[i][j][2] = 0.;
        EDGE[i][j][3] = 0.;
        g_mesh->edges[i][j].potential = potential;
        g_mesh->edges[i][j].n = 0.;
        g_mesh->edges[i][j].p = 0.;
      }
    }

// Then everything is ok in the contact definition...
    if(k!=2) {
        printf("CONTACT %s %g %g %s %g ---> Ok\n",
            pos,ipos,fpos,kind,potential);
    }
    else if(k==2){
     printf("CONTACT %s %g %g %s %g %g ",
            pos,ipos,fpos,kind,potential,dens);
     if(g_config->simulation_model==MCH || g_config->simulation_model==MEPH
        || g_config->simulation_model==MCEH || g_config->simulation_model==MEPEH) printf("%g ",denshole);
     printf("---> Ok\n");
    }
  }
  else if(strcmp(s,"QEP")==0){
   fscanf(fp,"%s",s);
   if(strcmp(s,"ON")==0) {
    g_config->quantum_flag = ON;
   }
   else if(strcmp(s,"OFF")==0) {
    g_config->quantum_flag = OFF;
   }
   else{
     printf("%s : QEP can be set on ON and OFF only.\n",progname);
     exit(0);
   }
   printf("QUANTUM EFFECTIVE POTENTIAL %s ---> Ok\n",s);
  }
  else if(strcmp(s,"NOQUANTUMEFFECTS")==0){
// WE KEEP IT ONLY FOR BACK COMPATIBILITY
    g_config->quantum_flag = OFF;
    printf("WARNING! This command is DEPRECATED! Use QEP instead!\n");
    printf("QUANTUM EFFECTIVE POTENTIAL = OFF --->Ok\n");
  }
  else if(strcmp(s,"QUANTUMEFFECTS")==0){
// WE KEEP IT ONLY FOR BACK COMPATIBILITY
    g_config->quantum_flag = ON;
    printf("WARNING! This command is DEPRECATED! Use QEP instead!\n");
    printf("QUANTUM EFFECTIVE POTENTIAL = ON --->Ok\n");
  }
  else if(strcmp(s,"MEDIA")==0){
    fscanf(fp,"%lf",&num);
    if(num<0){
      printf("%s: number of media is negative\n",progname);
      exit(EXIT_FAILURE);
    }
    g_config->avg_steps = (int)num;
    g_config->avg_alpha = 1. / (real)g_config->avg_steps;
    printf("MEDIA = %d ---> Ok\n",g_config->avg_steps);
  }
  else if(strcmp(s,"MAXIMINI")==0){
    g_config->max_min_output = ON;
    printf("MAXIMINI ---> Ok\n");
  }
  else if(strcmp(s,"NOMAXIMINI")==0){
    g_config->max_min_output = OFF;
    printf("NO MAXIMINI ---> Ok\n");
  }
  else if(strcmp(s,"SAVEEACHSTEP")==0){
    g_config->save_step_output = ON;
    printf("SAVE AT EACH TIME STEP ---> Ok\n");
  }
  else if(strcmp(s,"SCATTERING_OUTPUT") == 0) {
    g_config->scattering_output = ON;
    printf("OUTPUT SCATTERING RATES ---> Ok\n");
  }
  else if(strcmp(s,"FARADAY")==0){
// Faraday equation ON or OFF
   fscanf(fp,"%s",s);
   if(strcmp(s,"ON")==0) {
    g_config->faraday_flag = ON;
   }
   else if(strcmp(s,"OFF")==0) {
    g_config->faraday_flag = OFF;
   }
   else{
     printf("%s : command FARADAY accept ON or OFF.\n",progname);
     exit(0);
   }
   if(g_config->faraday_flag) printf("FARADAY = ON ---> Ok\n");
   else printf("FARADAY = OFF ---> Ok\n");
  }
  else if(strcmp(s,"CONSTANTMAGNETICFIELD")==0){
// specification of a constant magnetic field in a rectangular area of the device
    real xi,yi,xf,yf,value;
    fscanf(fp,"%lf %lf %lf %lf %lf",&xi,&yi,&xf,&yf,&value);
    if(xi>xf){
      printf("%s : %f > %f in CONSTANTMAGNETICFIELD.\n",progname,xi,xf);
      exit(0);
    }
    if(yi>yf){
      printf("%s : %f > %f in CONSTANTMAGNETICFIELD.\n",progname,yi,yf);
      exit(0);
    }
    for(int i=1;i<=g_mesh->nx+1;i++)
      for(int j=1;j<=g_mesh->ny+1;j++)
        if((i-0.5)*g_mesh->dx>=xi && (i-1.5)*g_mesh->dx<=xf
         &&(j-0.5)*g_mesh->dy>=yi && (j-1.5)*g_mesh->dy<=yf){
           g_mesh->nodes[i][j].magnetic_field = value;
        }
    printf("Constant Magnetic Field %f %f %f %f %f ---> Ok\n",xi,yi,xf,yf,value);
  }
  else if(strcmp(s,"OPTICALSCATTERING")==0){
// specify if the Optical phonons scateering has to be taken into account
   fscanf(fp,"%s",s);
   if(strcmp(s,"ON")==0){
    g_config->optical_phonon_scattering = ON;
    printf("OPTICAL PHONONS SCATTERING = ON ---> Ok\n");
   }
   else if(strcmp(s,"OFF")==0){
    g_config->optical_phonon_scattering = OFF;
    printf("OPTICAL PHONONS SCATTERING = OFF ---> Ok\n");
   }
   else {
    printf("%s : command OPTICALSCATTERING accept ON or OFF.\n",progname);
    exit(0);
   }
  }
  else if(strcmp(s,"ACOUSTICSCATTERING")==0){
// specify if the Acoustic phonons scattering has to be taken into account
   fscanf(fp,"%s",s);
   if(strcmp(s,"ON")==0){
    g_config->acoustic_phonon_scattering = ON;
    printf("ACOUSTIC PHONONS SCATTERING = ON ---> Ok\n");
   }
   else if(strcmp(s,"OFF")==0){
    g_config->acoustic_phonon_scattering = OFF;
    printf("ACOUSTIC PHONONS SCATTERING = OFF ---> Ok\n");
   }
   else {
    printf("%s : command ACOUSTINGSCATTERING accept ON or OFF.\n",progname);
    exit(0);
   }
  }
  else if(strcmp(s,"IMPURITYSCATTERING")==0){
// specify if the Acoustic phonons scattering has to be taken into account
   fscanf(fp,"%s",s);
   if(strcmp(s,"ON")==0){
    g_config->impurity_scattering = ON;
    printf("IMPURITY PHONONS SCATTERING = ON ---> Ok\n");
   }
   else if(strcmp(s,"OFF")==0){
    g_config->impurity_scattering = OFF;
    printf("IMPURITY PHONONS SCATTERING = OFF ---> Ok\n");
   }
   else {
    printf("%s : command IMPURITYSCATTERING accept ON or OFF.\n",progname);
    exit(0);
   }
  }
  else if(strcmp(s,"PIEZOELECTRIC")==0){
// specify if the Acoustic phonons scattering has to be taken into account
   fscanf(fp,"%s",s);
   if(strcmp(s,"ON")==0){
    g_config->piezoelectric_scattering = ON;
    printf("PIEZOELECTRIC SCATTERING = ON ---> Ok\n");
   }
   else if(strcmp(s,"OFF")==0){
    g_config->piezoelectric_scattering = OFF;
    printf("PIEZOELECTRIC SCATTERING = OFF ---> Ok\n");
   }
   else {
    printf("%s : command PIEZOELECTRIC accept ON or OFF.\n",progname);
    exit(0);
   }
  }
  else if(strcmp(s,"NEUTRALIMPURITY")==0){
   fscanf(fp,"%s",s);
   if(strcmp(s,"ON")==0){
    g_config->neutral_impurity_scattering = ON;
    printf("NEUTRAL IMPURITY SCATTERING = ON ---> Ok\n");
   }
   else if(strcmp(s,"OFF")==0){
    g_config->neutral_impurity_scattering = OFF;
    printf("NEUTRAL IMPURITY SCATTERING = OFF ---> Ok\n");
   }
   else {
    printf("%s : command NEUTRALIMPURITY accept ON or OFF.\n",progname);
    exit(0);
   }
  }
  else if(strcmp(s,"ELECTRONHOLE")==0){
   fscanf(fp,"%s",s);
   if(strcmp(s,"ON")==0){
    g_config->electron_hole_scattering = ON;
    printf("ELECTRON HOLE SCATTERING = ON ---> Ok\n");
   }
   else if(strcmp(s,"OFF")==0){
    g_config->electron_hole_scattering = OFF;
    printf("ELECTRON HOLE SCATTERING = OFF ---> Ok\n");
   }
   else {
    printf("%s : command ELECTRONHOLE accept ON or OFF.\n",progname);
    exit(0);
   }
  }
  else if(strcmp(s,"CONDUCTIONBAND")==0){
// Specify the conduction band model
// Possible choices are PARABOLIC, KANE, FULL
   fscanf(fp,"%s",s);
   if(strcmp(s,"PARABOLIC")==0){
    g_config->conduction_band = PARABOLIC;
    printf("CONDUCTION BAND = PARABOLIC ---> Ok\n");
   }
   else if(strcmp(s,"KANE")==0){
    g_config->conduction_band = KANE;
    printf("CONDUCTION BAND = KANE ---> Ok\n");
   }
   else if(strcmp(s,"FULL")==0){
    g_config->conduction_band = FULL;
    printf("CONDUCTION BAND = FULL BAND ---> Ok\n");
   }
   else {
    printf("%s : Unknown specified conduction band model.\n",progname);
    exit(0);
   }
  }
  else if(strcmp(s,"QEP_PARAMETERS")==0){
// Specify the parameters alpha and gamma for the
// quantum effective potential approximation
   real tmp;
   fscanf(fp,"%lf",&tmp);
   g_config->qep_alpha = tmp;
   fscanf(fp,"%lf",&tmp);
   g_config->qep_gamma = tmp;
   printf("QUAT. EFF. POT. PARAMETERS\nALPHA = %f --> Ok\nGAMMA = %f --> Ok\n",g_config->qep_alpha,g_config->qep_gamma);
  }
  else if(strcmp(s,"QEP_MODEL")==0){
// Specify the QEP model to be simulated
   fscanf(fp,"%s",s);
   if(strcmp(s,"CALIBRATED_BOHM")==0){
    g_config->qep_model = QEP_CALIBRATED_BOHM;
    printf("QEP_MODEL = CALIBRATED BOHM POTENTIAL ---> Ok\n");
   }
   else if(strcmp(s,"BOHM")==0){
    g_config->qep_model = QEP_BOHM;
    printf("QEP_MODEL = BOHM POTENTIAL ---> Ok\n");
   }
   else if(strcmp(s,"FULL")==0){
    g_config->qep_model = QEP_FULL;
    printf("QEP_MODEL = FULL EFFECTIVE POTENTIAL ---> Ok\n");
   }
   else if(strcmp(s,"DENSITY_GRADIENT")==0){
    g_config->qep_model = QEP_DENSITY_GRADIENT;
    printf("QEP_MODEL = DENSITY GRADIENT ---> Ok\n");
   } else {
    printf("Unknown specified quant. eff. potential!\n");
    exit(0);
   }
  }
  else if(strcmp(s,"SAVEMESH")==0){
   g_config->save_mesh = ON;
   printf("SAVE THE MESH --> Ok\n");
  }
    else if(strcmp(s, "PHOTON") == 0) {
        fscanf(fp, "%lf", &num);
        g_config->photon_energy = num;
        g_config->photoexcitation_flag = ON;
        printf("PHOTON ENERGY = %g ---> Ok\n", g_config->photon_energy);
    }
    else if(strcmp(s, "TRACKING") == 0) {
        int mod = 0;
        fscanf(fp, "%d", &mod);
        g_config->tracking_output = ON;
        g_config->tracking_mod = mod;
        printf("ELECTRON TRACKING = id %% %d ---> Ok\n", g_config->tracking_mod);
    }
    else if(strcmp(s, "EFIELD") == 0) {
        g_config->constant_efield_flag = ON;
        g_config->poisson_flag = OFF;
        printf("CONSTANT EFIELD = ON ---> Ok\n");
    }
    else if(strcmp(s, "POISSON") == 0) {
        fscanf(fp, "%s", s);
        if(strcmp(s, "ON") == 0) {
            g_config->poisson_flag = ON;
            printf("POISSON CALCULATION = ON ---> Ok\n");
        }
        else if(strcmp(s, "OFF") == 0) {
            g_config->poisson_flag = OFF;
            printf("POISSON CALCULATION = OFF ---> Ok\n");
        }
        else {
            printf("%s: command POISSON accept ON or OFF, given '%s'.\n", progname, s);
            exit(0);
        }
    }
    else if(strcmp(s, "THOMASFERMI") == 0) {
        fscanf(fp, "%s", s);
        if(strcmp(s, "ON") == 0) {
            g_config->thomas_fermi_screening = ON;
            printf("THOMAS FERMI SCREENING = ON ---> Ok\n");
        }
        else if(strcmp(s, "OFF") == 0) {
            g_config->thomas_fermi_screening = OFF;
            printf("THOMAS FERMI SCREENING = OFF ---> Ok\n");
        }
        else {
            printf("%s: command THOMASFERMI accept ON or OFF, given '%s'.\n", progname, s);
            exit(0);
        }
    }
    else if(strcmp(s, "TCAD") == 0) {
        char tcad[1024];
        fgets(tcad, sizeof(tcad), fp);
        char *filename = trim(tcad);
        FILE *input = fopen(filename, "r");
        if(!input) {
            printf("File '%s' does not exist or is unaccessible!\n", filename);
            exit(1);
        }
        printf("USING TCAD FILE '%s' ---> Ok\n", filename);
        int id = 0;
        double x  = 0.,
               Na = 0.,
               Nd = 0.,
                V = 0.,
                n = 0.,
                p = 0.,
                efieldX = 0.,
                efieldY = 0.;
        fgets(tcad, sizeof(tcad), input);
        for(int i = 1; i <= g_mesh->nx + 1; ++i){
            //             id x   Na  Nd  V   n   p   Ex  Ey
            //             1  2   3   4   5   6   7   8   9
            fscanf(input, "%d %lf %lf %lf %lf %lf %lf %lf %lf",
                   &id, &x, &Na, &Nd, &V, &n, &p, &efieldX, &efieldY);
            for(int j = 1; j <= g_mesh->ny + 1; ++j) {
                g_mesh->nodes[i][j].acceptor_conc = Na;
                g_mesh->nodes[i][j].h.density     = p;
                g_mesh->nodes[i][j].donor_conc    = Nd;
                g_mesh->nodes[i][j].e.density     = n;
                g_mesh->nodes[i][j].e.energy      = n * 1.5 * KB * g_config->lattice_temp;
                g_mesh->nodes[i][j].potential     = V;
                g_mesh->nodes[i][j].efield.x      = efieldX;
                g_mesh->nodes[i][j].efield.y      = 0.;
            }
        }
        fclose(input);
        g_config->tcad_data = ON;
    }
    else if(strcmp(s, "SURFACEBB") == 0) {
        double delV = 0.;
        int direction = 0;
        fscanf(fp, "%lf", &delV);
        fscanf(fp, "%s", s);

        if(strcmp(s, "LEFT")        == 0) { direction = direction_t.LEFT; }
        else if(strcmp(s, "RIGHT")  == 0) { direction = direction_t.RIGHT; }
        else if(strcmp(s, "TOP")    == 0) { direction = direction_t.TOP; }
        else if(strcmp(s, "BOTTOM") == 0) { direction = direction_t.BOTTOM; }

        g_config->surface_bb_flag = ON;
        g_config->surface_bb_direction = direction;
        g_config->surface_bb_delV = delV;

        printf("SURFACE BAND BENDING: %s %g eV ---> Ok\n", s, delV);
    }
// elseif(strcmp(s,"")==0){
 }while(!feof(fp));
// computation of the maximum doping density
 g_config->max_doping = 0.;
 for(int i=1;i<=g_mesh->nx+1;i++)
   for(int j=1;j<=g_mesh->ny+1;j++){
     if(g_config->max_doping<=g_mesh->nodes[i][j].donor_conc) {
        g_config->max_doping = g_mesh->nodes[i][j].donor_conc;
    }
   }
 g_config->carriers_per_superparticle = g_config->max_doping * g_mesh->dx * g_mesh->dy / g_config->particles_per_cell;
 printf("Max Doping: %.2e\n"
        "Carriers per Superparticle: %.2e\n",
        g_config->max_doping, g_config->carriers_per_superparticle);
 printf("=========================\n");
}

// =============================================
