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
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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

void
Read_Input_File(void)
{
 char s[180];
 double num,num0;
 int ini,fin;
 int i,j;
 int LXflag=0, LYflag=0;
 int transportflag=0;
 
// Thess are the default values
// i.e. if nothing is specified
// in the input file we have
// the following values.
// ============================
 OPTICALPHONONS=ON;  // Optical phonons scattering ON by default
 ACOUSTICPHONONS=ON; // Acoustic phonons scattering ON by default
 IMPURITYPHONONS=ON; // Impurity phonons scattering ON by default
 CONDUCTION_BAND=KANE;   // Kane approximation by default
 QEP_ALPHA=0.5; // alpha parameter for quantum effective potential
 QEP_GAMMA=1.0; // gamma parameter for quantum effective potential
 QEP_MODEL=QEP_BOHM; // Bohm potential by default (if QEP is ON)
 SAVE_MESH=OFF; // dont save the mesh by default
 XVAL[ALXINXSB]=XVAL[ALXIN1XSB]=XVAL[INXGA1XAS]=XVAL[INXAL1XAS]=XVAL[INXGAXXAS]=0.5;
 Model_Number=MCE; // here we choose the model
 nx=50; // default # of cells in x-direction
 ny=50; // default # of cells in y-direction
 TF=5.0e-12; // Default Final Time in seconds
 DT=0.001e-12; // Default constant time step
 TAUW=0.4e-12; // Default Energy relaxation time
 FARADAYFLAG=0; // no evolution of magnetic field taken into account as default (B = 0 or const.)
// standard doping concentration
 for(i=1;i<=nx+1;i++)
   for(j=1;j<=ny+1;j++){
      u2d[i][j][1]=N_D[i][j]=NI;
      h2d[i][j][1]=N_H[i][j]=NI;
      u2d[i][j][2]=u2d[i][j][3]=0.;
      h2d[i][j][2]=h2d[i][j][3]=0.;
      i_dom[i][j]=SILICON; // Silicon everywhere as default
      B[i][j]=0.; // no magnetic field as default
   }
 CIMP=1.e20; // Impurity concentration
 TL=300.; // Lattice temperature in Kelvin
 NP1=2500; // Number of particle in n+ cell
 MEDIA=500; // Number of step on which we compute the various mean average
 Quantum_Flag=0; // No quantum effects taken into account
 MAXIMINI=0; // No computation of max and min of variables during simulation
 SAVEALWAYS=0; // we don't save at each step by default
 File_Format=GNUPLOTFORMAT; // output file in mesh format
 leid_flag=0;
 SIO2_UP_FLAG=0; // No upper SiO2
 SIO2_DOWN_FLAG=0; // No lower SiO2
// =====================

// Reading the input file
// ======================
 printf("\
Processing the input file\n\
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
    else if(strcmp(s,"ALxINxSB")==0){
      real dum;
      type=ALXINXSB;
      fscanf(fp,"%lf",&dum);
      XVAL[type]=dum;
      if(XVAL[type]<0.){
        printf("%s : x negative for AlxInxSb!\n",progname);
        exit(0);
      }
      XVAL[type]=dum;
      printf("AlxInxSb X = %g\n",XVAL[type]);
    }
    else if(strcmp(s,"ALxIN1-xSB")==0){
      real dum;
      type=ALXIN1XSB;
      fscanf(fp,"%lf",&dum);
      XVAL[type]=dum;
      if(XVAL[type]<0.){
        printf("%s : x negative for AlxIn(1-x)Sb!\n",progname);
        exit(0);
      }
      XVAL[type]=dum;
      printf("AlxIn(1-x)Sb X = %g\n",XVAL[type]);
    }
    else if(strcmp(s,"INxGA1-xAS")==0){
      real dum;
      type=INXGA1XAS;
      fscanf(fp,"%lf",&dum);
      XVAL[type]=dum;
      if(XVAL[type]<0.){
        printf("%s : x negative for InxGa(1-x)As!\n",progname);
        exit(0);
      }
      XVAL[type]=dum;
      printf("InxGa(1-x)As X = %g\n",XVAL[type]);
    }
    else if(strcmp(s,"INxAL1-xAS")==0){
      real dum;
      type=INXAL1XAS;
      fscanf(fp,"%lf",&dum);
      XVAL[type]=dum;
      if(XVAL[type]<0.){
        printf("%s : x negative for InxAl(1-x)As!\n",progname);
        exit(0);
      }
      XVAL[type]=dum;
      printf("InxAl(1-x)As X = %g\n",XVAL[type]);
    }
    else if(strcmp(s,"INxGAxAS2")==0){
      real dum;
      type=INXGAXXAS;
      fscanf(fp,"%lf",&dum);
      XVAL[type]=dum;
      if(XVAL[type]<0.){
        printf("%s : x negative for InxGaxAs2!\n",progname);
        exit(0);
      }
      XVAL[type]=dum;
      printf("InxGaxAs2 X = %g\n",XVAL[type]);
    }
    else{
      printf("%s : unknown specified material!\n",progname);
      exit(0);
    }
    for(i=0;i<=nx+4;i++)
      for(j=0;j<=ny+4;j++){
        if((i-0.5)*dx>=xi && (i-1.5)*dx<=xf
         &&(j-0.5)*dy>=yi && (j-1.5)*dy<=yf){
           i_dom[i][j]=type;
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
        Model_Number=MCE;
        printf("TRANSPORT MC %s ---> Ok\n",s);
      }
      else if(strcmp(s,"HOLES")==0){
        Model_Number=MCH;
        printf("TRANSPORT %s ---> Ok\n",s);
        printf("%s: hole transport not yet implemented\n",progname);
        exit(EXIT_FAILURE);
      }
      else if(strcmp(s,"BIPOLAR")==0){
        Model_Number=MCEH;
        printf("TRANSPORT MC %s ---> Ok\n",s);
      }
      else{
        printf("\
%s: unknown particles in MC transport specification\n",progname);
        exit(EXIT_FAILURE);
      }
    }
    else if(strcmp(s,"MEP")==0){
      fscanf(fp,"%s",s);
      if(strcmp(s,"ELECTRONS")==0){
        Model_Number=MEPE;
        printf("TRANSPORT MEP %s ---> Ok\n",s);
      }
      else if(strcmp(s,"HOLES")==0){
        Model_Number=MEPH;
        printf("TRANSPORT %s ---> Ok\n",s);
      }
      else if(strcmp(s,"BIPOLAR")==0){
        Model_Number=MEPEH;
        printf("TRANSPORT %s ---> Ok\n",s);
      }
      else{
        printf("\
%s: unknown particles in MEP transport specification\n",progname);
        exit(EXIT_FAILURE);
      }
    }
    else{
        printf("\%s: unknown transport model specified\n",progname);
        exit(EXIT_FAILURE);
    }
  }
// Specify the number of cells in x direction
  else if(strcmp(s,"XSPATIALSTEP")==0){
    fscanf(fp,"%lf",&num);
    nx=(int) num;
    if(nx>NXM){
      printf("%s: too large x-spatial step\n",progname);
      exit(EXIT_FAILURE);
    }
    if(LXflag==0){
      printf("%s: you have to define the x-length first\n",progname);
      exit(EXIT_FAILURE);
    }
    dx=LX/nx; // length of cell in x-direction
    printf("XSPATIALSTEP = %d ---> Ok\n",nx);
  }
// Specify the number of cells in y direction
  else if(strcmp(s,"YSPATIALSTEP")==0){
    fscanf(fp,"%lf",&num);
    ny=(int) num;
    if(ny>NYM){
      printf("%s: too large y-spatial step\n",progname);
      exit(EXIT_FAILURE);
    }
    if(LYflag==0){
      printf("%s: you have to define the y-length first\n",progname);
      exit(EXIT_FAILURE);
    }
    dy=LY/ny; // length of cell in y-direction
    printf("YSPATIALSTEP = %d ---> Ok\n",ny);
  }
// specify the final time of simulation
  else if(strcmp(s,"FINALTIME")==0){
    fscanf(fp,"%lg",&num);
    TF=num;
    if(TF<=0.){
      printf("%s: not valid final time\n",progname);
    }
    printf("FINAL TIME = %lg ---> Ok\n",TF);
  }
// Setting the Silicon Oxyde interface
  else if(strcmp(s,"OXYDE")==0){
    fscanf(fp,"%s",s);
    if(strcmp(s,"UP")==0){
      SIO2_UP_FLAG=1;
      fscanf(fp,"%lf",&num);
      SIO2_INI[0]=num;
      if(num<0.){
        printf("%s: not valid upper SiO2 initial position\n",progname);
        exit(EXIT_FAILURE);
      }
      fscanf(fp,"%lf",&num);
      SIO2_FIN[0]=num;
      if(num<=SIO2_INI[0] || num>LX){
        printf("%s: not valid upper SiO2 final position\n",progname);
        exit(EXIT_FAILURE);
      }
      fscanf(fp,"%lf",&num);
      if(num<=0.){
        printf("%s: not valid upper SiO2 thickness\n",progname);
        exit(EXIT_FAILURE);
      }
      SIO2_THICKNESS[0]=num;
      fscanf(fp,"%lf",&num);
      SIO2_POT[0]=num;
      printf("OXYDE ---> UP %g %g %g %g\n",
              SIO2_INI[0],SIO2_FIN[0],SIO2_THICKNESS[0],SIO2_POT[0]);
    }
    else if(strcmp(s,"DOWN")==0){
      SIO2_DOWN_FLAG=1;
      fscanf(fp,"%lf",&num);
      SIO2_INI[1]=num;
      if(num<0.){
        printf("%s: not valid lower SiO2 initial position\n",progname);
        exit(EXIT_FAILURE);
      }
      fscanf(fp,"%lf",&num);
      SIO2_FIN[1]=num;
      if(num<=SIO2_INI[1] || num>LX){
        printf("%s: not valid lower SiO2 final position\n",progname);
        exit(EXIT_FAILURE);
      }
      fscanf(fp,"%lf",&num);
      if(num<=0.){
        printf("%s: not valid lower SiO2 thickness\n",progname);
        exit(EXIT_FAILURE);
      }
      SIO2_THICKNESS[1]=num;
      fscanf(fp,"%lf",&num);
      SIO2_POT[1]=num;
      printf("OXYDE ---> DOWN %g %g %g %g\n",
              SIO2_INI[1],SIO2_FIN[1],SIO2_THICKNESS[1],SIO2_POT[1]);
    }
    else {
      printf("%s: not valid final time\n",progname);
      exit(EXIT_FAILURE);
    }
    printf("FINAL TIME = %g ---> Ok\n",TF);
  }
// setting the energy relaxation time
  else if(strcmp(s,"TAUW")==0){
    fscanf(fp,"%lf",&num);
    TAUW=num;
    if(TAUW==0.){
      printf("%s: not valid energy relaxation time\n",progname);
      exit(EXIT_FAILURE);
    }
    printf("TAUW = %g ---> Ok\n",TAUW);
  }
// configuration of the time step
  else if(strcmp(s,"TIMESTEP")==0){
    fscanf(fp,"%lf",&num);
    DT=num;
    if(DT<=0.){
      printf("%s: not valid time step\n",progname);
      exit(EXIT_FAILURE);
    }
    printf("TIME STEP = %g ---> Ok\n",DT);
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
      printf("\
%s: fatal error in opening the density000.xyz input file\n",progname);
      exit(EXIT_FAILURE);
     }
// Open the initial electron energy file
    ep=fopen("energy_start.xyz","r");
// File Control, just in case the file does not exist...
     if(ep==NULL){
      printf("\
%s: fatal error in opening the energy000.xyz input file\n",progname);
      exit(EXIT_FAILURE);
     }
// Open the initial potential file
    pp=fopen("potential_start.xyz","r");
// File Control, just in case the file does not exist...
     if(pp==NULL){
      printf("\
%s: fatal error in opening the potential000.xyz input file\n",progname);
      exit(EXIT_FAILURE);
     }
// Load the initial data for electrons in case of Monte Carlo method
    if(Model_Number==MCE || Model_Number==MCEH){
      for(j=1;j<=ny+1;j++)
       for(i=1;i<=nx+1;i++){
        fscanf(dp,"%lf %lf %lf",&dum0,&dum1,&dum);
        u2d[i][j][1]=(real)(dum);
        fscanf(ep,"%lf %lf %lf",&dum0,&dum1,&dum);
        u2d[i][j][4]=(real)(dum*Q*u2d[i][j][1]);
        fscanf(pp,"%lf %lf %lf",&dum0,&dum1,&dum);
        PSI[i][j]=(real)(dum);
       }
    }
// Load the initial data for electrons in case of Simplified MEP model
    if(Model_Number==MEPE || Model_Number==MEPEH){
      for(j=1;j<=ny+1;j++)
       for(i=1;i<=nx+1;i++){
        fscanf(dp,"%lf %lf %lf",&dum0,&dum1,&dum);
        u2d[i+2][j+2][1]=(real)(dum);
        fscanf(ep,"%lf %lf %lf",&dum0,&dum1,&dum);
        u2d[i+2][j+2][4]=(real)(dum*Q*u2d[i+2][j+2][1]);
        fscanf(pp,"%lf %lf %lf",&dum0,&dum1,&dum);
        PSI[i][j]=(real)(dum);
       }
    }
    fclose(dp);
    fclose(ep);
    fclose(pp);
    leid_flag=1;
    printf("Electron initial data loaded ---> Ok\n");
  }
// choice of the GaAs impurity concentration
  else if(strcmp(s,"CIMP")==0){
    fscanf(fp,"%lf",&num);
    CIMP=num;
    if(CIMP<0.){
      printf("%s: not valid impurity concentration\n",progname);
      exit(EXIT_FAILURE);
    }
    printf("CIMP = %g ---> Ok\n",CIMP);
  }

// read the x length of the device
  else if(strcmp(s,"XLENGTH")==0){
    fscanf(fp,"%lf",&num);
    LX=num;
    if(LX<=0.){
      printf("%s: not valid x-length\n",progname);
      exit(EXIT_FAILURE);
    }
    LXflag=1;
    printf("LENGTH X = %g ---> Ok\n",LX);
  }
// read the y length of the device
  else if(strcmp(s,"YLENGTH")==0){
    fscanf(fp,"%lf",&num);
    LY=num;
    if(LY<=0.){
      printf("%s: not valid y-length\n",progname);
      exit(EXIT_FAILURE);
    }
    LYflag=1;
    printf("LENGTH Y = %g ---> Ok\n",LY);
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
    if(xmax>LX){
      printf("%s: not valid xmax value\n",progname);
      exit(EXIT_FAILURE);
    }
// read and check the ymax value
    fscanf(fp,"%lf",&num);
    ymax=num;
    if(ymax>LY){
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
    if(Model_Number==MCE || Model_Number==MCEH)
    for(i=1;i<=nx+1;i++)
      for(j=1;j<=ny+1;j++)
        if((i-0.5)*dx>=xmin && (i-1.5)*dx<=xmax
         &&(j-0.5)*dy>=ymin && (j-1.5)*dy<=ymax){
           u2d[i][j][1]=N_D[i][j]=conc;
        }
    if(Model_Number==MEPE || Model_Number==MEPEH || Model_Number==MEPH)
    for(i=1;i<=nx+1;i++)
      for(j=1;j<=ny+1;j++){
        if((i-0.5)*dx>=xmin && (i-1.5)*dx<=xmax
         && (j-0.5)*dy>=ymin && (j-1.5)*dy<=ymax){
           u2d[i+2][j+2][1]=N_D[i][j]=conc;
           u2d[i+2][j+2][4]=conc*1.5*KB*TL;
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
    if(xmax>LX){
      printf("%s: not valid xmax value\n",progname);
      exit(EXIT_FAILURE);
    }
// read and check the ymax value
    fscanf(fp,"%lf",&num);
    ymax=num;
    if(ymax>LY){
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
    if(Model_Number==MEPEH || Model_Number==MEPH || Model_Number==MCEH)
    for(i=1;i<=nx+1;i++)
      for(j=1;j<=ny+1;j++){
        if((i-0.5)*dx>=xmin && (i-1.5)*dx<=xmax
         && (j-0.5)*dy>=ymin && (j-1.5)*dy<=ymax){
           h2d[i+2][j+2][1]=N_H[i][j]=conc;
           h2d[i+2][j+2][4]=conc*1.5*KB*TL;
        }
      }
  }
// read the lattice temperature
  else if(strcmp(s,"LATTICETEMPERATURE")==0){
    fscanf(fp,"%lf",&num);
    TL=num;
    if(TL<=0.){
      printf("%s: not valid lattice temperature\n",progname);
      exit(EXIT_FAILURE);
    }
    printf("LATTICE TEMPERATURE = %g ---> Ok\n",TL);
  }
// read the statistical weight for electrons and holes
  else if(strcmp(s,"STATISTICALWEIGHT")==0){
    fscanf(fp,"%lf",&num);
    NP1=(int) num;
    printf("STATISTICAL WEIGHT = %d ---> Ok\n",NP1);
  }
// configuration of the output format
  else if(strcmp(s,"OUTPUTFORMAT")==0){
    fscanf(fp,"%s",s);
    if(strcmp(s,"GNUPLOT")==0){
     File_Format=GNUPLOTFORMAT;
     printf("OUTPUT FORMAT = GNUPLOT/XD3D\n");
    }
    else if(strcmp(s,"MESH")==0){
     File_Format=MESHFORMAT;
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
    real ipos,fpos,delt,dens,denshole;
    real potential;
    int i,j,k;
// read and check the "qualitative" position of the contact
    fscanf(fp,"%s",pos);
    if(strcmp(pos,"UP")!=0 && strcmp(pos,"DOWN")!=0
       && strcmp(pos,"LEFT")!=0 && strcmp(pos,"RIGHT")!=0){
      printf("%s: unknown position of contact\n",progname);
      exit(EXIT_FAILURE);
    }
    if(strcmp(pos,"DOWN")==0){
     i=0;
     delt=LX/nx;
    }
    if(strcmp(pos,"RIGHT")==0){
     i=1;
     delt=LY/ny;
    }
    if(strcmp(pos,"UP")==0){
     i=2;
     delt=LX/nx;
    }
    if(strcmp(pos,"LEFT")==0){
     i=3;
     delt=LY/ny;
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
// It can be : Insulator, Schottky or Ohmic.
    fscanf(fp,"%s",kind);
    if(strcmp(kind,"INSULATOR")!=0 && strcmp(kind,"SCHOTTKY")!=0
       && strcmp(kind,"OHMIC")!=0){
      printf("%s: specified physical contact unknown\n",progname);
      exit(EXIT_FAILURE);
    }
    if(strcmp(kind,"INSULATOR")==0) k=0;
    if(strcmp(kind,"SCHOTTKY")==0) k=1;
    if(strcmp(kind,"OHMIC")==0) k=2;
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
    if(k==2 && (Model_Number==MCH || Model_Number==MEPH
                   || Model_Number==MCEH || Model_Number==MEPEH)){
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
// k = 1
// ref is the applied potential reference
// k = 2
// ref is the density of electron reservoirs at the contact
// k = 3
// ref is the density of hole reservoirs at the contact
    ini=(int)(ipos/delt)+1;
    fin=(int)(fpos/delt)+2;
    for(j=ini;j<=fin;j++){
      EDGE[i][j][0]=k;
      if(k==0 || k==1){
        EDGE[i][j][1]=potential;
        EDGE[i][j][2]=0;
        if(k==1) EDGE[i][j][2]=NGATE;
        if(k==1 && (Model_Number==MCH || Model_Number==MEPH
                   || Model_Number==MCEH || Model_Number==MEPEH))
          EDGE[i][j][3]=NI*NI/NGATE;
      }
      else if(k==2){
        EDGE[i][j][1]=potential;
        EDGE[i][j][2]=dens;
        if(Model_Number==MCH || Model_Number==MEPH
          || Model_Number==MCEH || Model_Number==MEPEH)
         EDGE[i][j][3]=denshole;
      }
    }

// Then everything is ok in the contact definition...
    if(k!=2)printf("CONTACT %s %g %g %s %g ---> Ok\n",
            pos,ipos,fpos,kind,potential);
    else if(k==2){
     printf("CONTACT %s %g %g %s %g %g ",
            pos,ipos,fpos,kind,potential,dens);
     if(Model_Number==MCH || Model_Number==MEPH
        || Model_Number==MCEH || Model_Number==MEPEH) printf("%g ",denshole);
     printf("---> Ok\n");
    }
  }
  else if(strcmp(s,"QEP")==0){
   fscanf(fp,"%s",s);
   if(strcmp(s,"ON")==0) Quantum_Flag=1;
   else if(strcmp(s,"OFF")==0) Quantum_Flag=0;
   else{
     printf("%s : QEP can be set on ON and OFF only.\n",progname);
     exit(0);
   }
   printf("QUANTUM EFFECTIVE POTENTIAL %s ---> Ok\n",s);
  }
  else if(strcmp(s,"NOQUANTUMEFFECTS")==0){
// WE KEEP IT ONLY FOR BACK COMPATIBILITY
    Quantum_Flag=0;
    printf("WARNING! This command is DEPRECATED! Use QEP instead!\n");
    printf("QUANTUM EFFECTIVE POTENTIAL = OFF --->Ok\n");
  }
  else if(strcmp(s,"QUANTUMEFFECTS")==0){
// WE KEEP IT ONLY FOR BACK COMPATIBILITY
    Quantum_Flag=1;
    printf("WARNING! This command is DEPRECATED! Use QEP instead!\n");
    printf("QUANTUM EFFECTIVE POTENTIAL = ON --->Ok\n");
  }
  else if(strcmp(s,"MEDIA")==0){
    fscanf(fp,"%lf",&num);
    if(num<0){
      printf("%s: number of media is negative\n",progname);
      exit(EXIT_FAILURE);
    }
    MEDIA=(int) num;
    printf("MEDIA = %d ---> Ok\n",MEDIA);
  }
  else if(strcmp(s,"MAXIMINI")==0){
    MAXIMINI=1;
    printf("MAXIMINI ---> Ok\n");
  }
  else if(strcmp(s,"NOMAXIMINI")==0){
    MAXIMINI=0;
    printf("NO MAXIMINI ---> Ok\n");
  }
  else if(strcmp(s,"SAVEEACHSTEP")==0){
    SAVEALWAYS=1;
    printf("SAVE AT EACH TIME STEP ---> Ok\n");
  }
  else if(strcmp(s,"FARADAY")==0){
// Faraday equation ON or OFF
   fscanf(fp,"%s",s);
   if(strcmp(s,"ON")==0) FARADAYFLAG=1;
   else if(strcmp(s,"OFF")==0) FARADAYFLAG=0;
   else{
     printf("%s : command FARADAY accept ON or OFF.\n",progname);
     exit(0);
   }
   if(FARADAYFLAG) printf("FARADAY = ON ---> Ok\n");
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
    for(i=1;i<=nx+1;i++)
      for(j=1;j<=ny+1;j++)
        if((i-0.5)*dx>=xi && (i-1.5)*dx<=xf
         &&(j-0.5)*dy>=yi && (j-1.5)*dy<=yf){
           B[i][j]=value;
        }
    printf("Constant Magnetic Field %f %f %f %f %f ---> Ok\n",xi,yi,xf,yf,value);
  }
  else if(strcmp(s,"OPTICALSCATTERING")==0){
// specify if the Optical phonons scateering has to be taken into account
   fscanf(fp,"%s",s);
   if(strcmp(s,"ON")==0){
    OPTICALPHONONS=ON;
    printf("OPTICAL PHONONS SCATTERING = ON ---> Ok\n");
   }
   else if(strcmp(s,"OFF")==0){
    OPTICALPHONONS=OFF;
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
    ACOUSTICPHONONS=ON;
    printf("ACOUSTIC PHONONS SCATTERING = ON ---> Ok\n");
   }
   else if(strcmp(s,"OFF")==0){
    ACOUSTICPHONONS=OFF;
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
    IMPURITYPHONONS=ON;
    printf("IMPURITY PHONONS SCATTERING = ON ---> Ok\n");
   }
   else if(strcmp(s,"OFF")==0){
    IMPURITYPHONONS=OFF;
    printf("IMPURITY PHONONS SCATTERING = OFF ---> Ok\n");
   }
   else {
    printf("%s : command IMPURITYSCATTERING accept ON or OFF.\n",progname);
    exit(0);
   }
  }
  else if(strcmp(s,"CONDUCTIONBAND")==0){
// Specify the conduction band model
// Possible choices are PARABOLIC, KANE, FULL
   fscanf(fp,"%s",s);
   if(strcmp(s,"PARABOLIC")==0){
    CONDUCTION_BAND=PARABOLIC;
    printf("CONDUCTION BAND = PARABOLIC ---> Ok\n");
   }
   else if(strcmp(s,"KANE")==0){
    CONDUCTION_BAND=KANE;
    printf("CONDUCTION BAND = KANE ---> Ok\n");
   }
   else if(strcmp(s,"FULL")==0){
    CONDUCTION_BAND=FULL;
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
   QEP_ALPHA=tmp;
   fscanf(fp,"%lf",&tmp);
   QEP_GAMMA=tmp;
   printf("QUAT. EFF. POT. PARAMETERS\nALPHA = %f --> Ok\nGAMMA = %f --> Ok\n",QEP_ALPHA,QEP_GAMMA);
  }
  else if(strcmp(s,"QEP_MODEL")==0){
// Specify the QEP model to be simulated
   fscanf(fp,"%s",s);
   if(strcmp(s,"CALIBRATED_BOHM")==0){
    QEP_MODEL=QEP_CALIBRATED_BOHM;
    printf("QEP_MODEL = CALIBRATED BOHM POTENTIAL ---> Ok\n");
   }
   else if(strcmp(s,"BOHM")==0){
    QEP_MODEL=QEP_BOHM;
    printf("QEP_MODEL = BOHM POTENTIAL ---> Ok\n");
   }
   else if(strcmp(s,"FULL")==0){
    QEP_MODEL=QEP_FULL;
    printf("QEP_MODEL = FULL EFFECTIVE POTENTIAL ---> Ok\n");
   }
   else if(strcmp(s,"DENSITY_GRADIENT")==0){
    QEP_MODEL=QEP_DENSITY_GRADIENT;
    printf("QEP_MODEL = DENSITY GRADIENT ---> Ok\n");
   } else {
    printf("Unknown specified quant. eff. potential!\n");
    exit(0);
   }
  }
  else if(strcmp(s,"SAVEMESH")==0){
   SAVE_MESH=ON;
   printf("SAVE THE MESH --> Ok\n");
  }
// elseif(strcmp(s,"")==0){
 }while(!feof(fp));
// computation of the maximum doping density
 DDmax=0.;
 for(i=1;i<=nx+1;i++)
   for(j=1;j<=ny+1;j++){
     if(DDmax<=N_D[i][j]) DDmax=N_D[i][j];
   }
 printf("=========================\n");
}

// =============================================
