/* archimedes.c -- This belongs to GNU archimedes.

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


// =================================================================
// File Name : archimedes.c
// Version   : release 2.0.0
// Date of Creation : 07 Dec.2003, Paris, France, Jean Michel Sellier
// Last Revision : 15 Sept. 2011, Carry le Rouet, France, Jean Michel Sellier.
// System of Measure : M.K.S.C. System
// =================================================================

#include<getopt.h>
#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
#include<math.h>
#include<memory.h>
#include<time.h>
#ifdef	HAVE_STRING_H
#include<string.h>
#else
#include<strings.h>
#endif

// #include "rappture.h"

// Preprocessor Definitions
#define real double
#define ON 1
#define OFF 0
#define KANE 0
#define PARABOLIC 1
#define FULL 2
#define QEP_BOHM 0
#define QEP_CALIBRATED_BOHM 1
#define QEP_FULL 2
#define QEP_DENSITY_GRADIENT 3
#define MN3 4
#define NXM 308
#define NYM 308
#define DIME 1003
#define ITMAX 10000000
#define POISSONITMAX 1500
#define SMALL 1.e-5
#define VMAX 1000000
#define NPMAX 10000000 // maximum number of super-particles
#define MCE 0  // MCE stands for MC for electrons only
#define MCH 1  // MCH stands for MC for holes only
#define MCEH 2 // MCEH stands for MC both for electrons and holes
#define MEPE 3 // MEPE stands for MEP model for electrons only
#define MEPH 4 // MEPH stands for MEP model for holes only
#define MEPEH 5 // MEPEH stands for MEP model for electrons and holes
#define GNUPLOTFORMAT 0 // output file in GNUPLOT format
#define MESHFORMAT 1 // output file in Mesh format

// definition of the material reference table
#define NOAMTIA 16   // Number Of All Material Taken Into Account (excluding SiO2)
// ********************
#define iSIO2 -1     // reference number to Silicon Oxide
#define SILICON 0    // reference number to Silicon
#define GAAS 1       // reference number to GaAs
#define GERMANIUM 2  // reference number to Germanium
#define INSB 3       // reference number to InSb
#define ALSB 4       // reference number to AlSb
#define ALXINXSB 5   // reference number to Al_x In_x Sb
#define ALXIN1XSB 6  // reference number to Al_x In_(1-x) Sb
#define ALAS 7       // reference number to AlAs
#define ALP 8        // reference number to AlP
#define GAP 9        // reference number to GaP
#define GASB 10      // reference number to GaSb
#define INAS 11      // reference number to InAs
#define INP 12       // reference number to InP
#define INXGA1XAS 13 // reference number to In_x Ga_(1-x) As
#define INXAL1XAS 14 // reference number to In_x Al_(1-x) As
#define INXGAXXAS 15 // reference number to In_x Ga_(1-x) As (second zone)
// ********************

#define NUMSIO2 2 // maximum number of SiO2 interfaces

// ===============================
#include "constants.h"
#include "extrema.h"
#include "sign.h"
#include "mm.h"
#include "mm2.h"
// ===============================

// All integers here...
int NUM_VERT,NUM_EXAHEDRA,MEDIA,MAXIMINI;
int SAVEALWAYS;
int nx,ny;
int ISEED,NP1,INUM,IV;
int c,Model_Number,File_Format,Quantum_Flag;
int NG,NE;
int leid_flag;
int SIO2_UP_FLAG;
int SIO2_DOWN_FLAG;
int FARADAYFLAG;
int i_dom[NXM+1][NYM+1];
int NOVALLEY[NOAMTIA+1];
int ACOUSTICPHONONS;
int OPTICALPHONONS;
int IMPURITYPHONONS;
int CONDUCTION_BAND;
int SAVE_MESH;
int NODE_GEO[3][(NXM+1)*(NYM+1)];

// All "real"'s here...
real u2d[NXM+1][NYM+1][MN3+1];
real h2d[NXM+1][NYM+1][MN3+1];
real PSI[NXM+1][NYM+1],E[NXM+1][NYM+1][2];
real N_D[NXM+1][NYM+1],N_H[NXM+1][NYM+1];
real dx,dy;
real TEMPO=0.,TF;
real LX,LY;
real TL,DT;
real mstar2;
real BKTQ,QH;
real SMH[NOAMTIA+1][3];
real HHM[NOAMTIA+1][3];
real HM[NOAMTIA+1][3];
real GM[NOAMTIA+1];
real SWK[NOAMTIA+1][3][14][DIME+1];
real P[NPMAX+1][7];
real KX,KY,KZ,X,Y;
real TS,EPP,DDmax;
real EDGE[4][NXM+NYM+1][4];
real CIMP,QD2,TAUW;
real EPSRSIO2;
real bufx2d[NXM+1][NYM+1];
real bufy2d[NXM+1][NYM+1];
real ux2d[NXM+1][NYM+1][MN3+1];
real uy2d[NXM+1][NYM+1][MN3+1];
real f2d[NXM+1][NYM+1][MN3+1];
real g2d[NXM+1][NYM+1][MN3+1];
real fx2d[NXM+1][NYM+1][MN3+1];
real gy2d[NXM+1][NYM+1][MN3+1];
real c11[7],c12[7],c21[7],c22[7];
real u[7],f[7],g[7],cw[7];
real SIO2_INI[NUMSIO2],SIO2_FIN[NUMSIO2];
real SIO2_POT[NUMSIO2];
real SIO2_THICKNESS[NUMSIO2];
real SIO2[NUMSIO2][NXM+1][NYM+1];
real B[NXM+1][NYM+1];
real EPSR[NOAMTIA+1];
real EPF[NOAMTIA+1];
real MSTAR[NOAMTIA+1][6];
real alphaK[NOAMTIA+1][4];
real EG[NOAMTIA+1];
real HWO[NOAMTIA+1][6];
real DTK[NOAMTIA+1][6];
real ZF[NOAMTIA+1][6];
real RHO[NOAMTIA+1];
real DA[NOAMTIA+1];
real UL[NOAMTIA+1];
real EMIN[NOAMTIA+1][4];
real XVAL[NOAMTIA+1];
real LATTCONST[NOAMTIA+1];
real CB_FULL[NOAMTIA+1][11];
real QEP_ALPHA;
real QEP_GAMMA;
real QEP_MODEL;
real COORD[2][(NXM+1)*(NYM+1)];

// All structures here...
time_t binarytime;
struct tm *nowtm;
struct option longopts[] =
{
  { "version", no_argument, NULL, 'v' },
  { "help", no_argument, NULL, 'h' }
};
// All files here...
FILE *fp;
// All strings here...
static char *progname;

#include "mesher.h"
#include "poissonbcs.h"
#include "faradaybcs.h"
#include "media.h"
#include "saveoutput2dmeshformat.h"
#include "saveoutput2dgnuplot.h"
#include "saveoutput2dholegnuplot.h"
#include "saveoutput2dholemeshformat.h"
#include "saveoutputfiles.h"
#include "quantumeffectivepotential.h"
#include "electric_field.h"
#include "faraday.h"
#include "random.h"
#include "mcparameters.h"
#include "deviceconfig.h"
#include "particlecreation.h"
#include "drift.h"
#include "scattering.h"
#include "ensemblemontecarlo.h"
#include "charge.h"
#include "computecurrents.h"
#include "MEP_interpolation.h"
#include "electron_relaxation.h"
#include "HMEPbcs.h"
#include "ParabMEP2D.h"
#include "Hole_bcs.h"
#include "holemep2d.h"
#include "hole_relaxation.h"
#include "updating.h"
#include "readinputfile.h"
//#include "SaveRappture.h"

int 
main(int argc,char* argv[])
{

//  RpLibrary* lib=NULL;

  int optc;
  int h = 0, v = 0, lose = 0, z = 0;

//  lib=rpLibrary(argv[1]);

//  if(lib!=NULL) printf("Rappture Library loaded correctly.\n");
//  else{
//   printf("Unable to load the Rappture library.\n");
//   return(0);
//  }

  progname = argv[0];

  while ((optc = getopt_long (argc, argv, "hv", longopts, (int *) 0))
         != EOF)
    switch (optc)
      {
      case 'v':
        v = 1;
        break;
      case 'h':
        h = 1;
        break;
      default:
        lose = 1;
        break;
      }
  
  if (optind == argc - 1)
    z = 1;
  else if (lose || optind < argc)
    {
      /* Print error message and exit.  */
      if (optind < argc)
        printf("Too many arguments\n");
        printf("Try `%s --help' for more information.\n",progname);
      exit(1);
    }

  /* `help' should come first.  If `help' is requested, ignore the other
     options. */
  if (h)
    {
      /* Print help info and exit.  */
      /* TRANSLATORS: --help output 1
         no-wrap */
      printf("\
GNU archimedes, a simulator for submicron and nanoscale semiconductor devices.\nCopyright (C) 2004-2011 Jean Michel D. Sellier.\n");
      printf ("\n");
      /* TRANSLATORS: --help output 2
         no-wrap */
      printf ("\
Usage: %s [OPTION] file...\n",progname);

      printf ("\n");
      /* TRANSLATORS: --help output 3 : options 1/2
         no-wrap */
      printf("\
  -h, --help          display this help and exit\n\
  -v, --version       display version information and exit\n");

      printf ("\n");
      /* TRANSLATORS: --help output 5 (end)
         TRANSLATORS, please don't forget to add the contact address for
         your translation!
         no-wrap */
      printf ("\
Report bugs to jeanmichel.sellier@gmail.com or jsellier@purdue.edu\n");

      exit (0);
    }

  if (v)
    {
      /* Print version number.  */
      printf("archimedes - GNU archimedes 2.0.0\n");
      /* xgettext: no-wrap */
      printf("");
      printf("\
Copyright (C) %s Jean Michel D. Sellier.\n\n\
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A\n\
PARTICULAR PURPOSE.\n\
You may redistribute copies of GNU %s under the terms\n\
of the GNU General Public License.\n\
For more information about these matters, see the file named COPYING.\n",
              "2004 - 2011","Archimedes");
      exit (0);
    }
  else if (z){
// In case of filename specified
// =============================
     fp=fopen(argv[1],"r"); // here we open th input file...
// File Control, just in case the file does not exist...
     if(fp==NULL){
      printf("%s: fatal error in opening the input file %s\n",
             progname,argv[1]);
      exit(EXIT_FAILURE);
     }
// ========================
// Material constants here!
// ========================
// Silicon electrons in the X-valley
// Germanium electrons in the Gamma-valley
// GaAs electrons in the Gamma- and L-valley
// For the others, see comments below
     NOVALLEY[SILICON]=1;    // X-valley
     NOVALLEY[GERMANIUM]=1;  // G-valley
     NOVALLEY[GAAS]=2;       // G and L-valleys
     NOVALLEY[INSB]=1;       // G valley
     NOVALLEY[ALSB]=1;       // G-valley
     NOVALLEY[ALXINXSB]=1;   // G-valley
     NOVALLEY[ALXIN1XSB]=1;  // G-valley
     NOVALLEY[ALAS]=1;       // G-valley
     NOVALLEY[ALP]=1;        // G-valley
     NOVALLEY[GAP]=1;        // G-valley
     NOVALLEY[GASB]=1;       // G-valley
     NOVALLEY[INAS]=1;       // G-valley
     NOVALLEY[INP]=1;        // G-valley
     NOVALLEY[INXGA1XAS]=1;  // only G-valley
     NOVALLEY[INXAL1XAS]=1;  // G-valley
     NOVALLEY[INXGAXXAS]=1;  // only G-valley
// Dielectric constant for Silicon Oxide SiO2
     EPSRSIO2=3.9*EPS0;         // see http://en.wikipedia.org/wiki/Relative_permittivity
// Dielectric constant for Semiconducting materials
// STATIC
// ======
     EPSR[SILICON]=11.68;       // see http://en.wikipedia.org/wiki/Relative_permittivity
     EPSR[GERMANIUM]=16.2;      // see http://www.ioffe.ru/SVA/NSM/Semicond/Ge/basic.html
     EPSR[GAAS]=12.90;          // see http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/basic.html
     EPSR[INSB]=16.8;           // see http://www.ioffe.ru/SVA/NSM/Semicond/InSb/basic.html
     EPSR[ALSB]=12.04;          // Fischetti conversations
     EPSR[ALAS]=12.90-2.84;     // see http://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/basic.html (x=1.0)
     EPSR[ALP]=9.80;            // Fischetti conversations
     EPSR[GAP]=11.10;           // see http://www.ioffe.ru/SVA/NSM/Semicond/GaP/basic.html
     EPSR[GASB]=15.69;          // see http://www.ioffe.ru/SVA/NSM/Semicond/GaSb/basic.html
     EPSR[INAS]=15.15;          // see http://www.ioffe.ru/SVA/NSM/Semicond/InAs/basic.html
     EPSR[INP]=12.50;           // see http://www.ioffe.ru/SVA/NSM/Semicond/InP/basic.html 
// III-V semiconductor compounds high frequency dieletric constant
// HIGH FREQUENCY
// ==============
     EPF[GAAS]=10.89;           // see http://www.ioffe.ru/SVA/NSM/Semicond/GaAs/basic.html
     EPF[INSB]=15.68;           // Fischetti conversations
     EPF[ALSB]=9.88;            // Fiscehtti conversations
     EPF[ALAS]=10.89-2.73;      // see http://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/basic.html (x=1.0)
     EPF[ALP]=7.54;             // Fischetti conversations
     EPF[GAP]=9.11;             // see http://www.ioffe.ru/SVA/NSM/Semicond/GaP/basic.html
     EPF[GASB]=14.44;           // see http://www.ioffe.ru/SVA/NSM/Semicond/GaSb/basic.html
     EPF[INAS]=12.3;            // see http://www.ioffe.ru/SVA/NSM/Semicond/InAs/basic.html
     EPF[INP]=9.61;             // see http://www.ioffe.ru/SVA/NSM/Semicond/InP/basic.html
     int ii;
     for(ii=0;ii<NOAMTIA;ii++){
       int i;
       for(i=0;i<6;i++){
         HWO[ii][i]=0.;
         DTK[ii][i]=0.;
         ZF[ii][i]=0.;
       }
     }
// Optical phonon scattering energy (eV)
     HWO[SILICON][0]=0.0120;    // Sellier, Tomizawa
     HWO[SILICON][1]=0.0185;    // Sellier, Tomizawa
     HWO[SILICON][2]=0.0190;    // Sellier, Tomizawa
     HWO[SILICON][3]=0.0474;    // Sellier, Tomizawa
     HWO[SILICON][4]=0.0612;    // Sellier, Tomizawa
     HWO[SILICON][5]=0.0590;    // Sellier, Tomizawa
     HWO[GERMANIUM][0]=0.03704; // Fischetti
     HWO[GAAS][0]=0.03536;      // Fischetti
     HWO[INSB][0]=0.02404;      // Fischetti
     HWO[ALSB][0]=0.0360;       // Fischetti
     HWO[ALAS][0]=0.05009;      // Fischetti
     HWO[ALP][0]=0.06211;       // Fischetti
     HWO[GAP][0]=0.04523;       // Fischetti
     HWO[GASB][0]=0.02529;      // Fischetti
     HWO[INAS][0]=0.03008;      // Fischetti
     HWO[INP][0]=0.04240;       // Fischetti
// Optical coupling constants (eV/m)
     DTK[SILICON][0]=0.05e11;     // Jacoboni Reggiani
     DTK[SILICON][1]=0.08e11;     // Jacoboni Reggiani
     DTK[SILICON][2]=0.03e11;     // Jacoboni Reggiani
     DTK[SILICON][3]=0.20e11;     // Jacoboni Reggiani
     DTK[SILICON][4]=1.14e11;     // Jacoboni Reggiani
     DTK[SILICON][5]=0.20e11;     // Jacoboni Reggiani
     DTK[GERMANIUM][0]=0.0;       // Jacoboni Reggiani
     DTK[GERMANIUM][1]=0.079e11;  // Jacoboni Reggiani
     DTK[GERMANIUM][2]=0.0;       // Jacoboni Reggiani
     DTK[GERMANIUM][3]=0.0;       // Jacoboni Reggiani
     DTK[GERMANIUM][4]=0.95e11;   // Jacoboni Reggiani
     DTK[GERMANIUM][5]=0.0;       // Jacoboni Reggiani
     DTK[GAAS][0]=1.11e11;        // Sellier, Tomizawa
     DTK[INSB][0]=0.47e11; // see ???
     DTK[ALSB][0]=0.55e11; // see ???
     DTK[ALAS][0]=3.0e11;  // see ???
     DTK[ALP][0]=0.95e11;  // see ???
     DTK[GAP][0]=5.33e11;  // see ???
     DTK[GASB][0]=0.94e11; // see ???
     DTK[INAS][0]=3.59e11; // see ???
     DTK[INP][0]=2.46e11;  // see ???
// Optical phonon Z-factor
     ZF[SILICON][0]=1.;   // Sellier
     ZF[SILICON][1]=1.;   // Sellier
     ZF[SILICON][2]=4.;   // Sellier
     ZF[SILICON][3]=4.;   // Sellier
     ZF[SILICON][4]=1.;   // Sellier
     ZF[SILICON][5]=4.;   // Sellier
     ZF[GERMANIUM][0]=1.; // see ???
     ZF[GAAS][0]=1.;      // Sellier
     ZF[INSB][0]=1.;      // see ???
     ZF[ALSB][0]=1.;      // see ???
     ZF[ALAS][0]=1.;      // see ???
     ZF[ALP][0]=1.;       // see ???
     ZF[GAP][0]=1.;       // see ???
     ZF[GASB][0]=1.;      // see ???
     ZF[INAS][0]=1.;      // see ???
     ZF[INP][0]=1.;       // see ???
// Crystal Density (Kg/m^3)
     RHO[SILICON]=2.33e3;   // Fischetti conversations
     RHO[GERMANIUM]=5.32e3; // Fischetti conversations
     RHO[GAAS]=5.36e3;      // Fischetti conversations
     RHO[INSB]=5.78e3;      // Fischetti conversations
     RHO[ALSB]=4.26e3;      // Fischetti conversations
     RHO[ALAS]=3.76e3;      // Fischetti conversations
     RHO[ALP]=2.40e3;       // Fischetti conversations
     RHO[GAP]=4.14e3;       // Fischetti conversations
     RHO[GASB]=5.61e3;      // Fischetti conversations
     RHO[INAS]=5.67e3;      // Fischetti conversations
     RHO[INP]=4.81e3;       // Fischetti conversations
// Acoustic deformation potential (Joule)
     DA[SILICON]=9.*Q;      // Fischetti -- Jacoboni Reggiani
     DA[GERMANIUM]=9.*Q;    // Fischetti -- Jacoboni Reggiani
     DA[GAAS]=7.*Q;         // Fischetti - Gamma valley
     DA[INSB]=7.*Q;         // Fischetti
     DA[ALSB]=4.6*Q;        // Fischetti
     DA[ALAS]=9.3*Q;        // Fischetti
     DA[ALP]=9.3*Q;         // Fischetti
     DA[GAP]=7.4*Q;         // Fischetti
     DA[GASB]=9.*Q;         // Fischetti
     DA[INAS]=8.2*Q;        // Fischetti
     DA[INP]=6.2*Q;         // Fischetti
// Longitudinal sound velocity (m/sec)
     UL[SILICON]=9.18e3;   // Fischetti
     UL[GERMANIUM]=5.4e3;  // Fischetti
     UL[GAAS]=5.24e3;      // Fischetti
     UL[INSB]=3.41e3;      // Fischetti
     UL[ALSB]=4.25e3;      // Fischetti
     UL[ALAS]=5.65e3;      // Fischetti
     UL[ALP]=7.41e3;       // Fischetti
     UL[GAP]=5.85e3;       // Fischetti
     UL[GASB]=3.97e3;      // Fischetti
     UL[INAS]=4.28e3;      // Fischetti
     UL[INP]=5.13e3;       // Fischetti
// Band minimum energy
// first valley
     EMIN[SILICON][1]=0.0;     // Sellier, Fischetti, etc.
     EMIN[GERMANIUM][1]=0.173; 
     EMIN[GAAS][1]=0.0;        // Tomizawa
     EMIN[INSB][1]=0.0;
     EMIN[ALSB][1]=0.507;
     EMIN[ALAS][1]=0.767;
     EMIN[ALP][1]=1.237;
     EMIN[GAP][1]=0.496;
     EMIN[GASB][1]=0.0;
     EMIN[INAS][1]=0.0;
     EMIN[INP][1]=0.0;
// eventual second valley
     EMIN[GAAS][2]=0.323;

// Definition of effective mass for all materials in all valleys
     MSTAR[SILICON][1]=0.32;     // see Sellier, Tomizawa, etc.
     MSTAR[GAAS][1]=0.067;       // Gamma-valley -- see Tomizawa
     MSTAR[GAAS][2]=0.350;       // L-valley     -- see Tomizawa
     MSTAR[GERMANIUM][1]=0.12;   // Gamma valley -- see http://ecee.colorado.edu/~bart/book/effmass.htm#long
     MSTAR[INSB][1]=0.0135;      // Gamma-valley -- see Ram-Mohan
     MSTAR[ALSB][1]=0.14;        // Gamma-valley -- See Ram-Mohan
     MSTAR[ALAS][1]=0.149;       // Gamma-valley -- see Ram-Mohan J.App.Phys. Vol.89, Num.11
     MSTAR[ALP][1]=0.22;         // Gamma-valley -- see Ram-Mohan
     MSTAR[GAP][1]=0.13;         // Gamma-valley -- see Ram-Mohan J.App.Phys. Vol.89, Num.11
     MSTAR[GASB][1]=0.039;       // Gamma-valley -- see Ram-Mohan
     MSTAR[INAS][1]=0.026;       // Gamma-valley -- see Ram-Mohan J.App.Phys. Vol.89, Num.11
     MSTAR[INP][1]=0.0795;       // Gamma-valley -- see Ram-Mohan
// non-parabolicity coefficients
     alphaK[SILICON][1]=0.5;    // see Sellier, Tomizawa
     alphaK[GERMANIUM][1]=0.3;  // Gamma valley - Jacoboni Reggiani
// Lattice constants
     LATTCONST[GAAS]=565.35e-12;       // CODATA
     LATTCONST[SILICON]=543.102e-12;   // CODATA
     LATTCONST[GERMANIUM]=564.613e-12; // CODATA
     LATTCONST[ALP]=545.10e-12;        // CODATA
     LATTCONST[ALAS]=565.05e-12;       // CODATA
     LATTCONST[ALSB]=613.55e-12;       // CODATA
     LATTCONST[GAP]=545.12e-12;        // CODATA
     LATTCONST[GASB]=609.59e-12;       // CODATA
     LATTCONST[INP]=586.87e-12;        // CODATA
     LATTCONST[INAS]=605.83e-12;       // CODATA
     LATTCONST[INSB]=647.9e-12;        // CODATA
// sp3s* Conduction Band calculated parameters for full band simulations
#include "Silicon.h"
#include "GaAs.h"
#include "Germanium.h"
#include "AlP.h"
#include "AlAs.h"
#include "AlSb.h"
#include "GaP.h"
#include "GaSb.h"
#include "InP.h"
#include "InAs.h"
#include "InSb.h"
// ======================================================
// ======================================================

// We reset the some arrays
// ========================
     memset(&u2d,0,sizeof(u2d));
     memset(&E,0,sizeof(E));
     memset(&EDGE,0,sizeof(EDGE));
     memset(&SIO2,0,sizeof(SIO2));

// Read the geometrical and physical description of the MESFET
// ===========================================================
     Read_Input_File();
// ===========================================================

// Construction of the mesh for the electrostatic potential
// (to properly take into account the oxyde layers)
     mesher();

// III-V Semiconductor materials energy gap (depending on the lattice temperature)
     printf("\n");
     EG[SILICON]=1.21-3.333e-4*TL;
     printf("EG[SILICON]      = %g\n",EG[SILICON]);
     EG[GERMANIUM]=0.747-3.587e-4*TL;
     printf("EG[GERMANIUM]    = %g\n",EG[GERMANIUM]);
     EG[GAAS]=1.54-4.036e-4*TL;
     printf("EG[GAAS]         = %g\n",EG[GAAS]);
     EG[INSB]=0.2446-2.153e-4*TL;
     printf("EG[INSB]         = %g\n",EG[INSB]);
     EG[ALSB]=1.696-2.20e-4*TL;
     printf("EG[ALSB]         = %g\n",EG[ALSB]);
     EG[ALAS]=2.314-3.0e-4*TL;
     printf("EG[ALAS]         = %g\n",EG[ALAS]);
     EG[ALP]=2.51-3.333e-4*TL;
     printf("EG[ALP]         = %g\n",EG[ALP]);
     EG[GAP]=2.35-2.667e-4*TL;
     printf("EG[GAP]         = %g\n",EG[GAP]);
     EG[GASB]=0.81-3.667e-4*TL;
     printf("EG[GASB]         = %g\n",EG[GASB]);
     EG[INAS]=0.434-2.601e-4*TL;
     printf("EG[INAS]         = %g\n",EG[INAS]);
     EG[INP]=1.445-3.296e-4*TL;
     printf("EG[INP]         = %g\n",EG[INP]);
     printf("\n");

     if(CONDUCTION_BAND==KANE || CONDUCTION_BAND==PARABOLIC || CONDUCTION_BAND==FULL){
      // USED WHATEVER IS THE CONDUCTION BAND FOR THE INITIAL PSUEDO WAVE VECTOR
      // OF THE PSEUDO PARTICLES
// all the following non-parabolicity coefficients depend on lattice temperature
// non-parabolicity coefficient for GaAs in the GAMMA-valley
      alphaK[GAAS][1]=pow(1.-MSTAR[GAAS][1],2.)/(EG[GAAS]+EMIN[GAAS][1]);//expected value = 0.611
      printf("alphaK_gamma[GaAs] = %g\n",alphaK[GAAS][1]);
// non-parabolicity coefficient for GaAs in the L-valley
      alphaK[GAAS][2]=pow(1.-MSTAR[GAAS][2],2.)/(EG[GAAS]+EMIN[GAAS][2]);//expected value = 0.242;
      printf("alphaK_L[GaAs]     = %g\n",alphaK[GAAS][2]);
// non-parabolicity coefficient for InSb in the GAMMA-valley
      alphaK[INSB][1]=pow(1.-MSTAR[INSB][1],2.)/(EG[INSB]+EMIN[INSB][1]);//5.59;
      printf("alphaK_gamma[InSb] = %g\n",alphaK[INSB][1]);
// non-parabolicity coefficient for AlSb in the GAMMA-valley
      alphaK[ALSB][1]=pow(1.-MSTAR[ALSB][1],2.)/(EG[ALSB]+EMIN[ALSB][1]);//0.321;
      printf("alphaK_gamma[AlSb] = %g\n",alphaK[ALSB][1]);
// non-parabolicity coefficient for AlAs in the GAMMA-valley
      alphaK[ALAS][1]=pow(1.-MSTAR[ALAS][1],2.)/(EG[ALAS]+EMIN[ALAS][1]);
      printf("alphaK_gamma[AlAs] = %g\n",alphaK[ALAS][1]);
// non-parabolicity coefficient for AlP in the GAMMA-valley
      alphaK[ALP][1]=pow(1.-MSTAR[ALP][1],2.)/(EG[ALP]+EMIN[ALP][1]);
      printf("alphaK_gamma[AlP] = %g\n",alphaK[ALP][1]);
// non-parabolicity coefficient for GaP in the GAMMA-valley
      alphaK[GAP][1]=pow(1.-MSTAR[GAP][1],2.)/(EG[GAP]+EMIN[GAP][1]);
      printf("alphaK_gamma[GaP] = %g\n",alphaK[GAP][1]);
// non-parabolicity coefficient for GaSb in the GAMMA-valley
      alphaK[GASB][1]=pow(1.-MSTAR[GASB][1],2.)/(EG[GASB]+EMIN[GASB][1]);
      printf("alphaK_gamma[GaSb] = %g\n",alphaK[GASB][1]);
// non-parabolicity coefficient for InAs in the GAMMA-valley
      alphaK[INAS][1]=pow(1.-MSTAR[INAS][1],2.)/(EG[INAS]+EMIN[INAS][1]);
      printf("alphaK_gamma[InAs] = %g\n",alphaK[INAS][1]);
// non-parabolicity coefficient for InP in the GAMMA-valley
      alphaK[INP][1]=pow(1.-MSTAR[INP][1],2.)/(EG[INP]+EMIN[INP][1]);
      printf("alphaK_gamma[InP] = %g\n",alphaK[INP][1]);
     }

// Semiconductor compounds
// ***
// Relative dielectric constant for semiconductor compounds
     EPSR[ALXINXSB]=XVAL[ALXINXSB]*EPSR[ALSB]+XVAL[ALXINXSB]*EPSR[INSB];
     EPSR[ALXIN1XSB]=XVAL[ALXIN1XSB]*EPSR[ALSB]+(1.-XVAL[ALXIN1XSB])*EPSR[INSB];
     EPSR[INXGA1XAS]=XVAL[INXGA1XAS]*EPSR[INAS]+(1.-XVAL[INXGA1XAS])*EPSR[GAAS];
     EPSR[INXAL1XAS]=XVAL[INXAL1XAS]*EPSR[INAS]+(1.-XVAL[INXAL1XAS])*EPSR[ALAS];
     EPSR[INXGAXXAS]=XVAL[INXGAXXAS]*EPSR[INAS]+XVAL[INXGAXXAS]*EPSR[GAAS];
// semiconductor compounds high frequency dieletric constant
     EPF[ALXINXSB]=XVAL[ALXINXSB]*(EPF[ALSB]+EPF[INSB]);
     EPF[ALXIN1XSB]=XVAL[ALXIN1XSB]*EPF[ALSB]+(1.-XVAL[ALXIN1XSB])*EPF[INSB];
     EPF[INXGA1XAS]=XVAL[INXGA1XAS]*EPF[INAS]+(1.-XVAL[INXGA1XAS])*EPF[GAAS];
     EPF[INXAL1XAS]=XVAL[INXAL1XAS]*EPF[INAS]+(1.-XVAL[INXAL1XAS])*EPF[ALAS];
     EPF[INXGAXXAS]=XVAL[INXGAXXAS]*EPF[INAS]+XVAL[INXGAXXAS]*EPF[GAAS];
// semiconductor compounds optical phonon scattering energy (eV)
     HWO[ALXINXSB][0]=XVAL[ALXINXSB]*(HWO[ALSB][0]+HWO[INSB][0]);
     HWO[ALXIN1XSB][0]=XVAL[ALXIN1XSB]*HWO[ALSB][0]+(1.-XVAL[ALXIN1XSB])*HWO[INSB][0];
     HWO[INXGA1XAS][0]=XVAL[INXGA1XAS]*HWO[INAS][0]+(1.-XVAL[INXGA1XAS])*HWO[GAAS][0];
     HWO[INXAL1XAS][0]=XVAL[INXAL1XAS]*HWO[INAS][0]+(1.-XVAL[INXAL1XAS])*HWO[ALAS][0];
     HWO[INXGAXXAS][0]=XVAL[INXGAXXAS]*HWO[INAS][0]+XVAL[INXGAXXAS]*HWO[GAAS][0];
// semiconductor compounds optical coupling constants (eV/m)
     DTK[ALXINXSB][0]=XVAL[ALXINXSB]*(DTK[ALSB][0]+DTK[INSB][0]);
     DTK[ALXIN1XSB][0]=XVAL[ALXIN1XSB]*DTK[ALSB][0]+(1.-XVAL[ALXIN1XSB])*DTK[INSB][0];
     DTK[INXGA1XAS][0]=XVAL[INXGA1XAS]*DTK[INAS][0]+(1.-XVAL[INXGA1XAS])*DTK[GAAS][0];
     DTK[INXAL1XAS][0]=XVAL[INXAL1XAS]*DTK[INAS][0]+(1.-XVAL[INXAL1XAS])*DTK[ALAS][0];
     DTK[INXGAXXAS][0]=XVAL[INXGAXXAS]*DTK[INAS][0]+XVAL[INXGAXXAS]*DTK[GAAS][0];
// semiconductor compounds optical phonon Z-factor
     ZF[ALXINXSB][0]=XVAL[ALXINXSB]*(ZF[ALSB][0]+ZF[INSB][0]);
     ZF[ALXIN1XSB][0]=XVAL[ALXIN1XSB]*ZF[ALSB][0]+(1.-XVAL[ALXIN1XSB])*ZF[INSB][0];
     ZF[INXGA1XAS][0]=XVAL[INXGA1XAS]*ZF[INAS][0]+(1.-XVAL[INXGA1XAS])*ZF[GAAS][0];
     ZF[INXAL1XAS][0]=XVAL[INXAL1XAS]*ZF[INAS][0]+(1.-XVAL[INXAL1XAS])*ZF[ALAS][0];
     ZF[INXGAXXAS][0]=XVAL[INXGAXXAS]*ZF[INAS][0]+XVAL[INXGAXXAS]*ZF[GAAS][0];
// semiconductor compounds Crystal Density (Kg/m^3)
     RHO[ALXINXSB]=XVAL[ALXINXSB]*(RHO[ALSB]+RHO[INSB]);
     RHO[ALXIN1XSB]=XVAL[ALXIN1XSB]*RHO[ALSB]+(1.-XVAL[ALXIN1XSB])*RHO[INSB];
     RHO[INXGA1XAS]=XVAL[INXGA1XAS]*RHO[INAS]+(1.-XVAL[INXGA1XAS])*RHO[GAAS];
     RHO[INXAL1XAS]=XVAL[INXAL1XAS]*RHO[INAS]+(1.-XVAL[INXAL1XAS])*RHO[ALAS];
     RHO[INXGAXXAS]=XVAL[INXGAXXAS]*RHO[INAS]+XVAL[INXGAXXAS]*RHO[GAAS];
// semiconductor compounds acoustic deformation potential (Joule)
     DA[ALXINXSB]=XVAL[ALXINXSB]*(DA[ALSB]+DA[INSB]);
     DA[ALXIN1XSB]=XVAL[ALXIN1XSB]*DA[ALSB]+(1.-XVAL[ALXIN1XSB])*DA[INSB];
     DA[INXGA1XAS]=XVAL[INXGA1XAS]*DA[INAS]+(1.-XVAL[INXGA1XAS])*DA[GAAS];
     DA[INXAL1XAS]=XVAL[INXAL1XAS]*DA[INAS]+(1.-XVAL[INXAL1XAS])*DA[ALAS];
     DA[INXGAXXAS]=XVAL[INXGAXXAS]*DA[INAS]+XVAL[INXGAXXAS]*DA[GAAS];
// semiconductor compounds longitudinal sound velocity (m/sec)
     UL[ALXINXSB]=XVAL[ALXINXSB]*(UL[ALSB]+UL[INSB]);
     UL[ALXIN1XSB]=XVAL[ALXIN1XSB]*UL[ALSB]+(1.-XVAL[ALXIN1XSB])*UL[INSB];
     UL[INXGA1XAS]=XVAL[INXGA1XAS]*UL[INAS]+(1.-XVAL[INXGA1XAS])*UL[GAAS];
     UL[INXAL1XAS]=XVAL[INXAL1XAS]*UL[INAS]+(1.-XVAL[INXAL1XAS])*UL[ALAS];
     UL[INXGAXXAS]=XVAL[INXGAXXAS]*UL[INAS]+XVAL[INXGAXXAS]*UL[GAAS];
// semiconductor compounds energy gap
     EG[ALXINXSB]=XVAL[ALXINXSB]*(EG[ALSB]+EG[INSB]);
     EG[ALXIN1XSB]=XVAL[ALXIN1XSB]*EG[ALSB]+(1.-XVAL[ALXIN1XSB])*EG[INSB];
     EG[INXGA1XAS]=XVAL[INXGA1XAS]*EG[INAS]+(1.-XVAL[INXGA1XAS])*EG[GAAS];
     EG[INXAL1XAS]=XVAL[INXAL1XAS]*EG[INAS]+(1.-XVAL[INXAL1XAS])*EG[ALAS];
     EG[INXGAXXAS]=XVAL[INXGAXXAS]*EG[INAS]+XVAL[INXGAXXAS]*EG[GAAS];
// semiconductor compounds energy minimum of GAMMA-valley
     EMIN[ALXINXSB][1]=XVAL[ALXINXSB]*(EMIN[ALSB][1]+EMIN[INSB][1]);
     EMIN[ALXIN1XSB][1]=XVAL[ALXIN1XSB]*EMIN[ALSB][1]+(1.-XVAL[ALXIN1XSB])*EMIN[INSB][1];
     EMIN[INXGA1XAS][1]=XVAL[INXGA1XAS]*EMIN[INAS][1]+(1.-XVAL[INXGA1XAS])*EMIN[GAAS][1];
     EMIN[INXAL1XAS][1]=XVAL[INXAL1XAS]*EMIN[INAS][1]+(1.-XVAL[INXAL1XAS])*EMIN[ALAS][1];
     EMIN[INXGAXXAS][1]=XVAL[INXGAXXAS]*EMIN[INAS][1]+XVAL[INXGAXXAS]*EMIN[GAAS][1];
// semiconductor compounds energy minimum 0f L-valley
     EMIN[ALXINXSB][2]=XVAL[ALXINXSB]*(EMIN[ALSB][2]+EMIN[INSB][2]);
     EMIN[ALXIN1XSB][2]=XVAL[ALXIN1XSB]*EMIN[ALSB][2]+(1.-XVAL[ALXIN1XSB])*EMIN[INSB][2];
     EMIN[INXGA1XAS][2]=XVAL[INXGA1XAS]*EMIN[INAS][2]+(1.-XVAL[INXGA1XAS])*EMIN[GAAS][2];
     EMIN[INXAL1XAS][2]=XVAL[INXAL1XAS]*EMIN[INAS][2]+(1.-XVAL[INXAL1XAS])*EMIN[ALAS][2];
     EMIN[INXGAXXAS][2]=XVAL[INXGAXXAS]*EMIN[INAS][2]+XVAL[INXGAXXAS]*EMIN[GAAS][2];
// GAMMA-valley effective mass
     MSTAR[ALXINXSB][1]=XVAL[ALXINXSB]*(MSTAR[ALSB][1]+MSTAR[INSB][1]);
     MSTAR[ALXIN1XSB][1]=XVAL[ALXIN1XSB]*MSTAR[ALSB][1]+(1.-XVAL[ALXIN1XSB])*MSTAR[INSB][1];
     MSTAR[INXGA1XAS][1]=XVAL[INXGA1XAS]*MSTAR[INAS][1]+(1.-XVAL[INXGA1XAS])*MSTAR[GAAS][1];
     MSTAR[INXAL1XAS][1]=XVAL[INXAL1XAS]*MSTAR[INAS][1]+(1.-XVAL[INXAL1XAS])*MSTAR[ALAS][1];
     MSTAR[INXGAXXAS][1]=XVAL[INXGAXXAS]*MSTAR[INAS][1]+XVAL[INXGAXXAS]*MSTAR[GAAS][1];
// L-valley effective mass
     MSTAR[ALXINXSB][2]=XVAL[ALXINXSB]*(MSTAR[ALSB][2]+MSTAR[INSB][2]);
     MSTAR[ALXIN1XSB][2]=XVAL[ALXIN1XSB]*MSTAR[ALSB][2]+(1.-XVAL[ALXIN1XSB])*MSTAR[INSB][2];
     MSTAR[INXGA1XAS][2]=XVAL[INXGA1XAS]*MSTAR[INAS][2]+(1.-XVAL[INXGA1XAS])*MSTAR[GAAS][2];
     MSTAR[INXAL1XAS][2]=XVAL[INXAL1XAS]*MSTAR[INAS][2]+(1.-XVAL[INXAL1XAS])*MSTAR[ALAS][2];
     MSTAR[INXGAXXAS][2]=XVAL[INXGAXXAS]*MSTAR[INAS][2]+XVAL[INXGAXXAS]*MSTAR[GAAS][2];
// non-parabolicity coefficient for semiconductor compounds in the GAMMA-valley
     alphaK[ALXINXSB][1]=XVAL[ALXINXSB]*(alphaK[ALSB][1]+alphaK[INSB][1]);
     alphaK[ALXIN1XSB][1]=XVAL[ALXIN1XSB]*alphaK[ALSB][1]+(1.-XVAL[ALXIN1XSB])*alphaK[INSB][1];
     alphaK[INXGA1XAS][1]=XVAL[INXGA1XAS]*alphaK[INAS][1]+(1.-XVAL[INXGA1XAS])*alphaK[GAAS][1];
     alphaK[INXAL1XAS][1]=XVAL[INXAL1XAS]*alphaK[INAS][1]+(1.-XVAL[INXAL1XAS])*alphaK[ALAS][1];
     alphaK[INXGAXXAS][1]=XVAL[INXGAXXAS]*alphaK[INAS][1]+XVAL[INXGAXXAS]*alphaK[GAAS][1];
// non-parabolicity coefficient for Al_x In_(1-x) Sb in the L-valley
     alphaK[ALXINXSB][2]=XVAL[ALXINXSB]*(alphaK[ALSB][2]+alphaK[INSB][2]);
     alphaK[ALXIN1XSB][2]=XVAL[ALXIN1XSB]*alphaK[ALSB][2]+(1.-XVAL[ALXIN1XSB])*alphaK[INSB][2];
     alphaK[INXGA1XAS][2]=XVAL[INXGA1XAS]*alphaK[INAS][2]+(1.-XVAL[INXGA1XAS])*alphaK[GAAS][2];
     alphaK[INXAL1XAS][2]=XVAL[INXAL1XAS]*alphaK[INAS][2]+(1.-XVAL[INXAL1XAS])*alphaK[ALAS][2];
     alphaK[INXGAXXAS][2]=XVAL[INXGAXXAS]*alphaK[INAS][2]+XVAL[INXGAXXAS]*alphaK[GAAS][2];
// ***

// Closure of the input file
// =========================
     fclose(fp);

// Read all the coefficients for MEP simulation
// ============================================
     MEP_coefficients();

// Loading of the initial runtime
// ==============================
     binarytime=time(NULL);
     nowtm=localtime(&binarytime);

     printf("\n\nComputation Started at %s\n",asctime(nowtm));

// Boundary conditions for the model simulated
// ===========================================
     PoissonBCs();
     if(FARADAYFLAG) FaradayBCs();

// Initialisation for Monte Carlo
// ==============================
     if(Model_Number==MCE || Model_Number==MCEH){
      int i;
      for(i=0;i<NOAMTIA;i++) MCparameters(i);
      printf("\n");
      MCdevice_config();
     }
     printf("\n");
// HERE IS THE SIMULATION
// ======================
     for(c=1;c<=ITMAX;c++) updating(Model_Number);
// Here we save the outputs
// ========================
//     SaveRappture();
     SaveOutputFiles(File_Format,0);
     printf("Final Output has been saved\n");
    }
  else{
// No filename has been specified
// ==============================
   printf("%s: no input file\n",progname);
   exit(0);
// ########################################################
// ## NEW: if archimedes is called without a file input
// ##      it means that it is running under Rappture GUI
// ########################################################

//   int err=0;
//   const char* retstr=NULL;

//   rpGetString(lib,"input.number(materialx).current",&retstr);

//   rpResult(lib);
//   rpFreeLibrary(&lib);
//   exit(0);   
  }
    /* Print greeting message and exit. */
   binarytime=time(NULL);
   nowtm=localtime(&binarytime);
   printf("Computation Finished at %s\n",asctime(nowtm));
//   rpFreeLibrary(&lib);
   return(EXIT_SUCCESS); // Successfull exit
}

/* archimedes.c ends here */

// ***********************************************************************
// ************************************************************************
