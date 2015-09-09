/* deviceconfig.h -- This file is part of Archimedes release 0.1.2.
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements both the Monte Carlo method
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It also includes the quantum effects by means 
   of effective potential method. It is now able to simulate applied
   magnetic fields along with self consistent Faraday equation.

   Copyright (C) 2004, 2005, 2006, 2007 Jean Michel Sellier <sellier@dmi.unict.it>
 
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


// ######################################################
// Created on 05 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 13 Sep.2007, Siracusa, J.M.Sellier
// ######################################################

// Initial Configuration of the Particles (electrons)
// simulated in the devices.

void
MCdevice_config(void)
{
  long int n=0; 
  int i,j,np,m;
  real c1,c2,c3,c4,c5,c6,c7;

// Number of carriers per particle
 EPP=DDmax*dx*dy/NP1;

 for(i=1;i<=nx+1;i++)
   for(j=1;j<=ny+1;j++){
// np=number of particles in the (i,j)-th cell
     np=(int)(N_D[i][j]*dx*dy/EPP+0.5);
     if(leid_flag==1) np=(int)(u2d[i][j][1]*dx*dy/EPP+0.5);
     if((i==1) || (i==nx+1)) np/=2;
     if((j==1) || (j==ny+1)) np/=2;
     if(np>0){
      for(m=1;m<=np;m++){
        n++;
        if(n>NPMAX){
         printf("%s: too big number of particles\n",progname);
         exit(EXIT_FAILURE);
        }
// We assume that the particles are initially 
// at near thermal equilibrium

// In the case of two-valleys materials, 80% of the electrons are considered in the Gamma
// valley in the starting simulation time, while the other 20% are
// in the L-valley.
      if(leid_flag==0){
       IV=1;
       c1=log(rnd());
       if(NOVALLEY[i_dom[i][j]]==1) 
        c2=SMH[i_dom[i][j]][0]*sqrt(-1.5*BKTQ*c1*(1.-alphaK[i_dom[i][j]][1]*1.5*BKTQ*c1));
       if(NOVALLEY[i_dom[i][j]]==2){
// 80% of the created particles goes in the first valley
        IV=1;
        c2=SMH[i_dom[i][j]][IV]*sqrt(-1.5*BKTQ*c1*(1.-alphaK[i_dom[i][j]][1]*1.5*BKTQ*c1));
// 20% of the created particles goes in the second valley
        if(rnd()>0.8){
         IV=2;
         c2=SMH[i_dom[i][j]][IV]*sqrt(-1.5*BKTQ*c1*(1.-alphaK[i_dom[i][j]][2]*1.5*BKTQ*c1));
        }
       }
      }
// The following is in case of precendtly loaded initial data.
      if(leid_flag==1){
       IV=1;
// In this case c1 represents the mean electron energy
// loaded from precedent simulations and have nothing to 
// do with the lattice energy.
       c1=-u2d[i][j][4]/u2d[i][j][1]/Q;
       if(NOVALLEY[i_dom[i][j]]==1) 
        c2=SMH[i_dom[i][j]][0]*sqrt(-1.5*c1*(1.-alphaK[i_dom[i][j]][1]*1.5*c1));
       if(NOVALLEY[i_dom[i][j]]==2){
        IV=1;
        c2=SMH[i_dom[i][j]][IV]*sqrt(-1.5*c1*(1.-alphaK[i_dom[i][j]][1]*1.5*c1));
        if(rnd()>0.8){
         IV=2;
         c2=SMH[i_dom[i][j]][IV]*sqrt(-1.5*c1*(1.-alphaK[i_dom[i][j]][2]*1.5*c1));
        }
       }
      }
      c3=1.-2.*rnd();
      c4=sqrt(1.-c3*c3);
      c5=2.*PI*rnd();
      c6=sin(c5);
      c7=cos(c5);
      P[n][0]=IV;
      P[n][1]=c2*c3*c6;
      P[n][2]=c2*c4*c6;
      P[n][3]=c2*c7;
      P[n][4]=-log(rnd())/GM[i_dom[i][j]];
      P[n][5]=dx*(rnd()+(real)(i)-1.5);
      P[n][6]=dy*(rnd()+(real)(j)-1.5);
      if(i==1) P[n][5]=dx*0.5*rnd();
      if(j==1) P[n][6]=dy*0.5*rnd();
      if(i==nx+1) P[n][5]=LX-dx*0.5*rnd();
      if(j==ny+1) P[n][6]=LY-dy*0.5*rnd();
     }
    }
  }
 INUM = n;
 for(i=1;i<=nx+1;i++)
   for(j=1;j<=ny+1;j++){
     u2d[i][j][2]=0.;
     u2d[i][j][3]=0.;
     u2d[i][j][4]=0.;
   }
 printf("Initial Number of Electron Super-particles = %d\n", INUM);
}

// ==============================================================
