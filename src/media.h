/* media.h -- This file is part of Archimedes release 1.5.0
   Archimedes is a simulator for Submicron and Nanoscaled
   2D III-V Semiconductor Devices.

   Copyright (C) 2004-2011 Jean Michel Sellier
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

// ######################################################
// Created on 05 Oct.2004, Siracusa, J.M.Sellier
// Last modif. : 02 Sep. 2011, Carry le Rouet, J.M.Sellier
// ######################################################

void media(void)
{
 register int i,j,n;
 int iv;
 int cont[NXM+1][NYM+1];
 real xvel[NXM+1][NYM+1],yvel[NXM+1][NYM+1],ener[NXM+1][NYM+1];
 real x,y,xvelocity,yvelocity,ksquared,thesquareroot;
 real superparticle_energy;

 printf("Computation of macroscopic observables\n");

// resetting of the electronic density
// a simple way to avoid NaN propagation...
 for(i=0;i<=NXM;i++)
   for(j=0;j<=NYM;j++){
     cont[i][j]=0;
     xvel[i][j]=0.;
     yvel[i][j]=0.;
     ener[i][j]=0.;
   }
 for(n=1;n<=INUM;n++){
   int iaux;
   iv=(int)(P[n][0]);
   x=P[n][5];
   y=P[n][6];
   i=(int)(x/dx+1.5);
   j=(int)(y/dy+1.5);
   if(i<=1)  i=1;
   if(j<=1)  j=1;
   if(i>=nx+1) i=nx+1;
   if(j>=ny+1) j=ny+1;
   if(NOVALLEY[i_dom[i][j]]==1) iaux=0;
   if(NOVALLEY[i_dom[i][j]]==2) iaux=iv;
   if(CONDUCTION_BAND==FULL) iaux=0;

   ksquared=P[n][1]*P[n][1]+P[n][2]*P[n][2]+P[n][3]*P[n][3];

   if(CONDUCTION_BAND==PARABOLIC){
    superparticle_energy=HHM[i_dom[i][j]][0]*ksquared;
    xvelocity=P[n][1]*HM[i_dom[i][j]][iaux];
    yvelocity=P[n][2]*HM[i_dom[i][j]][iaux];
   }
   if(CONDUCTION_BAND==KANE){
    thesquareroot=sqrt(1.+4.*alphaK[i_dom[i][j]][iv]*HHM[i_dom[i][j]][iaux]*ksquared);
    superparticle_energy=(thesquareroot-1.)/(2.*alphaK[i_dom[i][j]][iv]);
    xvelocity=P[n][1]*HM[i_dom[i][j]][iaux]/thesquareroot;
    yvelocity=P[n][2]*HM[i_dom[i][j]][iaux]/thesquareroot;
   }
   if(CONDUCTION_BAND==FULL){
    real k,k2,k4;
    real dx,dy,d;
    k=sqrt(ksquared)*0.5/PI*1.e-12;
    // periodicity on reciprocal lattice
    k2=k*k;
    k4=k2*k2;
    superparticle_energy=CB_FULL[i_dom[i][j]][0]*k4*k4*k2
                        +CB_FULL[i_dom[i][j]][1]*k4*k4*k
                        +CB_FULL[i_dom[i][j]][2]*k4*k4
                        +CB_FULL[i_dom[i][j]][3]*k4*k2*k
                        +CB_FULL[i_dom[i][j]][4]*k4*k2
                        +CB_FULL[i_dom[i][j]][5]*k4*k
                        +CB_FULL[i_dom[i][j]][6]*k4
                        +CB_FULL[i_dom[i][j]][7]*k2*k
                        +CB_FULL[i_dom[i][j]][8]*k2
                        +CB_FULL[i_dom[i][j]][9]*k
                        +CB_FULL[i_dom[i][j]][10]; // in eV

    d=10.*CB_FULL[i_dom[i][j]][0]*k4*k4*k
     +9.*CB_FULL[i_dom[i][j]][1]*k4*k4
     +8.*CB_FULL[i_dom[i][j]][2]*k4*k2*k
     +7.*CB_FULL[i_dom[i][j]][3]*k4*k2
     +6.*CB_FULL[i_dom[i][j]][4]*k4*k
     +5.*CB_FULL[i_dom[i][j]][5]*k4
     +4.*CB_FULL[i_dom[i][j]][6]*k2*k
     +3.*CB_FULL[i_dom[i][j]][7]*k2
     +2.*CB_FULL[i_dom[i][j]][8]*k
     +CB_FULL[i_dom[i][j]][9];
    k*=1.e+12*2.*PI;
    d*=1.e-12*0.5/PI;
    xvelocity=QH*d*P[n][1]/k;
    yvelocity=QH*d*P[n][2]/k;
   }
// for the following two rows see
// pag.10 formula (1.19) of Tomizawa,
// "Numerical Simulation of Submicron Semiconductor
//  Devices", 1993, ARTECH HOUSE
// (this are due to the Kane non-parabolic energy band)
   i=(int)(x/dx+1.5);
   j=(int)(y/dy+1.5);
   if(i<1) i=1;
   if(j<1) j=1;
   if(i>=nx+1) i=nx+1;
   if(j>=ny+1) j=ny+1;
   cont[i][j]++;
   ener[i][j]+=superparticle_energy;
   ener[i][j]+=EMIN[i_dom[i][j]][iv];
   xvel[i][j]+=xvelocity;
   yvel[i][j]+=yvelocity;
 }
// Mean Value of the macroscopic variables
// =======================================
  for(i=1;i<=nx+1;i++)
    for(j=1;j<=ny+1;j++){
      if(cont[i][j]!=0){
          xvel[i][j]/=(real)cont[i][j];
          yvel[i][j]/=(real)cont[i][j];
          ener[i][j]/=(real)cont[i][j];
      }
    }
// =======================================
  for(i=1;i<=nx+1;i++)
    for(j=1;j<=ny+1;j++){
      u2d[i][j][2]+=xvel[i][j];
      u2d[i][j][3]+=yvel[i][j];
      u2d[i][j][4]+=ener[i][j];
    }
}

// ================================================
