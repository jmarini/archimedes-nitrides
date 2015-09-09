/* electric_field.h -- This file is part of GNU archimedes

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


// ######################################################
// Created on 10 Mar.2004, Siracusa, J.M.Sellier
// Last modif. : 01 Sep.2011, Carry le Rouet, J.M.Sellier
// ######################################################

// Computation of the electrostatic potential,
// i.e. resolution of the 2D Poisson equation,
// by means of the computation of the stationary 
// solution of a pseudo-transient Poisson equation.
// From version 0.1.0 on, the potential is added to the 
// the minimum energy of the semiconductor material
// in order to take into account heterostructures.
// For more information see the manual of 
// GNU Archimedes release 1.0.0.

void 
Electric_Field(void)
{
 register int i,j,k;

 real factor,kappa,deltat,rho;

 PoissonBCs();
// ===============================
 factor=0.9;
// calcolo del potenziale "stazionario"
 for(i=1;i<=POISSONITMAX;i++){
// Eventual Upper SiO2
// ###
   if(SIO2_UP_FLAG==1)
   {
    int sio2nx=(int)(fabs(SIO2_INI[0]-SIO2_FIN[0])/dx);
    int sio2ny=(int)(SIO2_THICKNESS[0]/dy);
    if(sio2nx==0) sio2nx++;
    if(sio2ny==0) sio2ny++;
// BCs for upper SiO2
// ***
// upper boundary
    for(j=1;j<=sio2nx+1;j++){
      SIO2[0][j][sio2ny+1]=SIO2_POT[0];
      SIO2[0][j][sio2ny+2]=SIO2_POT[0];
    }
// left and right boundaries
    for(j=1;j<=sio2ny+1;j++){
      SIO2[0][0][j]=SIO2[0][3][j];
      SIO2[0][1][j]=SIO2[0][2][j];
      SIO2[0][sio2nx+1][j]=SIO2[0][sio2nx-1][j];
      SIO2[0][sio2nx+2][j]=SIO2[0][sio2nx][j];
    }
// lower transmissive boundary
    for(j=1;j<=sio2nx+1;j++){
      SIO2[0][j][0]=u2d[j+(int)(SIO2_INI[0]/dx)][ny-1][0];
      SIO2[0][j][1]=u2d[j+(int)(SIO2_INI[0]/dx)][ny][0];
    }
// ***
    for(k=2;k<=sio2ny;k++)
     for(j=2;j<=sio2nx;j++){
       kappa=(EPSRSIO2)/Q;
       deltat=factor*0.5/kappa/(1./(dx*dx)+1./(dy*dy));
       SIO2[0][j][k]=SIO2[0][j][k]+deltat*kappa*
          ((SIO2[0][j+1][k]-2.0*SIO2[0][j][k]+SIO2[0][j-1][k])/(dx*dx)
          +(SIO2[0][j][k+1]-2.0*SIO2[0][j][k]+SIO2[0][j][k-1])/(dy*dy));
     }
   }
// ###
// Semiconducting material aprt
   PoissonBCs();
   for(j=0;j<=ny+2;j++)
     for(k=0;k<=nx+2;k++)
       PSI[k][j]=u2d[k][j][0];
   for(k=2;k<=ny;k++)
     for(j=2;j<=nx;j++){
       kappa=(EPSR[i_dom[k][j]]*EPS0)/Q;
       deltat=factor*0.5/kappa/(1./(dx*dx)+1./(dy*dy));
       rho=(u2d[j][k][1]-N_D[j][k]-h2d[j][k][1]+N_H[j][k]);
       u2d[j][k][0]=PSI[j][k]-deltat*rho+deltat*kappa*
             ((PSI[j+1][k]-2.0*PSI[j][k]+PSI[j-1][k])/(dx*dx)
             +(PSI[j][k+1]-2.0*PSI[j][k]+PSI[j][k-1])/(dy*dy));
     }
// Eventual lower SiO2
// ###
   if(SIO2_DOWN_FLAG==1)
   {
    int sio2nx=(int)(fabs(SIO2_INI[1]-SIO2_FIN[1])/dx);
    int sio2ny=(int)(SIO2_THICKNESS[1]/dy);
    if(sio2nx==0) sio2nx++;
    if(sio2ny==0) sio2ny++;
// BCs for lower SiO2
// ***
// lower boundary
    for(j=1;j<=sio2nx+1;j++){
      SIO2[1][j][0]=SIO2_POT[1];
      SIO2[1][j][1]=SIO2_POT[1];
    }
// left and right boundaries
    for(j=1;j<=sio2ny+1;j++){
      SIO2[1][0][j]=SIO2[1][3][j];
      SIO2[1][1][j]=SIO2[1][2][j];
      SIO2[1][sio2nx+1][j]=SIO2[1][sio2nx-1][j];
      SIO2[1][sio2nx+2][j]=SIO2[1][sio2nx][j];
    }
// upper transmissive boundary
    for(j=1;j<=sio2nx+1;j++){
      SIO2[1][j][sio2ny+1]=u2d[j+(int)(SIO2_INI[1]/dx)][1][0];
      SIO2[1][j][sio2ny+2]=u2d[j+(int)(SIO2_INI[1]/dx)][0][0];
    }
// ***
    for(k=2;k<=sio2ny;k++)
     for(j=2;j<=sio2nx;j++){
       SIO2[1][j][k]=SIO2[1][j][k]+deltat*kappa*
          ((SIO2[1][j+1][k]-2.0*SIO2[1][j][k]+SIO2[1][j-1][k])/(dx*dx)
          +(SIO2[1][j][k+1]-2.0*SIO2[1][j][k]+SIO2[1][j][k-1])/(dy*dy));
     }
   }
// ###
 }

// ===============================
 PoissonBCs();

// We save the classical potential
// and we substract the energy minimum of the
// semiconductor material in order to 
// take into account heterostructures
   for(j=0;j<=ny+1;j++)
     for(i=0;i<=nx+1;i++)
       u2d[i][j][0]-=EMIN[i_dom[i][j]][1];

 if(Quantum_Flag==1){
   printf("Calculation of Quantum Effective Potential\n");
// We take in account the Quantum Effects
   quantum_effective_potential();
   PoissonBCs();
  }

// Computation of the X-component of the electric Field
// ====================================================
  for(k=1;k<=ny+1;k++)
   for(j=2;j<=nx;j++)
    E[j][k][0]=-0.5*(u2d[j+1][k][0]-u2d[j-1][k][0])/dx;
  for(k=1;k<=ny+1;k++){
    E[1][k][0]=E[2][k][0];
    E[nx+1][k][0]=E[nx][k][0];
  }

// Computation of the Y-component of the electric Field
// ====================================================
  for(k=2;k<=ny;k++)
   for(j=1;j<=nx+1;j++)
    E[j][k][1]=-0.5*(u2d[j][k+1][0]-u2d[j][k-1][0])/dy;
  for(j=1;j<=nx+1;j++){
    E[j][1][1]=E[j][2][1];
    E[j][ny+1][1]=E[j][ny][1];
  }
}

// ==============================================================
