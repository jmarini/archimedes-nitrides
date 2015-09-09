/* HMEPbcs.h -- This file is part of GNU/Archimedes release 0.1.0.
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method and Hybrid MEP model
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
// Created on 26 Feb.2005, Siracusa, J.M.Sellier
// Last modif. : 09 Sep.2007, Siracusa, J.M.Sellier
// ######################################################

// Boundary Conditions for the Hybrid MEP model
// For more informations about this boundaries conditions
// see the manual of GNU/Archimedes release 0.0.3.

void
HMEPBCs(void)
{ 
 int i,j;
 real xvel,yvel;

// These are completely generic boundary conditions

// Bottom Edge
// ===========
   for(i=0;i<=nx+2;i++){
// INSULATOR without potential
     if(EDGE[0][i][0]==0){
// electron Density
       u2d[i+2][0][1]=u2d[i+2][5][1];
       u2d[i+2][1][1]=u2d[i+2][4][1];
       u2d[i+2][2][1]=u2d[i+2][3][1];
// x-component velocity
       u2d[i+2][0][2]=u2d[i+2][5][2];
       u2d[i+2][1][2]=u2d[i+2][4][2];
       u2d[i+2][2][2]=u2d[i+2][3][2];
// y-component velocity
       u2d[i+2][0][3]=u2d[i+2][5][3];
       u2d[i+2][1][3]=u2d[i+2][4][3];
       u2d[i+2][2][3]=u2d[i+2][3][3];
// electron energy
       u2d[i+2][0][4]=u2d[i+2][5][4];
       u2d[i+2][1][4]=u2d[i+2][4][4];
       u2d[i+2][2][4]=u2d[i+2][3][4];
     }
// SCHOTTKY or OHMIC
     if(EDGE[0][i][0]==1 || EDGE[0][i][0]==2){
       yvel=u2d[i+2][3][3]/u2d[i+2][3][1];
// density
       u2d[i+2][0][1]=EDGE[0][i][2];
       u2d[i+2][1][1]=EDGE[0][i][2];
       u2d[i+2][2][1]=EDGE[0][i][2];
// x-component velocity
       u2d[i+2][0][2]=0.;
       u2d[i+2][1][2]=0.;
       u2d[i+2][2][2]=0.;
// y-component velocity
       u2d[i+2][0][3]=EDGE[0][i][2]*yvel;
       u2d[i+2][1][3]=EDGE[0][i][2]*yvel;
       u2d[i+2][2][3]=EDGE[0][i][2]*yvel;
// energy
       u2d[i+2][0][4]=0.5*MSTAR[i_dom[i+2][3]][1]*M*yvel*yvel*EDGE[0][i][2]
                   +1.5*KB*TL*EDGE[0][i][2];
       u2d[i+2][1][4]=0.5*MSTAR[i_dom[i+2][3]][1]*M*yvel*yvel*EDGE[0][i][2]
                   +1.5*KB*TL*EDGE[0][i][2];
       u2d[i+2][2][4]=0.5*MSTAR[i_dom[i+2][3]][1]*M*yvel*yvel*EDGE[0][i][2]
                   +1.5*KB*TL*EDGE[0][i][2];
     }
   }
// Left Edge
// =========
   for(j=1;j<=ny+2;j++){
// INSULATOR without potential
     if(EDGE[3][j][0]==0){
// density
       u2d[0][j+2][1]=u2d[5][j+2][1];
       u2d[1][j+2][1]=u2d[4][j+2][1];
       u2d[2][j+2][1]=u2d[3][j+2][1];
// x-component velocity
       u2d[0][j][2]=u2d[5][j+2][2];
       u2d[1][j][2]=u2d[4][j+2][2];
       u2d[2][j][2]=u2d[3][j+2][2];
// y-component velocity
       u2d[0][j+2][3]=u2d[5][j+2][3];
       u2d[1][j+2][3]=u2d[4][j+2][3];
       u2d[2][j+2][3]=u2d[3][j+2][3];
// electron energy
       u2d[0][j+2][4]=u2d[5][j+2][4];
       u2d[1][j+2][4]=u2d[4][j+2][4];
       u2d[2][j+2][4]=u2d[3][j+2][4];
     }
// SCHOTTKY or OHMIC
     if(EDGE[3][j][0]==1 || EDGE[3][j][0]==2){
       xvel=u2d[3][j+2][2]/u2d[3][j+2][1];
// Density
       u2d[0][j+2][1]=EDGE[3][j][2];
       u2d[1][j+2][1]=EDGE[3][j][2];
       u2d[2][j+2][1]=EDGE[3][j][2];
// x-component velocity
       u2d[0][j+2][2]=EDGE[3][j][2]*xvel;
       u2d[1][j+2][2]=EDGE[3][j][2]*xvel;
       u2d[2][j+2][2]=EDGE[3][j][2]*xvel;
// y-component velocity
       u2d[0][j+2][3]=0.;
       u2d[1][j+2][3]=0.;
       u2d[2][j+2][3]=0.;
// energy
       u2d[0][j+2][4]=0.5*MSTAR[i_dom[3][j+2]][1]*M*xvel*xvel*EDGE[3][j][2]
                   +1.5*KB*TL*EDGE[3][j][2];
       u2d[1][j+2][4]=0.5*MSTAR[i_dom[3][j+2]][1]*M*xvel*xvel*EDGE[3][j][2]
                   +1.5*KB*TL*EDGE[3][j][2];
       u2d[2][j+2][4]=0.5*MSTAR[i_dom[3][j+2]][1]*M*xvel*xvel*EDGE[3][j][2]
                   +1.5*KB*TL*EDGE[3][j][2];
     }
   }
// Right Edge
// ==========
   for(j=1;j<=ny+2;j++){
// INSULATOR
     if(EDGE[1][j][0]==0){
// density
       u2d[nx+2][j+2][1]=u2d[nx+1][j+2][1];
       u2d[nx+3][j+2][1]=u2d[nx][j+2][1];
       u2d[nx+4][j+2][1]=u2d[nx-1][j+2][1];
// x-component velocity
       u2d[nx+2][j+2][2]=u2d[nx+1][j+2][2];
       u2d[nx+3][j+2][2]=u2d[nx][j+2][2];
       u2d[nx+4][j+2][2]=u2d[nx-1][j+2][2];
// y-component velocity
       u2d[nx+2][j+2][3]=u2d[nx+1][j+2][3];
       u2d[nx+3][j+2][3]=u2d[nx][j+2][3];
       u2d[nx+4][j+2][3]=u2d[nx-1][j+2][3];
// energy
       u2d[nx+2][j+2][4]=u2d[nx+1][j+2][4];
       u2d[nx+3][j+2][4]=u2d[nx][j+2][4];
       u2d[nx+4][j+2][4]=u2d[nx-1][j+2][4];
     }
// SCHOTTKY or OHMIC
     if(EDGE[1][j][0]==1 || EDGE[1][j][0]==2){
       xvel=u2d[nx-1][j+2][2]/u2d[nx-1][j+2][1];
// density
       u2d[nx+2][j+2][1]=EDGE[1][j][2];
       u2d[nx+3][j+2][1]=EDGE[1][j][2];
       u2d[nx+4][j+2][1]=EDGE[1][j][2];
// x-component velocity
       u2d[nx+2][j+2][2]=EDGE[1][j][2]*xvel;
       u2d[nx+3][j+2][2]=EDGE[1][j][2]*xvel;
       u2d[nx+4][j+2][2]=EDGE[1][j][2]*xvel;
// y-component velocity
       u2d[nx+2][j+2][3]=0.;
       u2d[nx+3][j+2][3]=0.;
       u2d[nx+4][j+2][3]=0.;
// energy
       u2d[nx+2][j+2][4]=0.5*MSTAR[i_dom[nx-1][j+2]][1]*M*xvel*xvel*EDGE[1][j][2]
                      +1.5*KB*TL*EDGE[1][j][2];
       u2d[nx+3][j+2][4]=0.5*MSTAR[i_dom[nx-1][j+2]][1]*M*xvel*xvel*EDGE[1][j][2]
                      +1.5*KB*TL*EDGE[1][j][2];
       u2d[nx+4][j+2][4]=0.5*MSTAR[i_dom[nx-1][j+2]][1]*M*xvel*xvel*EDGE[1][j][2]
                      +1.5*KB*TL*EDGE[1][j][2];
     }
   }
// Upper Edge
// ==========   
   for(i=1;i<=nx+2;i++){
// INSULATOR
     if(EDGE[2][i][0]==0){
// density
       u2d[i+2][ny+2][1]=u2d[i+2][ny+1][1];
       u2d[i+2][ny+3][1]=u2d[i+2][ny][1];
       u2d[i+2][ny+4][1]=u2d[i+2][ny-1][1];
// x-component velocity
       u2d[i+2][ny+2][2]=u2d[i+2][ny+1][2];
       u2d[i+2][ny+3][2]=u2d[i+2][ny][2];
       u2d[i+2][ny+4][2]=u2d[i+2][ny-1][2];
// y-component velocity
       u2d[i+2][ny+2][3]=u2d[i+2][ny+1][3];
       u2d[i+2][ny+3][3]=u2d[i+2][ny][3];
       u2d[i+2][ny+4][3]=u2d[i+2][ny-1][3];
// energy
       u2d[i+2][ny+2][4]=u2d[i+2][ny+1][4];
       u2d[i+2][ny+3][4]=u2d[i+2][ny][4];
       u2d[i+2][ny+4][4]=u2d[i+2][ny-1][4];
     }
// SCHOTTKY or OHMIC
     if(EDGE[2][i][0]==1 || EDGE[2][i][0]==2){
       yvel=u2d[i+2][ny-1][3]/u2d[i+2][ny-1][1];
// density
       u2d[i+2][ny+2][1]=EDGE[2][i][2];
       u2d[i+2][ny+3][1]=EDGE[2][i][2];
       u2d[i+2][ny+4][1]=EDGE[2][i][2];
// x-component velocity
       u2d[i+2][ny+2][2]=0.;
       u2d[i+2][ny+3][2]=0.;
       u2d[i+2][ny+4][2]=0.;
// y-component velocity
       u2d[i+2][ny+2][3]=EDGE[2][i][2]*yvel;
       u2d[i+2][ny+3][3]=EDGE[2][i][2]*yvel;
       u2d[i+2][ny+4][3]=EDGE[2][i][2]*yvel;
// energy
       u2d[i+2][ny+2][4]=0.5*MSTAR[i_dom[i+2][ny-1]][1]*M*yvel*yvel*EDGE[2][i][2]
                      +1.5*KB*TL*EDGE[2][i][2];
       u2d[i+2][ny+3][4]=0.5*MSTAR[i_dom[i+2][ny-1]][1]*M*yvel*yvel*EDGE[2][i][2]
                      +1.5*KB*TL*EDGE[2][i][2];
       u2d[i+2][ny+4][4]=0.5*MSTAR[i_dom[i+2][ny-1]][1]*M*yvel*yvel*EDGE[2][i][2]
                      +1.5*KB*TL*EDGE[2][i][2];
     }
  }
}

// ######################################################
