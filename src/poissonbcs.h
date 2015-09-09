/* poissonbcs.h -- This file is part of Archimedes release 0.0.5.
   Archimedes is a simulator for Submicron 2D Silicon/GaAs
   Devices. It implements the Monte Carlo method and a simplified
   MEP model for the simulation of the semiclassical Boltzmann
   equation for both electrons and holes. It also includes the
   quantum effects by means of effective potential method.

   Copyright (C) 2004, 2005, 2006, 2007 Jean Michel Sellier <sellier@dmi.unict.it>
 
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA.  */


// ######################################################
// Created on 24 Mar.2004, Siracusa, J.M.Sellier
// Last modif. : 07 Apr.2005, Siracusa, J.M.Sellier
// ######################################################

// Boundary Conditions for the non-stationary Poisson equation
// For more informations about this equation
// see the manual of GNU Archimedes release 0.0.5.

void
PoissonBCs(void)
{ 
 int i,j;

// These are completely generic boundary conditions

// Bottom Edge
// ===========
   for(i=1;i<=nx+1;i++){
// INSULATOR without potential
     if(EDGE[0][i][0]==0){
       u2d[i][0][0]=u2d[i][3][0];
       u2d[i][1][0]=u2d[i][2][0];
     }
// INSULATOR with potential
     if(EDGE[0][i][0]==0 && EDGE[0][i][1]!=0.0){
       u2d[i][0][0]=EDGE[0][i][1];
       u2d[i][1][0]=EDGE[0][i][1];
     }
// SCHOTTKY or OHMIC
     if(EDGE[0][i][0]==1 || EDGE[0][i][0]==2){
       u2d[i][0][0]=EDGE[0][i][1];
       u2d[i][1][0]=EDGE[0][i][1];
     }
   }
// Eventual lower SiO2 interface
// ###
   if(SIO2_DOWN_FLAG==1)
   {
    int sio2nx=(int)(fabs(SIO2_INI[1]-SIO2_FIN[1])/dx);
    int sio2ny=(int)(SIO2_THICKNESS[1]/dy);
    if(sio2nx==0) sio2nx++;
    if(sio2ny==0) sio2ny++;
    for(i=1;i<=sio2nx+1;i++){
       u2d[i+(int)(SIO2_INI[1]/dy)][0][0]=SIO2[1][i][sio2ny+2];
       u2d[i+(int)(SIO2_INI[1]/dy)][1][0]=SIO2[1][i][sio2ny+1];
    }
   }
// ###
// Left Edge
// =========
   for(j=1;j<=ny+1;j++){
// INSULATOR without potential
     if(EDGE[3][j][0]==0){
       u2d[0][j][0]=u2d[3][j][0];
       u2d[1][j][0]=u2d[2][j][0];
     }
// INSULATOR with potential
     if(EDGE[3][j][0]==0 && EDGE[3][j][1]!=0.0){
       u2d[0][j][0]=EDGE[3][j][1];
       u2d[1][j][0]=EDGE[3][j][1];
     }
// SCHOTTKY or OHMIC
     if(EDGE[3][j][0]==1 || EDGE[3][j][0]==2){
       u2d[0][j][0]=EDGE[3][j][1];
       u2d[1][j][0]=EDGE[3][j][1];
     }
   }
// Right Edge
// ==========
   for(j=1;j<=ny+1;j++){
// INSULATOR without potential
     if(EDGE[1][j][0]==0){
       u2d[nx+1][j][0]=u2d[nx-1][j][0];
       u2d[nx+2][j][0]=u2d[nx][j][0];
     }
// INSULATOR with potential
     if(EDGE[1][j][0]==0 && EDGE[1][j][1]!=0){
       u2d[nx+1][j][0]=EDGE[1][j][1];
       u2d[nx+2][j][0]=EDGE[1][j][1];
     }
// SCHOTTKY or OHMIC
     if(EDGE[1][j][0]==1 || EDGE[1][j][0]==2){
       u2d[nx+1][j][0]=EDGE[1][j][1];
       u2d[nx+2][j][0]=EDGE[1][j][1];
     }
   }
// Upper Edge
// ==========   
   for(i=1;i<=nx+1;i++){
// INSULATOR without potential
     if(EDGE[2][i][0]==0){
       u2d[i][ny+1][0]=u2d[i][ny][0];
       u2d[i][ny+2][0]=u2d[i][ny-1][0];
     }
// INSULATOR with potential
     if(EDGE[2][i][0]==0 && EDGE[2][i][1]!=0){
       u2d[i][ny+1][0]=EDGE[2][i][1];
       u2d[i][ny+2][0]=EDGE[2][i][1];
     }
// SCHOTTKY or OHMIC
     if(EDGE[2][i][0]==1 || EDGE[2][i][0]==2){
       u2d[i][ny+1][0]=EDGE[2][i][1];
       u2d[i][ny+2][0]=EDGE[2][i][1];
     }
  }
// Eventual upper SiO2 interface
// ###
   if(SIO2_UP_FLAG==1)
   {
    int sio2nx=(int)(fabs(SIO2_INI[0]-SIO2_FIN[0])/dx);
    int sio2ny=(int)(SIO2_THICKNESS[0]/dy);
    if(sio2nx==0) sio2nx++;
    if(sio2ny==0) sio2ny++;
    for(i=1;i<=sio2nx+1;i++){
       u2d[i+(int)(SIO2_INI[0]/dy)][ny+1][0]=SIO2[0][i][0];
       u2d[i+(int)(SIO2_INI[0]/dy)][ny+2][0]=SIO2[0][i][1];
    }
   }
// ###
}

// ######################################################
