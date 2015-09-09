/* Hole_bcs.h -- This file is part of GNU/Archimedes release 0.0.3.
   Archimedes is a simulator for Submicron 2D Silicon/GaAs
   Devices. It implements the Monte Carlo method and Hybrid MEP model
   for the simulation of the semiclassical Boltzmann equation for both
   holes and holes. It also includes the quantum effects by means 
   of effective potential method.

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
// Created on 04 Mar.2005, Siracusa, J.M.Sellier
// Last modif. : 05 Mar.2005, Siracusa, J.M.Sellier
// ######################################################

// Boundary Conditions for the Hybrid MEP model
// For more informations about this boundaries conditions
// see the manual of GNU/Archimedes release 0.0.3.

void
HoleHMEPBCs(void)
{ 
 int i,j;
 real xvel,yvel;

// These are completely generic boundary conditions

// Bottom Edge
// ===========
   for(i=0;i<=nx+2;i++){
// INSULATOR without potential
     if(EDGE[0][i][0]==0){
// hole Density
       h2d[i+2][0][1]=h2d[i+2][5][1];
       h2d[i+2][1][1]=h2d[i+2][4][1];
       h2d[i+2][2][1]=h2d[i+2][3][1];
// x-component velocity
       h2d[i+2][0][2]=h2d[i+2][5][2];
       h2d[i+2][1][2]=h2d[i+2][4][2];
       h2d[i+2][2][2]=h2d[i+2][3][2];
// y-component velocity
       h2d[i+2][0][3]=h2d[i+2][5][3];
       h2d[i+2][1][3]=h2d[i+2][4][3];
       h2d[i+2][2][3]=h2d[i+2][3][3];
// hole energy
       h2d[i+2][0][4]=h2d[i+2][5][4];
       h2d[i+2][1][4]=h2d[i+2][4][4];
       h2d[i+2][2][4]=h2d[i+2][3][4];
     }
// SCHOTTKY or OHMIC
     if(EDGE[0][i][0]==1 || EDGE[0][i][0]==2){
       yvel=h2d[i+2][3][3]/h2d[i+2][3][1];
// density
       h2d[i+2][0][1]=EDGE[0][i][3];
       h2d[i+2][1][1]=EDGE[0][i][3];
       h2d[i+2][2][1]=EDGE[0][i][3];
// x-component velocity
       h2d[i+2][0][2]=0.;
       h2d[i+2][1][2]=0.;
       h2d[i+2][2][2]=0.;
// y-component velocity
       h2d[i+2][0][3]=EDGE[0][i][3]*yvel;
       h2d[i+2][1][3]=EDGE[0][i][3]*yvel;
       h2d[i+2][2][3]=EDGE[0][i][3]*yvel;
// energy
       h2d[i+2][0][4]=0.5*mstarhole*M*yvel*yvel*EDGE[0][i][3]
                   +1.5*KB*TL*EDGE[0][i][3];
       h2d[i+2][1][4]=0.5*mstarhole*M*yvel*yvel*EDGE[0][i][3]
                   +1.5*KB*TL*EDGE[0][i][3];
       h2d[i+2][2][4]=0.5*mstarhole*M*yvel*yvel*EDGE[0][i][3]
                   +1.5*KB*TL*EDGE[0][i][3];
     }
   }
// Left Edge
// =========
   for(j=1;j<=ny+2;j++){
// INSULATOR without potential
     if(EDGE[3][j][0]==0){
// density
       h2d[0][j+2][1]=h2d[5][j+2][1];
       h2d[1][j+2][1]=h2d[4][j+2][1];
       h2d[2][j+2][1]=h2d[3][j+2][1];
// x-component velocity
       h2d[0][j][2]=h2d[5][j+2][2];
       h2d[1][j][2]=h2d[4][j+2][2];
       h2d[2][j][2]=h2d[3][j+2][2];
// y-component velocity
       h2d[0][j+2][3]=h2d[5][j+2][3];
       h2d[1][j+2][3]=h2d[4][j+2][3];
       h2d[2][j+2][3]=h2d[3][j+2][3];
// hole energy
       h2d[0][j+2][4]=h2d[5][j+2][4];
       h2d[1][j+2][4]=h2d[4][j+2][4];
       h2d[2][j+2][4]=h2d[3][j+2][4];
     }
// SCHOTTKY or OHMIC
     if(EDGE[3][j][0]==1 || EDGE[3][j][0]==2){
       xvel=h2d[3][j+2][2]/h2d[3][j+2][1];
// Density
       h2d[0][j+2][1]=EDGE[3][j][3];
       h2d[1][j+2][1]=EDGE[3][j][3];
       h2d[2][j+2][1]=EDGE[3][j][3];
// x-component velocity
       h2d[0][j+2][2]=EDGE[3][j][3]*xvel;
       h2d[1][j+2][2]=EDGE[3][j][3]*xvel;
       h2d[2][j+2][2]=EDGE[3][j][3]*xvel;
// y-component velocity
       h2d[0][j+2][3]=0.;
       h2d[1][j+2][3]=0.;
       h2d[2][j+2][3]=0.;
// energy
       h2d[0][j+2][4]=0.5*mstarhole*M*xvel*xvel*EDGE[3][j][3]
                   +1.5*KB*TL*EDGE[3][j][3];
       h2d[1][j+2][4]=0.5*mstarhole*M*xvel*xvel*EDGE[3][j][3]
                   +1.5*KB*TL*EDGE[3][j][3];
       h2d[2][j+2][4]=0.5*mstarhole*M*xvel*xvel*EDGE[3][j][3]
                   +1.5*KB*TL*EDGE[3][j][3];
     }
   }
// Right Edge
// ==========
   for(j=1;j<=ny+2;j++){
// INSULATOR
     if(EDGE[1][j][0]==0){
// density
       h2d[nx+2][j+2][1]=h2d[nx+1][j+2][1];
       h2d[nx+3][j+2][1]=h2d[nx][j+2][1];
       h2d[nx+4][j+2][1]=h2d[nx-1][j+2][1];
// x-component velocity
       h2d[nx+2][j+2][2]=h2d[nx+1][j+2][2];
       h2d[nx+3][j+2][2]=h2d[nx][j+2][2];
       h2d[nx+4][j+2][2]=h2d[nx-1][j+2][2];
// y-component velocity
       h2d[nx+2][j+2][3]=h2d[nx+1][j+2][3];
       h2d[nx+3][j+2][3]=h2d[nx][j+2][3];
       h2d[nx+4][j+2][3]=h2d[nx-1][j+2][3];
// energy
       h2d[nx+2][j+2][4]=h2d[nx+1][j+2][4];
       h2d[nx+3][j+2][4]=h2d[nx][j+2][4];
       h2d[nx+4][j+2][4]=h2d[nx-1][j+2][4];
     }
// SCHOTTKY or OHMIC
     if(EDGE[1][j][0]==1 || EDGE[1][j][0]==2){
       xvel=h2d[nx-1][j+2][2]/h2d[nx-1][j+2][1];
// density
       h2d[nx+2][j+2][1]=EDGE[1][j][3];
       h2d[nx+3][j+2][1]=EDGE[1][j][3];
       h2d[nx+4][j+2][1]=EDGE[1][j][3];
// x-component velocity
       h2d[nx+2][j+2][2]=EDGE[1][j][3]*xvel;
       h2d[nx+3][j+2][2]=EDGE[1][j][3]*xvel;
       h2d[nx+4][j+2][2]=EDGE[1][j][3]*xvel;
// y-component velocity
       h2d[nx+2][j+2][3]=0.;
       h2d[nx+3][j+2][3]=0.;
       h2d[nx+4][j+2][3]=0.;
// energy
       h2d[nx+2][j+2][4]=0.5*mstarhole*M*xvel*xvel*EDGE[1][j][3]
                      +1.5*KB*TL*EDGE[1][j][3];
       h2d[nx+3][j+2][4]=0.5*mstarhole*M*xvel*xvel*EDGE[1][j][3]
                      +1.5*KB*TL*EDGE[1][j][3];
       h2d[nx+4][j+2][4]=0.5*mstarhole*M*xvel*xvel*EDGE[1][j][3]
                      +1.5*KB*TL*EDGE[1][j][3];
     }
   }
// Upper Edge
// ==========   
   for(i=1;i<=nx+2;i++){
// INSULATOR
     if(EDGE[2][i][0]==0){
// density
       h2d[i+2][ny+2][1]=h2d[i+2][ny+1][1];
       h2d[i+2][ny+3][1]=h2d[i+2][ny][1];
       h2d[i+2][ny+4][1]=h2d[i+2][ny-1][1];
// x-component velocity
       h2d[i+2][ny+2][2]=h2d[i+2][ny+1][2];
       h2d[i+2][ny+3][2]=h2d[i+2][ny][2];
       h2d[i+2][ny+4][2]=h2d[i+2][ny-1][2];
// y-component velocity
       h2d[i+2][ny+2][3]=h2d[i+2][ny+1][3];
       h2d[i+2][ny+3][3]=h2d[i+2][ny][3];
       h2d[i+2][ny+4][3]=h2d[i+2][ny-1][3];
// energy
       h2d[i+2][ny+2][4]=h2d[i+2][ny+1][4];
       h2d[i+2][ny+3][4]=h2d[i+2][ny][4];
       h2d[i+2][ny+4][4]=h2d[i+2][ny-1][4];
     }
// SCHOTTKY or OHMIC
     if(EDGE[2][i][0]==1 || EDGE[2][i][0]==2){
       yvel=h2d[i+2][ny-1][3]/h2d[i+2][ny-1][1];
// density
       h2d[i+2][ny+2][1]=EDGE[2][i][3];
       h2d[i+2][ny+3][1]=EDGE[2][i][3];
       h2d[i+2][ny+4][1]=EDGE[2][i][3];
// x-component velocity
       h2d[i+2][ny+2][2]=0.;
       h2d[i+2][ny+3][2]=0.;
       h2d[i+2][ny+4][2]=0.;
// y-component velocity
       h2d[i+2][ny+2][3]=EDGE[2][i][3]*yvel;
       h2d[i+2][ny+3][3]=EDGE[2][i][3]*yvel;
       h2d[i+2][ny+4][3]=EDGE[2][i][3]*yvel;
// energy
       h2d[i+2][ny+2][4]=0.5*mstarhole*M*yvel*yvel*EDGE[2][i][3]
                      +1.5*KB*TL*EDGE[2][i][3];
       h2d[i+2][ny+3][4]=0.5*mstarhole*M*yvel*yvel*EDGE[2][i][3]
                      +1.5*KB*TL*EDGE[2][i][3];
       h2d[i+2][ny+4][4]=0.5*mstarhole*M*yvel*yvel*EDGE[2][i][3]
                      +1.5*KB*TL*EDGE[2][i][3];
     }
  }
}

// ######################################################
