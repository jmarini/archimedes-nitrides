/* charge.h -- This file is part of Archimedes release 0.1.2
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method
   for the simulation of the semiclassical Boltzmann equation for
   electrons. It also includes the quantum effects by means 
   of Bohm effective potential method. It is now able to simulate applied
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
// Created on 06 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 12 Sep.2007, Siracusa, J.M.Sellier
// ######################################################

// computation of the macroscopic variables
// (density, x momentum, y momentum, energy)

void
Charge(void)
{
 int i,j,n;
 real x,y;

// resetting of the electronic density
// a simple way to avoid NaN propagation...
 for(i=1;i<=nx+1;i++)
   for(j=1;j<=ny+1;j++)
     u2d[i][j][1]=0.;

//if(Material==SILICON || Material==GERMANIUM){
// well known "cloud in cell" method
  for(n=1;n<=INUM;n++){
    x=P[n][5]/dx;
    y=P[n][6]/dy;
    i=(int)(x+1.);
    j=(int)(y+1.);
// Cloud in cell method
    u2d[i][j][1]+=(1.-(x-(real)(i-1)))*(1.-(y-(real)(j-1)));
    if(i<=nx) 
      u2d[i+1][j][1]+=(1.-(1.-(x-(real)(i-1))))*(1.-(y-(real)(j-1)));
    if(j<=ny)
      u2d[i][j+1][1]+=(1.-(x-(real)(i-1)))*(1.-(1.-(y-(real)(j-1))));
    if(i<=nx && j<=ny)
     u2d[i+1][j+1][1]+=(1.-(1.-(x-(real)(i-1))))*(1.-(1.-(y-(real)(j-1))));
  }
// =======================================
  for(i=1;i<=nx+1;i++)
    for(j=1;j<=ny+1;j++){
      u2d[i][j][1]*=EPP/(dx*dy);
      if(i==1 || i==nx+1) u2d[i][j][1]*=2.;
      if(j==1 || j==ny+1) u2d[i][j][1]*=2.;
    }
// This trick is to avoid the strange oscillations
// in the ghost cells, when we plot the density.
// This does not influence the resolution of the Poisson
// equation, since this last is not computed on the
// ghost cells. 08-12 oct.2004, J.M.Sellier, Siracusa.
// bottom and upper edge contacts
  u2d[nx+1][ny+1][1]=u2d[nx][ny+1][1];
  u2d[1][ny+1][1]=u2d[1][ny][1];
// }
/*if(Material==GAAS){
// well known "cloud in cell" method
  for(n=1;n<=INUM;n++){
    x=P[n][5]/dx;
    y=P[n][6]/dy;
    i=(int)(x+1.);
    j=(int)(y+1.);
// Cloud in cell method
// Gamma-valley
    if(P[n][0]==1){
      DG[i][j]+=(1.-(x-(real)(i-1)))*(1.-(y-(real)(j-1)));
      if(i<=nx) 
        DG[i+1][j]+=(1.-(1.-(x-(real)(i-1))))*(1.-(y-(real)(j-1)));
      if(j<=ny)
        DG[i][j+1]+=(1.-(x-(real)(i-1)))*(1.-(1.-(y-(real)(j-1))));
      if(i<=nx && j<=ny)
       DG[i+1][j+1]+=(1.-(1.-(x-(real)(i-1))))*(1.-(1.-(y-(real)(j-1))));
    }
// L-valley
    if(P[n][0]==2){
      DL[i][j]+=(1.-(x-(real)(i-1)))*(1.-(y-(real)(j-1)));
      if(i<=nx) 
        DL[i+1][j]+=(1.-(1.-(x-(real)(i-1))))*(1.-(y-(real)(j-1)));
      if(j<=ny)
        DL[i][j+1]+=(1.-(x-(real)(i-1)))*(1.-(1.-(y-(real)(j-1))));
      if(i<=nx && j<=ny)
       DL[i+1][j+1]+=(1.-(1.-(x-(real)(i-1))))*(1.-(1.-(y-(real)(j-1))));
    }

  }
// =======================================
  for(i=1;i<=nx+1;i++)
    for(j=1;j<=ny+1;j++){
      DG[i][j]*=EPP/(dx*dy);
      DL[i][j]*=EPP/(dx*dy)/4.;
      u2d[i][j][1]+=DG[i][j]+4.*DL[i][j];
      if(i==1 || i==nx+1){
        u2d[i][j][1]*=2.;
        DG[i][j]*=2.;
        DL[i][j]*=2.;
      }
      if(j==1 || j==ny+1){
        u2d[i][j][1]*=2.;
        DG[i][j]*=2.;
        DL[i][j]*=2.;
      }
    }
// This trick is to avoid the strange oscillations
// in the ghost cells, when we plot the density.
// This does not influence the resolution of the Poisson
// equation, since this last is not computed on the
// ghost cells. 08-12 oct.2004, J.M.Sellier, Siracusa.
// bottom and upper edge contacts
  u2d[nx+1][ny+1][1]=u2d[nx][ny+1][1];
  u2d[1][ny+1][1]=u2d[1][ny][1];
 }*/
}

// ===========================================================
