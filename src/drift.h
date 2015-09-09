/* drift.h -- This file is part of Archimedes release 1.2.0
   Archimedes is a simulator for Submicron 2D III-V Semiconductor
   Devices. It implements the Monte Carlo method
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It also includes the quantum effects by means 
   of effective potential method. It is now able to simulate applied
   magnetic fields along with self consistent Faraday equation.

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
// Created on 06 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 31 Aug.2011, Carry le Rouet, J.M.Sellier
// ######################################################

// calculation of drift process

void
drift(real tau)
{
 int iaux;
 int i,j;
 real dkx,dky,hmt,ksquared;
 real vx,vy;

 if(IV==9) return;

 i=(int)(X/dx)+1;
 j=(int)(Y/dy)+1;
 if(i<=1) i=1;
 if(j<=1) j=1;
 if(i>=nx) i=nx;
 if(j>=ny) j=ny;

 if(NOVALLEY[i_dom[i][j]]==1) iaux=0;
 if(NOVALLEY[i_dom[i][j]]==2) iaux=IV;

// Electron drift process
// second order Runge-Kutta method
 hmt=HM[i_dom[i][j]][iaux]*tau;
 ksquared=KX*KX+KY*KY+KZ*KZ;
 if(CONDUCTION_BAND==KANE){
  real thesquareroot,gk;
  gk=HHM[i_dom[i][j]][iaux]*ksquared;
  thesquareroot=sqrt(1.+4.*alphaK[i_dom[i][j]][IV]*gk);
  vx=KX*HM[i_dom[i][j]][iaux]/thesquareroot;
  vy=KY*HM[i_dom[i][j]][iaux]/thesquareroot;
  dkx=-QH*(E[i][j][0]+vy*B[i][j])*tau;
  dky=-QH*(E[i][j][1]-vx*B[i][j])*tau;
  X+=hmt*(KX+0.5*dkx)/thesquareroot;
  Y+=hmt*(KY+0.5*dky)/thesquareroot;
  KX+=dkx;
  KY+=dky;
 }
 if(CONDUCTION_BAND==PARABOLIC){
  vx=KX*HM[i_dom[i][j]][iaux];
  vy=KY*HM[i_dom[i][j]][iaux];
  dkx=-QH*(E[i][j][0]+vy*B[i][j])*tau;
  dky=-QH*(E[i][j][1]-vx*B[i][j])*tau;
  X+=hmt*(KX+0.5*dkx);
  Y+=hmt*(KY+0.5*dky);
  KX+=dkx;
  KY+=dky;
 }
 if(CONDUCTION_BAND==FULL){
  real k4,k2,ks;
  real dx,dy,d;
  vx=KX*HM[i_dom[i][j]][iaux];
  vy=KY*HM[i_dom[i][j]][iaux];
  dkx=-QH*(E[i][j][0]+vy*B[i][j])*tau;
  dky=-QH*(E[i][j][1]-vx*B[i][j])*tau;
  k2=(KX+0.5*dkx)*(KX+0.5*dkx)+(KY+0.5*dky)*(KY+0.5*dky)+KZ*KZ;
  ks=sqrt(k2)*1.e-12*0.5/PI;
  k2=ks*ks;
  k4=k2*k2;
  d=10.*CB_FULL[i_dom[i][j]][0]*k4*k4*ks
    +9.*CB_FULL[i_dom[i][j]][1]*k4*k4
    +8.*CB_FULL[i_dom[i][j]][2]*k4*k2*ks
    +7.*CB_FULL[i_dom[i][j]][3]*k4*k2
    +6.*CB_FULL[i_dom[i][j]][4]*k4*ks
    +5.*CB_FULL[i_dom[i][j]][5]*k4
    +4.*CB_FULL[i_dom[i][j]][6]*k2*ks
    +3.*CB_FULL[i_dom[i][j]][7]*k2
    +2.*CB_FULL[i_dom[i][j]][8]*ks
    +CB_FULL[i_dom[i][j]][9];
  ks*=1.e+12*2.*PI;
  d*=1.e-12*0.5/PI;
  dx=QH*d*tau*(KX+0.5*dkx)/ks;
  dy=QH*d*tau*(KY+0.5*dky)/ks;
  KX+=dkx;
  KY+=dky;
  X+=dx;
  Y+=dy;
 }

// check if some particles are out of the device
 i=(int)(X/dx)+1;
 j=(int)(Y/dy)+1;
 if(i<=1) i=1;
 if(j<=1) j=1;
 if(i>=nx) i=nx;
 if(j>=ny) j=ny;

// Generic boundary conditions for the super-particles
// ===================================================

// left edge
// =========
// ---Insulator---
 if(X<=0. && EDGE[3][j][0]==0){
  X=-X;
  KX=-KX;
  return;
 }
// ---Schottky or ohmic contact---
 else if(X<=0. && (EDGE[3][j][0]==1 || EDGE[3][j][0]==2)){
   IV=9;
   return;
 }
// right edge
// ==========
// ---Insulator---
 if(X>=LX && EDGE[1][j][0]==0){
  X=LX-(X-LX);
  KX=-KX;
  return;
 }
// ---Schottky or ohmic contact---
 else if(X>=LX && (EDGE[1][j][0]==1 || EDGE[1][j][0]==2)){
   IV=9;
   return;
 }
// bottom edge
// ===========
// ---Insulator---
 if(Y<=0. && EDGE[0][i][0]==0){
  Y=-Y;
  KY=-KY;
  return;
 }
// ---Schottky or ohmic contact---
 else if(Y<=0. && (EDGE[0][i][0]==1 || EDGE[0][i][0]==2)){
   IV=9;
   return;
 }
// upper edge
// ==========
// ---Insulator---
 if(Y>=LY && EDGE[2][i][0]==0){
  Y=LY-(Y-LY);
  KY=-KY;
  return;
 }
// ---Schottky or ohmic contact---
 else if(Y>=LY && (EDGE[2][i][0]==1 || EDGE[2][i][0]==2)){
   IV=9;
   return;
 }
}

// ============================================================
