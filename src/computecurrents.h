/* computecurrents.h -- This file is part of Archimedes release 0.0.3.
   Archimedes is a simulator for Submicron 2D Silicon/GaAs
   Devices. It implements both the Monte Carlo method and Hybrid MEP model
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It also includes the quantum effects by means 
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
// Created on 18 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 05 mar.2005, Siracusa, J.M.Sellier
// ######################################################

// Computes and shows the currents on the contacts

void
Compute_Currents(void)
{
 register int i;
 int cn;
 real sum;

// Electron current computations
// =============================
 sum=0.;
 cn=1;
// compute currents on the contacts of the bottom edge
 for(i=1;i<=nx;i++){
   if(EDGE[0][i][0]==1 || EDGE[0][i][0]==2){
    if(Model_Number==MCE || Model_Number==MCEH)
     sum+=u2d[i][2][1]*u2d[i][2][3];
    else if(Model_Number==MEPE || Model_Number==MEPEH)
     sum+=u2d[i][3][3];
   }
   if((EDGE[0][i+1][0]==0 && sum!=0.)
    || (i==nx && sum!=0.)){
    if(Model_Number==MCE || Model_Number==MCEH)
     sum*=-Q*dx/MEDIA;
    else if(Model_Number==MEPE || Model_Number==MEPEH)
     sum*=-Q*dx;
    printf("\
Bottom Edge : Electron Current on contact #%d = %g (A/m)\n",cn,sum);
    cn++;
    sum=0.;
   }
 }
 sum=0.;
 cn=1;
// compute currents on the contacts of the right edge
 for(i=1;i<=ny;i++){
   if(EDGE[1][i][0]==1 || EDGE[1][i][0]==2){
    if(Model_Number==MCE || Model_Number==MCEH)
     sum+=u2d[nx-1][i][1]*u2d[nx-1][i][2];
    else if(Model_Number==MEPE || Model_Number==MEPEH)
     sum+=u2d[nx+3][i][2];
   }
   if((EDGE[1][i+1][0]==0 && sum!=0.)
    || (i==ny && sum!=0.)){
    if(Model_Number==MCE || Model_Number==MCEH)
     sum*=-Q*dy/MEDIA;
    else if(Model_Number==MEPE || Model_Number==MEPEH)
     sum*=-Q*dy;
    printf("\
Right Edge : Electron Current on contact #%d = %g (A/m)\n",cn,sum);
    cn++;
    sum=0.;
   }
 }
 sum=0.;
 cn=1;
// compute currents on the contacts of the upper edge
 for(i=1;i<=nx;i++){
   if(EDGE[2][i][0]==1 || EDGE[2][i][0]==2){
    if(Model_Number==MCE || Model_Number==MCEH)
     sum+=u2d[i][ny-1][1]*u2d[i][ny-1][3];
    else if(Model_Number==MEPE || Model_Number==MEPEH)
     sum+=u2d[i][ny+3][3];
   }
   if((EDGE[2][i+1][0]==0 && sum!=0.)
    || (i==nx && sum!=0.)){
    if(Model_Number==MCE || Model_Number==MCEH)
     sum*=-Q*dx/MEDIA;
    else if(Model_Number==MEPE || Model_Number==MEPEH)
     sum*=-Q*dx;
    printf("\
Upper Edge : Electron Current on contact #%d = %g (A/m)\n",cn,sum);
    cn++;
    sum=0.;
   }
 }
 sum=0.;
 cn=1;
// compute currents on the contacts of the left edge
 for(i=1;i<=ny;i++){
   if(EDGE[3][i][0]==1 || EDGE[3][i][0]==2){
    if(Model_Number==MCE || Model_Number==MCEH)
     sum+=u2d[2][i][1]*u2d[2][i][2];
    else if(Model_Number==MEPE || Model_Number==MEPEH)
     sum+=u2d[3][i][2];
   }
   if((EDGE[3][i+1][0]==0 && sum!=0.)
    || (i==ny && sum!=0.)){
    if(Model_Number==MCE || Model_Number==MCEH)
     sum*=-Q*dy/MEDIA;
    else if(Model_Number==MEPE || Model_Number==MEPEH)
     sum*=-Q*dy;
    printf("\
Left Edge : Electron Current on contact #%d = %g (A/m)\n",cn,sum);
    cn++;
    sum=0.;
   }
 }
// Hole current computations
// =========================
 if(Model_Number==MCEH || Model_Number==MCH
  || Model_Number==MEPEH || Model_Number==MEPH){
  sum=0.;
  cn=1;
// compute currents on the contacts of the bottom edge
  for(i=1;i<=nx;i++){
    if(EDGE[0][i][0]==1 || EDGE[0][i][0]==2) sum+=h2d[i][3][3];
    if((EDGE[0][i+1][0]==0 && sum!=0.)
     || (i==nx && sum!=0.)){
      sum*=-Q*dx;
     printf("\
Bottom Edge : Hole Current on contact #%d = %g (A/m)\n",cn,sum);
     cn++;
     sum=0.;
    }
  }
  sum=0.;
  cn=1;
// compute currents on the contacts of the right edge
  for(i=1;i<=ny;i++){
    if(EDGE[1][i][0]==1 || EDGE[1][i][0]==2) sum+=h2d[nx+3][i][2];
    if((EDGE[1][i+1][0]==0 && sum!=0.)
     || (i==ny && sum!=0.)){
      sum*=-Q*dy;
     printf("\
Right Edge : Hole Current on contact #%d = %g (A/m)\n",cn,sum);
     cn++;
     sum=0.;
    }
  }
  sum=0.;
  cn=1;
// compute currents on the contacts of the upper edge
  for(i=1;i<=nx;i++){
    if(EDGE[2][i][0]==1 || EDGE[2][i][0]==2) sum+=h2d[i][ny+3][3];
    if((EDGE[2][i+1][0]==0 && sum!=0.)
     || (i==nx && sum!=0.)){
      sum*=-Q*dx;
     printf("\
Upper Edge : Hole Current on contact #%d = %g (A/m)\n",cn,sum);
     cn++;
     sum=0.;
    }
  }
  sum=0.;
  cn=1;
// compute currents on the contacts of the left edge
  for(i=1;i<=ny;i++){
    if(EDGE[3][i][0]==1 || EDGE[3][i][0]==2) sum+=h2d[3][i][2];
    if((EDGE[3][i+1][0]==0 && sum!=0.)
     || (i==ny && sum!=0.)){
      sum*=-Q*dy;
     printf("\
Left Edge : Hole Current on contact #%d = %g (A/m)\n",cn,sum);
     cn++;
     sum=0.;
    }
  }
 }
}

// =============================================
