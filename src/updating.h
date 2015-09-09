/* updating.h -- This file is part of GNU archimedes.

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


// ###############################################################
// Created on 08 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 02 September 2011, Carry le Rouet, J.M.D. Sellier
// ###############################################################

// Updates the macroscopic variables
// Here we call all the subroutines we
// need for simulating the dynamics of the
// super-particles in the device.
// The solution is writted in the array named u2d.

void 
updating(int model)
{
 register int i,j;
 real maxi=0.,mini=0.;
// Computation of the electric field
// =================================
  Electric_Field();
  if(FARADAYFLAG) Faraday();
// Monte Carlo Simulation
// ======================
  if(model==MCE || model==MCEH){
   EMC();
   Charge();
   if((int)(TF/DT)-(int)(TEMPO/DT)<MEDIA){
     media();
   }
   if((TEMPO+DT)>=TF) DT=TF-TEMPO;
   TEMPO+=DT;
  }
// Electron MEP Simulation
// =======================
  if(model==MEPE || model==MEPEH){
   DT/=2.;
   ParabMEP2D(nx,ny,dx,dy,0.475,1.0);
   electron_relaxation_step();
   DT*=2.;
  }
// Hole MEP Simulation
// ===================
  if(model==MEPH || model==MEPEH){
   DT/=2.;
   Hole_MEP2D(nx,ny,dx,dy,0.475,1.0);
   Relaxation_Step_Hole();
   DT*=2.;
  }
// Output on some usefull informations about the simulation
  printf("%5d   TIME = %10.4g  (picosec)\n",c,TEMPO*1.e12);
  if(MAXIMINI==1){
// Compute the maximum and minimum of various macroscopic variables
// Max and Min of Potential
    maxi=u2d[8][8][0];
    mini=u2d[8][8][0];
    for(i=1;i<=nx+1;i++)
      for(j=1;j<=ny+1;j++){
       if(Model_Number==MCE || Model_Number==MCEH){
        if(u2d[i][j][0]>=maxi) maxi=u2d[i][j][0];
        if(u2d[i][j][0]<=mini) mini=u2d[i][j][0];
       }
       if(Model_Number==MEPE || Model_Number==MEPEH){
        if(u2d[i+2][j+2][0]>=maxi) maxi=u2d[i+2][j+2][0];
        if(u2d[i+2][j+2][0]<=mini) mini=u2d[i+2][j+2][0];
       }
      }
    printf("Max. Potential = %g\n",maxi);
    printf("Min. Potential = %g\n",mini);
// Max and Min of x-component of electric field
    maxi=E[8][8][0];
    mini=E[8][8][0];
    for(i=1;i<=nx+1;i++)
      for(j=1;j<=ny+1;j++){
        if(E[i][j][0]>=maxi) maxi=E[i][j][0];
        if(E[i][j][0]<=mini) mini=E[i][j][0];
      }
    printf("Max. x-elec.field = %g\n",maxi);
    printf("Min. x-elec.field = %g\n",mini);
// Max and Min of y-component of electric field
    maxi=E[8][8][1];
    mini=E[8][8][1];
    for(i=1;i<=nx+1;i++)
      for(j=1;j<=ny+1;j++){
        if(E[i][j][1]>=maxi) maxi=E[i][j][1];
        if(E[i][j][1]<=mini) mini=E[i][j][1];
      }
    printf("Max. y-elec.field = %g\n",maxi);
    printf("Min. y-elec.field = %g\n",mini);
// Max and Min of Density
    maxi=0.;
    mini=DDmax;
    for(i=1;i<=nx+1;i++)
      for(j=1;j<=ny+1;j++){
       if(Model_Number==MCE || Model_Number==MCEH){
        if(u2d[i][j][1]>=maxi) maxi=u2d[i][j][1];
        if(u2d[i][j][1]<=mini) mini=u2d[i][j][1];
       }
       if(Model_Number==MEPE || Model_Number==MEPEH){
        if(u2d[i+2][j+2][1]>=maxi) maxi=u2d[i+2][j+2][1];
        if(u2d[i+2][j+2][1]<=mini) mini=u2d[i+2][j+2][1];
       }
      }
    printf("Max. Density = %g\n",maxi);
    printf("Min. Density = %g\n",mini);
  }

// Here we save at each step if this option has been choosed
  if(SAVEALWAYS==1){
    SaveOutputFiles(File_Format,c);
    printf("Output number %d has been saved\n",c);
  }
  if(fabs(TEMPO-TF)/fabs(TF)<SMALL){
// Compute the various currents on the various defined contacts
    Compute_Currents();
    c=ITMAX+10;
  }
}

// ========================================
