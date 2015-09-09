/* saveoutput2dgnuplot.h -- This file is part of Archimedes release 0.1.0.
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
// Created on 10 Sep.2004, Siracusa (Italy), J.M.Sellier
// Last modif. : 10 Sep.2007, Siracusa (Italy), J.M.Sellier
// ######################################################

// Here we save all the macroscopic variables 
// related to the electrons.
// The format of saving is compatible with GNUPLOT
// and is as follows (constituted of 3 columns):
// x1 y1 u(x1,y1)
// x1 y2 u(x1,y2)
// ...
// xn yn u(xn,yn)
// Actually this file format is completly compatible
// with the "xyz" format of xd3d package. It is higly
// recomended that you use this program in order to
// obtain a very good scientifical visualisation of
// the solutions. In this way the post-processing become
// very easy and clear, so let use this program! It
// is free and you download it from :
// www.gnu.org ---> Free Directory ---> Science ---> Scientific visualisation

// INPUT : je (integer). It is the integer which will be the index of
// the solution. In this way you can easily create a succession of
// solutions like this:
// density000.xyz, density001.xyz, density002.xyz, ..., density010.xyz
// and so on, in order to create movies or animated gif.

void
SaveOutput2DGNUPLOT(int je)
{
 char s[150];
 FILE *fp;
 FILE *up;
 FILE *vp;
 FILE *lp;
 FILE *lxp;
 FILE *lyp;
 FILE *ep;
 FILE *qp;
 FILE *fM;
 register int i,j;

//if(Material==SILICON || Material==GERMANIUM){
 if(je<=9){
   sprintf(s,"density00%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"x_velocity00%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity00%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"energy00%d.xyz",je);
   ep=fopen(s,"w");
   sprintf(s,"potential00%d.xyz",je);
   lp=fopen(s,"w");
   sprintf(s,"magnetic_field00%d.xyz",je);
   fM=fopen(s,"w");
   sprintf(s,"x_Efield00%d.xyz",je);
   lxp=fopen(s,"w");
   sprintf(s,"y_Efield00%d.xyz",je);
   lyp=fopen(s,"w");
   sprintf(s,"quantum_potential00%d.xyz",je);
   qp=fopen(s,"w");
 }
 else if(je>9 && je<=99){
   sprintf(s,"density0%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"x_velocity0%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity0%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"energy0%d.xyz",je);
   ep=fopen(s,"w");
   sprintf(s,"potential0%d.xyz",je);
   lp=fopen(s,"w");
   sprintf(s,"magnetic_field0%d.xyz",je);
   fM=fopen(s,"w");
   sprintf(s,"x_Efield0%d.xyz",je);
   lxp=fopen(s,"w");
   sprintf(s,"y_Efield0%d.xyz",je);
   lyp=fopen(s,"w");
   sprintf(s,"quantum_potential0%d.xyz",je);
   qp=fopen(s,"w");
 }
 else{
   sprintf(s,"density%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"x_velocity%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"energy%d.xyz",je);
   ep=fopen(s,"w");
   sprintf(s,"potential%d.xyz",je);
   lp=fopen(s,"w");
   sprintf(s,"magnetic_field%d.xyz",je);
   fM=fopen(s,"w");
   sprintf(s,"x_Efield%d.xyz",je);
   lxp=fopen(s,"w");
   sprintf(s,"y_Efield%d.xyz",je);
   lyp=fopen(s,"w");
   sprintf(s,"quantum_potential%d.xyz",je);
   qp=fopen(s,"w");
 }
// Save the Monte Carlo results
 if(Model_Number==MCE || Model_Number==MCEH){
// Electron Density Output
// =======================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(fp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,u2d[i][j][1]);
      fprintf(fp,"\n");
  }
// X-component of electronic velocity output
// =========================================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(up,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i][j][2]/MEDIA);
      fprintf(up,"\n");
  }
// Y-component of electronic velocity output
// =========================================
  for (j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(vp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i][j][3]/MEDIA);
      fprintf(vp,"\n");
  }
// Electrostatic Potential
// =======================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,PSI[i][j]);
      fprintf(lp,"\n");
  }
// Magnetic Field
// ==============
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(fM,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,B[i][j]);
      fprintf(fM,"\n");
  }
// X-component of electric field
// =============================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lxp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][0]);
      fprintf(lxp,"\n");
  }
// Y-component of electric field
// =============================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lyp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][1]);
      fprintf(lyp,"\n");
  }
// Electron Energy (in eV)
// ===============
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(ep,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i][j][4]/MEDIA);
      fprintf(ep,"\n");
  }
// Quantum Effective Potential
// ===========================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(qp,"%g %g %g\n",
          1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,u2d[i][j][0]);
      fprintf(qp,"\n");
  }

// Closure of output files
// =======================
   fclose(fp);
   fclose(up);
   fclose(vp);
   fclose(lp);
   fclose(fM);
   fclose(lxp);
   fclose(lyp);
   fclose(ep);
   fclose(qp);
  }
// Save the Hybrid MEP results
 if(Model_Number==MEPE || Model_Number==MEPEH){
// Electron Density Output
// =======================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(fp,"%g %g %g\n",1.e6*(i-1.)*dx,
              1.e6*(j-1.)*dy,u2d[i+2][j+2][1]);
      fprintf(fp,"\n");
  }
// X-component of electronic velocity output
// =========================================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(up,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i+2][j+2][2]/u2d[i+2][j+2][1]);
      fprintf(up,"\n");
  }
// Y-component of electronic velocity output
// =========================================
  for (j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(vp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i+2][j+2][3]/u2d[i+2][j+2][1]);
      fprintf(vp,"\n");
  }
// Electrostatic Potential
// =======================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,PSI[i][j]);
      fprintf(lp,"\n");
  }
// X-component of electric field
// =============================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lxp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][0]);
      fprintf(lxp,"\n");
  }
// Y-component of electric field
// =============================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(lyp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][1]);
      fprintf(lyp,"\n");
  }
// Electron Energy (in eV)
// ===============
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(ep,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
              u2d[i+2][j+2][4]/u2d[i+2][j+2][1]/Q);
      fprintf(ep,"\n");
  }
// Quantum Effective Potential
// ===========================
  for(j=1;j<=ny+1;j++){
    for(i=1;i<=nx+1;i++)
      fprintf(qp,"%g %g %g\n",
          1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,u2d[i][j][0]);
      fprintf(qp,"\n");
  }

// Closure of output files
// =======================
   fclose(fp);
   fclose(up);
   fclose(vp);
   fclose(lp);
   fclose(lxp);
   fclose(lyp);
   fclose(ep);
   fclose(qp);
  }
// }
// ****************************************************************
/*if(Material==GAAS){
 if(je<=9){
   sprintf(s,"density00%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"densityGAMMA_valley00%d.xyz",je);
   fG=fopen(s,"w");
   sprintf(s,"densityL_valley00%d.xyz",je);
   fL=fopen(s,"w");
   sprintf(s,"x_velocity00%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity00%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"energy00%d.xyz",je);
   ep=fopen(s,"w");
   sprintf(s,"potential00%d.xyz",je);
   lp=fopen(s,"w");
   sprintf(s,"x_Efield00%d.xyz",je);
   lxp=fopen(s,"w");
   sprintf(s,"y_Efield00%d.xyz",je);
   lyp=fopen(s,"w");
   sprintf(s,"quantum_potential00%d.xyz",je);
   qp=fopen(s,"w");
 }
 else if(je>9 && je<=99){
   sprintf(s,"density0%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"densityGAMMA_valley0%d.xyz",je);
   fG=fopen(s,"w");
   sprintf(s,"densityL_valley0%d.xyz",je);
   fL=fopen(s,"w");
   sprintf(s,"x_velocity0%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity0%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"energy0%d.xyz",je);
   ep=fopen(s,"w");
   sprintf(s,"potential0%d.xyz",je);
   lp=fopen(s,"w");
   sprintf(s,"x_Efield0%d.xyz",je);
   lxp=fopen(s,"w");
   sprintf(s,"y_Efield0%d.xyz",je);
   lyp=fopen(s,"w");
   sprintf(s,"quantum_potential0%d.xyz",je);
   qp=fopen(s,"w");
 }
 else{
   sprintf(s,"density%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"densityGAMMA_valley%d.xyz",je);
   fG=fopen(s,"w");
   sprintf(s,"densityL_valley%d.xyz",je);
   fL=fopen(s,"w");
   sprintf(s,"x_velocity%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"energy%d.xyz",je);
   ep=fopen(s,"w");
   sprintf(s,"potential%d.xyz",je);
   lp=fopen(s,"w");
   sprintf(s,"x_Efield%d.xyz",je);
   lxp=fopen(s,"w");
   sprintf(s,"y_Efield%d.xyz",je);
   lyp=fopen(s,"w");
   sprintf(s,"quantum_potential%d.xyz",je);
   qp=fopen(s,"w");
 }

// Electron Density Output
// =======================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(fp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,u2d[i][j][1]);
     fprintf(fp,"\n");
 }
// Electron Density GAMMA-valley Output
// ====================================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(fG,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,DG[i][j]);
     fprintf(fG,"\n");
 }
// Electron Density L-valley Output
// ================================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(fL,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,DL[i][j]);
     fprintf(fL,"\n");
 }
// X-component of electronic velocity output
// =========================================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(up,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
             u2d[i][j][2]/MEDIA);
     fprintf(up,"\n");
 }
// Y-component of electronic velocity output
// =========================================
 for (j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(vp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
             u2d[i][j][3]/MEDIA);
     fprintf(vp,"\n");
 }
// Electrostatic Potential
// =======================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(lp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,PSI[i][j]);
     fprintf(lp,"\n");
 }
// X-component of electric field
// =============================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(lxp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][0]);
     fprintf(lxp,"\n");
 }
// Y-component of electric field
// =============================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(lyp,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,E[i][j][1]);
     fprintf(lyp,"\n");
 }
// Electron Energy (in eV)
// ===============
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(ep,"%g %g %g\n",1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,
             u2d[i][j][4]/MEDIA);
     fprintf(ep,"\n");
 }
// Quantum Effective Potential
// ===========================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(qp,"%g %g %g\n",
         1.e6*(i-1.)*dx,1.e6*(j-1.)*dy,u2d[i][j][0]);
     fprintf(qp,"\n");
 }

// Closure of output files
// =======================
 fclose(fp);
 fclose(up);
 fclose(vp);
 fclose(lp);
 fclose(lxp);
 fclose(lyp);
 fclose(ep);
 fclose(qp);
 fclose(fG);
 fclose(fL);
 }*/
}

// ============================================================
