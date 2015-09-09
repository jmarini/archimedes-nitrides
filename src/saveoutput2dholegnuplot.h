/* saveoutput2dholegnuplot.h --
   This file is part of Archimedes release 0.0.3.
   Archimedes is a simulator for Submicron 2D Silicon/GaAs
   Devices. It implements the Monte Carlo method and MEP model
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
// Created on 10 Sep.2004, Siracusa (Italy), J.M.Sellier
// Last modif. : 05 Mar.2005, Siracusa (Italy), J.M.Sellier
// ######################################################

// Here we save all the macroscopic variables 
// related to the holes.
// The format of saving is compatible with GNUPLOT
// and is as follows (constituted of 3 columns):
// x1 y1 h(x1,y1)
// x1 y2 h(x1,y2)
// ...
// xn yn h(xn,yn)

void
SaveOutput2DHoleGNUPLOT(int je)
{
 char s[50];
 FILE *fp;
 FILE *up;
 FILE *vp;
 FILE *lp;
 FILE *lxp;
 FILE *lyp;
 FILE *ep;
 FILE *qp;
 register int i,j;

 if(je<=9){
   sprintf(s,"hole_density00%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"hole_x_velocity00%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"hole_y_velocity00%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"hole_energy00%d.xyz",je);
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
   sprintf(s,"hole_density0%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"hole_x_velocity0%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"hole_y_velocity0%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"hole_energy0%d.xyz",je);
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
   sprintf(s,"hole_density%d.xyz",je);
   fp=fopen(s,"w");
   sprintf(s,"hole_x_velocity%d.xyz",je);
   up=fopen(s,"w");
   sprintf(s,"hole_y_velocity%d.xyz",je);
   vp=fopen(s,"w");
   sprintf(s,"hole_energy%d.xyz",je);
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
     fprintf(fp,"%g %g %g\n",
             1.e6*(i-0.5)*dx,1.e6*(j-0.5)*dy,h2d[i+2][j+2][1]);
     fprintf(fp,"\n");
 }
// X-component of electronic velocity output
// =========================================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(up,"%g %g %g\n",
          1.e6*(i-0.5)*dx,1.e6*(j-0.5)*dy,
          h2d[i+2][j+2][2]/h2d[i+2][j+2][1]);
     fprintf(up,"\n");
 }
// Y-component of electronic velocity output
// =========================================
 for (j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(vp,"%g %g %g\n",
          1.e6*(i-0.5)*dx,1.e6*(j-0.5)*dy,
          h2d[i+2][j+2][3]/h2d[i+2][j+2][1]);
     fprintf(vp,"\n");
 }
// Electrostatic Potential
// =======================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(lp,"%g %g %g\n",1.e6*(i-0.5)*dx,1.e6*(j-0.5)*dy,PSI[i][j]);
     fprintf(lp,"\n");
 }
// X-component of electric field
// =============================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(lxp,"%g %g %g\n",1.e6*(i-0.5)*dx,1.e6*(j-0.5)*dy,E[i][j][0]);
     fprintf(lxp,"\n");
 }
// Y-component of electric field
// =============================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(lyp,"%g %g %g\n",1.e6*(i-0.5)*dx,1.e6*(j-0.5)*dy,E[i][j][1]);
     fprintf(lyp,"\n");
 }
// Electron Energy (in eV)
// ===============
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(ep,"%g %g %g\n",
         1.e6*(i-0.5)*dx,1.e6*(j-0.5)*dy,
         h2d[i+2][j+2][4]/h2d[i+2][j+2][1]/Q);
     fprintf(ep,"\n");
 }
// Quantum Effective Potential
// ===========================
 for(j=1;j<=ny+1;j++){
   for(i=1;i<=nx+1;i++)
     fprintf(qp,"%g %g %g\n",
         1.e6*(i-0.5)*dx,1.e6*(j-0.5)*dy,u2d[i][j][0]);
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

// ============================================================
