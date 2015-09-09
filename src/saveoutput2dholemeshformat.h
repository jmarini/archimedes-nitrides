/* saveoutput2dholemeshformat.h --
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
// Created on 30 Apr.2004, Siracusa, J.M.Sellier
// Last modif. : 05 Mar.2005, Siracusa, J.M.Sellier
// ######################################################

// Here we save all the macroscopic variables 
// related to the holes.
// The output for such a variable is constituted by
// two files, one denominated "variable.BB"
// in which we find the solution on the vertex
// of the mesh, and one denominated "variable.mesh"
// which describes the topology and geometry
// of the mesh. This is a very diffused format
// in the scientifical computing comunity.
// For more information about this format see 
// the manual of GNU archimedes release 0.0.3

void
SaveOutput2DHole_MeshFormat(int je)
{
  char s[50];
  FILE *fp;
  FILE *up;
  FILE *vp;
  FILE *ep;
  FILE *mp;
  FILE *pop;
  FILE *Xp;
  FILE *Yp;
  FILE *Qp;
  int i,j;
  int size,dim,nbmet,nbval,type;

  printf("SaveOutput3DHole_Mesh_Format\n\n");

// Definition of the number of vertex and hexaedra
// ===============================================
 NUM_VERT=nx*ny;
 printf("\nNUM_VERT=%d\n",NUM_VERT);
 NUM_EXAHEDRA=(nx-1)*(ny-1);
 printf("\nNUM_EXAHEDRA=%d\n",NUM_EXAHEDRA);
// =========================================


// Open the output files
// =====================
 if(je<10){
   sprintf(s,"hole_density00%d.BB",je);
   fp=fopen(s,"w");
   sprintf(s,"hole_x_velocity00%d.BB",je);
   up=fopen(s,"w");
   sprintf(s,"hole_y_velocity00%d.BB",je);
   vp=fopen(s,"w");
   sprintf(s,"hole_energy00%d.BB",je);
   ep=fopen(s,"w");
   sprintf(s,"mesh00%d.mesh",je);
   mp=fopen(s,"w");
   sprintf(s,"potential00%d.BB",je);
   pop=fopen(s,"w");
   sprintf(s,"x_Efield00%d.BB",je);
   Xp=fopen(s,"w");
   sprintf(s,"y_Efield00%d.BB",je);
   Yp=fopen(s,"w");
   sprintf(s,"quant_poten00%d.BB",je);
   Qp=fopen(s,"w");
 }
 else if(je>=10 && je<100){
   sprintf(s,"hole_density0%d.BB",je);
   fp=fopen(s,"w");
   sprintf(s,"hole_x_velocity0%d.BB",je);
   up=fopen(s,"w");
   sprintf(s,"hole_y_velocity0%d.BB",je);
   vp=fopen(s,"w");
   sprintf(s,"hole_energy0%d.BB",je);
   ep=fopen(s,"w");
   sprintf(s,"mesh0%d.mesh",je);
   mp=fopen(s,"w");
   sprintf(s,"potential0%d.BB",je);
   pop=fopen(s,"w");
   sprintf(s,"x_Efield0%d.BB",je);
   Xp=fopen(s,"w");
   sprintf(s,"y_Efield0%d.BB",je);
   Yp=fopen(s,"w");
   sprintf(s,"quant_poten0%d.BB",je);
   Qp=fopen(s,"w");
 }
 else{
   sprintf(s,"hole_density%d.BB",je);
   fp=fopen(s,"w");
   sprintf(s,"hole_x_velocity%d.BB",je);
   up=fopen(s,"w");
   sprintf(s,"hole_y_velocity%d.BB",je);
   vp=fopen(s,"w");
   sprintf(s,"hole_energy%d.BB",je);
   ep=fopen(s,"w");
   sprintf(s,"mesh%d.mesh",je);
   mp=fopen(s,"w");
   sprintf(s,"potential%d.BB",je);
   pop=fopen(s,"w");
   sprintf(s,"x_Efield%d.BB",je);
   Xp=fopen(s,"w");
   sprintf(s,"y_Efield%d.BB",je);
   Yp=fopen(s,"w");
   sprintf(s,"quant_poten%d.BB",je);
   Qp=fopen(s,"w");
 }

 printf("Files opened...\n");

 size=NUM_VERT;
 dim=2; // dimension of the space
 nbmet=1; // number of related fields (in this case a scalar value)
 nbval=NUM_VERT; // number of information attached to the vertex
 type=2; // solution in the vertices

// Save the geometry information of the mesh
// =========================================
 fprintf(fp,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(up,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(vp,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(ep,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(pop,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(Xp,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(Yp,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(Qp,"%d %d %d %d\n",dim,nbmet,nbval,type);
 
// Save the macroscopic variables
// ==============================
   for(j=1;j<=ny;j++)
    for(i=1;i<=nx;i++)
       fprintf(fp,"%g\n",h2d[i+2][j+2][1]);
   for(j=1;j<=ny;j++)
    for(i=1;i<=nx;i++)
       fprintf(up,"%g\n",h2d[i+2][j+2][2]/h2d[i+2][j+2][1]);
   for(j=1;j<=ny;j++)
    for(i=1;i<=nx;i++)
       fprintf(vp,"%g\n",h2d[i+2][j+2][3]/h2d[i+2][j+2][1]);
   for(j=1;j<=ny;j++)
    for(i=1;i<=nx;i++)
       fprintf(ep,"%g\n",h2d[i+2][j+2][4]/h2d[i+2][j+2][1]/Q); //energy (eV)
   for(j=1;j<=ny;j++)
     for(i=1;i<=nx;i++)
       fprintf(Qp,"%g\n",u2d[i][j][0]);
   for(j=1;j<=ny;j++)
     for(i=1;i<=nx;i++)
       fprintf(pop,"%g\n",PSI[i][j]);
    for(i=1;i<=nx;i++)
       fprintf(Xp,"%g\n",E[i][j][0]);
   for(j=1;j<=ny;j++)
    for(i=1;i<=nx;i++)
       fprintf(Yp,"%g\n",E[i][j][1]);

// Save the output file "variable.mesh"
// ====================================
 fprintf(mp,"MeshVersionFormatted 1\n");
 fprintf(mp,"Dimension %d\n",dim);
 fprintf(mp,"Vertices\n%d\n",NUM_VERT);
 for(j=1;j<=ny;j++)
   for(i=1;i<=nx;i++)
       fprintf(mp,"%g %g 0\n",(i-0.5)*dx,(j-0.5)*dy);
 fprintf(mp,"Quadrilaterals\n%d\n",6*NUM_EXAHEDRA);
 for(j=1;j<=(ny-1);j++)
  for(i=1;i<=(nx-1);i++)
     fprintf(mp,"%d %d %d %d 0\n",
         i+(j-1)*nx,i+1+(j-1)*nx,
         i+1+j*nx,i+j*nx); // 1 2 4 3
// =========================================

// Closure of files
    fclose(fp);
    fclose(up);
    fclose(vp);
    fclose(ep);
    fclose(mp);
    fclose(pop);
    fclose(Xp);
    fclose(Yp);
    fclose(Qp);
 printf("End of SaveOutput3DHole_Mesh_Format\n");
}

// =============================================================

