/* saveoutput2dmeshformat.h -- This file is part of Archimedes release 0.1.0.
   Archimedes is a simulator for Submicron 2D Silicon/GaAs/Germanium
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
// Created on 24 Apr.2004, Siracusa, J.M.Sellier
// Last modif. : 10 Sep.2007, Siracusa, J.M.Sellier
// ######################################################

// Here we save all the macroscopic variables 
// related to the electrons.
// The output for such a variable is constituted by
// two files, one denominated "variable.BB"
// in which we find the solution on the vertex
// of the mesh, and one denominated "variable.mesh"
// which describes the topology and geometry
// of the mesh. This is a very diffused format
// in the scientifical computing comunity (specially at INRIA).
// For more information about this format see 
// the manual of GNU archimedes release 0.0.1

// INPUT : je (integer). It is the integer which will be the index of
// the solution. In this way you can easily create a succession of
// solutions like this:
// density000.BB, density001.BB, density002.BB, ..., density010.BB
// mesh000.mesh, mesh001.mesh, mesh002.mesh, ..., mesh010.mesh
// and so on, in order to create movies or animated gif.

void
SaveOutput2D_MeshFormat(int je)
{
  char s[50];
  FILE *fp;
  FILE *up;
  FILE *vp;
  FILE *ep;
  FILE *pop;
  FILE *mp;
  FILE *Xp;
  FILE *Yp;
  FILE *Qp;
  FILE *fM;
  int i,j;
  int size,dim,nbmet,nbval,type;

  printf("SaveOutput2D_MeshFormat\n\n");

 NUM_VERT=(nx+1)*(ny+1);
 printf("\nNUM_VERT=%d\n",NUM_VERT);
 NUM_EXAHEDRA=nx*ny;
 printf("\nNUM_EXAHEDRA=%d\n",NUM_EXAHEDRA);

//if(Material==SILICON || Material==GERMANIUM){
// Open the files related to electrons
// ===================================
if(je<10){
   sprintf(s,"density00%d.BB",je);
   fp=fopen(s,"w");
   sprintf(s,"x_velocity00%d.BB",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity00%d.BB",je);
   vp=fopen(s,"w");
   sprintf(s,"potential00%d.BB",je);
   pop=fopen(s,"w");
   sprintf(s,"magnetic_field00%d.BB",je);
   fM=fopen(s,"w");
   sprintf(s,"energy00%d.BB",je);
   ep=fopen(s,"w");
   sprintf(s,"mesh00%d.mesh",je);
   mp=fopen(s,"w");
   sprintf(s,"x_Efield00%d.BB",je);
   Xp=fopen(s,"w");
   sprintf(s,"y_Efield00%d.BB",je);
   Yp=fopen(s,"w");
   sprintf(s,"quant_poten00%d.BB",je);
   Qp=fopen(s,"w");
}
else if(je>=10 && je<100){
   sprintf(s,"density0%d.BB",je);
   fp=fopen(s,"w");
   sprintf(s,"x_velocity0%d.BB",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity0%d.BB",je);
   vp=fopen(s,"w");
   sprintf(s,"potential0%d.BB",je);
   pop=fopen(s,"w");
   sprintf(s,"magnetic_field0%d.BB",je);
   fM=fopen(s,"w");
   sprintf(s,"energy0%d.BB",je);
   ep=fopen(s,"w");
   sprintf(s,"mesh0%d.mesh",je);
   mp=fopen(s,"w");
   sprintf(s,"x_Efield0%d.BB",je);
   Xp=fopen(s,"w");
   sprintf(s,"y_Efield0%d.BB",je);
   Yp=fopen(s,"w");
   sprintf(s,"quant_poten0%d.BB",je);
   Qp=fopen(s,"w");
}
else{
   sprintf(s,"density%d.BB",je);
   fp=fopen(s,"w");
   sprintf(s,"x_velocity%d.BB",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity%d.BB",je);
   vp=fopen(s,"w");
   sprintf(s,"potential%d.BB",je);
   pop=fopen(s,"w");
   sprintf(s,"magnetic_field%d.BB",je);
   fM=fopen(s,"w");
   sprintf(s,"energy%d.BB",je);
   ep=fopen(s,"w");
   sprintf(s,"mesh%d.mesh",je);
   mp=fopen(s,"w");
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

 fprintf(fp,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(up,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(vp,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(pop,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(fM,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(ep,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(Xp,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(Yp,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(Qp,"%d %d %d %d\n",dim,nbmet,nbval,type);

// Save the Monte Carlo results 
 if(Model_Number==MCE || Model_Number==MCEH){
// Save the macroscopic variables related to electrons
// ===================================================
   for(j=1;j<=ny+1;j++)
     for(i=1;i<=nx+1;i++)
       fprintf(Qp,"%g\n",u2d[i][j][0]);
   for(j=1;j<=ny+1;j++)
     for(i=1;i<=nx+1;i++)
       fprintf(pop,"%g\n",PSI[i][j]);
   for(j=1;j<=ny+1;j++)
     for(i=1;i<=nx+1;i++)
       fprintf(fM,"%g\n",B[i][j]);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(fp,"%g\n",u2d[i][j][1]);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(up,"%g\n",u2d[i][j][2]/MEDIA);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(vp,"%g\n",u2d[i][j][3]/MEDIA);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(ep,"%g\n",u2d[i][j][4]/MEDIA); //energy (eV)
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(Xp,"%g\n",E[i][j][0]);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(Yp,"%g\n",E[i][j][1]);

// Save the file "variable.mesh"
// =============================
 fprintf(mp,"MeshVersionFormatted 1\n");
 fprintf(mp,"Dimension %d\n",dim);
 fprintf(mp,"Vertices\n%d\n",NUM_VERT);
 for(j=1;j<=ny+1;j++)
   for(i=1;i<=nx+1;i++)
       fprintf(mp,"%g %g 0\n",(i-1.)*dx,(j-1.)*dy);
 fprintf(mp,"Quadrilaterals\n%d\n",6*NUM_EXAHEDRA);
 for(j=1;j<=ny;j++)
  for(i=1;i<=nx;i++)
     fprintf(mp,"%d %d %d %d 0\n",
         i+(j-1)*(nx+1),i+1+(j-1)*(nx+1),
         i+1+j*(nx+1),i+j*(nx+1)); // 1 2 4 3
 }
// Save the Hybrid MEP results 
 if(Model_Number==MEPE || Model_Number==MEPEH){
// Save the macroscopic variables related to electrons
// ===================================================
   for(j=1;j<=ny+1;j++)
     for(i=1;i<=nx+1;i++)
       fprintf(Qp,"%g\n",u2d[i][j][0]);
   for(j=1;j<=ny+1;j++)
     for(i=1;i<=nx+1;i++)
       fprintf(pop,"%g\n",PSI[i][j]);
   for(j=1;j<=ny+1;j++)
     for(i=1;i<=nx+1;i++)
       fprintf(fM,"%g\n",B[i][j]);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(fp,"%g\n",u2d[i+2][j+2][1]);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(up,"%g\n",u2d[i+2][j+2][2]/u2d[i+2][j+2][1]);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(vp,"%g\n",u2d[i+2][j+2][3]/u2d[i+2][j+2][1]);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(ep,"%g\n",u2d[i+2][j+2][4]/u2d[i+2][j+2][1]/Q); //energy (eV)
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(Xp,"%g\n",E[i][j][0]);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(Yp,"%g\n",E[i][j][1]);

// Save the file "variable.mesh"
// =============================
 fprintf(mp,"MeshVersionFormatted 1\n");
 fprintf(mp,"Dimension %d\n",dim);
 fprintf(mp,"Vertices\n%d\n",NUM_VERT);
 for(j=1;j<=ny+1;j++)
   for(i=1;i<=nx+1;i++)
       fprintf(mp,"%g %g 0\n",(i-1.)*dx,(j-1.)*dy);
 fprintf(mp,"Quadrilaterals\n%d\n",6*NUM_EXAHEDRA);
 for(j=1;j<=ny;j++)
  for(i=1;i<=nx;i++)
     fprintf(mp,"%d %d %d %d 0\n",
         i+(j-1)*(nx+1),i+1+(j-1)*(nx+1),
         i+1+j*(nx+1),i+j*(nx+1)); // 1 2 4 3
 }
// Closure of files
// ================
    fclose(fp);
    fclose(up);
    fclose(vp);
    fclose(pop);
    fclose(fM);
    fclose(ep);
    fclose(mp);
    fclose(Xp);
    fclose(Yp);
    fclose(Qp);
//}
// **************************************************
/*if(Material==GAAS){
// Open the files related to electrons
// ===================================
if(je<10){
   sprintf(s,"density00%d.BB",je);
   fp=fopen(s,"w");
   sprintf(s,"densityGAMMA_valley00%d.BB",je);
   fG=fopen(s,"w");
   sprintf(s,"densityL_valley00%d.BB",je);
   fL=fopen(s,"w");
   sprintf(s,"x_velocity00%d.BB",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity00%d.BB",je);
   vp=fopen(s,"w");
   sprintf(s,"potential00%d.BB",je);
   pop=fopen(s,"w");
   sprintf(s,"energy00%d.BB",je);
   ep=fopen(s,"w");
   sprintf(s,"mesh00%d.mesh",je);
   mp=fopen(s,"w");
   sprintf(s,"x_Efield00%d.BB",je);
   Xp=fopen(s,"w");
   sprintf(s,"y_Efield00%d.BB",je);
   Yp=fopen(s,"w");
   sprintf(s,"quant_poten00%d.BB",je);
   Qp=fopen(s,"w");
}
else if(je>=10 && je<100){
   sprintf(s,"density0%d.BB",je);
   fp=fopen(s,"w");
   sprintf(s,"densityGAMMA_valley0%d.BB",je);
   fG=fopen(s,"w");
   sprintf(s,"densityL_valley0%d.BB",je);
   fL=fopen(s,"w");
   sprintf(s,"x_velocity0%d.BB",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity0%d.BB",je);
   vp=fopen(s,"w");
   sprintf(s,"potential0%d.BB",je);
   pop=fopen(s,"w");
   sprintf(s,"energy0%d.BB",je);
   ep=fopen(s,"w");
   sprintf(s,"mesh0%d.mesh",je);
   mp=fopen(s,"w");
   sprintf(s,"x_Efield0%d.BB",je);
   Xp=fopen(s,"w");
   sprintf(s,"y_Efield0%d.BB",je);
   Yp=fopen(s,"w");
   sprintf(s,"quant_poten0%d.BB",je);
   Qp=fopen(s,"w");
}
else{
   sprintf(s,"density%d.BB",je);
   fp=fopen(s,"w");
   sprintf(s,"densityGAMMA_valley%d.BB",je);
   fG=fopen(s,"w");
   sprintf(s,"densityL_valley%d.BB",je);
   fL=fopen(s,"w");
   sprintf(s,"x_velocity%d.BB",je);
   up=fopen(s,"w");
   sprintf(s,"y_velocity%d.BB",je);
   vp=fopen(s,"w");
   sprintf(s,"potential%d.BB",je);
   pop=fopen(s,"w");
   sprintf(s,"energy%d.BB",je);
   ep=fopen(s,"w");
   sprintf(s,"mesh%d.mesh",je);
   mp=fopen(s,"w");
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

 fprintf(fp,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(fG,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(fL,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(up,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(vp,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(pop,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(ep,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(Xp,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(Yp,"%d %d %d %d\n",dim,nbmet,nbval,type);
 fprintf(Qp,"%d %d %d %d\n",dim,nbmet,nbval,type);
 
// Save the macroscopic variables related to electrons
// ===================================================
   for(j=1;j<=ny+1;j++)
     for(i=1;i<=nx+1;i++)
       fprintf(Qp,"%g\n",u2d[i][j][0]);
   for(j=1;j<=ny+1;j++)
     for(i=1;i<=nx+1;i++)
       fprintf(pop,"%g\n",PSI[i][j]);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(fp,"%g\n",u2d[i][j][1]);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(fG,"%g\n",DG[i][j]);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(fL,"%g\n",DL[i][j]);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(up,"%g\n",u2d[i][j][2]/MEDIA);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(vp,"%g\n",u2d[i][j][3]/MEDIA);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(ep,"%g\n",u2d[i][j][4]/MEDIA); //energy (eV)
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(Xp,"%g\n",E[i][j][0]);
   for(j=1;j<=ny+1;j++)
    for(i=1;i<=nx+1;i++)
       fprintf(Yp,"%g\n",E[i][j][1]);

// Save the file "variable.mesh"
// =============================
 fprintf(mp,"MeshVersionFormatted 1\n");
 fprintf(mp,"Dimension %d\n",dim);
 fprintf(mp,"Vertices\n%d\n",NUM_VERT);
 for(j=1;j<=ny+1;j++)
   for(i=1;i<=nx+1;i++)
       fprintf(mp,"%g %g 0\n",(i-1.)*dx,(j-1.)*dy);
 fprintf(mp,"Quadrilaterals\n%d\n",6*NUM_EXAHEDRA);
 for(j=1;j<=ny;j++)
  for(i=1;i<=nx;i++)
     fprintf(mp,"%d %d %d %d 0\n",
         i+(j-1)*(nx+1),i+1+(j-1)*(nx+1),
         i+1+j*(nx+1),i+j*(nx+1)); // 1 2 4 3

// Closure of files
// ================
    fclose(fp);
    fclose(fG);
    fclose(fL);
    fclose(up);
    fclose(vp);
    fclose(pop);
    fclose(ep);
    fclose(mp);
    fclose(Xp);
    fclose(Yp);
    fclose(Qp);
}        */
 printf("End of SaveOutput2D_MeshFormat\n");
}

// =============================================================

