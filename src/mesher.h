/* mesher.h -- This file is part of GNU archimedes.
   Archimedes is a simulator for Submicron and Nanoscaled
   2D III-V Semiconductor Devices.

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
// Created on 05 sep.2011, Carry le Rouet, France J.M.Sellier
// Last modif. : 05 sep.2011, Carry le Rouet, J.M.Sellier
// ######################################################

// Construction of the mesh for finite element method
// applied to Poisson eauqtion. This is done to properly
// take into account the oxyde layers and the dynamic of
// the electrostatic potential.

void mesher(void){
 int i,j;
 FILE *fp;

 // number of mesh nodes
 NG=(nx+1)*(ny+1);
 // definition of coordinates
 for(i=1;i<=(nx+1);i++)
  for(j=1;j<=(ny+1);j++){
   COORD[0][(j-1)*(nx+1)+i]=(i-1)*dx; // X-ccordinate
   COORD[1][(j-1)*(nx+1)+i]=(j-1)*dy; // Y-coordinate
  }
 // number of triangles
 NE=2*nx*ny;
 for(i=1;i<=nx;i++)
  for(j=1;j<=ny;j++){
   // bottom up triangle
   NODE_GEO[0][(j-1)*nx+i]=(j-1)*(nx+1)+i;
   NODE_GEO[1][(j-1)*nx+i]=(j-1)*(nx+1)+i+1;
   NODE_GEO[2][(j-1)*nx+i]=j*(nx+1)+i+1;
   // top down triangle
   NODE_GEO[0][(j-1)*nx+i+nx*ny]=(j-1)*(nx+1)+i;
   NODE_GEO[1][(j-1)*nx+i+nx*ny]=j*(nx+1)+i+1;
   NODE_GEO[2][(j-1)*nx+i+nx*ny]=j*(nx+1)+i;
  }
 // saves the mesh in mesh format (INRIA standard)
 if(SAVE_MESH==ON){
  fp=fopen("device.mesh","w");
  if(fp!=NULL){
   fprintf(fp,"MeshVersionFormatted 1\n\n");
   fprintf(fp,"Dimension\n2\n\n");
   fprintf(fp,"Vertices\n");
   fprintf(fp,"%d\n",NG);
   for(i=1;i<=NG;i++){
    fprintf(fp,"%g %g 0\n",COORD[0][i],COORD[1][i]);
   }
   fprintf(fp,"\n");
   fprintf(fp,"Triangles\n");
   fprintf(fp,"%d\n",NE);
   for(i=1;i<=NE;i++){
    fprintf(fp,"%d %d %d 0\n",NODE_GEO[0][i],NODE_GEO[1][i],NODE_GEO[2][i]);
   }
  } else {
   printf("an error occured while trying to open device.mesh!\n");
   exit(0);
  }
  fclose(fp);
 }
}
