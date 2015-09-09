/* holemep2d.c -- This file is part of GNU/Archimedes 0.0.3
   This code is a simulator for Submicron 2D Silicon/GaAs
   Devices. It implements the Monte Carlo method and Hybrid MEP model
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


// =================================================================
// File Name : ParabMEP2D.c
// Version   : release 0.0.3
// Date of Creation : 30 Mar.2004, Siracusa, Italy, Jean Michel Sellier
// Last Revision : 05 Mar.2005, Siracusa, Italy, Jean Michel Sellier.
// System of Measure : M.K.S.C. System
// =================================================================

// This part is for the simulation of the ballistic part
// of the Hybrid MEP model for holes. This is a parabolic model
// which is based on the Maximum Entropy Principle.
// For more informations and references see the manual.

// See the paper:
// "Parabolic Hydrodynamical Model for Bipolar Semiconductor Devices
//  and Low Field Hole Mobility", G.Mascali, V.Romano, J.M.Sellier

void
Hole_MEP2D(int nx,int ny,real dx,
           real dy,real cfl,real theta)
{
/*
************************************************************
* INPUT : nx, ny :   # of cells in x-, y-direction
*         dx, dy :   step size  in x-, y-direction
*         cfl : CFL #;
*         tf : final time to stop iteration; 
*         to : current time
*         h2d  : conservative variables at time to
*              entries of "h" needed : h[3:nx+3][3:ny+3][4]
* OUTPUT : h2d : conservative variable after 2-stage time iteration
* REMARK : Reset parameters "NXM","NYM" to adjust dimension of
*          arrays, Modify boundary conditions below when necessary.
* CAUTION : h2d[i][j][*] when io=1 is on the cell with solid-line
*           i.e. meaning h2d[i+1/2][j+1/2][*]
************************************************************
*/
 const int ND=2;

 register int m,i,j;
 int io,NXE,NYE;
 real den,xmt,ymt;
 real eng,vex,vey;
 real dtodx2,dtody2,dt;

 NXE=nx+ND;
 NYE=ny+ND;

// reset everything
 memset(bufx2d,0,sizeof(&bufx2d));
 memset(bufy2d,0,sizeof(&bufy2d));
 memset(ux2d,0,sizeof(&ux2d));
 memset(uy2d,0,sizeof(&uy2d));
 memset(f2d,0,sizeof(&f2d));
 memset(g2d,0,sizeof(&g2d));
 memset(fx2d,0,sizeof(&fx2d));
 memset(gy2d,0,sizeof(&gy2d));

// Start a 2-stage time iteration
 for(io=0;io<=1;io++){
   HoleHMEPBCs();
// Compute Numerical Slopes for u in x-, y-directions
// (denoted as "ux", "uy", resp.)
   for(m=1;m<=MN3;m++){
     for(j=ND-2;j<=NYE+2;j++)
       for(i=ND-2;i<=NXE+2;i++){
         bufx2d[i][j]=h2d[i+1][j][m]-h2d[i][j][m];
         bufy2d[i][j]=h2d[i][j+1][m]-h2d[i][j][m];
       }
    for(j=ND;j<=NYE+1;j++)
      for(i=ND;i<=NXE+1;i++){
        ux2d[i][j][m]=MM2(theta,bufx2d[i-1][j],bufx2d[i][j]);
        uy2d[i][j][m]=MM2(theta,bufy2d[i][j-1],bufy2d[i][j]);
      }
   }
// Compute the fluxes f and g
    for(j=ND-2;j<=NYE+3;j++)
      for(i=ND-2;i<=NXE+3;i++){
        den=h2d[i][j][1];
        xmt=h2d[i][j][2];
        ymt=h2d[i][j][3];
        eng=h2d[i][j][4];
        vex=xmt/den;
        vey=ymt/den;
        f2d[i][j][1]=xmt;
        f2d[i][j][2]=2./(3.*mstarhole*M)*eng;
        f2d[i][j][3]=0.;
        f2d[i][j][4]=4./3.*eng*vex;
        g2d[i][j][1]=ymt;
        g2d[i][j][2]=0.0;
        g2d[i][j][3]=2./(3.*mstarhole*M)*eng;
        g2d[i][j][4]=4./3.*eng*vey;
      }
// Compute numerical slopes for f, g in x-, y-direction
// (denoted as "fx", "gy", resp.)
   for(m=1;m<=MN3;m++){
     for(j=ND-2;j<=NYE+2;j++)
       for(i=ND-2;i<=NXE+2;i++){
         bufx2d[i][j]=f2d[i+1][j][m]-f2d[i][j][m];
         bufy2d[i][j]=g2d[i][j+1][m]-g2d[i][j][m];
       }
    for(j=ND;j<=NYE+1;j++)
      for(i=ND;i<=NXE+1;i++){
        fx2d[i][j][m]=MM2(theta,bufx2d[i-1][j],bufx2d[i][j]);
        gy2d[i][j][m]=MM2(theta,bufy2d[i][j-1],bufy2d[i][j]);
      }
   }
// Compute the flux values of f, g and h at the center of the four faces
   for(j=ND;j<=NYE+1;j++)
     for(i=ND;i<=NXE+1;i++){
       den=h2d[i][j][1]-dtodx2*fx2d[i][j][1]-dtody2*gy2d[i][j][1];
       xmt=h2d[i][j][2]-dtodx2*fx2d[i][j][2]-dtody2*gy2d[i][j][2];
       ymt=h2d[i][j][3]-dtodx2*fx2d[i][j][3]-dtody2*gy2d[i][j][3];
       eng=h2d[i][j][4]-dtodx2*fx2d[i][j][4]-dtody2*gy2d[i][j][4];
       vex=xmt/den;
       vey=ymt/den;
       f2d[i][j][1]=xmt;
       f2d[i][j][2]=2./(3.*mstarhole*M)*eng;
       f2d[i][j][3]=0.0;
       f2d[i][j][4]=4./3.*eng*vex;
       g2d[i][j][1]=ymt;
       g2d[i][j][2]=0.0;
       g2d[i][j][3]=2./(3.*mstarhole*M)*eng;
       g2d[i][j][4]=4./3.*eng*vey;
     }
// Compute time step size
   if(io==0){
     dt=DT;
     if((TEMPO+2.*dt)>=TF) dt=0.5*(TF-TEMPO);
   }
   dtodx2=0.5*dt/dx;
   dtody2=0.5*dt/dy;
// Compute the values of "u2d" at the next time level
   for(m=1;m<=MN3;m++){
     for(j=ND+1-io;j<=NYE-io;j++)
       for(i=ND+1-io;i<=NXE-io;i++)
         bufx2d[i][j]=0.25*(h2d[i][j][m]+h2d[i+1][j][m]
           +h2d[i][j+1][m]+h2d[i+1][j+1][m])
           +0.0625*(ux2d[i][j][m]-ux2d[i+1][j][m]
           +ux2d[i][j+1][m]-ux2d[i+1][j+1][m]
           +uy2d[i][j][m]+uy2d[i+1][j][m]
           -uy2d[i][j+1][m]-uy2d[i+1][j+1][m])
           +dtodx2*(f2d[i][j][m]-f2d[i+1][j][m]
           +f2d[i][j+1][m]-f2d[i+1][j+1][m])
           +dtody2*(g2d[i][j][m]+g2d[i+1][j][m]
           -g2d[i][j+1][m]-g2d[i+1][j+1][m]);
     for(j=ND+1;j<=NYE;j++)
       for(i=ND+1;i<=NXE;i++)
         h2d[i][j][m]=bufx2d[i-io][j-io];
   }
 }
// End of the 2-stage time iteration
}

// ######################################################

