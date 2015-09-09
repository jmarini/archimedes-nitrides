/* ParabMEP2D.c -- This file is part of GNU/Archimedes 0.1.0
   This code is a simulator for Submicron 2D III-V semiconductor
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


// =================================================================
// Date of Creation : 19 Mar.2004, Siracusa, Italy, Jean Michel Sellier
// Last Revision : 10 Sep.2007, Siracusa, Italy, Jean Michel Sellier.
// System of Measure : M.K.S.C. System
// =================================================================

// This part is for the simulation of the ballistic part
// of the Hybrid MEP model. This is a parabolic model
// which is based on the Maximum Entropy Principle.
// For more informations and references see the manual.

// Resolution of the CONVECTION STEP for the 2D MEP model
// (Maximum Entropy Principle), Parabolic case, for electrons
// by means of Jiang-Tadmor method.
// For more information see the following papers :
// For the numerical part:
// "Non-oscillatory Central Schemes for Multidimensional
//  Hyperbolic Conservation Laws", 
//  Guang-Shan Jiang, Eitan Tadmor, 30 Sep.1996.
// For the modeling part:
// "Two dimensional MESFET simulation of transients and
//  steady state with kinetic based hydrodynamical models",
//  A.M.Anile, S.F.Liotta, G.Mascali, S.Rinaudo.
// "Non parabolic band transport in semiconductors:
//  closure of the moment equations",
// A.M.Anile, V.Romano.
// "Non-parabolic band hydrodynamical model of silicon
//  semiconductors and simulation of electron devices",
//  V.Romano, Math.Meth.Appl.Sci. 2001, 24; 439-471.


void ParabMEP2D(int nx,int ny,real dx,real dy,real cfl,real theta)
{
/*
************************************************************
* INPUT : nx, ny :   # of cells in x-, y-direction
*         dx, dy :   step size  in x-, y-direction
*         cfl : CFL #;
*         tf : final time to stop iteration; 
*         to : current time
*         u2d  : conservative variables at time to
*              entries of "u" needed : u[3:nx+3][3:ny+3][4]
* OUTPUT : u2d : conservative variable after 2-stage time iteration
* REMARK : Reset parameters "NXM","NYM" to adjust dimension of
*          arrays, Modify boundary conditions below when necessary.
* CAUTION : u[i][j][*] when io=1 is on the cell with solid-line
*           i.e. meaning u[i+1/2][j+1/2][*]
************************************************************
*/
 const int ND=2;

 register int m,i,j;
 real den;
 real xmt;
 real ymt;
 real eng;
 real vex;
 real vey;
 real dt;
 real dtodx2;
 real dtody2;
 int io,NXE,NYE;

 NXE=nx+ND;
 NYE=ny+ND;

// Start a 2-stage time iteration
 for(io=0;io<=1;io++){
   HMEPBCs();
// Compute Numerical Slopes for u in x-, y-directions
// (denoted as "ux", "uy", resp.)
   for(m=1;m<=MN3;m++){
     for(j=ND-2;j<=NYE+2;j++)
       for(i=ND-2;i<=NXE+2;i++){
         bufx2d[i][j]=u2d[i+1][j][m]-u2d[i][j][m];
         bufy2d[i][j]=u2d[i][j+1][m]-u2d[i][j][m];
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
       den=u2d[i][j][1];
       xmt=u2d[i][j][2];
       ymt=u2d[i][j][3];
       eng=u2d[i][j][4];
       vex=xmt/den;
       vey=ymt/den;
       f2d[i][j][1]=xmt;
       f2d[i][j][2]=2./(3.*MSTAR[i_dom[i][j]][1]*M)*eng;
       f2d[i][j][3]=0.;
       f2d[i][j][4]=4./3.*eng*vex;
       g2d[i][j][1]=ymt;
       g2d[i][j][2]=0.0;
       g2d[i][j][3]=2./(3.*MSTAR[i_dom[i][j]][1]*M)*eng;
       g2d[i][j][4]=4./3.*eng*vey;
     }
// Compute numerical slopes for f, g and h in x-, y-direction
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
// Compute the flux values of f, g at the center of the four faces
   for(j=ND;j<=NYE+1;j++)
     for(i=ND;i<=NXE+1;i++){
       den=u2d[i][j][1]-dtodx2*fx2d[i][j][1]-dtody2*gy2d[i][j][1];
       xmt=u2d[i][j][2]-dtodx2*fx2d[i][j][2]-dtody2*gy2d[i][j][2];
       ymt=u2d[i][j][3]-dtodx2*fx2d[i][j][3]-dtody2*gy2d[i][j][3];
       eng=u2d[i][j][4]-dtodx2*fx2d[i][j][4]-dtody2*gy2d[i][j][4];
       vex=xmt/den;
       vey=ymt/den;
       f2d[i][j][1]=xmt;
       f2d[i][j][2]=2./(3.*MSTAR[i_dom[i][j]][1]*M)*eng;
       f2d[i][j][3]=0.0;
       f2d[i][j][4]=4./3.*eng*vex;
       g2d[i][j][1]=ymt;
       g2d[i][j][2]=0.0;
       g2d[i][j][3]=2./(3.*MSTAR[i_dom[i][j]][1]*M)*eng;
       g2d[i][j][4]=4./3.*eng*vey;
     }
// Compute time step size
   if(io==0){
     dt=DT;
     if((TEMPO+2.*dt)>=TF) dt=0.5*(TF-TEMPO);
   }
   dtodx2=0.5*dt/dx;
   dtody2=0.5*dt/dy;
   TEMPO+=dt;
// Compute the values of "u2d" at the next time level
   for(m=1;m<=MN3;m++){
     for(j=ND+1-io;j<=NYE-io;j++)
       for(i=ND+1-io;i<=NXE-io;i++){
         bufx2d[i][j]=0.25*(u2d[i][j][m]+u2d[i+1][j][m]
            +u2d[i][j+1][m]+u2d[i+1][j+1][m])
            +0.0625*(ux2d[i][j][m]-ux2d[i+1][j][m]
            +ux2d[i][j+1][m]-ux2d[i+1][j+1][m]
            +uy2d[i][j][m]+uy2d[i+1][j][m]
            -uy2d[i][j+1][m]-uy2d[i+1][j+1][m])
            +dtodx2*(f2d[i][j][m]-f2d[i+1][j][m]
            +f2d[i][j+1][m]-f2d[i+1][j+1][m])
            +dtody2*(g2d[i][j][m]+g2d[i+1][j][m]
            -g2d[i][j+1][m]-g2d[i+1][j+1][m]);
       }
    for(j=ND+1;j<=NYE;j++)
      for(i=ND+1;i<=NXE;i++)
        u2d[i][j][m]=bufx2d[i-io][j-io];
   }
 }
// End of the 2-stage time iteration
}

// ######################################################

