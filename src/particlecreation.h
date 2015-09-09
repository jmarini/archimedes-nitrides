/* particlecreation.h -- This file is part of Archimedes release 1.2.0
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method
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
// Created on 06 sep.2004, Siracusa, Italy, J.M.Sellier
// Last modif. : 26 Aug.2011, Carry le Rouet, France, J.M.Sellier
// ######################################################

// create a new super-particle on the edge
// edge = 0 Bottom edge
// edge = 1 Right edge
// edge = 2 Upper edge
// edge = 3 Left edge

inline
void
creation(int i,real t,int edge)
{
 int iaux;
 int iv;
 int ii,j;
 real c1,c2,c3,
      c4,c5,c6,c7;

// We assume that the particles are initially 
// at near thermal equilibrium

// creation of the particle position vector r=(X,Y)
 X=dx*(rnd()+((real) i)-1.5);
 Y=dy*(rnd()+((real) i)-1.5);
 if((edge==0 || edge==2) && i==1) X=dx*0.5*rnd();
 if((edge==0 || edge==2) && i==nx) X=LX-dx*0.5*rnd();
 if(edge==0) Y=dy*0.5*rnd();
 if(edge==2) Y=LY-dy*0.5*rnd();
 if((edge==1 || edge==3) && i==1) Y=dy*0.5*rnd();
 if((edge==1 || edge==3) && i==ny) Y=LY-dy*0.5*rnd();
 if(edge==1) X=LX-dx*0.5*rnd();
 if(edge==3) X=dx*0.5*rnd();

// creation of the particle pseudo-wave vector k=(KX,KY,KZ)
// in the (i,j)-th cell
 ii=(int)(X/dx+1.5);
 j=(int)(Y/dy+1.5);
 if(ii<=1) ii=1;
 if(j<=1) j=1;
 if(ii>=nx+1) ii=nx+1;
 if(j>=ny+1) j=ny+1;
 if(NOVALLEY[i_dom[ii][j]]==1){
   iv=1;
   iaux=0;
 }
 if(NOVALLEY[i_dom[ii][j]]==2){
   iv=iaux=1;
// 20% of the created particles belongs to the L-valley
   if(rnd()>=0.8) iv=iaux=2;
 }
 c1=log(rnd());
 c2=SMH[i_dom[ii][j]][iaux]*sqrt(-1.5*BKTQ*c1*(1.-alphaK[i_dom[ii][j]][iv]*1.5*BKTQ*c1));
 c3=rnd();
 c4=sqrt(1.-c3*c3);
 c5=2.*PI*rnd();
 c6=sin(c5);
 c7=cos(c5);
 KX=c2*c3;
 KY=c2*c4*c6;
 KZ=c2*c4*c7;
 IV=iv;
 TS=t-log(rnd())/GM[i_dom[ii][j]];
 if(edge==2) KY=-KY;
 if(edge==1) KX=-KX;
}

// =================================================
