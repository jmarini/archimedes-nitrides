/* hole_relaxation.c -- This file is part of GNU/Archimedes 0.0.3
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
// File Name : electron_relaxation.c
// Version   : release 0.0.3
// Date of Creation : 30 Aug.2004, Siracusa, Italy, Jean Michel Sellier
// Last Revision : 05 Mar.2005, Siracusa, Italy, Jean Michel Sellier.
// System of Measure : M.K.S.C. System
// =================================================================

// This part is for the simulation of the relaxation part
// of the Hybrid MEP model for holes. This is a parabolic model
// which is based on the Maximum Entropy Principle.
// For more informations and references see the manual.


void Relaxation_Step_Hole(void)
{
 int ND=2;
 register int i,j,c;
 real taup,tauw,ktaup,ktauw,t;

// Euler Relaxation step
for(c=1;c<=50;c++)
 for(i=ND;i<=nx+ND;i++)
   for(j=ND;j<=ny+ND;j++){
     t=(2./3.)*(h2d[i][j][4]/h2d[i][j][1])/KB;
     ktaup=M*mstarhole*MIU0hole*TL/Q;
     ktauw=MIU0hole*TL*KB/(Q*VShole*VShole);
     taup=ktaup/t;
     tauw=0.5*ktaup/t+1.5*ktauw*t/(t+TL);
     h2d[i][j][4]+=-DT/50.*(-Q*(h2d[i][j][2]*E[i-1][j-1][0]
                  +h2d[i][j][3]*E[i-1][j-1][1])
                  +h2d[i][j][1]*1.5*KB*(t-TL)/tauw);
     h2d[i][j][2]+=-DT/50.*(-h2d[i][j][1]*Q*E[i-1][j-1][0]/(M*mstarhole)
                  +h2d[i][j][2]/taup);
     h2d[i][j][3]+=-DT/50.*(-h2d[i][j][1]*Q*E[i-1][j-1][1]/(M*mstarhole)
                  +h2d[i][j][3]/taup);
   }
}
// ================================================
