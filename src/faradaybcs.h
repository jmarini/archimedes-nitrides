/* faradaybcs.h -- This file is part of Archimedes release 0.0.8.
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method and a simplified
   MEP model for the simulation of the semiclassical Boltzmann
   equation for both electrons and holes. It also includes the
   quantum effects by means of effective potential method.
   It is now able to simulate applied magnetic fields along with
   self consistent Faraday equation.

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
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA.  */


// ######################################################
// Created on 18 Aug.2007, Siracusa, J.M.Sellier
// Last modif. : 18 Aug.2007, Siracusa, J.M.Sellier
// ######################################################

// Boundary Conditions for the Faraday equation
// For more informations about this equation
// see the manual of GNU Archimedes release 0.0.8.

void
FaradayBCs(void)
{
 int i,j;
// Bottom Edge
// ===========
   for(i=1;i<=nx+1;i++){
     B[i][0]=B[i][3];
     B[i][1]=B[i][2];
   }
// Left Edge
// =========
   for(j=1;j<=ny+1;j++){
     B[0][j]=B[3][j];
     B[1][j]=B[2][j];
   }
// Right Edge
// ==========
   for(j=1;j<=ny+1;j++){
     B[nx+1][j]=B[nx-1][j];
     B[nx+2][j]=B[nx][j];
   }
// Upper Edge
// ==========   
   for(i=1;i<=nx+1;i++){
     B[i][ny+1]=B[i][ny];
     B[i][ny+2]=B[i][ny-1];
   }
}
