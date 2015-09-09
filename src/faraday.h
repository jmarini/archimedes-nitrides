/* faraday.h -- This file is part of Archimedes release 0.0.8.
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method and a simplfied
   version of the MEP model (Maximum Entropy Principle model)
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
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA.  */


// ######################################################
// Created on 18 Aug.2007, Siracusa, J.M.Sellier
// Last modif. : 18 Apr.2007, Siracusa, J.M.Sellier
// ######################################################

// Computation of the electrostatic potential,
// i.e. resolution of the 2D Faraday equation.

void
Faraday(void)
{
 int i,j;
// Lax-Friedrichs-Sellier method (first-order in space)
 for(i=2;i<=nx;i++)
  for(j=2;j<=ny;j++)
   B[i][j]=0.25*(B[i+1][j]+B[i][j+1]+B[i-1][j]+B[i][j-1])
          -DT*(0.5*(E[i+1][j][1]-E[i-1][j][1])/dx
	      -0.5*(E[i][j+1][0]-E[i][j-1][0])/dy);
}
