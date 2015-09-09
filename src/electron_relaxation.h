/* electron_relaxation.c -- This file is part of GNU/Archimedes 0.1.0
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
// Date of Creation : 12 Aug.2004, Siracusa, Italy, Jean Michel Sellier
// Last Revision : 09 Sep.2007, Siracusa, Italy, Jean Michel Sellier.
// System of Measure : M.K.S.C. System
// =================================================================

// This part is for the simulation of the relaxation part
// of the Hybrid MEP model. This is a parabolic model
// which is based on the Maximum Entropy Principle.
// For more informations and references see the manual.

void electron_relaxation_step(void)
{
 int ND=2;
 register int i,j,c;
 real taup,t;
 real ktaup,ktauw;
 
// Source Term for the Hybrid MEP model
// ====================================
// We take into account all the relevant scattering
// effects in Silicon. For this purpose, we
// use the same relaxation scattering model used
// in BBW (Baccarani et al.) model.

// This is a simple explicit Euler step.
// DT/=2.;
 for(c=1;c<=50;c++)
  for(i=ND;i<=nx+ND;i++)
   for(j=ND;j<=ny+ND;j++){
     t=(2./3.)*(u2d[i][j][4]/u2d[i][j][1])/KB;
     ktaup=M*MSTAR[i_dom[i][j]][1]*MIU0*TL/Q;
     ktauw=MIU0*TL*KB/(Q*VS*VS);
     taup=ktaup/t;
     u2d[i][j][4]+=-DT/50.*(Q*(u2d[i][j][2]*E[i-1][j-1][0]
                  +u2d[i][j][3]*E[i-1][j-1][1])
                  +u2d[i][j][1]*1.5*KB*(t-TL)/tauwi(1.5*KB*t));
     u2d[i][j][2]+=-DT/50.*(u2d[i][j][1]*Q*E[i-1][j-1][0]/(M*MSTAR[i_dom[i][j]][1])
                  +u2d[i][j][2]/taup);
     u2d[i][j][3]+=-DT/50.*(u2d[i][j][1]*Q*E[i-1][j-1][1]/(M*MSTAR[i_dom[i][j]][1])
                  +u2d[i][j][3]/taup);
   }
// DT*=2.;
}
// ================================================
