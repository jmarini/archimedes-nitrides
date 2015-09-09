/* MEP_interpolation.c -- This file is part of GNU/Archimedes 0.1.0
   This code is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method and MEP model
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
// Date of Creation : 05 May 2004, Siracusa, Italy, Jean Michel Sellier
// Last Revision : 09 Sep.2007, Siracusa, Italy, Jean Michel Sellier.
// System of Measure : M.K.S.C. System
// =================================================================

// In this file, one finds all the interpolation function
// for the non-parabolic (Kane) MEP model for electrons.

void
MEP_coefficients(void)
{
// All the coefficients below are needed in order
// to approximate all the scatterings and fluxes
// present in MEP model.

// c11 coefficients
 c11[1]=7.4985705e-016;
 c11[2]=-3.6122279e-015;
 c11[3]=3.9232487e-014;
 c11[4]=9.1529095e-014;
 c11[5]=-4.4414692e-015;
 c11[6]=-1.7323853e-016;
// c12 coefficients
 c12[1]=-1.3646041e-012;
 c12[2]=8.1197009e-012;
 c12[3]=-1.7070415e-011;
 c12[4]=1.5010866e-011;
 c12[5]=-5.2241407e-012;
 c12[6]=4.2032095e-013;
// c21 coefficients
 c21[1]=2.1940086e12;
 c21[2]=-1.4192452e13;
 c21[3]=1.5400896e14;
 c21[4]=2.2927873e14;
 c21[5]=2.3734694e12;
 c21[6]=1.2666909e11;
// c22 coefficients
 c22[1]=-1.1216571e13;
 c22[2]=6.4708042e13;
 c22[3]=-1.2457850e14;
 c22[4]=-1.1229164e14;
 c22[5]=-3.3463063e14;
 c22[6]=-4.0441071e12;
// U coefficients
 u[1]=1.2031722e-2;
 u[2]=-8.2472023e-2;
 u[3]=2.2075893e-1;
 u[4]=-3.1204242e-1;
 u[5]=6.3963512e-1;
 u[6]=7.4669684e-4;
// F coefficients
 f[1]=-2.2156802e-2;
 f[2]=1.4648397e-1;
 f[3]=-3.6560368e-1;
 f[4]=4.3888833e-1;
 f[5]=6.7589719e-2;
 f[6]=-1.7639096e-3;
// G coefficients
 g[1]=8.5792683e-2;
 g[2]=-5.7152259e-1;
 g[3]=1.4517340;
 g[4]=-1.8382828;
 g[5]=1.4079332;
 g[6]=6.7206574e-3;
// CW coefficients
 cw[1]=-7.2429491e10;
 cw[2]=5.1342003e11;
 cw[3]=-1.3655677e12;
 cw[4]=3.6907365e11;
 cw[5]=-3.0788470e12;
 cw[6]=9.9708156e10;
}

// ****************************

real
c11i(real t)
{
// 05 May 2004, Siracusa, J.M. Sellier
// input : energy in Joule
// output : c11 in kg/sec
 real W=t/Q; // energy in eV
 return(c11[1]*pow(W,5.)+c11[2]*pow(W,4.)
	 +c11[3]*pow(W,3.)+c11[4]*pow(W,2.)+c11[5]*W+c11[6])*1.e-3; 
}

// ****************************

real
c12i(real t)
{
// 05 May 2004, Siracusa, J.M. Sellier
// input : energy in Joule
// output : c12 in kg/(sec*Joule)
 real W=t/Q; // energy in eV
 return(c12[1]*pow(W,5.)+c12[2]*pow(W,4.)
	 +c12[3]*pow(W,3.)+c12[4]*pow(W,2.)+c12[5]*W+c12[6])*1.e-3/Q; 
}

// ****************************

real
c21i(real t)
{
// 05 May 2004, Siracusa, J.M. Sellier
// input : energy in Joule
// output : c21 in Joule/sec
 double W=t/Q; // energy in eV
 return(c21[1]*pow(W,5.)+c21[2]*pow(W,4.)
	 +c21[3]*pow(W,3.)+c21[4]*pow(W,2.)+c21[5]*W+c21[6])*Q; 
}

// *****************************

real
c22i(real t)
{
// 05 May 2004, Siracusa, J.M. Sellier
// input : energy in Joule
// output : c22 in 1/sec
 real W=t/Q; // energy in eV
 return c22[1]*pow(W,5.)+c22[2]*pow(W,4.)
	 +c22[3]*pow(W,3.)+c22[4]*pow(W,2.)+c22[5]*W+c22[6]; 
}

// *****************************

real
cwi(real t)
{
// Non-parabolic MEP production of the energy flux balance equation
// 06 May 2004, Siracusa, J.M. Sellier
// input : energy in Joule
// output : cw in Joule/sec
 real W=t/Q; // energy in eV
 return (cw[1]*pow(W,5.)+cw[2]*pow(W,4.)
	 +cw[3]*pow(W,3.)+cw[4]*pow(W,2.)+cw[5]*W+cw[6])*Q; 
}

// *****************************

real
tauwi(real t)
{
// Non-parabolic MEP energy relaxation time
// 06 May 2004, Siracusa, J.M. Sellier
// input : energy in Joule
// output : tauw in sec
 real W=t; // energy in Joule
 real W0=1.5*KB*TL; // lattice energy in Joule
 if(-(W-W0)/cwi(W)!=0.) return -(W-W0)/cwi(W);
 else return TAUW;
}

// ******************************

real
Ui(real t)
{
// Non-parabolic MEP flux of crystal momentum
// 05 May 2004, Siracusa, J.M. Sellier
// input : energy in Joule
// output : U in Joule
 real W=t/Q; // energy in eV
 return(u[1]*pow(W,5.)+u[2]*pow(W,4.)
	 +u[3]*pow(W,3.)+u[4]*pow(W,2.)+u[5]*W+u[6])*Q; 
}

// *******************************

real
Fi(real t)
{
// Non-parabolic MEP flux of the energy flux
// 05 May 2004, Siracusa, J.M. Sellier
// input : energy in Joule
// output : F in (Joule^2)/kg
 real W=t/Q; // energy in eV
 return(f[1]*pow(W,5.)+f[2]*pow(W,4.)
	 +f[3]*pow(W,3.)+f[4]*pow(W,2.)+f[5]*W+f[6])*Q*Q/MSTAR[SILICON][1]; 
}

// ********************************

real
Gi(real t)
{
// 05 May 2004, Siracusa, J.M. Sellier
// input : energy in Joule
// output : G in Joule/kg
 real W=t/Q; // energy in eV
 return(g[1]*pow(W,5.)+g[2]*pow(W,4.)
	 +g[3]*pow(W,3.)+g[4]*pow(W,2.)+g[5]*W+g[6])*Q/MSTAR[SILICON][1];
}

// *********************************
