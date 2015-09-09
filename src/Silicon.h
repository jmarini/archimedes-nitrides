/* Silicon.h -- This file is part of Archimedes release 1.2.0
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
// Created on 22 aug.2011, Carry le Rouet, France J.M.Sellier
// Last modif. : 31 Aug.2011, Carry le Rouet, J.M.Sellier
// ######################################################

// conduction band calculated using the model : sp3s* (see Vogl)
// data interpolated by a 10-th order polynomial

// Parabolic conduction band
// X-valley [1 0 0] direction
CB_FULL[SILICON][0]=-4.1432e+21;
CB_FULL[SILICON][1]=+2.4419e+20;
CB_FULL[SILICON][2]=-6.1106e+18;
CB_FULL[SILICON][3]=+8.4609e+16;
CB_FULL[SILICON][4]=-7.0825e+14;
CB_FULL[SILICON][5]=+3.6800e+12;
CB_FULL[SILICON][6]=-1.1741e+10;
CB_FULL[SILICON][7]=+2.1981e+07;
CB_FULL[SILICON][8]=+4.6785e+06;
CB_FULL[SILICON][9]=+9.3312e+00;
CB_FULL[SILICON][10]=0.0; //-8.8397e-04; 
