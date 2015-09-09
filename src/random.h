/* random.h -- This file is part of Archimedes release 0.0.1.
   Archimedes is a simulator for Submicron 2D Silicon
   Devices. It implements the Monte Carlo method
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


// ######################################################
// Created on 05 sep.2004, Siracusa, J.M.Sellier
// Last modif. : 22 Sep.2004, Siracusa, J.M.Sellier
// ######################################################

// a simple (but well working...) generator of random numbers

inline
real
rnd(void){
  ISEED = fmod(1027.*ISEED, 1048576.);
  return (ISEED/1048576.);
}

// =========================================
