/* vec.h -- This file is part of Archimedes release 1.2.0.
   Archimedes is a simulator for Submicron 2D III-V semiconductor
   Devices. It implements the Monte Carlo method
   for the simulation of the semiclassical Boltzmann equation for both
   electrons and holes. It includes some quantum effects by means
   of effective potential method. It is also able to simulate applied
   magnetic fields along with self consistent Faraday equation.

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


#ifndef ARCHIMEDES_VEC_H
#define ARCHIMEDES_VEC_H


// struct to hold real xy data
typedef struct Vec2 {
    double x;
    double y;
} Vec2;


// struct to hold real xyz data
typedef struct Vec3 {
    double x;
    double y;
    double z;
} Vec3;


// struct to hold node index data
typedef struct Index {
    int i;
    int j;
} Index;


// struct to hold size information
typedef struct Dimensions {
    double width;
    double height;
} Dimensions;

#endif
