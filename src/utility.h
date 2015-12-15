/* utility.h -- This file is part of Archimedes release 1.2.0.
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


// ######################################################
// Created on 21 oct.2015, J. Marini
// Last modif. : 21 oct.2015, J. Marini
// ######################################################

#ifndef ARCHIMEDES_UTILITY_H
#define ARCHIMEDES_UTILITY_H


inline particle_t creation(int, real, int);


struct {
    const int BOTTOM;
    const int RIGHT;
    const int TOP;
    const int LEFT;
} direction_t = { .BOTTOM=0, .RIGHT=1, .TOP=2, .LEFT=3 };
#define NUMBER_DIRECTIONS 4


struct {
    const int VERTICAL;
    const int HORIZONTAL;
} axis_t = { .VERTICAL=0, .HORIZONTAL=1 };


struct {
    const int INSULATOR;
    const int SCHOTTKY;
    const int OHMIC;
} boundary_t = { .INSULATOR=0, .SCHOTTKY=1, .OHMIC=2 };


inline int mc_boundary_type(int direction, int index) {
    return EDGE[direction][index][0];
}


inline int mc_is_boundary_insulator(int direction, int index) {
    return mc_boundary_type(direction, index) == boundary_t.INSULATOR;
}


inline int mc_is_boundary_contact(int direction, int index) {
    return mc_boundary_type(direction, index) == boundary_t.SCHOTTKY
        || mc_boundary_type(direction, index) == boundary_t.OHMIC;
}


inline int mc_parallel_axis(int direction) {
    return direction % 2;
}


inline int mc_perpendicular_axis(int direction) {
    return !mc_parallel_axis(direction);
}


inline int mc_is_direction_horizontal(int direction) {
    return mc_parallel_axis(direction) == axis_t.HORIZONTAL;
}


inline int mc_is_direction_vertical(int direction) {
    return mc_parallel_axis(direction) == axis_t.VERTICAL;
}


inline int mc_is_index_in_bounds_direction(int i, int j, int direction) {
    if(direction == direction_t.BOTTOM) { return j  > 1; }
    if(direction == direction_t.RIGHT)  { return i < nx + 1; }
    if(direction == direction_t.TOP)    { return j < ny + 1; }
    if(direction == direction_t.LEFT)   { return i > 1; }
    return 0;
}


inline int mc_is_index_in_bounds_axis(int i, int j, int axis) {
    if(axis == axis_t.VERTICAL) {
        return mc_is_index_in_bounds_direction(i, j, direction_t.TOP)
            && mc_is_index_in_bounds_direction(i, j, direction_t.BOTTOM);
    }
    if(axis == axis_t.HORIZONTAL) {
        return mc_is_index_in_bounds_direction(i, j, direction_t.LEFT)
            && mc_is_index_in_bounds_direction(i, j, direction_t.RIGHT);
    }
    return 0;
}


void mc_check_particle_leaving(particle_t *particle, int i, int j, int direction, int npt[][NUMBER_DIRECTIONS]) {
    if(!mc_does_particle_exist(particle)) { return; }

    int axis = mc_parallel_axis(direction);
    int index = (axis == axis_t.VERTICAL ? i : j);

    if(!mc_is_index_in_bounds_direction(i, j, direction)
            && mc_is_boundary_contact(direction, index)) {
        mc_remove_particle(particle);
        if(npt[index][direction] < (NP1/2) && mc_is_index_in_bounds_axis(i, j, !axis)) {
            npt[index][direction]++;
            particle->valley = 1;
        }
        else if(npt[index][direction] < (NP1/4) && !mc_is_index_in_bounds_axis(i, j, !axis)) {
            npt[index][direction]++;
            particle->valley = 1;
        }
    }
}


inline void mc_particle_coords(particle_t *particle, int *i, int *j) {
    // uses globabl variables dx & dy
    *i = (int)(particle->x / dx) + 1;
    if(*i < 1) { *i = 1; }
    if(*i > nx ) { *i = 1; }

    *j = (int)(particle->y / dy) + 1;
    if(*j < 1) { *j = 1; }
    if(*j > ny ) { *j = 1; }
}


#endif
