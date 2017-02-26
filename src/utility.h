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

#include "particle.h"
#include "mesh.h"


struct {
    int BOTTOM;
    int RIGHT;
    int TOP;
    int LEFT;
} direction_t = {.BOTTOM=0, .RIGHT=1, .TOP=2, .LEFT=3};


struct {
    int INSULATOR;
    int SCHOTTKY;
    int OHMIC;
    int VACUUM;
} boundary_t = {.INSULATOR=0, .SCHOTTKY=1, .OHMIC=2, .VACUUM=3};


inline int mc_is_boundary_insulator(int direction, int index) {
    return g_mesh->edges[direction][index].boundary == boundary_t.INSULATOR;
}


inline int mc_is_boundary_schottky(int direction, int index) {
    return g_mesh->edges[direction][index].boundary == boundary_t.SCHOTTKY;
}


inline int mc_is_boundary_ohmic(int direction, int index) {
    return g_mesh->edges[direction][index].boundary == boundary_t.OHMIC;
}


inline int mc_is_boundary_vacuum(int direction, int index) {
    return g_mesh->edges[direction][index].boundary == boundary_t.VACUUM;
}


inline int mc_is_boundary_contact(int direction, int index) {
    return mc_is_boundary_schottky(direction, index)
        || mc_is_boundary_ohmic(direction, index);
}


char* mc_band_model_name(int model) {
    switch(model) {
        case PARABOLIC: return "Parabolic";
        case KANE: return "Kane";
        case FULL: return "Full";
        default: return "Unknown Band Model";
    }
}


char *trim(char *s) {
    char *start = s;
    char *end = NULL;

    // skip spaces at start
    while(*start && isspace(*start)) {
        ++start;
    }

    char *i = start;
    // iterate over the rest remebering last non-whitespace
    while(*i) {
        if(!isspace(*(i++))) {
            end = i;
        }
    }

    // white the terminating zero after last non-whitespace
    *end = 0;

    return start;
}


#endif
