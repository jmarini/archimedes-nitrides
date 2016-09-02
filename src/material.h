/* particle.h -- This file is part of Archimedes release 1.2.0.
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
// Created on 21 apr.2016, J. Marini
// Last modif. : 21 apr.2016, J. Marini
// ######################################################

#ifndef ARCHIMEDES_MATERIAL_H
#define ARCHIMEDES_MATERIAL_H

#include "global_defines.h"


typedef struct {
    real emin;
    real mstar;
    real alpha;
} Band_Info;


typedef struct {
    int num_valleys;
    real Eg;
    Band_Info cb[6];
    Band_Info vb[6];

    int  zscatter[6][6];

    real eps_static;
    real eps_hf;
    real hwo[6];
    real dtk[6];
    int  zf[6];
    real rho;
    real da;
    real ul;
    real lattice_const;
    real kav;

} Material;


extern Material g_materials[NOAMTIA];


#endif
