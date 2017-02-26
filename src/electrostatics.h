/* electrostatics.h -- This file is part of GNU archimedes

   Archimedes is a simulator for Submicron and Nanoscaled
   2D III-V Semiconductor Devices.

   Copyright (C) 2004-2011 Jean Michel D. Sellier
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


#ifndef ARCHIMEDES_ELECTROSTATICS_H
#define ARCHIMEDES_ELECTROSTATICS_H


#include "mesh.h"


int electrostatics(Mesh *mesh);
int poisson(Mesh *mesh);
int faraday(Mesh *mesh);

int calculate_potential(Mesh *mesh);
int quantum_effective_potential(Mesh *mesh);

int electric_field(Mesh *mesh);
int poisson_boundary_conditions(Mesh *mesh);

int magnetic_field(Mesh *mesh);
int faraday_boundary_conditions(Mesh *mesh);


#endif
